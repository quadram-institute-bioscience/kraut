from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pandas as pd

from kraut.models.kraken_data import KrakenNode, KrakenReport


FEATURE_RANKS = ("D", "K", "P", "C", "O", "F", "G", "S")
TAXONOMY_COLUMNS = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAXONOMY_PREFIXES = {
    "Domain": "d__",
    "Phylum": "p__",
    "Class": "c__",
    "Order": "o__",
    "Family": "f__",
    "Genus": "g__",
    "Species": "s__",
}
RANK_TO_TAXONOMY_COLUMN = {
    "D": "Domain",
    "K": "Domain",
    "P": "Phylum",
    "C": "Class",
    "O": "Order",
    "F": "Family",
    "G": "Genus",
    "S": "Species",
}


@dataclass
class MAExportResult:
    counts_file: Path
    taxonomy_file: Path
    metadata_file: Path
    tree_file: Path
    warnings: list[str]


@dataclass
class TaxonInfo:
    tax_id: int
    name: str
    rank_code: str
    parent_id: int | None
    counts: list[int]
    children: list[int] = field(default_factory=list)
    feature_id: str | None = None


@dataclass
class TaxonIndex:
    taxa: dict[int, TaxonInfo]
    feature_order: list[int]
    sample_names: list[str]


def export_microbiome_analyst(
    input_files: Iterable[Path],
    outdir: Path,
    metadata_file: Path | None = None,
    metadata_sample_col: str | None = None,
    pseudo_col: str = "Random_label",
    metric: str = "LVL",
) -> MAExportResult:
    """Export Kraken reports as MicrobiomeAnalyst CSV inputs."""
    paths = [Path(path) for path in input_files]
    if not paths:
        raise ValueError("At least one input Kraken report is required")

    for path in paths:
        if not path.exists():
            raise ValueError(f"Input file does not exist: {path}")

    sample_names = [_sample_name(path) for path in paths]
    _validate_unique_sample_names(sample_names)
    count_attr = _count_attribute(metric)

    reports = [KrakenReport.from_file(str(path)) for path in paths]
    taxon_index = _collect_taxa(reports, sample_names, count_attr)
    _assign_feature_ids(taxon_index)

    if not taxon_index.feature_order:
        raise ValueError("No exportable taxa found in the input reports")

    warnings: list[str] = []
    counts_df = _counts_dataframe(taxon_index)
    taxonomy_df = _taxonomy_dataframe(taxon_index)
    metadata_df, metadata_warnings = _metadata_dataframe(
        sample_names,
        metadata_file,
        metadata_sample_col,
        pseudo_col,
    )
    warnings.extend(metadata_warnings)
    tree = _newick_tree(taxon_index)

    outdir.mkdir(parents=True, exist_ok=True)
    counts_file = outdir / "counts.csv"
    taxonomy_file = outdir / "taxonomy.csv"
    metadata_output_file = outdir / "metadata.csv"
    tree_file = outdir / "tree.nwk"

    counts_df.to_csv(counts_file, index=False)
    taxonomy_df.to_csv(taxonomy_file, index=False)
    metadata_df.to_csv(metadata_output_file, index=False)
    tree_file.write_text(tree + "\n")

    return MAExportResult(
        counts_file=counts_file,
        taxonomy_file=taxonomy_file,
        metadata_file=metadata_output_file,
        tree_file=tree_file,
        warnings=warnings,
    )


def generate_pseudo_labels(sample_count: int) -> tuple[list[str], list[str]]:
    """Generate deterministic group labels without singleton groups when possible."""
    if sample_count < 0:
        raise ValueError("sample_count must be >= 0")
    if sample_count == 0:
        return [], []
    if sample_count == 1:
        return ["A"], ["Only one sample was exported; pseudo labels contain a singleton"]
    if sample_count < 4:
        return ["A"] * sample_count, []
    return ["A" if idx % 2 == 0 else "B" for idx in range(sample_count)], []


def _sample_name(path: Path) -> str:
    sample_name = path.stem
    if sample_name.endswith(".krep"):
        sample_name = Path(sample_name).stem
    return sample_name


def _validate_unique_sample_names(sample_names: list[str]) -> None:
    seen = set()
    duplicates = []
    for sample_name in sample_names:
        if sample_name in seen and sample_name not in duplicates:
            duplicates.append(sample_name)
        seen.add(sample_name)

    if duplicates:
        names = ", ".join(duplicates)
        raise ValueError(f"Duplicate sample name(s) from input files: {names}")


def _count_attribute(metric: str) -> str:
    metric = metric.upper()
    if metric == "LVL":
        return "taxon_counts"
    if metric == "TOT":
        return "clade_counts"
    raise ValueError("--metric must be one of: LVL, TOT")


def _collect_taxa(
    reports: list[KrakenReport],
    sample_names: list[str],
    count_attr: str,
) -> TaxonIndex:
    taxon_index = TaxonIndex(taxa={}, feature_order=[], sample_names=sample_names)
    sample_count = len(sample_names)

    for sample_idx, report in enumerate(reports):
        if report.root is None:
            continue

        for node in _iter_tree(report.root):
            parent_id = node.parent.tax_id if node.parent is not None else None
            info = _add_or_validate_taxon(taxon_index, node, parent_id, sample_count)
            if parent_id is not None:
                _add_child(taxon_index, parent_id, node.tax_id)

            if _is_feature_rank(node.rank_code):
                if node.tax_id not in taxon_index.feature_order:
                    taxon_index.feature_order.append(node.tax_id)
                info.counts[sample_idx] = int(getattr(node, count_attr))

    return taxon_index


def _iter_tree(node: KrakenNode):
    yield node
    for child in node.children:
        yield from _iter_tree(child)


def _add_or_validate_taxon(
    taxon_index: TaxonIndex,
    node: KrakenNode,
    parent_id: int | None,
    sample_count: int,
) -> TaxonInfo:
    rank_code = node.rank_code.upper()
    existing = taxon_index.taxa.get(node.tax_id)
    if existing is None:
        info = TaxonInfo(
            tax_id=node.tax_id,
            name=node.name,
            rank_code=rank_code,
            parent_id=parent_id,
            counts=[0] * sample_count,
        )
        taxon_index.taxa[node.tax_id] = info
        return info

    conflicts = []
    if existing.name != node.name:
        conflicts.append(f"name '{existing.name}' vs '{node.name}'")
    if existing.rank_code != rank_code:
        conflicts.append(f"rank '{existing.rank_code}' vs '{rank_code}'")
    if existing.parent_id != parent_id:
        conflicts.append(f"parent '{existing.parent_id}' vs '{parent_id}'")
    if conflicts:
        details = "; ".join(conflicts)
        raise ValueError(f"TaxID {node.tax_id} has conflicting taxonomy: {details}")

    return existing


def _add_child(taxon_index: TaxonIndex, parent_id: int, child_id: int) -> None:
    parent = taxon_index.taxa.get(parent_id)
    if parent is None:
        raise ValueError(f"TaxID {child_id} has missing parent TaxID {parent_id}")
    if child_id not in parent.children:
        parent.children.append(child_id)


def _is_feature_rank(rank_code: str) -> bool:
    return rank_code.upper() in FEATURE_RANKS


def _assign_feature_ids(taxon_index: TaxonIndex) -> None:
    for idx, tax_id in enumerate(taxon_index.feature_order, start=1):
        taxon_index.taxa[tax_id].feature_id = f"Feat_{idx}"


def _feature_taxa(taxon_index: TaxonIndex) -> list[TaxonInfo]:
    return [taxon_index.taxa[tax_id] for tax_id in taxon_index.feature_order]


def _counts_dataframe(taxon_index: TaxonIndex) -> pd.DataFrame:
    rows = []
    for taxon in _feature_taxa(taxon_index):
        row = {"#NAME": taxon.feature_id}
        for sample_name, count in zip(taxon_index.sample_names, taxon.counts):
            row[sample_name] = count
        rows.append(row)
    return pd.DataFrame(rows, columns=["#NAME"] + taxon_index.sample_names)


def _taxonomy_dataframe(taxon_index: TaxonIndex) -> pd.DataFrame:
    rows = []
    for taxon in _feature_taxa(taxon_index):
        values = _taxonomy_values(taxon, taxon_index)
        row = {"#TAXONOMY": taxon.feature_id}
        row.update(values)
        rows.append(row)
    return pd.DataFrame(rows, columns=["#TAXONOMY"] + list(TAXONOMY_COLUMNS))


def _taxonomy_values(taxon: TaxonInfo, taxon_index: TaxonIndex) -> dict[str, str]:
    values = {column: TAXONOMY_PREFIXES[column] for column in TAXONOMY_COLUMNS}
    domain_from_d = None
    domain_from_k = None

    for lineage_taxon in _lineage(taxon, taxon_index):
        rank = lineage_taxon.rank_code
        column = RANK_TO_TAXONOMY_COLUMN.get(rank)
        if column is None:
            continue

        if rank == "D":
            domain_from_d = lineage_taxon.name
        elif rank == "K":
            domain_from_k = lineage_taxon.name
        else:
            values[column] = f"{TAXONOMY_PREFIXES[column]}{lineage_taxon.name}"

    if taxon.rank_code == "K":
        domain_name = taxon.name
    else:
        domain_name = domain_from_d if domain_from_d is not None else domain_from_k
    if domain_name is not None:
        values["Domain"] = f"{TAXONOMY_PREFIXES['Domain']}{domain_name}"

    return values


def _lineage(taxon: TaxonInfo, taxon_index: TaxonIndex) -> list[TaxonInfo]:
    lineage = []
    current: TaxonInfo | None = taxon
    seen = set()
    while current is not None:
        if current.tax_id in seen:
            raise ValueError(f"Taxonomy cycle detected at TaxID {current.tax_id}")
        seen.add(current.tax_id)
        lineage.append(current)
        if current.parent_id is None:
            current = None
        else:
            current = taxon_index.taxa.get(current.parent_id)
    return list(reversed(lineage))


def _metadata_dataframe(
    sample_names: list[str],
    metadata_file: Path | None,
    metadata_sample_col: str | None,
    pseudo_col: str,
) -> tuple[pd.DataFrame, list[str]]:
    if pseudo_col == "#NAME":
        raise ValueError("--pseudo-col must not be #NAME")

    if metadata_file is None:
        metadata_df = pd.DataFrame({"#NAME": sample_names})
    else:
        if metadata_sample_col is None:
            raise ValueError("--metadata-sample-col is required when --metadata is supplied")
        metadata_df = _read_metadata(metadata_file, metadata_sample_col, sample_names)

    labels, warnings = generate_pseudo_labels(len(sample_names))
    metadata_df[pseudo_col] = labels
    columns = ["#NAME"] + [column for column in metadata_df.columns if column != "#NAME"]
    return metadata_df[columns], warnings


def _read_metadata(
    metadata_file: Path,
    metadata_sample_col: str,
    sample_names: list[str],
) -> pd.DataFrame:
    if not metadata_file.exists():
        raise ValueError(f"Metadata file does not exist: {metadata_file}")

    try:
        metadata = pd.read_csv(
            metadata_file,
            sep=None,
            engine="python",
            dtype=str,
            keep_default_na=False,
        )
    except Exception as exc:
        raise ValueError(f"Could not read metadata file '{metadata_file}': {exc}") from exc

    if metadata_sample_col not in metadata.columns:
        raise ValueError(f"Metadata sample column not found: {metadata_sample_col}")

    if "#NAME" in metadata.columns and metadata_sample_col != "#NAME":
        raise ValueError("Metadata contains a #NAME column that conflicts with the output header")

    metadata = metadata.copy()
    metadata[metadata_sample_col] = metadata[metadata_sample_col].astype(str)
    duplicate_rows = metadata[
        metadata[metadata_sample_col].isin(sample_names)
        & metadata[metadata_sample_col].duplicated(keep=False)
    ]
    if not duplicate_rows.empty:
        duplicates = sorted(duplicate_rows[metadata_sample_col].unique())
        names = ", ".join(duplicates)
        raise ValueError(f"Metadata contains duplicate row(s) for sample(s): {names}")

    metadata_samples = set(metadata[metadata_sample_col])
    missing = [sample_name for sample_name in sample_names if sample_name not in metadata_samples]
    if missing:
        names = ", ".join(missing)
        raise ValueError(f"Metadata is missing input sample(s): {names}")

    selected = (
        metadata.set_index(metadata_sample_col, drop=False)
        .loc[sample_names]
        .reset_index(drop=True)
    )
    other_columns = [column for column in metadata.columns if column != metadata_sample_col]
    selected = selected[[metadata_sample_col] + other_columns].copy()
    selected.rename(columns={metadata_sample_col: "#NAME"}, inplace=True)
    return selected


def _newick_tree(taxon_index: TaxonIndex) -> str:
    has_feature_cache: dict[int, bool] = {}

    def has_feature(tax_id: int) -> bool:
        if tax_id in has_feature_cache:
            return has_feature_cache[tax_id]
        taxon = taxon_index.taxa[tax_id]
        result = taxon.feature_id is not None or any(has_feature(child) for child in taxon.children)
        has_feature_cache[tax_id] = result
        return result

    root_ids = _newick_root_ids(taxon_index, has_feature)
    rendered = [_render_newick_subtree(taxon_index, tax_id, has_feature) for tax_id in root_ids]
    rendered = [subtree for subtree in rendered if subtree]

    if not rendered:
        raise ValueError("Cannot build tree without exported features")
    if len(rendered) == 1:
        return f"{rendered[0]};"
    return f"({','.join(rendered)})Root;"


def _newick_root_ids(taxon_index: TaxonIndex, has_feature) -> list[int]:
    if 1 in taxon_index.taxa and has_feature(1):
        return [1]

    roots = []
    for tax_id, taxon in taxon_index.taxa.items():
        if taxon.parent_id not in taxon_index.taxa and has_feature(tax_id):
            roots.append(tax_id)
    return roots


def _render_newick_subtree(
    taxon_index: TaxonIndex,
    tax_id: int,
    has_feature,
) -> str:
    taxon = taxon_index.taxa[tax_id]
    child_parts = [
        _render_newick_subtree(taxon_index, child_id, has_feature)
        for child_id in taxon.children
        if has_feature(child_id)
    ]
    child_parts = [part for part in child_parts if part]

    if taxon.feature_id is not None:
        feature_leaf = f"{taxon.feature_id}:1"
        if not child_parts:
            return feature_leaf
        child_parts.insert(0, feature_leaf)

    if child_parts:
        return f"({','.join(child_parts)})Tax_{taxon.tax_id}:1"

    return ""
