from collections import OrderedDict
from pathlib import Path
from typing import List, Optional

import typer

RANK_LETTERS = ("k", "p", "c", "o", "f", "g", "s", "t")
RANK_WORDS = {
    "kingdom": "k",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
    "sgb": "t",
    "strain": "t",
}
UNCLASSIFIED = "UNCLASSIFIED"


def run(
    input_files: List[Path] = typer.Argument(..., help="Input MetaPhlAn profile files"),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output TSV file (default: stdout)",
    ),
    level: str = typer.Option(
        "S",
        "--level",
        "-l",
        help="Taxonomic level: k/p/c/o/f/g/s/t or full rank name",
    ),
    keep_taxid: bool = typer.Option(
        False,
        "--keep-taxid",
        help="Add an NCBI_tax_id column to the output",
    ),
    short_names: bool = typer.Option(
        False,
        "--short-names",
        help="Use terminal clade names only (genus plus species for species level)",
    ),
    drop_unclassified: bool = typer.Option(
        False,
        "--drop-unclassified",
        help="Drop the UNCLASSIFIED row",
    ),
    normalise: bool = typer.Option(
        False,
        "--normalise",
        help="Scale each sample column so retained rows sum to 100",
    ),
) -> None:
    """
    Merge multiple MetaPhlAn profiles into a single abundance table.
    """
    try:
        result = make_mpa_table(
            input_files=input_files,
            level=level,
            keep_taxid=keep_taxid,
            short_names=short_names,
            drop_unclassified=drop_unclassified,
            normalise=normalise,
        )
    except (OSError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if output_file:
        output_file.write_text(result.table)
    else:
        typer.echo(result.table, nl=False)

    typer.echo(
        (
            f"Merged {len(result.sample_names)} samples, {result.taxon_count} taxa "
            f"at level '{result.level}'."
        ),
        err=True,
    )


class MpaTableResult:
    def __init__(
        self,
        table: str,
        sample_names: list[str],
        taxon_count: int,
        level: str,
    ) -> None:
        self.table = table
        self.sample_names = sample_names
        self.taxon_count = taxon_count
        self.level = level


def make_mpa_table(
    input_files: list[Path],
    level: str,
    keep_taxid: bool = False,
    short_names: bool = False,
    drop_unclassified: bool = False,
    normalise: bool = False,
) -> MpaTableResult:
    resolved_level = resolve_level(level)
    if not input_files:
        raise ValueError("provide at least one MetaPhlAn profile file")

    sample_names = [_sample_name(path) for path in input_files]
    _check_duplicate_sample_names(sample_names, input_files)

    table: OrderedDict[str, dict] = OrderedDict()
    any_level_rows = False

    for sample, path in zip(sample_names, input_files):
        for clade, tax_id, abundance, rank in _parse_profile(path, resolved_level):
            if rank is None and drop_unclassified:
                continue
            if rank is not None:
                any_level_rows = True

            entry = table.setdefault(clade, {"tax_id": tax_id, "abund": {}})
            if not entry["tax_id"]:
                entry["tax_id"] = tax_id
            entry["abund"][sample] = abundance

    if not any_level_rows:
        raise ValueError(
            f"no taxa found at level {level!r} in any input file. "
            "Check that profiles were generated deep enough for this rank."
        )

    if normalise:
        _normalise_table(table, sample_names)

    ordered_clades = _sort_by_mean_abundance(table, sample_names)
    return MpaTableResult(
        table=_format_table(
            table=table,
            ordered_clades=ordered_clades,
            sample_names=sample_names,
            level=resolved_level,
            keep_taxid=keep_taxid,
            short_names=short_names,
        ),
        sample_names=sample_names,
        taxon_count=len(ordered_clades),
        level=resolved_level,
    )


def resolve_level(value: str) -> str:
    resolved = value.strip().lower()
    if resolved in RANK_LETTERS:
        return resolved
    if resolved in RANK_WORDS:
        return RANK_WORDS[resolved]

    valid_letters = ", ".join(RANK_LETTERS)
    valid_words = ", ".join(RANK_WORDS)
    raise ValueError(
        f"unrecognised taxonomic level {value!r}. "
        f"Valid levels: {valid_letters} or full names: {valid_words}"
    )


def _sample_name(path: Path) -> str:
    return path.stem


def _check_duplicate_sample_names(
    sample_names: list[str],
    input_files: list[Path],
) -> None:
    seen: dict[str, Path] = {}
    for sample_name, path in zip(sample_names, input_files):
        if sample_name in seen:
            raise ValueError(
                f"duplicate sample name {sample_name!r} from {seen[sample_name]!s} "
                f"and {path!s}. Rename one of the inputs."
            )
        seen[sample_name] = path


def _parse_profile(path: Path, level: str) -> list[tuple[str, str, float, Optional[str]]]:
    rows: list[tuple[str, str, float, Optional[str]]] = []
    header: Optional[dict[str, int]] = None

    with path.open() as handle:
        for line in handle:
            if line.startswith("#clade_name"):
                header = _parse_header(line)
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                raise ValueError(
                    f"{path}: no '#clade_name' header line found before data"
                )
            if "relative_abundance" not in header:
                raise ValueError(f"{path}: no 'relative_abundance' column in header")

            fields = line.rstrip("\n").split("\t")
            clade = fields[0]
            abundance_index = header["relative_abundance"]
            try:
                abundance = float(fields[abundance_index])
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"{path}: cannot read relative_abundance for {clade!r}"
                ) from exc

            tax_id = ""
            tax_id_index = header.get("NCBI_tax_id")
            if tax_id_index is not None and tax_id_index < len(fields):
                tax_id = fields[tax_id_index]

            if clade == UNCLASSIFIED:
                rows.append((UNCLASSIFIED, tax_id, abundance, None))
            elif _clade_rank(clade) == level:
                rows.append((clade, tax_id, abundance, level))

    return rows


def _parse_header(line: str) -> dict[str, int]:
    columns = line.lstrip("#").rstrip("\n").split("\t")
    return {name: index for index, name in enumerate(columns)}


def _clade_rank(clade: str) -> Optional[str]:
    last_token = clade.split("|")[-1]
    if "__" not in last_token:
        return None
    return last_token.split("__", 1)[0]


def _sort_by_mean_abundance(
    table: OrderedDict[str, dict],
    sample_names: list[str],
) -> list[str]:
    def mean_abundance(clade: str) -> float:
        abundances = table[clade]["abund"]
        return sum(abundances.get(sample, 0.0) for sample in sample_names) / len(
            sample_names
        )

    return sorted(table, key=mean_abundance, reverse=True)


def _normalise_table(
    table: OrderedDict[str, dict],
    sample_names: list[str],
) -> None:
    for sample in sample_names:
        total = sum(entry["abund"].get(sample, 0.0) for entry in table.values())
        if total == 0:
            continue

        for entry in table.values():
            entry["abund"][sample] = entry["abund"].get(sample, 0.0) / total * 100


def _format_table(
    table: OrderedDict[str, dict],
    ordered_clades: list[str],
    sample_names: list[str],
    level: str,
    keep_taxid: bool,
    short_names: bool,
) -> str:
    lines: list[str] = []
    header = ["clade_name"]
    if keep_taxid:
        header.append("NCBI_tax_id")
    header.extend(sample_names)
    lines.append("\t".join(header))

    for clade in ordered_clades:
        entry = table[clade]
        cells = [_row_label(clade, level, short_names)]
        if keep_taxid:
            cells.append(entry["tax_id"])
        cells.extend(f"{entry['abund'].get(sample, 0.0):g}" for sample in sample_names)
        lines.append("\t".join(cells))

    return "\n".join(lines) + "\n"


def _row_label(clade: str, level: str, short_names: bool) -> str:
    if clade == UNCLASSIFIED or not short_names:
        return clade

    tokens = clade.split("|")
    if level == "s":
        genus = next((token for token in tokens if token.startswith("g__")), None)
        species = tokens[-1]
        return f"{genus}|{species}" if genus else species

    return tokens[-1]
