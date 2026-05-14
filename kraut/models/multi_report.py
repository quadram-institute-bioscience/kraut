from typing import Dict, List

import pandas as pd

from .kraken_data import KrakenNode, KrakenReport


class MultiKrakenReport:
    def __init__(self, samples: List[str] = None):
        self.samples = samples if samples else []
        self._validate_unique_samples(self.samples)
        self.data: Dict[int, Dict] = {}
        self.tax_id_to_rank: Dict[int, str] = {}
        self.tax_id_to_name: Dict[int, str] = {}

    def add_report(self, report: KrakenReport, sample_name: str):
        if sample_name in self.samples:
            raise ValueError(f"Duplicate sample name: {sample_name}")

        nodes_to_visit = self._report_nodes(report)
        self._validate_report_taxa(nodes_to_visit)

        self.samples.append(sample_name)
        sample_idx = len(self.samples) - 1

        for node in nodes_to_visit:
            parent_tax_id = self._parent_tax_id(node)
            if node.tax_id not in self.data:
                self.data[node.tax_id] = {
                    "name": node.name,
                    "rank": node.rank_code,
                    "parent_tax_id": parent_tax_id,
                    "clade_counts": {},
                    "taxon_counts": {},
                }
                self.tax_id_to_rank[node.tax_id] = node.rank_code
                self.tax_id_to_name[node.tax_id] = node.name

            self.data[node.tax_id]["clade_counts"][sample_idx] = node.clade_counts
            self.data[node.tax_id]["taxon_counts"][sample_idx] = node.taxon_counts

    def _report_nodes(self, report: KrakenReport) -> List[KrakenNode]:
        nodes_to_visit = []
        if report.unclassified:
            nodes_to_visit.append(report.unclassified)
        if report.root:
            nodes_to_visit.append(report.root)
            self._collect_nodes(report.root, nodes_to_visit)
        return nodes_to_visit

    def _validate_report_taxa(self, nodes: List[KrakenNode]) -> None:
        seen_in_report: Dict[int, tuple[str, str, int | None]] = {}

        for node in nodes:
            signature = (node.name, node.rank_code, self._parent_tax_id(node))
            existing_signature = seen_in_report.get(node.tax_id)
            if existing_signature is not None:
                if existing_signature != signature:
                    raise ValueError(
                        f"TaxID {node.tax_id} appears more than once with conflicting taxonomy"
                    )
                raise ValueError(
                    f"TaxID {node.tax_id} appears more than once in a report"
                )
            seen_in_report[node.tax_id] = signature

            existing = self.data.get(node.tax_id)
            if existing is not None:
                self._validate_existing_taxon(node, existing)

    def _validate_existing_taxon(self, node: KrakenNode, existing: Dict) -> None:
        conflicts = []
        if existing["name"] != node.name:
            conflicts.append(f"name '{existing['name']}' vs '{node.name}'")
        if existing["rank"] != node.rank_code:
            conflicts.append(f"rank '{existing['rank']}' vs '{node.rank_code}'")

        existing_parent = existing.get("parent_tax_id")
        parent_tax_id = self._parent_tax_id(node)
        if existing_parent != parent_tax_id:
            conflicts.append(f"parent TaxID {existing_parent} vs {parent_tax_id}")

        if conflicts:
            details = ", ".join(conflicts)
            raise ValueError(f"TaxID {node.tax_id} has conflicting taxonomy: {details}")

    @staticmethod
    def _validate_unique_samples(sample_names: List[str]) -> None:
        seen = set()
        duplicates = []
        for sample_name in sample_names:
            if sample_name in seen and sample_name not in duplicates:
                duplicates.append(sample_name)
            seen.add(sample_name)

        if duplicates:
            names = ", ".join(duplicates)
            raise ValueError(f"Duplicate sample name(s): {names}")

    @staticmethod
    def _parent_tax_id(node: KrakenNode) -> int | None:
        if node.parent is None:
            return None
        return node.parent.tax_id

    def _collect_nodes(self, node: KrakenNode, list_acc: List[KrakenNode]):
        for child in node.children:
            list_acc.append(child)
            self._collect_nodes(child, list_acc)

    def to_dataframe(
        self,
        metric: str = "TOT",
        level: str = "S",
        use_taxid: bool = False,
        rank_prefix: bool = False,
        add_lineage: bool = False,
    ) -> pd.DataFrame:
        """
        metric: 'TOT' (Clade counts), 'LVL' (Taxon counts), 'PERCENTAGE' (Clade %)
        level: Rank code or 'ALL'
        """
        if use_taxid and add_lineage:
            raise ValueError("--taxid cannot be used with --add-lineage")

        rows = []
        label_tax_ids: Dict[str, List[int]] = {}

        for tax_id, info in self.data.items():
            rank = info["rank"]
            if level != "ALL" and rank != level:
                continue

            key = self._taxon_key(tax_id, use_taxid, rank_prefix, add_lineage)
            label_tax_ids.setdefault(key, []).append(tax_id)

            row = {"#Taxon": key}
            source_dict = (
                info["taxon_counts"] if metric == "LVL" else info["clade_counts"]
            )

            for i, sample in enumerate(self.samples):
                row[sample] = source_dict.get(i, 0)

            rows.append(row)

        self._validate_unique_taxon_labels(label_tax_ids)
        return pd.DataFrame(rows)

    def _taxon_key(
        self,
        tax_id: int,
        use_taxid: bool,
        rank_prefix: bool,
        add_lineage: bool,
    ) -> str:
        if use_taxid:
            return str(tax_id)

        info = self.data[tax_id]
        if add_lineage:
            return self._lineage_key(tax_id)

        name = info["name"].strip()
        if rank_prefix:
            return f"{self._rank_prefix(info['rank'])}{name}"

        return name

    def _lineage_key(self, tax_id: int) -> str:
        info = self.data[tax_id]
        if tax_id == 0 or info["rank"] == "U":
            return f"u__{info['name'].strip()}"

        lineage = []
        seen = set()
        current_tax_id: int | None = tax_id

        while current_tax_id is not None:
            if current_tax_id in seen:
                raise ValueError(f"Taxonomy cycle detected at TaxID {current_tax_id}")
            seen.add(current_tax_id)

            current = self.data.get(current_tax_id)
            if current is None:
                break

            rank = current["rank"]
            if rank not in {"R", "U"}:
                lineage.append(f"{self._rank_prefix(rank)}{current['name'].strip()}")
            current_tax_id = current.get("parent_tax_id")

        lineage.reverse()
        if not lineage:
            return info["name"].strip()
        return ",".join(lineage)

    @staticmethod
    def _rank_prefix(rank: str) -> str:
        prefix_rank = "k" if rank in {"D", "K"} else rank.lower()
        return f"{prefix_rank}__"

    @staticmethod
    def _validate_unique_taxon_labels(label_tax_ids: Dict[str, List[int]]) -> None:
        duplicates = {
            label: tax_ids
            for label, tax_ids in label_tax_ids.items()
            if len(set(tax_ids)) > 1
        }
        if not duplicates:
            return

        details = []
        for label, tax_ids in duplicates.items():
            tax_id_list = ", ".join(str(tax_id) for tax_id in sorted(set(tax_ids)))
            details.append(f"{label} (TaxIDs: {tax_id_list})")
        duplicate_list = "; ".join(details)
        raise ValueError(
            "Duplicate taxon label(s) would be written to #Taxon: "
            f"{duplicate_list}. Use --taxid or --add-lineage to disambiguate."
        )

    def _unclassified_row(
        self,
        metric: str,
        use_taxid: bool,
        rank_prefix: bool,
        add_lineage: bool,
    ) -> dict | None:
        """Build a single unclassified row dict, or None if not present."""
        if 0 not in self.data:
            return None

        source_dict = (
            self.data[0]["taxon_counts"]
            if metric == "LVL"
            else self.data[0]["clade_counts"]
        )
        row = {"#Taxon": self._taxon_key(0, use_taxid, rank_prefix, add_lineage)}
        for i, sample in enumerate(self.samples):
            row[sample] = source_dict.get(i, 0)
        return row

    def to_tsv(
        self,
        metric: str = "TOT",
        level: str = "S",
        no_unclassified: bool = False,
        use_taxid: bool = False,
        rank_prefix: bool = False,
        min_perc: float = 0.0,
        add_lineage: bool = False,
    ) -> str:
        if metric == "COUNTS":
            metric = "TOT"

        df = self.to_dataframe(metric, level, use_taxid, rank_prefix, add_lineage)

        if df.empty:
            return ""

        if no_unclassified and level == "ALL" and 0 in self.data:
            unclassified_key = self._taxon_key(0, use_taxid, rank_prefix, add_lineage)
            df = df[df["#Taxon"] != unclassified_key]
        elif not no_unclassified and level != "ALL":
            unclassified_row = self._unclassified_row(
                metric,
                use_taxid,
                rank_prefix,
                add_lineage,
            )
            if unclassified_row:
                df = pd.concat(
                    [pd.DataFrame([unclassified_row]), df], ignore_index=True
                )

        if metric == "PERCENTAGE":
            sample_cols = self.samples
            for col in sample_cols:
                total = df[col].sum()
                if total > 0:
                    df[col] = (df[col] / total) * 100

        if min_perc > 0.0:
            sample_cols = self.samples
            col_totals = {col: df[col].sum() for col in sample_cols}

            def row_max_perc(row):
                return max(
                    (row[col] / col_totals[col] * 100 if col_totals[col] > 0 else 0.0)
                    for col in sample_cols
                )

            mask = df.apply(row_max_perc, axis=1) >= min_perc
            df = df[mask]

        return df.to_csv(sep="\t", index=False)
