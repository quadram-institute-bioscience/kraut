from typing import List, Dict
import pandas as pd
from .kraken_data import KrakenReport, KrakenNode

class MultiKrakenReport:
    def __init__(self, samples: List[str] = None):
        self.samples = samples if samples else []
        # Structure: {tax_id: {'name': name, 'rank': rank, 'counts': [count_sample_1, count_sample_2, ...]}}
        self.data: Dict[int, Dict] = {}
        self.tax_id_to_rank: Dict[int, str] = {}
        self.tax_id_to_name: Dict[int, str] = {}
        # We need to maintain tree relationships for aggregation if needed, 
        # but for simple table output, flat structure with TaxID key is enough 
        # provided we have name and rank.
        # Add taxon_counts support
    
    def add_report(self, report: KrakenReport, sample_name: str):
        self.samples.append(sample_name)
        sample_idx = len(self.samples) - 1
        
        nodes_to_visit = []
        if report.unclassified:
            nodes_to_visit.append(report.unclassified)
        if report.root:
            nodes_to_visit.append(report.root)
            self._collect_nodes(report.root, nodes_to_visit)
            
        for node in nodes_to_visit:
            if node.tax_id not in self.data:
                self.data[node.tax_id] = {
                    'name': node.name, 
                    'rank': node.rank_code,
                    'clade_counts': {}, # keys: sample_idx
                    'taxon_counts': {}  # keys: sample_idx
                }
                self.tax_id_to_rank[node.tax_id] = node.rank_code
                self.tax_id_to_name[node.tax_id] = node.name
            
            # Store both counts
            self.data[node.tax_id]['clade_counts'][sample_idx] = node.clade_counts 
            self.data[node.tax_id]['taxon_counts'][sample_idx] = node.taxon_counts
            
    def _collect_nodes(self, node: KrakenNode, list_acc: List[KrakenNode]):
        for child in node.children:
            list_acc.append(child)
            self._collect_nodes(child, list_acc)

    def to_dataframe(self, metric: str = 'TOT', level: str = 'S', 
                     use_taxid: bool = False, rank_prefix: bool = False) -> pd.DataFrame:
        """
        metric: 'TOT' (Clade counts), 'LVL' (Taxon counts), 'PERCENTAGE' (Clade %)
        level: Rank code or 'ALL'
        """
        rows = []
        
        for tax_id, info in self.data.items():
            rank = info['rank']
            if level != 'ALL' and rank != level:
                continue
            
            # Key Name
            if use_taxid:
                key = str(tax_id)
            else:
                name = info['name'].strip()
                if rank_prefix:
                    # Map rank code to prefix
                    # Simple mapping: S->s__, G->g__ etc.
                    prefix = f"{rank.lower()}__"
                    # Determine display name
                    # If species and not taxid: "Genus species" usually already in name for Kraken?
                    # Yes, Kraken name is usually full scientific name. 
                    # If rank is species, print "Genus species"? It IS "Genus species". 
                    # Prompt says: "If the rank is "S" please print "Genus species", else only print the rank (eg. Prevotella)"
                    # Kraken data already provides full scientific name in `name`.
                    # So we just prepend prefix.
                    key = f"{prefix}{name}"
                else:
                    # Default: Taxon name.
                    # Prompt: "If the rank is "S" please print "Genus species", else only print the rank (eg. Prevotella)"
                    # Kraken names ARE "Genus species" (e.g. "Segatella copri").
                    # But for higher ranks? "Prevotella".
                    # So default behavior matches Kraken names.
                    key = name

            row = {'#Taxon': key} # Prompt wants #Taxon as header
            
            # Retrieve counts
            # Metric logic
            # TOT = Clade counts
            # LVL = Taxon counts
            # PERCENTAGE = Clade counts / Total (handled after)
            
            source_dict = info['taxon_counts'] if metric == 'LVL' else info['clade_counts']
            
            current_counts = []
            for i, sample in enumerate(self.samples):
                count = source_dict.get(i, 0)
                row[sample] = count
                current_counts.append(count)
                
            rows.append(row)
            
        df = pd.DataFrame(rows)
        
        return df

    def _unclassified_row(self, metric: str, use_taxid: bool, rank_prefix: bool) -> dict | None:
        """Build a single unclassified row dict, or None if not present."""
        unclassified_entry = self.data.get(0)
        if unclassified_entry is None:
            return None

        if use_taxid:
            key = "0"
        else:
            name = unclassified_entry['name'].strip()
            key = f"u__{name}" if rank_prefix else name

        source_dict = unclassified_entry['taxon_counts'] if metric == 'LVL' else unclassified_entry['clade_counts']
        row = {'#Taxon': key}
        for i, sample in enumerate(self.samples):
            row[sample] = source_dict.get(i, 0)
        return row

    def to_tsv(self, metric: str = 'TOT', level: str = 'S',
               no_unclassified: bool = False,
               use_taxid: bool = False, rank_prefix: bool = False,
               min_perc: float = 0.0) -> str:

        # Map old 'COUNTS' to 'TOT' for backward compatibility defaults if needed,
        # or just assume caller uses new constants.
        if metric == 'COUNTS':
            metric = 'TOT'

        df = self.to_dataframe(metric, level, use_taxid, rank_prefix)

        if df.empty:
            return ""

        if not no_unclassified:
            unclassified_row = self._unclassified_row(metric, use_taxid, rank_prefix)
            if unclassified_row:
                df = pd.concat([pd.DataFrame([unclassified_row]), df], ignore_index=True)

        if metric == 'PERCENTAGE':
            # Normalize sample columns
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

        return df.to_csv(sep='\t', index=False)
