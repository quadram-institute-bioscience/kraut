# make-table

Generates a highly customizable abundance table from multiple reports, with support for different metrics and formatting options.

## Syntax

```bash
kraut make-table [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One or more Kraken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output` | `-o` | PATH | Output file (default: stdout). |
| `--metric` | `-m` | TEXT | Metric: `LVL` (taxon-specific), `TOT` (cumulative), `PERCENTAGE` (default: `TOT`). |
| `--rank` | `-r` | TEXT | Taxonomic rank (K, P, C, O, F, G, S, ALL) (default: `S`). |
| `--rank-prefix` | `-p` | | Add rank prefix to names (e.g., `s__Species`). |
| `--taxid` | | | Use TaxID instead of taxon name as the identifier. |
| `--no-unclassified` | | | Exclude unclassified reads. |
| `--min-perc` | | FLOAT | Remove rows where no sample reaches this % abundance (default: 0.0). |

## Examples

### Creating a standard species table with cumulative counts
```bash
kraut make-table data/*.krep -o species_table.tsv
```

This is the standard command for merging multiple Kraken reports into one comparison table.

### Creating a genus table with rank prefixes and percentages
```bash
kraut make-table data/*.krep -r G -p -m PERCENTAGE -o genus_table.tsv
```

### Filtering low-abundance taxa
Keep only taxa that reach at least 0.5% abundance in at least one sample:
```bash
kraut make-table data/*.krep --min-perc 0.5 -o filtered_table.tsv
```

---
[← Back to Commands](../commands.md)
