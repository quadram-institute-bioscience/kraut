# merge-reports

Merges multiple Kraken reports into a single table at a specific taxonomic level.

## Syntax

```bash
kraut merge-reports [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One or more Kraken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output` | `-o` | PATH | Output file (default: stdout). |
| `--metric` | `-m` | TEXT | Metric to report: `COUNTS` or `PERCENTAGES` (default: `COUNTS`). |
| `--level` | `-l` | TEXT | Taxonomic level to report (default: `S`). |
| `--no-unclassified` | | | Exclude unclassified reads from the output. |

## Examples

### Merging multiple reports into a counts table
```bash
kraut merge-reports sample1.krep sample2.krep sample3.krep -o merged_counts.tsv
```

### Merging at Genus level with percentages
```bash
kraut merge-reports reports/*.krep -l G -m PERCENTAGES -o genus_percentages.tsv
```

---
[← Back to Commands](../commands.md)
