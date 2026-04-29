# plot-multi

Plots taxonomic composition across multiple Kraken reports. Supports stacked bar charts and bubble charts.

## Syntax

```bash
kraut plot-multi [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One or more Kraken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output` | `-o` | PATH | **Required**. Output plot file. |
| `--rank` | `-r` | TEXT | Taxonomic rank to plot (default: `S`). |
| `--metric` | `-m` | TEXT | Metric: `TOT` or `LVL` (default: `TOT`). |
| `--kind` | | TEXT | Plot kind: `stacked` or `bubble` (default: `stacked`). |
| `--min-perc` | | FLOAT | Fold taxa below this % abundance into "Others" (default: 1.0). |
| `--top-taxa` | | INTEGER | Keep top N taxa and fold others into "Others" (default: 0). |
| `--no-unclassified` | | | Exclude unclassified reads. |
| `--title` | | TEXT | Plot title. |

## Examples

### Stacked bar chart of multiple samples
```bash
kraut plot-multi reports/*.krep -o stacked_bars.png
```

### Interactive bubble chart
```bash
kraut plot-multi reports/*.krep --kind bubble -o bubbles.html
```

### Focusing on top 10 Genera
```bash
kraut plot-multi reports/*.krep -r G --top-taxa 10 -o top_genera.png
```

---
[← Back to Commands](../commands.md)
