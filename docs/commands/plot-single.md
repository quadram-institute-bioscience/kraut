# plot-single

Plots the taxonomic composition for a single Kraken report. Supports both static (PNG, SVG, PDF) and interactive (HTML) outputs.

## Syntax

```bash
kraut plot-single [OPTIONS] INPUT_FILE
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILE` | PATH | **Required**. Input Kraken report file. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output` | `-o` | PATH | **Required**. Output plot file. Extension determines format. |
| `--rank` | `-r` | TEXT | Taxonomic rank to plot (default: `S`). |
| `--metric` | `-m` | TEXT | Metric: `TOT` or `LVL` (default: `TOT`). |
| `--min-perc` | | FLOAT | Fold taxa below this % abundance into "Others" (default: 1.0). |
| `--top-taxa` | | INTEGER | Keep top N taxa and fold others into "Others" (default: 0). |
| `--no-unclassified` | | | Exclude unclassified reads from the plot. |
| `--title` | | TEXT | Plot title. |

## Examples

### Generating a static composition plot
```bash
kraut plot-single sample.krep -o composition.png
```

### Generating an interactive HTML plot at Phylum level
```bash
kraut plot-single sample.krep -r P -o interactive.html
```

### Customizing appearance
```bash
kraut plot-single sample.krep --min-perc 0.5 --title "Sample A1 Composition" -o plot.svg
```

---
[← Back to Commands](../commands.md)
