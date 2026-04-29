# alpha

Calculates alpha diversity metrics for Kraken or Bracken report files and optionally generates a visualization.

## Syntax

```bash
kraut alpha [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One or more Kraken or Bracken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output-table` | `-o` | PATH | Output alpha diversity table (default: `alpha.tsv`). |
| `--plot` | `-p` | PATH | Optional plot file (.html, .png, .svg, .pdf). |
| `--rank` | `-r` | TEXT | Taxonomic rank to use (default: `S`). |
| `--metric` | `-m` | TEXT | Abundance metric: `TOT` (cumulative) or `LVL` (taxon-specific) (default: `TOT`). |
| `--metrics` | | TEXT | Alpha metric preset or comma-separated list (default: `core`). |
| `--add-metrics` | | TEXT | Comma-separated additional alpha metrics to append. |
| `--include-unclassified` | | | Include unclassified reads as a feature in calculations. |
| `--min-perc` | | FLOAT | Remove taxa where no sample reaches this % abundance (default: 0.0). |

## Examples

### Basic alpha diversity calculation
```bash
kraut alpha reports/*.krep -o diversity_stats.tsv
```

### Calculating specific metrics and plotting
Calculate Shannon and Simpson indices and save a PNG plot:
```bash
kraut alpha reports/*.krep --metrics "shannon,simpson" -p alpha_plot.png
```

### Using Bracken reports at Genus level
```bash
kraut alpha bracken_outputs/*.brep -r G -o genus_alpha.tsv
```
