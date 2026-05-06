# dendrogram

Plots a hierarchical clustering dendrogram from a wide abundance table or from multiple Kraken/Bracken report files. Distances use the same methods as `kraut beta`.

## Syntax

```bash
kraut dendrogram [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One wide abundance table, or two or more Kraken/Bracken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output` | `-o` | PATH | Output dendrogram plot (`.html`, `.png`, `.svg`, `.pdf`). |
| `--distance` | `-d` | TEXT | Distance method: `braycurtis`, `aitchison`, or `jaccard` (default: `braycurtis`). |
| `--clustering` | `-c` | TEXT | Linkage method: `ward`, `average`, `single`, or `complete` (default: `ward`). |
| `--rank` | `-r` | TEXT | Taxonomic rank to use when reading report files (default: `S`). |
| `--abundance-metric` | | TEXT | Report abundance metric: `TOT` cumulative counts or `LVL` taxon-specific counts (default: `TOT`). |
| `--include-unclassified` | | | Include unclassified reads as a feature. |
| `--pseudocount` | | FLOAT | Pseudocount added before CLR for Aitchison distance (default: 1). |
| `--min-feature-count` | | FLOAT | Keep taxa with total abundance at least this value. |
| `--min-samples` | | INTEGER | Keep taxa detected in at least this many samples (default: 1). |
| `--presence-threshold` | | FLOAT | Presence cutoff for Jaccard; a feature is present when value is greater than this (default: 0). |
| `--metadata` | | PATH | Optional sample metadata CSV or TSV file; the first column must match sample names. |
| `--color-by` | | TEXT | Metadata column used to color sample labels. |

## Clustering Methods

- `ward`: treats clustering as a variance minimization problem.
- `average`: distance between clusters is the average distance between all points in the two clusters.
- `single`: distance between clusters is the nearest pair of points.
- `complete`: distance between clusters is the furthest pair of points.

## Examples

```bash
kraut dendrogram reports/*.tsv \
  --distance braycurtis \
  --clustering ward \
  --output dendrogram.html
```

```bash
kraut dendrogram reports/*.tsv \
  --metadata metadata.tsv \
  --color-by Treatment \
  --output dendrogram.svg
```

---
[← Back to Commands](../commands.md)
