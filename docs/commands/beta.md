# beta

Calculates beta diversity distance matrices from a wide abundance table or from multiple Kraken/Bracken report files. It also accepts combined Bracken output tables with `*_num` columns. Optional plots can be written as HTML or static image files.

## Syntax

```bash
kraut beta [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One wide abundance table, or two or more Kraken/Bracken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--metric`, `--method` | `-m` | TEXT | Beta diversity metric: `braycurtis`, `aitchison`, or `jaccard` (default: `braycurtis`). |
| `--output` | `-o` | PATH | Output distance matrix TSV file (default: stdout). |
| `--plot` | `-p` | PATH | Optional distance heatmap (`.html`, `.png`, `.svg`, `.pdf`). |
| `--pca` | | PATH | Optional 2D ordination plot (`.html`, `.png`, `.svg`, `.pdf`). Uses PCA for Aitchison/CLR data and PCoA for Bray-Curtis or Jaccard distances. |
| `--rank` | `-r` | TEXT | Taxonomic rank to use when reading report files (default: `S`). |
| `--abundance-metric` | | TEXT | Report abundance metric: `TOT` cumulative counts or `LVL` taxon-specific counts (default: `TOT`). |
| `--include-unclassified` | | | Include unclassified reads as a feature. |
| `--pseudocount` | | FLOAT | Pseudocount added before CLR for Aitchison distance (default: 1). |
| `--min-feature-count` | | FLOAT | Keep taxa with total abundance at least this value. |
| `--min-samples` | | INTEGER | Keep taxa detected in at least this many samples (default: 1). |
| `--presence-threshold` | | FLOAT | Presence cutoff for Jaccard; a feature is present when value is greater than this (default: 0). |

## Metrics

### Bray-Curtis

The default metric. Use this for abundance-weighted ecological differences between samples.

```bash
kraut beta table.tsv --metric braycurtis -o braycurtis.tsv
```

### Aitchison

Use this for compositional analysis of taxon-to-taxon relative proportions. Kraut adds a pseudocount before CLR transformation and warns when the table is very sparse, because Aitchison distances are sensitive to zero handling.

```bash
kraut beta table.tsv \
  --metric aitchison \
  --pseudocount 1 \
  --min-feature-count 10 \
  --min-samples 2 \
  -o aitchison.tsv
```

### Jaccard

Use this for presence/absence comparisons. The default presence rule is abundance greater than zero. A stricter threshold can reduce low-count Kraken2 noise.

```bash
kraut beta table.tsv \
  --metric jaccard \
  --presence-threshold 5 \
  -o jaccard.tsv
```

## Plots

`--plot` writes a heatmap of the distance matrix. `--pca` writes a 2D ordination plot.

```bash
kraut beta reports/*.tsv \
  --metric braycurtis \
  --output beta.tsv \
  --plot beta_heatmap.html \
  --pca beta_pca.html
```

---
[← Back to Commands](../commands.md)
