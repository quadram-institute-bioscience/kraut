# kraut

[![Run Pytest](https://github.com/quadram-institute-bioscience/kraut/actions/workflows/pytest.yml/badge.svg)](https://github.com/quadram-institute-bioscience/kraut/actions/workflows/pytest.yml)
[![GitHub Release](https://img.shields.io/github/v/release/quadram-institute-bioscience/kraut)](https://github.com/quadram-institute-bioscience/kraut/releases)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/krautils)](https://pypi.org/project/krautils/)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/kraut)](https://bioconda.github.io/recipes/kraut/README.html)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/m/quadram-institute-bioscience/kraut)](https://pypi.org/project/krautils/)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)


![Kraut logo](https://raw.githubusercontent.com/quadram-institute-bioscience/kraut/main/docs/kraut.png)


A python package for parsing, merging, and analyzing Kraken2 output files.

### Installation

From PyPI:

```bash
pip install krautils
```

From a local checkout:

```bash
pip install .
```

For developer's dependencies:
```bash
pip install .[dev]
pytest
```

### Command

`kraut` parses, filters, merges, tabulates, splits, and plots Kraken2 reports.

- `single-report`: filter or reformat one Kraken report.
- `make-table`: build a configurable multi-sample abundance table.
- `make-mpa-table`: merge MetaPhlAn profiles into a multi-sample abundance table.
- `table-summary`: summarize a table produced by `make-table`.
- `ma-export`: export reports for MicrobiomeAnalyst as counts, taxonomy, metadata, and tree files.
- `alpha`: calculate alpha diversity from Kraken or Bracken reports.
- `beta`: calculate beta diversity distance matrices and heatmap/PCA plots.
- `dendrogram`: plot hierarchical clustering dendrograms from beta distances.
- `split-combine-table`: split a KrakenTools combined table into ALL/LVL tables.
- `plot-single`: plot one sample as an HTML or static composition chart.
- `plot-multi`: plot multiple samples as an HTML or static stacked/bubble chart.

```bash
kraut alpha reports/*.tsv -o alpha.tsv -p alpha.html --metrics core --add-metrics chao1,ace
kraut beta reports/*.tsv -o beta.tsv --plot beta.html --pca beta_pca.html
kraut dendrogram reports/*.tsv -o dendrogram.html --distance braycurtis --clustering ward
kraut ma-export reports/*.tsv -o microbiomeanalyst --metadata metadata.tsv --metadata-sample-col Sample
kraut make-mpa-table profiles/*_profile.tsv -o metaphlan_species.tsv --keep-taxid --drop-unclassified --normalise
```

`ma-export` writes `counts.csv`, `taxonomy.csv`, `metadata.csv`, and `tree.nwk`.
It uses taxon-specific counts by default (`--metric LVL`); use `--metric TOT`
for cumulative clade counts.


Example plots:


![Kraut Multi](https://raw.githubusercontent.com/quadram-institute-bioscience/kraut/main/docs/multi.png)


![Kraut Single](https://raw.githubusercontent.com/quadram-institute-bioscience/kraut/main/docs/single.png)

### License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
