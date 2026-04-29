# kraut

![GitHub Release](https://img.shields.io/github/v/release/quadram-institute-bioscience/kraut)
![PyPI - Status](https://img.shields.io/pypi/status/krautils)
![PyPI - Downloads](https://img.shields.io/pypi/dm/krautils)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/quadram-institute-bioscience/kraut)
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
- `merge-reports`: merge reports into a simple comparison table.
- `make-table`: build a configurable multi-sample abundance table.
- `alpha`: calculate alpha diversity from Kraken or Bracken reports.
- `beta`: calculate beta diversity distance matrices and heatmap/PCA plots.
- `split-combine-table`: split a KrakenTools combined table into ALL/LVL tables.
- `plot-single`: plot one sample as an HTML or static composition chart.
- `plot-multi`: plot multiple samples as an HTML or static stacked/bubble chart.

```bash
kraut alpha reports/*.tsv -o alpha.tsv -p alpha.html --metrics core --add-metrics chao1,ace
kraut beta reports/*.tsv -o beta.tsv --plot beta.html --pca beta_pca.html
```


Example plots:


![Kraut Multi](https://raw.githubusercontent.com/quadram-institute-bioscience/kraut/main/docs/multi.png)


![Kraut Single](https://raw.githubusercontent.com/quadram-institute-bioscience/kraut/main/docs/single.png)

### License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
