# kraut

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
- `split-combine-table`: split a KrakenTools combined table into ALL/LVL tables.
- `plot-single`: plot one sample as an HTML or static composition chart.
- `plot-multi`: plot multiple samples as an HTML or static stacked/bubble chart.

```bash
kraut alpha reports/*.tsv -o alpha.tsv -p alpha.html --metrics core --add-metrics chao1,ace
```


Example plots:


![Kraut Multi](https://raw.githubusercontent.com/quadram-institute-bioscience/kraut/main/docs/multi.png)


![Kraut Single](https://raw.githubusercontent.com/quadram-institute-bioscience/kraut/main/docs/single.png)
