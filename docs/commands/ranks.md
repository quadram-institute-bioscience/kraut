# ranks

Profiles the reads assigned at each canonical taxonomic rank (`U`, `R`, `D`, `K`, `P`, `C`, `O`, `F`, `G`, `S`) for one or more Kraken or Bracken reports.

## Syntax

```bash
kraut ranks [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One or more Kraken or Bracken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output` | `-o` | PATH | Output TSV file (default: stdout). |
| `--counts` | | | Print raw read counts instead of percentages. |
| `--plot` | `-p` | PATH | Write a stacked rank-composition bar chart (`.html`, `.png`, `.pdf`, `.svg`). |

## Examples

### Rank composition table
```bash
kraut ranks reports/*.krep -o ranks.tsv
```

### Read counts and a plot
```bash
kraut ranks reports/*.krep --counts -p ranks_plot.png
```

---
[← Back to Commands](../commands.md)