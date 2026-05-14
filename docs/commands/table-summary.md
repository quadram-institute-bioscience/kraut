# table-summary

Summarizes a wide abundance table produced by `kraut make-table`.

## Syntax

```bash
kraut table-summary [OPTIONS] INPUT_FILE
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILE` | PATH | **Required**. A table produced by `kraut make-table`. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--top-taxa` | | INTEGER | Number of top taxa to report (default: `10`). |
| `--output` | `-o` | PATH | Save the summary report to a text file (default: stdout). |

## Report Contents

The report includes:

- Table size as samples x taxa.
- Total counts per sample when the table contains counts.
- Top taxa across the table, calculated after normalising each sample so larger samples do not dominate the ranking.

Relative abundance tables skip the total counts section.

## Examples

```bash
kraut make-table reports/*.tsv -o species_table.tsv
kraut table-summary species_table.tsv --top-taxa 20 -o species_summary.txt
```

---
[← Back to Commands](../commands.md)
