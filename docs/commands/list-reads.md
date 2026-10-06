# list-reads

Lists read names from a Kraken raw output file that match a taxon selection.

## Syntax

```bash
kraut list-reads [OPTIONS] RAW_OUTPUT
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `RAW_OUTPUT` | PATH | **Required**. A Kraken raw output TSV. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--taxon` | `-t` | INTEGER | TaxID to match; repeatable. |
| `--children` | | | Include each selected taxon and its descendants (requires `--report`). |
| `--report` | `-r` | PATH | Kraken report TSV; required with `--children`. |
| `--unclassified` | | | Include unclassified reads. |
| `--invert` | `-v` | | Output reads that do not match the selection. |
| `--output` | `-o` | PATH | Output file (default: stdout). |

## Examples

### Reads assigned to one or more taxa
```bash
kraut list-reads sample.out --taxon 562 --taxon 208964
```

### Reads for a taxon and its descendants
```bash
kraut list-reads sample.out --taxon 2 --children --report sample.krep
```

---
[← Back to Commands](../commands.md)