# ma-export

Exports Kraken reports as MicrobiomeAnalyst-compatible files: `counts.csv`, `taxonomy.csv`, `metadata.csv`, and `tree.nwk`.

## Syntax

```bash
kraut ma-export [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One or more Kraken report files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--outdir` | `-o` | PATH | **Required**. Output directory. |
| `--metadata` | | PATH | Optional metadata CSV or TSV file. |
| `--metadata-sample-col` | | TEXT | Column in the metadata matching input sample names. |
| `--pseudo-col` | | TEXT | Metadata column to create or overwrite with pseudo group labels (default: `Random_label`). |
| `--metric` | `-m` | TEXT | Count metric: `LVL` (taxon counts) or `TOT` (clade counts) (default: `LVL`). |

## Example

```bash
kraut ma-export reports/*.krep -o ma_output
```

---
[← Back to Commands](../commands.md)