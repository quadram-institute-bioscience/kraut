# split-combine-table

Splits a KrakenTools combined table into two separate tables: one with all taxonomic information (`ALL`) and one filtered by a specific rank (`LVL`), while cleaning up headers for better compatibility with downstream tools.

## Syntax

```bash
kraut split-combine-table [OPTIONS]
```

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--input` | `-i` | PATH | **Required**. Input table in KrakenTools combined format. |
| `--output` | `-o` | TEXT | **Required**. Output basename (will create `_all.tsv` and `_lvl.tsv`). |
| `--rank` | `-r` | TEXT | Filter by rank code (e.g., `S`). |
| `--taxid` | | | Use TaxID as the first column instead of the taxon name. |

## Examples

### Splitting a combined table at Species level
```bash
kraut split-combine-table -i combined_table.txt -o split_table -r S
```
This will generate `split_table_all.tsv` and `split_table_lvl.tsv`.

### Using TaxIDs for compatibility
```bash
kraut split-combine-table -i combined_table.txt -o taxid_output --taxid
```

---
[← Back to Commands](../commands.md)
