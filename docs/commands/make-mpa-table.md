# make-mpa-table

Merges multiple MetaPhlAn profile files into one wide abundance table. Rows are taxa, columns are samples, and missing taxa are filled with `0`.

## Syntax

```bash
kraut make-mpa-table [OPTIONS] INPUT_FILES...
```

### Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `INPUT_FILES...` | PATH | **Required**. One or more MetaPhlAn profile files. |

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--output` | `-o` | PATH | Output TSV file (default: stdout). |
| `--level` | `-l` | TEXT | Taxonomic level: `k`, `p`, `c`, `o`, `f`, `g`, `s`, `t`, or a full rank name (default: `S`). |
| `--keep-taxid` | | | Add an `NCBI_tax_id` column. |
| `--short-names` | | | Use terminal clade names only. At species level, include genus plus species. |
| `--drop-unclassified` | | | Drop the `UNCLASSIFIED` row. |
| `--normalise` | | | Scale each sample column so retained rows sum to 100. |

Rows are sorted by mean abundance across samples in descending order. If `--normalise` is used, sorting uses the normalised abundances. Sample names are inferred from the input filename with the final extension removed.

## Examples

### Creating a species-level MetaPhlAn table

```bash
kraut make-mpa-table profiles/*_profile.tsv -o metaphlan_species.tsv
```

### Adding taxids and shorter labels

```bash
kraut make-mpa-table profiles/*.tsv -l species --keep-taxid --short-names -o metaphlan_species.tsv
```

### Creating a genus table without unclassified abundance

```bash
kraut make-mpa-table profiles/*.tsv -l genus --drop-unclassified -o metaphlan_genus.tsv
```

### Dropping unclassified and recalibrating to 100

```bash
kraut make-mpa-table profiles/*.tsv --drop-unclassified --normalise -o metaphlan_species.tsv
```

---
[Back to Commands](../commands.md)
