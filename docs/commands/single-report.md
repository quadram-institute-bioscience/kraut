# single-report

Parses a Kraken report and prints it in text format, optionally filtering by abundance or taxonomic level.

## Syntax

```bash
kraut single-report [OPTIONS]
```

### Options

| Option | Short | Type | Description |
| :--- | :--- | :--- | :--- |
| `--input` | `-i` | PATH | **Required**. Input Kraken report file (KREP). |
| `--output` | `-o` | PATH | Output file (default: stdout). |
| `--min-fract` | `-m` | FLOAT | Minimum fraction of reads to keep a taxon (default: 0.0). |
| `--min-count` | `-c` | INTEGER | Minimum count of reads to keep a taxon (default: 0). |
| `--min-level` | `-l` | TEXT | Minimum level to keep (K, P, C, O, F, G, S). |
| `--max-level` | `-L` | TEXT | Maximum level to keep (e.g., do not print below Species). |

## Examples

### Basic conversion
```bash
kraut single-report -i sample.krep -o output.txt
```

### Filtering by abundance
Keep only taxa with at least 1% abundance and at least 100 reads:
```bash
kraut single-report -i sample.krep -m 0.01 -c 100
```

### Filtering by rank
Only show results from Phylum to Genus:
```bash
kraut single-report -i sample.krep -l P -L G
```

---
[← Back to Commands](../commands.md)
