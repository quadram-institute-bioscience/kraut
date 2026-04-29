from pathlib import Path
from typing import List, Optional

import typer

from kraut.models.kraken_data import KrakenReport
from kraut.models.ranks import CANONICAL_RANKS, rank_read_counts


RANK_ROWS = ("U",) + CANONICAL_RANKS


def _direct_default(value, default):
    if isinstance(value, typer.models.OptionInfo):
        return default
    return value


def run(
    input_files: List[Path] = typer.Argument(
        ...,
        help="Input Kraken or Bracken report files",
    ),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output TSV file (default: stdout)",
    ),
    counts: bool = typer.Option(
        False,
        "--counts",
        help="Print raw read counts instead of percentages",
    ),
):
    """Profile reads assigned at each canonical taxonomic rank."""
    output_file = _direct_default(output_file, None)
    counts = _direct_default(counts, False)

    missing = [input_file for input_file in input_files if not input_file.exists()]
    if missing:
        typer.echo(f"Error: Input file {missing[0]} does not exist.", err=True)
        raise typer.Exit(code=1)

    sample_names = [_sample_name_from_path(input_file) for input_file in input_files]
    sample_counts = []
    for input_file in input_files:
        report = KrakenReport.from_file(str(input_file))
        sample_counts.append(
            rank_read_counts(
                report,
                mode="exact",
                include_unclassified=True,
                include_root=True,
            )
        )

    table = _rank_table(sample_names, sample_counts, as_counts=counts)
    if output_file:
        output_file.write_text(table)
    else:
        typer.echo(table, nl=False)


def _rank_table(
    sample_names: list[str],
    sample_counts: list[dict[str, int]],
    as_counts: bool,
) -> str:
    rows = [["#Rank", *sample_names]]
    totals = [sum(counts.values()) for counts in sample_counts]

    for rank in RANK_ROWS:
        row = [rank]
        for counts, total in zip(sample_counts, totals):
            value = counts.get(rank, 0)
            if as_counts:
                row.append(str(value))
            elif total > 0:
                row.append(str(value / total * 100))
            else:
                row.append("0.0")
        rows.append(row)

    return "\n".join("\t".join(row) for row in rows) + "\n"


def _sample_name_from_path(path: Path) -> str:
    sample_name = path.stem
    if sample_name.endswith(".krep"):
        sample_name = Path(sample_name).stem
    return sample_name
