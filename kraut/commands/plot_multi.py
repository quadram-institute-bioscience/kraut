from pathlib import Path
from typing import List, Optional

import typer

from kraut.plotting import (
    build_multi_report,
    composition_dataframe,
    render_multi_composition,
)


def _direct_default(value, default):
    if isinstance(value, typer.models.OptionInfo):
        return default
    return value


def run(
    input_files: List[Path] = typer.Argument(..., help="Input Kraken report files"),
    output_file: Path = typer.Option(..., "--output", "-o", help="Output plot file"),
    rank: str = typer.Option("S", "--rank", "-r", help="Taxonomic rank to plot"),
    metric: str = typer.Option("TOT", "--metric", "-m", help="Metric: TOT or LVL"),
    min_perc: float = typer.Option(
        1.0,
        "--min-perc",
        help="Fold taxa below this percent abundance into Others",
    ),
    top_taxa: int = typer.Option(
        0,
        "--top-taxa",
        help="Keep the top N classified taxa and fold the rest into Others",
    ),
    no_unclassified: bool = typer.Option(
        False,
        "--no-unclassified",
        help="Exclude unclassified reads from the plot",
    ),
    title: Optional[str] = typer.Option(None, "--title", help="Plot title"),
    width: float = typer.Option(9.0, "--width", help="Static plot width in inches"),
    height: float = typer.Option(5.5, "--height", help="Static plot height in inches"),
    dpi: int = typer.Option(300, "--dpi", help="Static plot DPI"),
    kind: str = typer.Option(
        "stacked",
        "--kind",
        help="Plot kind: stacked or bubble",
    ),
):
    """Plot taxonomy composition across multiple Kraken reports."""
    missing = [input_file for input_file in input_files if not input_file.exists()]
    if missing:
        typer.echo(f"Error: Input file {missing[0]} does not exist.", err=True)
        raise typer.Exit(code=1)

    rank = _direct_default(rank, "S")
    metric = _direct_default(metric, "TOT")
    min_perc = _direct_default(min_perc, 1.0)
    top_taxa = _direct_default(top_taxa, 0)
    no_unclassified = _direct_default(no_unclassified, False)
    title = _direct_default(title, None)
    width = _direct_default(width, 9.0)
    height = _direct_default(height, 5.5)
    dpi = _direct_default(dpi, 300)
    kind = _direct_default(kind, "stacked")

    try:
        multi_report = build_multi_report(input_files)
        df = composition_dataframe(
            multi_report,
            rank=rank,
            metric=metric,
            min_perc=min_perc,
            top_taxa=top_taxa,
            no_unclassified=no_unclassified,
        )
        render_multi_composition(
            df,
            output_file,
            title=title,
            width=width,
            height=height,
            dpi=dpi,
            kind=kind,
        )
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
