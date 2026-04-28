from pathlib import Path
from typing import Optional

import typer

from kraut.plotting import (
    build_multi_report,
    composition_dataframe,
    render_single_composition,
)


def _direct_default(value, default):
    if isinstance(value, typer.models.OptionInfo):
        return default
    return value


def run(
    input_file: Path = typer.Argument(..., help="Input Kraken report file"),
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
    width: float = typer.Option(7.0, "--width", help="Static plot width in inches"),
    height: float = typer.Option(5.0, "--height", help="Static plot height in inches"),
    dpi: int = typer.Option(300, "--dpi", help="Static plot DPI"),
):
    """Plot taxonomy composition for a single Kraken report."""
    if not input_file.exists():
        typer.echo(f"Error: Input file {input_file} does not exist.", err=True)
        raise typer.Exit(code=1)

    rank = _direct_default(rank, "S")
    metric = _direct_default(metric, "TOT")
    min_perc = _direct_default(min_perc, 1.0)
    top_taxa = _direct_default(top_taxa, 0)
    no_unclassified = _direct_default(no_unclassified, False)
    title = _direct_default(title, None)
    width = _direct_default(width, 7.0)
    height = _direct_default(height, 5.0)
    dpi = _direct_default(dpi, 300)

    try:
        multi_report = build_multi_report([input_file])
        df = composition_dataframe(
            multi_report,
            rank=rank,
            metric=metric,
            min_perc=min_perc,
            top_taxa=top_taxa,
            no_unclassified=no_unclassified,
        )
        render_single_composition(
            df,
            output_file,
            title=title,
            width=width,
            height=height,
            dpi=dpi,
        )
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
