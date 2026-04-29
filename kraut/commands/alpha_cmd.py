from pathlib import Path
from typing import List, Optional

import typer

from kraut.alpha_diversity import (
    abundance_dataframe,
    calculate_alpha_diversity,
    render_alpha_diversity,
    select_alpha_metrics,
)
from kraut.plotting import build_multi_report


def _direct_default(value, default):
    if isinstance(value, typer.models.OptionInfo):
        return default
    return value


def run(
    input_files: List[Path] = typer.Argument(
        ...,
        help="Input Kraken or Bracken report files",
    ),
    output_table: Path = typer.Option(
        Path("alpha.tsv"),
        "--output-table",
        "-o",
        help="Output alpha diversity table",
    ),
    plot_file: Optional[Path] = typer.Option(
        None,
        "--plot",
        "-p",
        help="Optional alpha diversity plot file (.html, .png, .svg, .pdf)",
    ),
    rank: str = typer.Option("S", "--rank", "-r", help="Taxonomic rank to use"),
    metric: str = typer.Option(
        "TOT",
        "--metric",
        "-m",
        help="Abundance metric: TOT (cumulative) or LVL (taxon-specific)",
    ),
    metrics: str = typer.Option(
        "core",
        "--metrics",
        help="Alpha metric preset or comma-separated metric list",
    ),
    add_metrics: Optional[str] = typer.Option(
        None,
        "--add-metrics",
        help="Comma-separated additional alpha metrics to append",
    ),
    include_unclassified: bool = typer.Option(
        False,
        "--include-unclassified",
        help="Include unclassified reads as a feature",
    ),
    min_perc: float = typer.Option(
        0.0,
        "--min-perc",
        help="Remove taxa where no sample reaches this percent abundance",
    ),
    title: Optional[str] = typer.Option(None, "--title", help="Plot title"),
    width: float = typer.Option(9.0, "--width", help="Static plot width in inches"),
    height: float = typer.Option(
        2.2,
        "--height",
        help="Static plot height in inches per metric",
    ),
    dpi: int = typer.Option(300, "--dpi", help="Static plot DPI"),
):
    """Calculate alpha diversity for Kraken or Bracken report files."""
    output_table = _direct_default(output_table, Path("alpha.tsv"))
    plot_file = _direct_default(plot_file, None)
    rank = _direct_default(rank, "S")
    metric = _direct_default(metric, "TOT")
    metrics = _direct_default(metrics, "core")
    add_metrics = _direct_default(add_metrics, None)
    include_unclassified = _direct_default(include_unclassified, False)
    min_perc = _direct_default(min_perc, 0.0)
    title = _direct_default(title, None)
    width = _direct_default(width, 9.0)
    height = _direct_default(height, 2.2)
    dpi = _direct_default(dpi, 300)

    missing = [input_file for input_file in input_files if not input_file.exists()]
    if missing:
        typer.echo(f"Error: Input file {missing[0]} does not exist.", err=True)
        raise typer.Exit(code=1)

    try:
        multi_report = build_multi_report(input_files)
        counts_df = abundance_dataframe(
            multi_report,
            rank=rank,
            metric=metric,
            include_unclassified=include_unclassified,
            min_perc=min_perc,
        )
        metric_names = select_alpha_metrics(metrics, add_metrics)
        alpha_df = calculate_alpha_diversity(counts_df, metric_names)
        alpha_df.to_csv(output_table, sep="\t", index=False)

        if plot_file:
            render_alpha_diversity(
                alpha_df,
                plot_file,
                title=title,
                width=width,
                height=height,
                dpi=dpi,
            )
    except (RuntimeError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
