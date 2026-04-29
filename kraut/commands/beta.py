from pathlib import Path
from typing import List, Optional

import typer

from kraut.models.beta import (
    abundance_dataframe,
    calculate_beta_diversity,
    distance_matrix_table,
    looks_like_abundance_table,
    read_abundance_table,
    render_beta_heatmap,
    render_beta_ordination,
)
from kraut.plotting import build_multi_report


def _direct_default(value, default):
    if isinstance(value, typer.models.OptionInfo):
        return default
    return value


def run(
    input_files: List[Path] = typer.Argument(
        ...,
        help="Input wide abundance table, or Kraken/Bracken report files",
    ),
    metric: str = typer.Option(
        "braycurtis",
        "--metric",
        "--method",
        "-m",
        help="Beta diversity metric: braycurtis, aitchison, or jaccard",
    ),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output distance matrix TSV file (default: stdout)",
    ),
    plot_file: Optional[Path] = typer.Option(
        None,
        "--plot",
        "-p",
        help="Optional beta distance heatmap (.html, .png, .pdf, .svg)",
    ),
    pca_file: Optional[Path] = typer.Option(
        None,
        "--pca",
        help="Optional 2D ordination plot (.html, .png, .pdf, .svg)",
    ),
    rank: str = typer.Option(
        "S",
        "--rank",
        "-r",
        help="Taxonomic rank to use when reading report files",
    ),
    abundance_metric: str = typer.Option(
        "TOT",
        "--abundance-metric",
        help="Report abundance metric: TOT (cumulative) or LVL (taxon-specific)",
    ),
    include_unclassified: bool = typer.Option(
        False,
        "--include-unclassified",
        help="Include unclassified reads as a feature",
    ),
    pseudocount: float = typer.Option(
        1.0,
        "--pseudocount",
        help="Pseudocount added before CLR for Aitchison distance",
    ),
    min_feature_count: float = typer.Option(
        0.0,
        "--min-feature-count",
        help="Keep taxa with total abundance at least this value",
    ),
    min_samples: int = typer.Option(
        1,
        "--min-samples",
        help="Keep taxa detected in at least this many samples",
    ),
    presence_threshold: float = typer.Option(
        0.0,
        "--presence-threshold",
        help=(
            "Presence cutoff for Jaccard; feature is present when value is "
            "greater than this"
        ),
    ),
    title: Optional[str] = typer.Option(None, "--title", help="Plot title"),
    width: float = typer.Option(6.0, "--width", help="Static plot width in inches"),
    height: float = typer.Option(5.0, "--height", help="Static plot height in inches"),
    dpi: int = typer.Option(300, "--dpi", help="Static plot DPI"),
):
    """Calculate beta diversity distances from Kraken or Bracken abundances."""
    metric = _direct_default(metric, "braycurtis")
    output_file = _direct_default(output_file, None)
    plot_file = _direct_default(plot_file, None)
    pca_file = _direct_default(pca_file, None)
    rank = _direct_default(rank, "S")
    abundance_metric = _direct_default(abundance_metric, "TOT")
    include_unclassified = _direct_default(include_unclassified, False)
    pseudocount = _direct_default(pseudocount, 1.0)
    min_feature_count = _direct_default(min_feature_count, 0.0)
    min_samples = _direct_default(min_samples, 1)
    presence_threshold = _direct_default(presence_threshold, 0.0)
    title = _direct_default(title, None)
    width = _direct_default(width, 6.0)
    height = _direct_default(height, 5.0)
    dpi = _direct_default(dpi, 300)

    missing = [input_file for input_file in input_files if not input_file.exists()]
    if missing:
        typer.echo(f"Error: Input file {missing[0]} does not exist.", err=True)
        raise typer.Exit(code=1)

    try:
        abundance_df = _load_abundance_dataframe(
            input_files,
            rank=rank,
            abundance_metric=abundance_metric,
            include_unclassified=include_unclassified,
        )
        result = calculate_beta_diversity(
            abundance_df,
            metric=metric,
            pseudocount=pseudocount,
            min_feature_count=min_feature_count,
            min_samples=min_samples,
            presence_threshold=presence_threshold,
        )

        for warning in result.warnings:
            typer.echo(f"Warning: {warning}", err=True)

        if plot_file:
            render_beta_heatmap(
                result.distance_df,
                plot_file,
                title=title,
                width=width,
                height=height,
                dpi=dpi,
            )
        if pca_file:
            render_beta_ordination(
                result.ordination_df,
                pca_file,
                title=title,
                width=width,
                height=height,
                dpi=dpi,
                kind=result.ordination_kind,
            )

        table = distance_matrix_table(result.distance_df).to_csv(sep="\t", index=False)
        if output_file:
            output_file.write_text(table)
        else:
            typer.echo(table, nl=False)
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc


def _load_abundance_dataframe(
    input_files: List[Path],
    rank: str,
    abundance_metric: str,
    include_unclassified: bool,
):
    if len(input_files) == 1 and looks_like_abundance_table(input_files[0]):
        return read_abundance_table(input_files[0], rank=rank)

    multi_report = build_multi_report(input_files)
    return abundance_dataframe(
        multi_report,
        rank=rank,
        abundance_metric=abundance_metric,
        include_unclassified=include_unclassified,
    )
