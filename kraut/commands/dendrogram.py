from pathlib import Path
from typing import List, Optional

import typer

from kraut.dendrogram import metadata_color_annotations, render_dendrogram
from kraut.models.beta import (
    abundance_dataframe,
    calculate_beta_diversity,
    looks_like_abundance_table,
    read_abundance_table,
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
    output_file: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output dendrogram plot file (.html, .png, .pdf, .svg)",
    ),
    distance: str = typer.Option(
        "braycurtis",
        "--distance",
        "-d",
        help="Distance method: braycurtis, aitchison, or jaccard",
    ),
    clustering: str = typer.Option(
        "ward",
        "--clustering",
        "-c",
        help="Clustering method: ward, average, single, or complete",
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
    metadata_file: Optional[Path] = typer.Option(
        None,
        "--metadata",
        help="Optional sample metadata CSV or TSV file; first column is sample ID",
    ),
    color_by: Optional[str] = typer.Option(
        None,
        "--color-by",
        help="Metadata column used to color sample labels",
    ),
    title: Optional[str] = typer.Option(None, "--title", help="Plot title"),
    width: float = typer.Option(8.0, "--width", help="Static plot width in inches"),
    height: float = typer.Option(5.0, "--height", help="Static plot height in inches"),
    dpi: int = typer.Option(300, "--dpi", help="Static plot DPI"),
):
    """Plot a hierarchical clustering dendrogram from Kraken or Bracken reports."""
    output_file = _direct_default(output_file, None)
    distance = _direct_default(distance, "braycurtis")
    clustering = _direct_default(clustering, "ward")
    rank = _direct_default(rank, "S")
    abundance_metric = _direct_default(abundance_metric, "TOT")
    include_unclassified = _direct_default(include_unclassified, False)
    pseudocount = _direct_default(pseudocount, 1.0)
    min_feature_count = _direct_default(min_feature_count, 0.0)
    min_samples = _direct_default(min_samples, 1)
    presence_threshold = _direct_default(presence_threshold, 0.0)
    metadata_file = _direct_default(metadata_file, None)
    color_by = _direct_default(color_by, None)
    title = _direct_default(title, None)
    width = _direct_default(width, 8.0)
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
            metric=distance,
            pseudocount=pseudocount,
            min_feature_count=min_feature_count,
            min_samples=min_samples,
            presence_threshold=presence_threshold,
        )

        for warning in result.warnings:
            typer.echo(f"Warning: {warning}", err=True)

        colors = metadata_color_annotations(
            result.distance_df.index.astype(str).tolist(),
            metadata_file,
            color_by,
        )
        render_dendrogram(
            result.distance_df,
            output_file,
            clustering=clustering,
            metadata_colors=colors,
            title=title,
            width=width,
            height=height,
            dpi=dpi,
        )
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
