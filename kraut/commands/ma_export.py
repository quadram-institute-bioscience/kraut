from pathlib import Path
from typing import List, Optional

import typer

from kraut.models.ma_export import export_microbiome_analyst


def run(
    input_files: List[Path] = typer.Argument(..., help="Input Kraken report files"),
    outdir: Path = typer.Option(
        ...,
        "--outdir",
        "-o",
        help="Output directory for MicrobiomeAnalyst files",
    ),
    metadata: Optional[Path] = typer.Option(
        None,
        "--metadata",
        help="Optional metadata CSV or TSV file",
    ),
    metadata_sample_col: Optional[str] = typer.Option(
        None,
        "--metadata-sample-col",
        help="Column in metadata matching input sample names",
    ),
    pseudo_col: str = typer.Option(
        "Random_label",
        "--pseudo-col",
        help="Metadata column to create or overwrite with pseudo group labels",
    ),
    metric: str = typer.Option(
        "LVL",
        "--metric",
        "-m",
        help="Count metric: LVL for taxon counts, TOT for clade counts",
    ),
):
    """
    Export Kraken reports as MicrobiomeAnalyst counts, taxonomy, metadata, and tree files.
    """
    try:
        result = export_microbiome_analyst(
            input_files=input_files,
            outdir=outdir,
            metadata_file=metadata,
            metadata_sample_col=metadata_sample_col,
            pseudo_col=pseudo_col,
            metric=metric,
        )
    except (OSError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    for warning in result.warnings:
        typer.echo(f"Warning: {warning}", err=True)

    typer.echo(f"Created {result.counts_file}")
    typer.echo(f"Created {result.taxonomy_file}")
    typer.echo(f"Created {result.metadata_file}")
    typer.echo(f"Created {result.tree_file}")
