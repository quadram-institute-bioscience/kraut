from pathlib import Path
from typing import Optional

import pandas as pd
import typer


def _direct_default(value, default):
    if isinstance(value, typer.models.OptionInfo):
        return default
    return value


def run(
    input_file: Path = typer.Argument(
        ...,
        help="Input table produced by `kraut make-table`",
    ),
    top_taxa: int = typer.Option(
        10,
        "--top-taxa",
        min=1,
        help="Number of top taxa to report",
    ),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output report text file (default: stdout)",
    ),
):
    """Summarize a wide abundance table produced by make-table."""
    top_taxa = _direct_default(top_taxa, 10)
    output_file = _direct_default(output_file, None)

    if top_taxa < 1:
        typer.echo("Error: --top-taxa must be at least 1.", err=True)
        raise typer.Exit(code=1)

    if not input_file.exists():
        typer.echo(f"Error: Input file {input_file} does not exist.", err=True)
        raise typer.Exit(code=1)

    try:
        summary = summarize_table(input_file, top_taxa=top_taxa)
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if output_file:
        output_file.write_text(summary)
    else:
        typer.echo(summary, nl=False)


def summarize_table(input_file: Path, top_taxa: int = 10) -> str:
    try:
        df = pd.read_csv(input_file, sep="\t")
    except Exception as exc:
        raise ValueError(f"Could not read table '{input_file}': {exc}") from exc

    if "#Taxon" not in df.columns:
        raise ValueError("Input table must contain a #Taxon column")

    sample_cols = [col for col in df.columns if col != "#Taxon"]
    if not sample_cols:
        raise ValueError("Input table must contain at least one sample column")

    numeric_df = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    if numeric_df.isna().any().any():
        raise ValueError("Sample columns must contain only numeric values")

    is_count_table = _is_count_table(numeric_df)
    taxa_count = len(df)
    sample_count = len(sample_cols)

    lines = [
        "Table summary",
        f"Input table: {input_file}",
        f"Table size: {sample_count} samples x {taxa_count} taxa",
        f"Abundance type: {'counts' if is_count_table else 'relative'}",
    ]

    if is_count_table:
        lines.append("")
        lines.append("Total counts per sample:")
        totals = numeric_df.sum(axis=0)
        for sample in sample_cols:
            lines.append(f"{sample}: {int(totals[sample])}")

    lines.append("")
    top_rows = _top_normalised_taxa(df["#Taxon"], numeric_df, top_taxa)
    top_n = len(top_rows)
    lines.append(f"Top {top_n} taxa across the table (mean normalised abundance):")
    if top_rows.empty:
        lines.append("No taxa available.")
    else:
        for rank, row in enumerate(top_rows.itertuples(index=False), start=1):
            lines.append(f"{rank}. {row.taxon}: {row.score:.2f}%")

    return "\n".join(lines) + "\n"


def _is_count_table(df: pd.DataFrame) -> bool:
    return all(pd.api.types.is_integer_dtype(df[col]) for col in df.columns)


def _top_normalised_taxa(
    taxa: pd.Series,
    abundance_df: pd.DataFrame,
    top_taxa: int,
) -> pd.DataFrame:
    sample_totals = abundance_df.sum(axis=0)
    denominators = sample_totals.where(sample_totals != 0)
    normalised = abundance_df.div(denominators, axis=1).fillna(0.0)

    ranking = pd.DataFrame(
        {
            "taxon": taxa.astype(str),
            "score": normalised.mean(axis=1) * 100,
        }
    )
    ranking = ranking.sort_values(
        ["score", "taxon"],
        ascending=[False, True],
        kind="mergesort",
    )
    return ranking.head(top_taxa).reset_index(drop=True)
