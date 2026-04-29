from pathlib import Path
from typing import List, Optional

import typer

from kraut.models.kraken_data import KrakenReport
from kraut.models.ranks import CANONICAL_RANKS, rank_read_counts
from kraut.plotting import PALETTE, output_kind


RANK_ROWS = ("U",) + CANONICAL_RANKS
RANK_COLORS = {
    "U": "#6C757D",
    "R": "#D0D0D0",
    **{
        rank: PALETTE[idx % len(PALETTE)]
        for idx, rank in enumerate(CANONICAL_RANKS[1:])
    },
}


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
    plot_file: Optional[Path] = typer.Option(
        None,
        "--plot",
        "-p",
        help="Write a stacked rank-composition bar chart (.html, .png, .pdf, .svg)",
    ),
):
    """Profile reads assigned at each canonical taxonomic rank."""
    output_file = _direct_default(output_file, None)
    counts = _direct_default(counts, False)
    plot_file = _direct_default(plot_file, None)

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

    try:
        table = _rank_table(sample_names, sample_counts, as_counts=counts)
        if plot_file:
            _render_rank_plot(sample_names, sample_counts, plot_file)

        if output_file:
            output_file.write_text(table)
        else:
            typer.echo(table, nl=False)
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc


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


def _rank_percentages(
    sample_counts: list[dict[str, int]],
) -> dict[str, list[float]]:
    totals = [sum(counts.values()) for counts in sample_counts]
    percentages = {}
    for rank in RANK_ROWS:
        percentages[rank] = [
            counts.get(rank, 0) / total * 100 if total > 0 else 0.0
            for counts, total in zip(sample_counts, totals)
        ]
    return percentages


def _render_rank_plot(
    sample_names: list[str],
    sample_counts: list[dict[str, int]],
    output_file: Path,
) -> None:
    kind = output_kind(output_file)
    percentages = _rank_percentages(sample_counts)
    if kind == "html":
        _render_rank_html(sample_names, percentages, output_file)
    else:
        _render_rank_static(sample_names, percentages, output_file)


def _render_rank_html(
    sample_names: list[str],
    percentages: dict[str, list[float]],
    output_file: Path,
) -> None:
    import plotly.graph_objects as go

    fig = go.Figure()
    for rank in RANK_ROWS:
        fig.add_bar(
            name=rank,
            x=sample_names,
            y=percentages[rank],
            marker_color=RANK_COLORS[rank],
            hovertemplate="%{x}<br>%{y:.2f}%<extra>" + rank + "</extra>",
        )
    fig.update_layout(
        title="Rank composition",
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Composition (%)",
        yaxis_range=[0, 100],
    )
    fig.write_html(output_file)


def _render_rank_static(
    sample_names: list[str],
    percentages: dict[str, list[float]],
    output_file: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    width = max(7.0, len(sample_names) * 0.6)
    fig, ax = plt.subplots(figsize=(width, 5.0))
    bottoms = [0.0] * len(sample_names)
    for rank in RANK_ROWS:
        values = percentages[rank]
        ax.bar(
            sample_names,
            values,
            bottom=bottoms,
            label=rank,
            color=RANK_COLORS[rank],
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]

    ax.set_ylim(0, 100)
    ax.set_ylabel("Composition (%)")
    ax.set_xlabel("Sample")
    ax.set_title("Rank composition")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _sample_name_from_path(path: Path) -> str:
    sample_name = path.stem
    if sample_name.endswith(".krep"):
        sample_name = Path(sample_name).stem
    return sample_name
