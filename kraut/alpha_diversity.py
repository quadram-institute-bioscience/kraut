from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from kraut.models.multi_report import MultiKrakenReport
from kraut.plotting import PALETTE, output_kind


CORE_METRICS = [
    "observed_features",
    "shannon",
    "simpson",
    "inv_simpson",
    "pielou_e",
    "dominance",
    "goods_coverage",
]

METRIC_PRESETS = {
    "core": CORE_METRICS,
    "richness": [
        "observed_features",
        "chao1",
        "ace",
        "margalef",
        "menhinick",
        "singles",
        "doubles",
    ],
    "diversity": [
        "shannon",
        "simpson",
        "inv_simpson",
        "enspie",
        "brillouin_d",
        "fisher_alpha",
    ],
    "evenness": [
        "pielou_e",
        "simpson_e",
        "heip_e",
        "mcintosh_e",
    ],
    "dominance": [
        "dominance",
        "berger_parker_d",
        "simpson_d",
        "mcintosh_d",
        "strong",
    ],
    "coverage": [
        "goods_coverage",
        "robbins",
    ],
}
METRIC_PRESETS["extended"] = list(
    dict.fromkeys(
        CORE_METRICS
        + METRIC_PRESETS["richness"]
        + METRIC_PRESETS["diversity"]
        + METRIC_PRESETS["evenness"]
        + METRIC_PRESETS["dominance"]
        + METRIC_PRESETS["coverage"]
    )
)


def select_alpha_metrics(
    metrics: str = "core",
    add_metrics: Optional[str] = None,
) -> list[str]:
    """Resolve an alpha metric preset plus optional comma-separated additions."""
    selected = _resolve_metric_argument(metrics)
    selected.extend(_split_metric_list(add_metrics))
    selected = [metric.lower() for metric in selected]
    return list(dict.fromkeys(selected))


def abundance_dataframe(
    multi_report: MultiKrakenReport,
    rank: str = "S",
    metric: str = "TOT",
    include_unclassified: bool = False,
    min_perc: float = 0.0,
) -> pd.DataFrame:
    """Build the taxon-count matrix used as input to alpha diversity metrics."""
    if min_perc < 0:
        raise ValueError("--min-perc must be >= 0")

    metric = metric.upper()
    if metric == "COUNTS":
        metric = "TOT"
    if metric not in {"TOT", "LVL"}:
        raise ValueError("--metric must be one of: TOT, LVL")

    rank = rank.upper()
    source_key = "taxon_counts" if metric == "LVL" else "clade_counts"
    sample_cols = list(multi_report.samples)
    rows = []

    for tax_id, info in multi_report.data.items():
        is_unclassified = tax_id == 0
        if is_unclassified:
            if not include_unclassified:
                continue
        elif rank != "ALL" and info["rank"].upper() != rank:
            continue

        row = {"#Taxon": str(tax_id)}
        counts = info[source_key]
        for sample_idx, sample in enumerate(sample_cols):
            row[sample] = counts.get(sample_idx, 0)
        rows.append(row)

    df = pd.DataFrame(rows, columns=["#Taxon"] + sample_cols)
    if df.empty:
        return df

    for sample in sample_cols:
        df[sample] = pd.to_numeric(df[sample], errors="coerce").fillna(0).astype(int)

    if min_perc > 0:
        df = _filter_by_min_perc(df, sample_cols, min_perc)

    return df


def calculate_alpha_diversity(
    abundance_df: pd.DataFrame,
    metrics: Iterable[str],
    validate: bool = True,
) -> pd.DataFrame:
    """Calculate one or more scikit-bio alpha diversity metrics."""
    sample_cols = [col for col in abundance_df.columns if col != "#Taxon"]
    if abundance_df.empty or not sample_cols:
        raise ValueError("No taxa found for the requested alpha diversity options")

    counts = abundance_df.set_index("#Taxon")[sample_cols].T
    empty_samples = counts.sum(axis=1)
    empty_samples = empty_samples[empty_samples <= 0]
    if not empty_samples.empty:
        names = ", ".join(empty_samples.index)
        raise ValueError(f"No non-zero counts found for sample(s): {names}")

    metric_names = list(metrics)
    _validate_alpha_metrics(metric_names)
    alpha_diversity, _ = _load_skbio_diversity()

    result = pd.DataFrame(index=metric_names, columns=sample_cols, dtype=object)
    for metric_name in metric_names:
        try:
            series = alpha_diversity(
                metric_name,
                counts,
                ids=sample_cols,
                validate=validate,
            )
        except Exception as exc:
            raise ValueError(
                f"Could not calculate alpha metric '{metric_name}': {exc}"
            ) from exc
        result.loc[metric_name, sample_cols] = series.reindex(sample_cols).tolist()

    result.insert(0, "#Metric", result.index)
    return result.reset_index(drop=True)


def render_alpha_diversity(
    alpha_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str] = None,
    width: float = 9.0,
    height: float = 2.2,
    dpi: int = 300,
) -> None:
    """Render alpha diversity values to an HTML or static bar plot."""
    plot_df = _numeric_alpha_dataframe(alpha_df)
    if plot_df.empty:
        raise ValueError("No numeric alpha diversity metrics found to plot")

    if output_kind(output_file) == "html":
        _render_alpha_html(plot_df, output_file, title)
    else:
        _render_alpha_static(plot_df, output_file, title, width, height, dpi)


def _resolve_metric_argument(metrics: str) -> list[str]:
    values = _split_metric_list(metrics)
    if not values:
        return list(CORE_METRICS)
    if len(values) == 1 and values[0].lower() in METRIC_PRESETS:
        return list(METRIC_PRESETS[values[0].lower()])
    return values


def _split_metric_list(value: Optional[str]) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def _filter_by_min_perc(
    df: pd.DataFrame,
    sample_cols: list[str],
    min_perc: float,
) -> pd.DataFrame:
    col_totals = {sample: df[sample].sum() for sample in sample_cols}

    def row_max_perc(row) -> float:
        return max(
            (
                row[sample] / col_totals[sample] * 100
                if col_totals[sample] > 0
                else 0.0
            )
            for sample in sample_cols
        )

    return df[df.apply(row_max_perc, axis=1) >= min_perc].reset_index(drop=True)


def _validate_alpha_metrics(metric_names: list[str]) -> None:
    _, get_alpha_diversity_metrics = _load_skbio_diversity()
    available = set(get_alpha_diversity_metrics())
    unknown = [metric for metric in metric_names if metric not in available]
    if unknown:
        supported = ", ".join(sorted(available))
        requested = ", ".join(unknown)
        raise ValueError(
            f"Unknown alpha metric(s): {requested}. Supported: {supported}"
        )


def _load_skbio_diversity():
    try:
        from skbio.diversity import alpha_diversity, get_alpha_diversity_metrics
    except ImportError as exc:
        raise RuntimeError(
            "scikit-bio is required for `kraut alpha`; "
            "install krautils with dependencies."
        ) from exc
    return alpha_diversity, get_alpha_diversity_metrics


def _numeric_alpha_dataframe(alpha_df: pd.DataFrame) -> pd.DataFrame:
    if "#Metric" not in alpha_df.columns:
        raise ValueError("Alpha diversity table must contain a #Metric column")
    plot_df = alpha_df.set_index("#Metric").apply(pd.to_numeric, errors="coerce")
    return plot_df.dropna(axis=0, how="any")


def _render_alpha_static(
    plot_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
    width: float,
    height: float,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = plot_df.index.tolist()
    samples = plot_df.columns.tolist()
    fig_height = max(height * len(metrics), 2.4)
    fig, axes = plt.subplots(
        len(metrics),
        1,
        figsize=(width, fig_height),
        squeeze=False,
        sharex=True,
    )

    for idx, metric_name in enumerate(metrics):
        ax = axes[idx][0]
        ax.bar(
            samples,
            plot_df.loc[metric_name].astype(float),
            color=PALETTE[idx % len(PALETTE)],
        )
        ax.set_ylabel(metric_name)
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.8)
        ax.set_axisbelow(True)

    axes[0][0].set_title(title or "Alpha diversity")
    axes[-1][0].set_xlabel("Sample")
    axes[-1][0].tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _render_alpha_html(
    plot_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    metrics = plot_df.index.tolist()
    samples = plot_df.columns.tolist()
    fig = make_subplots(
        rows=len(metrics),
        cols=1,
        shared_xaxes=True,
        subplot_titles=metrics,
        vertical_spacing=min(0.08, 0.4 / max(len(metrics), 1)),
    )

    for idx, metric_name in enumerate(metrics, start=1):
        values = plot_df.loc[metric_name].astype(float).tolist()
        fig.add_trace(
            go.Bar(
                x=samples,
                y=values,
                marker_color=PALETTE[(idx - 1) % len(PALETTE)],
                hovertemplate=(
                    "Sample: %{x}<br>"
                    f"{metric_name}: "
                    "%{y}<extra></extra>"
                ),
                showlegend=False,
            ),
            row=idx,
            col=1,
        )
        fig.update_yaxes(title_text=metric_name, row=idx, col=1)

    fig.update_layout(
        title=title or "Alpha diversity",
        height=max(260 * len(metrics), 360),
    )
    fig.update_xaxes(title_text="Sample", row=len(metrics), col=1)
    fig.write_html(output_file)
