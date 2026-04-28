from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

from kraut.models.kraken_data import KrakenReport
from kraut.models.multi_report import MultiKrakenReport


STATIC_SUFFIXES = {".png", ".pdf", ".svg"}
HTML_SUFFIXES = {".html"}
UNCLASSIFIED_LABEL = "unclassified"
OTHERS_LABEL = "Others"

PALETTE = [
    "#4E79A7",
    "#F28E2B",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#E377C2",
    "#7F7F7F",
    "#BCBD22",
    "#17BECF",
]


def sample_name_from_path(path: Path) -> str:
    sample_name = path.stem
    if sample_name.endswith(".krep"):
        sample_name = Path(sample_name).stem
    return sample_name


def output_kind(output_file: Path) -> str:
    suffix = output_file.suffix.lower()
    if suffix in HTML_SUFFIXES:
        return "html"
    if suffix in STATIC_SUFFIXES:
        return "static"
    supported = ", ".join(sorted(HTML_SUFFIXES | STATIC_SUFFIXES))
    raise ValueError(f"Unsupported output extension '{suffix}'. Supported: {supported}")


def build_multi_report(input_files: List[Path]) -> MultiKrakenReport:
    multi_report = MultiKrakenReport()
    for input_file in input_files:
        report = KrakenReport.from_file(str(input_file))
        multi_report.add_report(report, sample_name_from_path(input_file))
    return multi_report


def composition_dataframe(
    multi_report: MultiKrakenReport,
    rank: str = "S",
    metric: str = "TOT",
    min_perc: float = 1.0,
    top_taxa: int = 0,
    no_unclassified: bool = False,
) -> pd.DataFrame:
    if min_perc < 0:
        raise ValueError("--min-perc must be >= 0")
    if top_taxa < 0:
        raise ValueError("--top-taxa must be >= 0")

    metric = metric.upper()
    if metric == "COUNTS":
        metric = "TOT"
    if metric not in {"TOT", "LVL"}:
        raise ValueError("--metric must be one of: TOT, LVL")

    rows = _raw_rows(multi_report, rank.upper(), metric, no_unclassified)
    sample_cols = list(multi_report.samples)
    if not rows:
        return pd.DataFrame(columns=["#Taxon"] + sample_cols)

    df = pd.DataFrame(rows)
    for sample in sample_cols:
        total = df[sample].sum()
        if total > 0:
            df[sample] = df[sample] / total * 100
        else:
            df[sample] = 0.0

    classified = df[~df["_is_unclassified"]].copy()
    unclassified = df[df["_is_unclassified"]].copy()
    selected = _select_classified_taxa(classified, sample_cols, min_perc, top_taxa)
    selected_ids = set(selected["_tax_id"])
    others = classified[~classified["_tax_id"].isin(selected_ids)]

    output_rows = []
    if not no_unclassified and not unclassified.empty:
        output_rows.extend(_strip_internal_columns(unclassified).to_dict("records"))

    if not selected.empty:
        output_rows.extend(_strip_internal_columns(selected).to_dict("records"))

    if not others.empty:
        others_row = {"#Taxon": OTHERS_LABEL}
        for sample in sample_cols:
            others_row[sample] = others[sample].sum()
        if any(others_row[sample] > 0 for sample in sample_cols):
            output_rows.append(others_row)

    result = pd.DataFrame(output_rows, columns=["#Taxon"] + sample_cols)
    if result.empty:
        return pd.DataFrame(columns=["#Taxon"] + sample_cols)
    return result


def render_single_composition(
    df: pd.DataFrame,
    output_file: Path,
    title: Optional[str] = None,
    width: float = 7.0,
    height: float = 5.0,
    dpi: int = 300,
) -> None:
    _require_data(df)
    kind = output_kind(output_file)
    sample = _sample_columns(df)[0]
    if kind == "html":
        _render_single_html(df, sample, output_file, title)
    else:
        _render_single_static(df, sample, output_file, title, width, height, dpi)


def render_multi_composition(
    df: pd.DataFrame,
    output_file: Path,
    title: Optional[str] = None,
    width: float = 9.0,
    height: float = 5.5,
    dpi: int = 300,
) -> None:
    _require_data(df)
    kind = output_kind(output_file)
    if kind == "html":
        _render_multi_html(df, output_file, title)
    else:
        _render_multi_static(df, output_file, title, width, height, dpi)


def _raw_rows(
    multi_report: MultiKrakenReport,
    rank: str,
    metric: str,
    no_unclassified: bool,
) -> List[dict]:
    source_key = "taxon_counts" if metric == "LVL" else "clade_counts"
    rows = []

    for tax_id, info in multi_report.data.items():
        is_unclassified = tax_id == 0
        if is_unclassified:
            if no_unclassified:
                continue
        elif info["rank"].upper() != rank:
            continue

        row = {
            "#Taxon": info["name"].strip(),
            "_tax_id": tax_id,
            "_is_unclassified": is_unclassified,
        }
        counts = info[source_key]
        for sample_idx, sample in enumerate(multi_report.samples):
            row[sample] = counts.get(sample_idx, 0)
        rows.append(row)

    return rows


def _select_classified_taxa(
    classified: pd.DataFrame,
    sample_cols: List[str],
    min_perc: float,
    top_taxa: int,
) -> pd.DataFrame:
    if classified.empty:
        return classified

    selected = classified.copy()
    selected["_total"] = selected[sample_cols].sum(axis=1)
    selected["_max"] = selected[sample_cols].max(axis=1)
    selected = selected[selected["_max"] >= min_perc]
    selected = selected.sort_values(["_total", "#Taxon"], ascending=[False, True])

    if top_taxa > 0:
        selected = selected.head(top_taxa)

    return selected.drop(columns=["_total", "_max"])


def _strip_internal_columns(df: pd.DataFrame) -> pd.DataFrame:
    internal_cols = [col for col in df.columns if col.startswith("_")]
    return df.drop(columns=internal_cols)


def _require_data(df: pd.DataFrame) -> None:
    if df.empty or not _sample_columns(df):
        raise ValueError("No plottable taxa found for the requested options")
    if df[_sample_columns(df)].sum().sum() <= 0:
        raise ValueError("No non-zero abundance values found for the requested options")


def _sample_columns(df: pd.DataFrame) -> List[str]:
    return [col for col in df.columns if col != "#Taxon"]


def _category_colors(categories: List[str]) -> List[str]:
    colors = []
    palette_idx = 0
    for category in categories:
        if category == UNCLASSIFIED_LABEL:
            colors.append("#6C757D")
        elif category == OTHERS_LABEL:
            colors.append("#D0D0D0")
        else:
            colors.append(PALETTE[palette_idx % len(PALETTE)])
            palette_idx += 1
    return colors


def _single_values(df: pd.DataFrame, sample: str) -> Tuple[List[str], List[float]]:
    labels = df["#Taxon"].tolist()
    values = df[sample].astype(float).tolist()
    pairs = [(label, value) for label, value in zip(labels, values) if value > 0]
    return [label for label, _ in pairs], [value for _, value in pairs]


def _render_single_static(
    df: pd.DataFrame,
    sample: str,
    output_file: Path,
    title: Optional[str],
    width: float,
    height: float,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels, values = _single_values(df, sample)
    colors = _category_colors(labels)
    fig, ax = plt.subplots(figsize=(width, height))
    wedges, _ = ax.pie(
        values,
        colors=colors,
        startangle=90,
        counterclock=False,
        wedgeprops={"width": 0.45, "edgecolor": "white"},
    )
    legend_labels = [f"{label} ({value:.1f}%)" for label, value in zip(labels, values)]
    ax.legend(
        wedges,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
    )
    ax.text(0, 0, sample, ha="center", va="center", fontsize=11)
    ax.set_title(title or sample)
    ax.set(aspect="equal")
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _render_multi_static(
    df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
    width: float,
    height: float,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sample_cols = _sample_columns(df)
    categories = df["#Taxon"].tolist()
    colors = _category_colors(categories)
    fig, ax = plt.subplots(figsize=(width, height))
    bottoms = [0.0] * len(sample_cols)

    for category, color in zip(categories, colors):
        values = df.loc[df["#Taxon"] == category, sample_cols].iloc[0].astype(float).tolist()
        ax.bar(sample_cols, values, bottom=bottoms, label=category, color=color)
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]

    ax.set_ylim(0, 100)
    ax.set_ylabel("Composition (%)")
    ax.set_xlabel("Sample")
    ax.set_title(title or "Taxonomy composition")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _render_single_html(
    df: pd.DataFrame,
    sample: str,
    output_file: Path,
    title: Optional[str],
) -> None:
    import plotly.graph_objects as go

    labels, values = _single_values(df, sample)
    fig = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.45,
                marker={"colors": _category_colors(labels)},
                hovertemplate="%{label}<br>%{percent}<extra></extra>",
            )
        ]
    )
    fig.update_layout(title=title or sample)
    fig.write_html(output_file)


def _render_multi_html(
    df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
) -> None:
    import plotly.graph_objects as go

    sample_cols = _sample_columns(df)
    categories = df["#Taxon"].tolist()
    colors = _category_colors(categories)
    fig = go.Figure()
    for category, color in zip(categories, colors):
        values = (
            df.loc[df["#Taxon"] == category, sample_cols]
            .iloc[0]
            .astype(float)
            .tolist()
        )
        fig.add_bar(
            name=category,
            x=sample_cols,
            y=values,
            marker_color=color,
            hovertemplate="%{x}<br>%{y:.2f}%<extra>" + category + "</extra>",
        )

    fig.update_layout(
        title=title or "Taxonomy composition",
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Composition (%)",
        yaxis_range=[0, 100],
    )
    fig.write_html(output_file)
