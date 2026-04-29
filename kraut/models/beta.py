from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from kraut.models.multi_report import MultiKrakenReport
from kraut.plotting import PALETTE, output_kind

BETA_METRICS = {"braycurtis", "aitchison", "jaccard"}
SPARSE_WARNING_THRESHOLD = 0.8


@dataclass
class BetaDiversityResult:
    distance_df: pd.DataFrame
    ordination_df: pd.DataFrame
    ordination_kind: str
    warnings: list[str]


def abundance_dataframe(
    multi_report: MultiKrakenReport,
    rank: str = "S",
    abundance_metric: str = "TOT",
    include_unclassified: bool = False,
) -> pd.DataFrame:
    """Build a feature-by-sample abundance table from parsed Kraken reports."""
    abundance_metric = abundance_metric.upper()
    if abundance_metric == "COUNTS":
        abundance_metric = "TOT"
    if abundance_metric not in {"TOT", "LVL"}:
        raise ValueError("--abundance-metric must be one of: TOT, LVL")

    rank = rank.upper()
    source_key = "taxon_counts" if abundance_metric == "LVL" else "clade_counts"
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

    return pd.DataFrame(rows, columns=["#Taxon"] + sample_cols)


def read_abundance_table(input_file: Path, rank: Optional[str] = None) -> pd.DataFrame:
    """Read a wide feature table such as the output of `kraut make-table`."""
    try:
        df = pd.read_csv(input_file, sep="\t")
    except Exception as exc:
        raise ValueError(
            f"Could not read abundance table '{input_file}': {exc}"
        ) from exc

    if df.empty:
        raise ValueError("Abundance table is empty")

    if _is_combined_bracken_table(df.columns):
        return _combined_bracken_abundance_dataframe(df, rank)

    feature_col = "#Taxon" if "#Taxon" in df.columns else df.columns[0]
    numeric_cols = []
    for col in df.columns:
        if col == feature_col:
            continue
        values = pd.to_numeric(df[col], errors="coerce")
        if values.notna().all():
            numeric_cols.append(col)

    if len(numeric_cols) < 2:
        raise ValueError(
            "Abundance table must contain at least two numeric sample columns"
        )

    result = df[[feature_col] + numeric_cols].copy()
    result.rename(columns={feature_col: "#Taxon"}, inplace=True)
    result["#Taxon"] = result["#Taxon"].astype(str)
    for col in numeric_cols:
        result[col] = pd.to_numeric(result[col], errors="raise")
    return result


def looks_like_abundance_table(input_file: Path) -> bool:
    """Return True when a single input file looks like a wide sample table."""
    try:
        with input_file.open() as handle:
            first_line = handle.readline()
    except OSError:
        return False

    columns = first_line.rstrip("\n").split("\t")
    if not columns:
        return False
    return (
        "#Taxon" in columns
        or columns[0] in {"Tax", "taxon", "Taxon"}
        or _is_combined_bracken_table(columns)
    )


def _is_combined_bracken_table(columns) -> bool:
    column_set = set(columns)
    return {"name", "taxonomy_id", "taxonomy_lvl"}.issubset(column_set) and len(
        [col for col in columns if col.endswith("_num")]
    ) >= 2


def _combined_bracken_abundance_dataframe(
    df: pd.DataFrame,
    rank: Optional[str],
) -> pd.DataFrame:
    selected = df.copy()
    if rank and rank.upper() != "ALL":
        selected = selected[
            selected["taxonomy_lvl"].astype(str).str.upper() == rank.upper()
        ]

    num_cols = [col for col in selected.columns if col.endswith("_num")]
    if len(num_cols) < 2:
        raise ValueError(
            "Combined Bracken table must contain at least two *_num columns"
        )

    result = pd.DataFrame({"#Taxon": selected["taxonomy_id"].astype(str)})
    for col in num_cols:
        result[_bracken_sample_name(col)] = pd.to_numeric(
            selected[col],
            errors="raise",
        )
    return result


def _bracken_sample_name(column: str) -> str:
    sample = column.removesuffix("_num")
    return sample.removesuffix(".brout")


def calculate_beta_diversity(
    abundance_df: pd.DataFrame,
    metric: str = "braycurtis",
    pseudocount: float = 1.0,
    min_feature_count: float = 0.0,
    min_samples: int = 1,
    presence_threshold: float = 0.0,
) -> BetaDiversityResult:
    """Calculate a beta-diversity distance matrix and 2D ordination."""
    metric = metric.lower()
    if metric not in BETA_METRICS:
        supported = ", ".join(sorted(BETA_METRICS))
        raise ValueError(f"--metric must be one of: {supported}")

    matrix = _prepare_abundance_matrix(
        abundance_df,
        min_feature_count=min_feature_count,
        min_samples=min_samples,
    )
    warnings = []

    if metric == "braycurtis":
        distance_df = _braycurtis_distance(matrix)
        ordination_df = _pcoa(distance_df)
        ordination_kind = "PCoA"
    elif metric == "aitchison":
        if pseudocount <= 0:
            raise ValueError("--pseudocount must be > 0")
        zero_fraction = float((matrix == 0).sum().sum()) / float(matrix.size)
        if zero_fraction >= SPARSE_WARNING_THRESHOLD:
            warnings.append(
                "Aitchison input is very sparse; distances are sensitive to "
                "zero handling and the chosen pseudocount."
            )
        clr_matrix = _clr_transform(matrix, pseudocount)
        distance_df = _euclidean_distance(clr_matrix)
        ordination_df = _pca(clr_matrix)
        ordination_kind = "PCA"
    else:
        if presence_threshold < 0:
            raise ValueError("--presence-threshold must be >= 0")
        binary = matrix > presence_threshold
        _require_non_empty_samples(binary.astype(int), "detected features")
        warnings.append(
            "Jaccard distances use presence/absence only and are highly "
            "sensitive to low-count false positives."
        )
        distance_df = _jaccard_distance(binary)
        ordination_df = _pcoa(distance_df)
        ordination_kind = "PCoA"

    return BetaDiversityResult(
        distance_df=distance_df,
        ordination_df=ordination_df,
        ordination_kind=ordination_kind,
        warnings=warnings,
    )


def distance_matrix_table(distance_df: pd.DataFrame) -> pd.DataFrame:
    """Format a square distance DataFrame for TSV output."""
    result = distance_df.copy()
    result.insert(0, "#Sample", result.index)
    return result.reset_index(drop=True)


def render_beta_heatmap(
    distance_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str] = None,
    width: float = 6.0,
    height: float = 5.0,
    dpi: int = 300,
) -> None:
    """Render a beta-diversity distance matrix as a heatmap."""
    if output_kind(output_file) == "html":
        _render_heatmap_html(distance_df, output_file, title)
    else:
        _render_heatmap_static(distance_df, output_file, title, width, height, dpi)


def render_beta_ordination(
    ordination_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str] = None,
    width: float = 6.0,
    height: float = 5.0,
    dpi: int = 300,
    kind: str = "PCA",
) -> None:
    """Render a 2D PCA/PCoA sample ordination scatter plot."""
    _require_ordination_columns(ordination_df)
    if output_kind(output_file) == "html":
        _render_ordination_html(ordination_df, output_file, title, kind)
    else:
        _render_ordination_static(
            ordination_df,
            output_file,
            title,
            width,
            height,
            dpi,
            kind,
        )


def _prepare_abundance_matrix(
    abundance_df: pd.DataFrame,
    min_feature_count: float,
    min_samples: int,
) -> pd.DataFrame:
    if min_feature_count < 0:
        raise ValueError("--min-feature-count must be >= 0")
    if min_samples < 1:
        raise ValueError("--min-samples must be >= 1")
    if "#Taxon" not in abundance_df.columns:
        raise ValueError("Abundance table must contain a #Taxon column")

    sample_cols = [col for col in abundance_df.columns if col != "#Taxon"]
    if len(sample_cols) < 2:
        raise ValueError("Beta diversity requires at least two samples")
    if abundance_df.empty:
        raise ValueError("No taxa found for the requested beta diversity options")

    values = abundance_df[sample_cols].apply(pd.to_numeric, errors="coerce")
    invalid_cols = [col for col in sample_cols if values[col].isna().any()]
    if invalid_cols:
        cols = ", ".join(invalid_cols)
        raise ValueError(f"Sample columns must be numeric: {cols}")
    if (values < 0).any().any():
        raise ValueError("Beta diversity requires non-negative counts or abundances")

    keep = pd.Series(True, index=values.index)
    if min_feature_count > 0:
        keep &= values.sum(axis=1) >= min_feature_count
    keep &= (values > 0).sum(axis=1) >= min_samples
    values = values.loc[keep]

    if values.empty:
        raise ValueError("No taxa remain after beta diversity filtering")

    matrix = values.T.astype(float)
    _require_non_empty_samples(matrix, "non-zero counts")
    return matrix


def _require_non_empty_samples(matrix: pd.DataFrame, description: str) -> None:
    empty_samples = matrix.sum(axis=1)
    empty_samples = empty_samples[empty_samples <= 0]
    if not empty_samples.empty:
        names = ", ".join(empty_samples.index.astype(str))
        raise ValueError(f"No {description} found for sample(s): {names}")


def _braycurtis_distance(matrix: pd.DataFrame) -> pd.DataFrame:
    samples = matrix.index.tolist()
    values = matrix.to_numpy(dtype=float)
    distances = np.zeros((len(samples), len(samples)), dtype=float)

    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            numerator = np.abs(values[i] - values[j]).sum()
            denominator = (values[i] + values[j]).sum()
            distance = 0.0 if denominator == 0 else numerator / denominator
            distances[i, j] = distance
            distances[j, i] = distance

    return pd.DataFrame(distances, index=samples, columns=samples)


def _clr_transform(matrix: pd.DataFrame, pseudocount: float) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float) + pseudocount
    log_values = np.log(values)
    clr_values = log_values - log_values.mean(axis=1, keepdims=True)
    return pd.DataFrame(clr_values, index=matrix.index, columns=matrix.columns)


def _euclidean_distance(matrix: pd.DataFrame) -> pd.DataFrame:
    samples = matrix.index.tolist()
    values = matrix.to_numpy(dtype=float)
    diffs = values[:, np.newaxis, :] - values[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diffs * diffs, axis=2))
    return pd.DataFrame(distances, index=samples, columns=samples)


def _jaccard_distance(binary: pd.DataFrame) -> pd.DataFrame:
    samples = binary.index.tolist()
    values = binary.to_numpy(dtype=bool)
    distances = np.zeros((len(samples), len(samples)), dtype=float)

    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            intersection = np.logical_and(values[i], values[j]).sum()
            union = np.logical_or(values[i], values[j]).sum()
            distance = 0.0 if union == 0 else 1.0 - intersection / union
            distances[i, j] = distance
            distances[j, i] = distance

    return pd.DataFrame(distances, index=samples, columns=samples)


def _pca(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    centered = values - values.mean(axis=0, keepdims=True)
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    components = centered @ vt.T
    variances = singular_values**2
    ratios = _variance_ratios(variances)
    return _ordination_dataframe(matrix.index.tolist(), components, ratios, "PC")


def _pcoa(distance_df: pd.DataFrame) -> pd.DataFrame:
    distances = distance_df.to_numpy(dtype=float)
    n_samples = distances.shape[0]
    centering = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    gram = -0.5 * centering @ (distances**2) @ centering
    eigvals, eigvecs = np.linalg.eigh(gram)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = np.where(eigvals > 0, eigvals, 0.0)
    coordinates = eigvecs * np.sqrt(positive)
    ratios = _variance_ratios(positive)
    return _ordination_dataframe(distance_df.index.tolist(), coordinates, ratios, "PCo")


def _variance_ratios(variances: np.ndarray) -> list[float]:
    total = variances[variances > 0].sum()
    if total <= 0:
        return [0.0, 0.0]
    ratios = (variances / total).tolist()
    return ratios + [0.0] * max(0, 2 - len(ratios))


def _ordination_dataframe(
    samples: list[str],
    coordinates: np.ndarray,
    ratios: list[float],
    prefix: str,
) -> pd.DataFrame:
    axis1 = coordinates[:, 0] if coordinates.shape[1] >= 1 else np.zeros(len(samples))
    axis2 = coordinates[:, 1] if coordinates.shape[1] >= 2 else np.zeros(len(samples))
    ratio1 = ratios[0] if ratios else 0.0
    ratio2 = ratios[1] if len(ratios) > 1 else 0.0
    return pd.DataFrame(
        {
            "#Sample": samples,
            "Axis1": axis1,
            "Axis2": axis2,
            "Axis1Explained": ratio1,
            "Axis2Explained": ratio2,
            "Axis1Label": f"{prefix}1 ({ratio1 * 100:.1f}%)",
            "Axis2Label": f"{prefix}2 ({ratio2 * 100:.1f}%)",
        }
    )


def _require_ordination_columns(ordination_df: pd.DataFrame) -> None:
    required = {"#Sample", "Axis1", "Axis2", "Axis1Label", "Axis2Label"}
    missing = required - set(ordination_df.columns)
    if missing:
        missing_names = ", ".join(sorted(missing))
        raise ValueError(f"Ordination table missing column(s): {missing_names}")


def _render_heatmap_static(
    distance_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
    width: float,
    height: float,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    samples = distance_df.index.tolist()
    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(distance_df.to_numpy(dtype=float), cmap="viridis", vmin=0)
    ax.set_xticks(range(len(samples)), labels=samples, rotation=45, ha="right")
    ax.set_yticks(range(len(samples)), labels=samples)
    ax.set_title(title or "Beta diversity distance")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Sample")
    fig.colorbar(im, ax=ax, label="Distance")
    fig.tight_layout()
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _render_heatmap_html(
    distance_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
) -> None:
    import plotly.graph_objects as go

    samples = distance_df.index.tolist()
    fig = go.Figure(
        data=go.Heatmap(
            z=distance_df.to_numpy(dtype=float),
            x=samples,
            y=samples,
            colorscale="Viridis",
            colorbar={"title": "Distance"},
            hovertemplate=(
                "Sample 1: %{y}<br>"
                "Sample 2: %{x}<br>"
                "Distance: %{z:.4f}<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        title=title or "Beta diversity distance",
        xaxis_title="Sample",
        yaxis_title="Sample",
    )
    fig.write_html(output_file)


def _render_ordination_static(
    ordination_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
    width: float,
    height: float,
    dpi: int,
    kind: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(width, height))
    colors = [PALETTE[idx % len(PALETTE)] for idx in range(len(ordination_df))]
    ax.scatter(ordination_df["Axis1"], ordination_df["Axis2"], s=70, color=colors)
    for _, row in ordination_df.iterrows():
        ax.annotate(
            row["#Sample"],
            (row["Axis1"], row["Axis2"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=9,
        )
    ax.axhline(0, color="#D0D0D0", linewidth=0.8)
    ax.axvline(0, color="#D0D0D0", linewidth=0.8)
    ax.set_xlabel(ordination_df["Axis1Label"].iloc[0])
    ax.set_ylabel(ordination_df["Axis2Label"].iloc[0])
    ax.set_title(title or f"Beta diversity {kind}")
    ax.grid(color="#EEEEEE", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _render_ordination_html(
    ordination_df: pd.DataFrame,
    output_file: Path,
    title: Optional[str],
    kind: str,
) -> None:
    import plotly.express as px

    fig = px.scatter(
        ordination_df,
        x="Axis1",
        y="Axis2",
        text="#Sample",
        hover_name="#Sample",
        color="#Sample",
        color_discrete_sequence=PALETTE,
        labels={
            "Axis1": ordination_df["Axis1Label"].iloc[0],
            "Axis2": ordination_df["Axis2Label"].iloc[0],
        },
        title=title or f"Beta diversity {kind}",
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(showlegend=False)
    fig.write_html(output_file)
