from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from kraut.plotting import PALETTE, output_kind


CLUSTERING_METHODS = {"ward", "average", "single", "complete"}
CLUSTERING_ALIASES = {
    "ward": "ward",
    "ward-linkage": "ward",
    "average": "average",
    "average-linkage": "average",
    "single": "single",
    "single-linkage": "single",
    "complete": "complete",
    "complete-linkage": "complete",
}
MISSING_METADATA_LABEL = "NA"


@dataclass(frozen=True)
class ClusterMerge:
    left: int
    right: int
    distance: float
    count: int


@dataclass(frozen=True)
class ClusterTree:
    labels: list[str]
    merges: list[ClusterMerge]
    root: int


@dataclass(frozen=True)
class DendrogramLayout:
    leaf_order: list[str]
    segments: list[tuple[list[float], list[float]]]
    max_distance: float


@dataclass(frozen=True)
class MetadataColors:
    column: str
    sample_to_color: dict[str, str]
    sample_to_value: dict[str, str]
    value_to_color: dict[str, str]


def cluster_distance_matrix(
    distance_df: pd.DataFrame,
    method: str = "ward",
) -> ClusterTree:
    """Cluster samples from a square distance matrix."""
    method = _normalize_clustering_method(method)

    labels, values = _validated_distance_matrix(distance_df)
    n_samples = len(labels)
    active = set(range(n_samples))
    sizes = {idx: 1 for idx in active}
    min_leaf = {idx: idx for idx in active}
    distances = {
        _pair_key(i, j): float(values[i, j])
        for i in range(n_samples)
        for j in range(i + 1, n_samples)
    }
    merges: list[ClusterMerge] = []
    next_cluster_id = n_samples

    while len(active) > 1:
        left, right = _closest_pair(distances, min_leaf)
        if min_leaf[left] > min_leaf[right]:
            left, right = right, left

        merge_distance = distances[_pair_key(left, right)]
        merged_id = next_cluster_id
        next_cluster_id += 1
        merged_size = sizes[left] + sizes[right]
        merges.append(
            ClusterMerge(
                left=left,
                right=right,
                distance=float(merge_distance),
                count=merged_size,
            )
        )

        remaining = sorted(active - {left, right})
        for cluster_id in remaining:
            distances[_pair_key(merged_id, cluster_id)] = _updated_distance(
                method,
                left,
                right,
                cluster_id,
                distances,
                sizes,
                merge_distance,
            )

        distances = {
            pair: distance
            for pair, distance in distances.items()
            if left not in pair and right not in pair
        }
        active.remove(left)
        active.remove(right)
        active.add(merged_id)
        sizes[merged_id] = merged_size
        min_leaf[merged_id] = min(min_leaf[left], min_leaf[right])

    return ClusterTree(labels=labels, merges=merges, root=next(iter(active)))


def build_dendrogram_layout(tree: ClusterTree) -> DendrogramLayout:
    """Convert a clustered tree into line segments for plotting."""
    n_leaves = len(tree.labels)
    merge_by_id = {
        n_leaves + merge_idx: merge for merge_idx, merge in enumerate(tree.merges)
    }
    leaf_ids = _ordered_leaf_ids(tree.root, merge_by_id, n_leaves)
    leaf_order = [tree.labels[leaf_id] for leaf_id in leaf_ids]
    leaf_positions = {leaf_id: float(idx) for idx, leaf_id in enumerate(leaf_ids)}
    segments: list[tuple[list[float], list[float]]] = []

    def visit(cluster_id: int) -> tuple[float, float]:
        if cluster_id < n_leaves:
            return leaf_positions[cluster_id], 0.0

        merge = merge_by_id[cluster_id]
        left_x, left_y = visit(merge.left)
        right_x, right_y = visit(merge.right)
        height = float(merge.distance)
        segments.append(
            (
                [left_x, left_x, right_x, right_x],
                [left_y, height, height, right_y],
            )
        )
        return (left_x + right_x) / 2.0, height

    _, max_distance = visit(tree.root)
    return DendrogramLayout(
        leaf_order=leaf_order,
        segments=segments,
        max_distance=float(max_distance),
    )


def metadata_color_annotations(
    samples: list[str],
    metadata_file: Optional[Path],
    color_by: Optional[str],
) -> Optional[MetadataColors]:
    """Return per-sample colors from a metadata column."""
    if metadata_file is None and color_by is None:
        return None
    if metadata_file is None:
        raise ValueError("--metadata is required when --color-by is supplied")
    if color_by is None:
        raise ValueError("--color-by is required when --metadata is supplied")
    if not metadata_file.exists():
        raise ValueError(f"Metadata file does not exist: {metadata_file}")

    df, sample_column = _read_metadata_file(metadata_file)
    if color_by not in df.columns:
        raise ValueError(f"--color-by column not found: {color_by}")

    df = df.set_index(sample_column, drop=False)
    missing = [sample for sample in samples if sample not in df.index]
    if missing:
        missing_names = ", ".join(missing)
        raise ValueError(f"Metadata is missing sample(s): {missing_names}")

    sample_to_value = {
        sample: _metadata_value_label(df.loc[sample, color_by]) for sample in samples
    }
    values = list(dict.fromkeys(sample_to_value[sample] for sample in samples))
    value_to_color = _value_colors(values)
    sample_to_color = {
        sample: value_to_color[sample_to_value[sample]] for sample in samples
    }
    return MetadataColors(
        column=color_by,
        sample_to_color=sample_to_color,
        sample_to_value=sample_to_value,
        value_to_color=value_to_color,
    )


def render_dendrogram(
    distance_df: pd.DataFrame,
    output_file: Path,
    clustering: str = "ward",
    metadata_colors: Optional[MetadataColors] = None,
    title: Optional[str] = None,
    width: float = 8.0,
    height: float = 5.0,
    dpi: int = 300,
) -> None:
    """Render a sample dendrogram from a distance matrix."""
    clustering = _normalize_clustering_method(clustering)
    tree = cluster_distance_matrix(distance_df, clustering)
    layout = build_dendrogram_layout(tree)
    if output_kind(output_file) == "html":
        _render_dendrogram_html(layout, output_file, clustering, metadata_colors, title)
    else:
        _render_dendrogram_static(
            layout,
            output_file,
            clustering,
            metadata_colors,
            title,
            width,
            height,
            dpi,
        )


def _validated_distance_matrix(distance_df: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    if distance_df.empty:
        raise ValueError("Distance matrix is empty")
    if distance_df.shape[0] != distance_df.shape[1]:
        raise ValueError("Distance matrix must be square")
    if len(distance_df.index) < 2:
        raise ValueError("Dendrogram requires at least two samples")

    labels = [str(label) for label in distance_df.index]
    columns = [str(column) for column in distance_df.columns]
    if labels != columns:
        raise ValueError("Distance matrix rows and columns must have the same order")

    values = distance_df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Distance matrix must contain only finite numeric values")
    if (values < -1e-12).any():
        raise ValueError("Distance matrix cannot contain negative distances")
    values = np.where(values < 0, 0.0, values)

    if not np.allclose(values, values.T, atol=1e-9, rtol=1e-9):
        raise ValueError("Distance matrix must be symmetric")
    if not np.allclose(np.diag(values), 0.0, atol=1e-9, rtol=1e-9):
        raise ValueError("Distance matrix diagonal must be zero")

    return labels, values


def _normalize_clustering_method(method: str) -> str:
    key = method.lower().replace("_", "-").replace(" ", "-")
    if key in CLUSTERING_ALIASES:
        return CLUSTERING_ALIASES[key]

    supported = ", ".join(sorted(CLUSTERING_METHODS))
    raise ValueError(f"--clustering must be one of: {supported}")


def _closest_pair(
    distances: dict[tuple[int, int], float],
    min_leaf: dict[int, int],
) -> tuple[int, int]:
    return min(
        distances,
        key=lambda pair: (
            distances[pair],
            min(min_leaf[pair[0]], min_leaf[pair[1]]),
            max(min_leaf[pair[0]], min_leaf[pair[1]]),
            pair[0],
            pair[1],
        ),
    )


def _updated_distance(
    method: str,
    left: int,
    right: int,
    other: int,
    distances: dict[tuple[int, int], float],
    sizes: dict[int, int],
    merge_distance: float,
) -> float:
    left_distance = distances[_pair_key(left, other)]
    right_distance = distances[_pair_key(right, other)]

    if method == "single":
        return min(left_distance, right_distance)
    if method == "complete":
        return max(left_distance, right_distance)
    if method == "average":
        return (
            sizes[left] * left_distance + sizes[right] * right_distance
        ) / (sizes[left] + sizes[right])

    total_size = sizes[left] + sizes[right] + sizes[other]
    ward_squared = (
        ((sizes[other] + sizes[left]) / total_size) * left_distance**2
        + ((sizes[other] + sizes[right]) / total_size) * right_distance**2
        - (sizes[other] / total_size) * merge_distance**2
    )
    return float(np.sqrt(max(ward_squared, 0.0)))


def _pair_key(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _ordered_leaf_ids(
    cluster_id: int,
    merge_by_id: dict[int, ClusterMerge],
    n_leaves: int,
) -> list[int]:
    if cluster_id < n_leaves:
        return [cluster_id]

    merge = merge_by_id[cluster_id]
    left = _ordered_leaf_ids(merge.left, merge_by_id, n_leaves)
    right = _ordered_leaf_ids(merge.right, merge_by_id, n_leaves)
    if min(left) > min(right):
        left, right = right, left
    return left + right


def _metadata_value_label(value) -> str:
    if pd.isna(value):
        return MISSING_METADATA_LABEL
    if isinstance(value, pd.Timestamp):
        return value.date().isoformat()
    text = str(value).strip()
    return text if text else MISSING_METADATA_LABEL


def _read_metadata_file(path: Path) -> tuple[pd.DataFrame, str]:
    delimiter = _detect_delimiter(path)
    try:
        data = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False)
    except Exception as exc:
        raise ValueError(f"Could not read metadata file '{path}': {exc}") from exc

    if len(data.columns) == 0:
        raise ValueError("Metadata must contain a header row")

    data = data.copy()
    data.columns = [str(column) for column in data.columns]
    sample_column = data.columns[0]
    sample_ids = data[sample_column].astype(str).str.strip()
    if sample_ids.eq("").any():
        raise ValueError("Metadata sample column contains empty values")
    if sample_ids.duplicated().any():
        duplicates = sorted(sample_ids[sample_ids.duplicated()].unique())
        duplicate_names = ", ".join(duplicates)
        raise ValueError(
            f"Metadata sample column contains duplicate values: {duplicate_names}"
        )
    data[sample_column] = sample_ids
    return data, sample_column


def _detect_delimiter(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return "\t"
    if suffix == ".csv":
        return ","

    sample = path.read_text()[:4096]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        return "\t" if "\t" in first_line else ","


def _value_colors(values: list[str]) -> dict[str, str]:
    colors = {}
    palette_idx = 0
    for value in values:
        if value == MISSING_METADATA_LABEL:
            colors[value] = "#6C757D"
        else:
            colors[value] = PALETTE[palette_idx % len(PALETTE)]
            palette_idx += 1
    return colors


def _render_dendrogram_static(
    layout: DendrogramLayout,
    output_file: Path,
    clustering: str,
    metadata_colors: Optional[MetadataColors],
    title: Optional[str],
    width: float,
    height: float,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(width, height))
    for xs, ys in layout.segments:
        ax.plot(xs, ys, color="#343A40", linewidth=1.6)

    x_positions = list(range(len(layout.leaf_order)))
    leaf_colors = _leaf_colors(layout.leaf_order, metadata_colors)
    ax.scatter(x_positions, [0.0] * len(x_positions), color=leaf_colors, zorder=3)
    ax.set_xticks(x_positions, labels=layout.leaf_order, rotation=45, ha="right")
    for tick_label in ax.get_xticklabels():
        sample = tick_label.get_text()
        tick_label.set_color(_sample_color(sample, metadata_colors))

    ax.set_ylabel("Distance")
    ax.set_title(title or f"Dendrogram ({clustering} linkage)")
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_ylim(bottom=0, top=max(layout.max_distance * 1.08, 1.0))

    if metadata_colors:
        handles = [
            Patch(facecolor=color, edgecolor="none", label=value)
            for value, color in metadata_colors.value_to_color.items()
        ]
        ax.legend(
            handles=handles,
            title=metadata_colors.column,
            loc="best",
            frameon=False,
        )

    fig.tight_layout()
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _render_dendrogram_html(
    layout: DendrogramLayout,
    output_file: Path,
    clustering: str,
    metadata_colors: Optional[MetadataColors],
    title: Optional[str],
) -> None:
    import plotly.graph_objects as go

    line_x: list[float | None] = []
    line_y: list[float | None] = []
    for xs, ys in layout.segments:
        line_x.extend(xs + [None])
        line_y.extend(ys + [None])

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=line_x,
            y=line_y,
            mode="lines",
            line={"color": "#343A40", "width": 2},
            hoverinfo="skip",
            showlegend=False,
        )
    )

    if metadata_colors:
        for value, color in metadata_colors.value_to_color.items():
            samples = [
                sample
                for sample in layout.leaf_order
                if metadata_colors.sample_to_value[sample] == value
            ]
            positions = [layout.leaf_order.index(sample) for sample in samples]
            fig.add_trace(
                go.Scatter(
                    x=positions,
                    y=[0.0] * len(samples),
                    text=samples,
                    mode="markers+text",
                    name=value,
                    marker={"color": color, "size": 9},
                    textposition="bottom center",
                    hovertemplate=(
                        "Sample: %{text}<br>"
                        f"{metadata_colors.column}: {value}<extra></extra>"
                    ),
                )
            )
    else:
        fig.add_trace(
            go.Scatter(
                x=list(range(len(layout.leaf_order))),
                y=[0.0] * len(layout.leaf_order),
                text=layout.leaf_order,
                mode="markers+text",
                marker={"color": PALETTE[0], "size": 9},
                textposition="bottom center",
                hovertemplate="Sample: %{text}<extra></extra>",
                showlegend=False,
            )
        )

    fig.update_layout(
        title=title or f"Dendrogram ({clustering} linkage)",
        xaxis={
            "tickmode": "array",
            "tickvals": list(range(len(layout.leaf_order))),
            "ticktext": layout.leaf_order,
            "title": "Sample",
        },
        yaxis={"title": "Distance", "rangemode": "tozero"},
        margin={"l": 70, "r": 30, "t": 70, "b": 110},
    )
    fig.write_html(output_file)


def _leaf_colors(
    leaf_order: list[str],
    metadata_colors: Optional[MetadataColors],
) -> list[str]:
    return [_sample_color(sample, metadata_colors) for sample in leaf_order]


def _sample_color(sample: str, metadata_colors: Optional[MetadataColors]) -> str:
    if metadata_colors is None:
        return PALETTE[0]
    return metadata_colors.sample_to_color[sample]
