import pandas as pd
import pytest

from kraut.dendrogram import (
    build_dendrogram_layout,
    cluster_distance_matrix,
    metadata_color_annotations,
)


def distance_table():
    return pd.DataFrame(
        [
            [0.0, 0.1, 0.8],
            [0.1, 0.0, 0.7],
            [0.8, 0.7, 0.0],
        ],
        index=["alpha", "beta", "gamma"],
        columns=["alpha", "beta", "gamma"],
    )


def test_average_linkage_clusters_nearest_samples_first():
    tree = cluster_distance_matrix(distance_table(), method="average")

    assert len(tree.merges) == 2
    assert tree.merges[0].left == 0
    assert tree.merges[0].right == 1
    assert tree.merges[0].distance == pytest.approx(0.1)
    assert tree.merges[1].distance == pytest.approx(0.75)


def test_dendrogram_layout_preserves_leaf_labels():
    tree = cluster_distance_matrix(distance_table(), method="ward")
    layout = build_dendrogram_layout(tree)

    assert layout.leaf_order == ["alpha", "beta", "gamma"]
    assert len(layout.segments) == 2
    assert layout.max_distance > 0


def test_cluster_distance_matrix_rejects_invalid_clustering_method():
    with pytest.raises(ValueError, match="--clustering"):
        cluster_distance_matrix(distance_table(), method="centroid")


def test_cluster_distance_matrix_accepts_linkage_aliases():
    tree = cluster_distance_matrix(distance_table(), method="complete linkage")

    assert tree.merges[0].distance == pytest.approx(0.1)


def test_metadata_color_annotations_uses_requested_column(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample\tgroup\nalpha\tA\nbeta\tB\ngamma\tA\n")

    colors = metadata_color_annotations(
        ["alpha", "beta", "gamma"],
        metadata,
        "group",
    )

    assert colors.column == "group"
    assert colors.sample_to_value == {"alpha": "A", "beta": "B", "gamma": "A"}
    assert colors.sample_to_color["alpha"] == colors.sample_to_color["gamma"]
    assert colors.sample_to_color["alpha"] != colors.sample_to_color["beta"]
