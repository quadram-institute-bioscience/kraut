import math

import pandas as pd
import pytest

from kraut.models.beta import (
    calculate_beta_diversity,
    distance_matrix_table,
    read_abundance_table,
)


def abundance_table():
    return pd.DataFrame(
        {
            "#Taxon": ["t1", "t2"],
            "alpha": [10, 0],
            "beta": [0, 10],
            "gamma": [10, 10],
        }
    )


def test_braycurtis_distance_matrix():
    result = calculate_beta_diversity(abundance_table(), metric="braycurtis")

    assert result.distance_df.loc["alpha", "alpha"] == pytest.approx(0.0)
    assert result.distance_df.loc["alpha", "beta"] == pytest.approx(1.0)
    assert result.distance_df.loc["alpha", "gamma"] == pytest.approx(1 / 3)
    assert result.ordination_kind == "PCoA"

    output = distance_matrix_table(result.distance_df)
    assert list(output.columns) == ["#Sample", "alpha", "beta", "gamma"]


def test_jaccard_distance_matrix_warns_about_detection_noise():
    result = calculate_beta_diversity(abundance_table(), metric="jaccard")

    assert result.distance_df.loc["alpha", "beta"] == pytest.approx(1.0)
    assert result.distance_df.loc["alpha", "gamma"] == pytest.approx(0.5)
    assert result.ordination_kind == "PCoA"
    assert result.warnings
    assert "false positives" in result.warnings[0]


def test_jaccard_presence_threshold_can_remove_low_counts():
    df = pd.DataFrame(
        {
            "#Taxon": ["t1", "t2"],
            "alpha": [5, 0],
            "beta": [1, 8],
        }
    )

    result = calculate_beta_diversity(
        df,
        metric="jaccard",
        presence_threshold=4,
    )

    assert result.distance_df.loc["alpha", "beta"] == pytest.approx(1.0)


def test_aitchison_uses_clr_euclidean_distance():
    result = calculate_beta_diversity(
        abundance_table(),
        metric="aitchison",
        pseudocount=1,
    )

    expected_alpha_beta = math.sqrt(2) * math.log(11)
    expected_alpha_gamma = math.log(11) / math.sqrt(2)
    assert result.distance_df.loc["alpha", "beta"] == pytest.approx(expected_alpha_beta)
    assert result.distance_df.loc["alpha", "gamma"] == pytest.approx(
        expected_alpha_gamma
    )
    assert result.ordination_kind == "PCA"


def test_filtering_removes_rare_taxa_before_beta_diversity():
    df = pd.DataFrame(
        {
            "#Taxon": ["rare", "shared"],
            "alpha": [1, 10],
            "beta": [0, 20],
        }
    )

    result = calculate_beta_diversity(
        df,
        metric="braycurtis",
        min_feature_count=10,
        min_samples=2,
    )

    assert list(result.distance_df.index) == ["alpha", "beta"]
    assert result.distance_df.loc["alpha", "beta"] == pytest.approx(1 / 3)


def test_beta_diversity_rejects_negative_counts():
    df = pd.DataFrame({"#Taxon": ["t1"], "alpha": [1], "beta": [-1]})

    with pytest.raises(ValueError, match="non-negative"):
        calculate_beta_diversity(df)


def test_read_abundance_table_selects_numeric_sample_columns(tmp_path):
    table = tmp_path / "table.tsv"
    table.write_text(
        "Tax\talpha\tbeta\t#Taxon\n"
        "Tax1\t10\t0\tEscherichia coli\n"
        "Tax2\t0\t5\tSalmonella enterica\n"
    )

    df = read_abundance_table(table)

    assert list(df.columns) == ["#Taxon", "alpha", "beta"]
    assert list(df["#Taxon"]) == ["Escherichia coli", "Salmonella enterica"]


def test_read_abundance_table_handles_combined_bracken_outputs(tmp_path):
    table = tmp_path / "combined_bracken.tsv"
    table.write_text(
        "name\ttaxonomy_id\ttaxonomy_lvl\talpha.brout_num\talpha.brout_frac"
        "\tbeta.brout_num\tbeta.brout_frac\n"
        "Escherichia coli\t562\tS\t10\t0.1\t0\t0.0\n"
        "Escherichia\t561\tG\t10\t0.1\t5\t0.05\n"
        "Salmonella enterica\t28901\tS\t0\t0.0\t20\t0.2\n"
    )

    df = read_abundance_table(table, rank="S")

    assert list(df.columns) == ["#Taxon", "alpha", "beta"]
    assert list(df["#Taxon"]) == ["562", "28901"]
    assert list(df["alpha"]) == [10, 0]
    assert list(df["beta"]) == [0, 20]
