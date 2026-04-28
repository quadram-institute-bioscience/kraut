import pytest

from kraut.models.kraken_data import KrakenNode, KrakenReport
from kraut.models.multi_report import MultiKrakenReport
from kraut.plotting import composition_dataframe


def build_species_report(unclassified_count, species_counts):
    report = KrakenReport()

    classified_total = sum(species_counts.values())
    total = unclassified_count + classified_total
    unclassified = KrakenNode(
        unclassified_count / total * 100,
        unclassified_count,
        unclassified_count,
        "U",
        0,
        "unclassified",
    )
    root = KrakenNode(100.0, classified_total, 0, "R", 1, "root")
    bacteria = KrakenNode(
        classified_total / total * 100,
        classified_total,
        0,
        "D",
        2,
        "Bacteria",
        depth=1,
    )
    root.add_child(bacteria)

    nodes = {0: unclassified, 1: root, 2: bacteria}
    tax_ids = {
        "Escherichia coli": 562,
        "Salmonella enterica": 28901,
        "Lactobacillus crispatus": 47770,
    }
    for name, count in species_counts.items():
        node = KrakenNode(
            count / total * 100,
            count,
            count,
            "S",
            tax_ids[name],
            name,
            depth=2,
        )
        bacteria.add_child(node)
        nodes[node.tax_id] = node

    report.unclassified = unclassified
    report.root = root
    report.nodes = nodes
    return report


@pytest.fixture
def plotting_multi_report():
    multi_report = MultiKrakenReport()
    multi_report.add_report(
        build_species_report(
            10,
            {
                "Escherichia coli": 65,
                "Salmonella enterica": 20,
                "Lactobacillus crispatus": 5,
            },
        ),
        "alpha",
    )
    multi_report.add_report(
        build_species_report(
            10,
            {
                "Escherichia coli": 5,
                "Salmonella enterica": 70,
                "Lactobacillus crispatus": 15,
            },
        ),
        "beta",
    )
    return multi_report


def test_composition_dataframe_normalizes_each_sample(plotting_multi_report):
    df = composition_dataframe(plotting_multi_report, min_perc=0.0)

    assert df["alpha"].sum() == pytest.approx(100.0)
    assert df["beta"].sum() == pytest.approx(100.0)
    assert list(df["#Taxon"]) == [
        "unclassified",
        "Salmonella enterica",
        "Escherichia coli",
        "Lactobacillus crispatus",
    ]


def test_composition_dataframe_folds_low_abundance_taxa(plotting_multi_report):
    df = composition_dataframe(plotting_multi_report, min_perc=16.0)

    assert list(df["#Taxon"]) == [
        "unclassified",
        "Salmonella enterica",
        "Escherichia coli",
        "Others",
    ]
    others = df[df["#Taxon"] == "Others"].iloc[0]
    assert others["alpha"] == pytest.approx(5.0)
    assert others["beta"] == pytest.approx(15.0)


def test_composition_dataframe_top_taxa_uses_total_abundance(plotting_multi_report):
    df = composition_dataframe(plotting_multi_report, min_perc=0.0, top_taxa=1)

    assert list(df["#Taxon"]) == ["unclassified", "Salmonella enterica", "Others"]
    salmonella = df[df["#Taxon"] == "Salmonella enterica"].iloc[0]
    others = df[df["#Taxon"] == "Others"].iloc[0]
    assert salmonella["alpha"] == pytest.approx(20.0)
    assert salmonella["beta"] == pytest.approx(70.0)
    assert others["alpha"] == pytest.approx(70.0)
    assert others["beta"] == pytest.approx(20.0)


def test_composition_dataframe_no_unclassified_removes_and_renormalizes(
    plotting_multi_report,
):
    df = composition_dataframe(
        plotting_multi_report,
        min_perc=0.0,
        no_unclassified=True,
    )

    assert "unclassified" not in df["#Taxon"].values
    assert df["alpha"].sum() == pytest.approx(100.0)
    ecoli = df[df["#Taxon"] == "Escherichia coli"].iloc[0]
    assert ecoli["alpha"] == pytest.approx(65 / 90 * 100)


def test_composition_dataframe_rejects_invalid_options(plotting_multi_report):
    with pytest.raises(ValueError, match="--metric"):
        composition_dataframe(plotting_multi_report, metric="PERCENTAGE")

    with pytest.raises(ValueError, match="--min-perc"):
        composition_dataframe(plotting_multi_report, min_perc=-1)

    with pytest.raises(ValueError, match="--top-taxa"):
        composition_dataframe(plotting_multi_report, top_taxa=-1)
