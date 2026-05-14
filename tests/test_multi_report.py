from io import StringIO

import pandas as pd
import pytest

from kraut.models.kraken_data import KrakenNode, KrakenReport
from kraut.models.multi_report import MultiKrakenReport


def build_report(
    unclassified_counts,
    ecoli_counts,
    salmonella_counts=None,
):
    report = KrakenReport()

    unclassified = KrakenNode(0.0, *unclassified_counts, "U", 0, "unclassified")
    root = KrakenNode(100.0, 100, 0, "R", 1, "root")
    bacteria = KrakenNode(90.0, 90, 0, "D", 2, "Bacteria", depth=1)
    escherichia = KrakenNode(
        80.0,
        ecoli_counts[0],
        0,
        "G",
        561,
        "Escherichia",
        depth=2,
    )
    ecoli = KrakenNode(70.0, *ecoli_counts, "S", 562, "Escherichia coli", depth=3)

    root.add_child(bacteria)
    bacteria.add_child(escherichia)
    escherichia.add_child(ecoli)

    nodes = {
        0: unclassified,
        1: root,
        2: bacteria,
        561: escherichia,
        562: ecoli,
    }

    if salmonella_counts is not None:
        salmonella = KrakenNode(
            12.0,
            *salmonella_counts,
            "S",
            28901,
            "Salmonella enterica",
            depth=2,
        )
        bacteria.add_child(salmonella)
        nodes[28901] = salmonella

    report.unclassified = unclassified
    report.root = root
    report.nodes = nodes
    return report


def make_multi_report():
    multi_report = MultiKrakenReport()
    multi_report.add_report(build_report((10, 10), (70, 30)), "alpha")
    multi_report.add_report(
        build_report((20, 20), (8, 5), salmonella_counts=(12, 11)),
        "beta",
    )
    return multi_report


def build_lineage_report():
    report = KrakenReport()

    unclassified = KrakenNode(0.0, 0, 0, "U", 0, "unclassified")
    root = KrakenNode(100.0, 100, 0, "R", 1, "root")
    bacteria = KrakenNode(95.0, 95, 0, "D", 2, "Bacteria", depth=1)
    phylum = KrakenNode(90.0, 90, 0, "P", 201174, "Actinobacteria", depth=2)
    klass = KrakenNode(85.0, 85, 0, "C", 1760, "Actinobacteria", depth=3)
    order = KrakenNode(80.0, 80, 0, "O", 85004, "Bifidobacteriales", depth=4)
    family = KrakenNode(75.0, 75, 0, "F", 31953, "Bifidobacteriaceae", depth=5)
    genus = KrakenNode(70.0, 70, 10, "G", 1678, "Bifidobacterium", depth=6)
    species = KrakenNode(
        60.0,
        60,
        30,
        "S",
        216816,
        "Bifidobacterium longum",
        depth=7,
    )

    root.add_child(bacteria)
    bacteria.add_child(phylum)
    phylum.add_child(klass)
    klass.add_child(order)
    order.add_child(family)
    family.add_child(genus)
    genus.add_child(species)

    report.unclassified = unclassified
    report.root = root
    report.nodes = {
        node.tax_id: node
        for node in [
            unclassified,
            root,
            bacteria,
            phylum,
            klass,
            order,
            family,
            genus,
            species,
        ]
    }
    return report


def test_tot_and_lvl_metrics_use_the_expected_kraken_counts():
    multi_report = make_multi_report()

    tot_df = multi_report.to_dataframe(metric="TOT", level="S")
    lvl_df = multi_report.to_dataframe(metric="LVL", level="S")

    ecoli_tot = tot_df[tot_df["#Taxon"] == "Escherichia coli"].iloc[0]
    ecoli_lvl = lvl_df[lvl_df["#Taxon"] == "Escherichia coli"].iloc[0]
    salmonella_tot = tot_df[tot_df["#Taxon"] == "Salmonella enterica"].iloc[0]

    assert ecoli_tot["alpha"] == 70
    assert ecoli_tot["beta"] == 8
    assert ecoli_lvl["alpha"] == 30
    assert ecoli_lvl["beta"] == 5
    assert salmonella_tot["alpha"] == 0
    assert salmonella_tot["beta"] == 12


def test_dataframe_can_use_rank_prefixes_or_taxids_as_taxon_keys():
    multi_report = make_multi_report()

    prefixed = multi_report.to_dataframe(metric="TOT", level="S", rank_prefix=True)
    taxids = multi_report.to_dataframe(metric="TOT", level="S", use_taxid=True)

    assert "s__Escherichia coli" in prefixed["#Taxon"].values
    assert "562" in taxids["#Taxon"].values


def test_tsv_includes_unclassified_first_unless_excluded():
    multi_report = make_multi_report()

    with_unclassified = multi_report.to_tsv(metric="TOT", level="S")
    without_unclassified = multi_report.to_tsv(
        metric="TOT",
        level="S",
        no_unclassified=True,
    )

    assert with_unclassified.splitlines()[1].startswith("unclassified\t10\t20")
    assert "unclassified" not in without_unclassified


def test_percentage_metric_normalizes_selected_rows_per_sample():
    multi_report = make_multi_report()

    df = pd.read_csv(
        StringIO(multi_report.to_tsv(metric="PERCENTAGE", level="S")),
        sep="\t",
    )

    unclassified = df[df["#Taxon"] == "unclassified"].iloc[0]
    ecoli = df[df["#Taxon"] == "Escherichia coli"].iloc[0]
    salmonella = df[df["#Taxon"] == "Salmonella enterica"].iloc[0]

    assert unclassified["alpha"] == 12.5
    assert ecoli["alpha"] == 87.5
    assert salmonella["alpha"] == 0.0
    assert unclassified["beta"] == 50.0
    assert ecoli["beta"] == 20.0
    assert salmonella["beta"] == 30.0


def test_min_percentage_filters_rows_that_never_reach_threshold():
    multi_report = make_multi_report()

    df = pd.read_csv(
        StringIO(multi_report.to_tsv(metric="TOT", level="S", min_perc=40.0)),
        sep="\t",
    )

    assert set(df["#Taxon"]) == {"unclassified", "Escherichia coli"}


def test_duplicate_taxon_labels_are_rejected():
    report = build_report((10, 10), (70, 30), salmonella_counts=(12, 11))
    report.nodes[28901].name = "Escherichia coli"

    multi_report = MultiKrakenReport()
    multi_report.add_report(report, "alpha")

    with pytest.raises(ValueError, match="Duplicate taxon label"):
        multi_report.to_dataframe(metric="TOT", level="S")


def test_duplicate_sample_names_are_rejected():
    multi_report = MultiKrakenReport()
    multi_report.add_report(build_report((10, 10), (70, 30)), "alpha")

    with pytest.raises(ValueError, match="Duplicate sample name"):
        multi_report.add_report(build_report((20, 20), (8, 5)), "alpha")


def test_conflicting_taxonomy_for_existing_taxid_is_rejected():
    multi_report = MultiKrakenReport()
    multi_report.add_report(build_report((10, 10), (70, 30)), "alpha")

    conflicting_report = build_report((20, 20), (8, 5))
    conflicting_report.nodes[562].name = "Shigella flexneri"

    with pytest.raises(ValueError, match="TaxID 562 has conflicting taxonomy"):
        multi_report.add_report(conflicting_report, "beta")


def test_dataframe_can_use_lineage_labels_up_to_requested_rank():
    multi_report = MultiKrakenReport()
    multi_report.add_report(build_lineage_report(), "alpha")

    df = multi_report.to_dataframe(metric="TOT", level="G", add_lineage=True)

    assert df.to_dict("records") == [
        {
            "#Taxon": (
                "k__Bacteria,p__Actinobacteria,c__Actinobacteria,"
                "o__Bifidobacteriales,f__Bifidobacteriaceae,g__Bifidobacterium"
            ),
            "alpha": 70,
        }
    ]


def test_lineage_labels_are_incompatible_with_taxid_labels():
    multi_report = make_multi_report()

    with pytest.raises(ValueError, match="--taxid cannot be used with --add-lineage"):
        multi_report.to_dataframe(
            metric="TOT",
            level="S",
            use_taxid=True,
            add_lineage=True,
        )
