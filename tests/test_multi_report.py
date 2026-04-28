from io import StringIO

import pandas as pd

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
    escherichia = KrakenNode(80.0, ecoli_counts[0], 0, "G", 561, "Escherichia", depth=2)
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
