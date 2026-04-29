from kraut.alpha_diversity import (
    CORE_METRICS,
    abundance_dataframe,
    select_alpha_metrics,
)
from kraut.models.kraken_data import KrakenNode, KrakenReport
from kraut.models.multi_report import MultiKrakenReport


def build_report():
    report = KrakenReport()
    report.unclassified = KrakenNode(10.0, 10, 10, "U", 0, "unclassified")
    report.root = KrakenNode(90.0, 90, 0, "R", 1, "root")
    bacteria = KrakenNode(90.0, 90, 0, "D", 2, "Bacteria", depth=1)
    ecoli = KrakenNode(70.0, 70, 70, "S", 562, "Escherichia coli", depth=2)
    salmonella = KrakenNode(
        20.0,
        20,
        20,
        "S",
        28901,
        "Salmonella enterica",
        depth=2,
    )
    report.root.add_child(bacteria)
    bacteria.add_child(ecoli)
    bacteria.add_child(salmonella)
    return report


def test_select_alpha_metrics_expands_preset_and_additions():
    assert select_alpha_metrics("core", "chao1,ace,chao1") == CORE_METRICS + [
        "chao1",
        "ace",
    ]


def test_abundance_dataframe_excludes_unclassified_by_default():
    multi_report = MultiKrakenReport()
    multi_report.add_report(build_report(), "alpha")

    df = abundance_dataframe(multi_report, rank="S", metric="TOT")

    assert list(df["#Taxon"]) == ["562", "28901"]
    assert list(df["alpha"]) == [70, 20]


def test_abundance_dataframe_can_include_unclassified():
    multi_report = MultiKrakenReport()
    multi_report.add_report(build_report(), "alpha")

    df = abundance_dataframe(
        multi_report,
        rank="S",
        metric="TOT",
        include_unclassified=True,
    )

    assert list(df["#Taxon"]) == ["0", "562", "28901"]
    assert list(df["alpha"]) == [10, 70, 20]
