from kraut.models.kraken_data import KrakenNode, KrakenReport
from kraut.models.ranks import canonical_rank, lowest_rank, rank_read_counts


def make_report(*nodes):
    report = KrakenReport()
    report.nodes = {node.tax_id: node for node in nodes}
    report.unclassified = report.nodes.get(0)
    report.root = report.nodes.get(1)
    return report


def test_canonical_rank_normalizes_suffixes():
    assert canonical_rank("S") == "S"
    assert canonical_rank("S1") == "S"
    assert canonical_rank("s2") == "S"
    assert canonical_rank("G1") == "G"
    assert canonical_rank("R1") == "R"
    assert canonical_rank("U") == "U"
    assert canonical_rank("-") is None
    assert canonical_rank("") is None


def test_rank_read_counts_can_count_cumulative_or_exact_assignments():
    unclassified = KrakenNode(5.0, 5, 5, "U", 0, "unclassified")
    root = KrakenNode(95.0, 95, 0, "R", 1, "root")
    bacteria = KrakenNode(90.0, 90, 3, "D", 2, "Bacteria", depth=1)
    pseudomonadota = KrakenNode(70.0, 70, 7, "P", 1224, "Pseudomonadota", depth=2)
    escherichia = KrakenNode(50.0, 50, 5, "G", 561, "Escherichia", depth=3)
    ecoli = KrakenNode(40.0, 40, 31, "S", 562, "Escherichia coli", depth=4)
    ealbertii = KrakenNode(
        10.0,
        10,
        10,
        "S",
        61645,
        "Escherichia albertii",
        depth=4,
    )
    archaea = KrakenNode(5.0, 5, 5, "D", 2157, "Archaea", depth=1)

    root.add_child(bacteria)
    root.add_child(archaea)
    bacteria.add_child(pseudomonadota)
    pseudomonadota.add_child(escherichia)
    escherichia.add_child(ecoli)
    escherichia.add_child(ealbertii)
    report = make_report(
        unclassified,
        root,
        bacteria,
        pseudomonadota,
        escherichia,
        ecoli,
        ealbertii,
        archaea,
    )

    assert rank_read_counts(report, mode="cumulative") == {
        "D": 95,
        "K": 0,
        "P": 70,
        "C": 0,
        "O": 0,
        "F": 0,
        "G": 50,
        "S": 50,
    }
    assert rank_read_counts(report, mode="exact") == {
        "D": 8,
        "K": 0,
        "P": 7,
        "C": 0,
        "O": 0,
        "F": 0,
        "G": 5,
        "S": 41,
    }
    assert rank_read_counts(
        report,
        include_unclassified=True,
        include_root=True,
    ) == {
        "U": 5,
        "R": 95,
        "D": 95,
        "K": 0,
        "P": 70,
        "C": 0,
        "O": 0,
        "F": 0,
        "G": 50,
        "S": 50,
    }


def test_suffix_ranks_are_collapsed_without_cumulative_double_counting():
    root = KrakenNode(100.0, 100, 0, "R", 1, "root")
    bacteria = KrakenNode(100.0, 100, 0, "D", 2, "Bacteria", depth=1)
    escherichia = KrakenNode(80.0, 80, 10, "G", 561, "Escherichia", depth=2)
    unclassified_genus = KrakenNode(
        30.0,
        30,
        6,
        "G1",
        2608889,
        "unclassified Escherichia",
        depth=3,
    )
    ecoli = KrakenNode(20.0, 20, 5, "S", 562, "Escherichia coli", depth=4)
    ecoli_strain = KrakenNode(
        15.0,
        15,
        15,
        "S1",
        83333,
        "Escherichia coli K-12",
        depth=5,
    )

    root.add_child(bacteria)
    bacteria.add_child(escherichia)
    escherichia.add_child(unclassified_genus)
    unclassified_genus.add_child(ecoli)
    ecoli.add_child(ecoli_strain)
    report = make_report(
        root,
        bacteria,
        escherichia,
        unclassified_genus,
        ecoli,
        ecoli_strain,
    )

    cumulative = rank_read_counts(report, mode="cumulative")
    exact = rank_read_counts(report, mode="exact")

    assert cumulative["G"] == 80
    assert cumulative["S"] == 20
    assert exact["G"] == 16
    assert exact["S"] == 20


def test_cumulative_species_counts_include_species_suffix_without_species_parent():
    root = KrakenNode(100.0, 100, 0, "R", 1, "root")
    bacteria = KrakenNode(100.0, 100, 0, "D", 2, "Bacteria", depth=1)
    escherichia = KrakenNode(100.0, 100, 0, "G", 561, "Escherichia", depth=2)
    ecoli_strain = KrakenNode(
        15.0,
        15,
        15,
        "S1",
        83333,
        "Escherichia coli K-12",
        depth=3,
    )

    root.add_child(bacteria)
    bacteria.add_child(escherichia)
    escherichia.add_child(ecoli_strain)
    report = make_report(root, bacteria, escherichia, ecoli_strain)

    assert rank_read_counts(report, mode="cumulative")["S"] == 15


def test_lowest_rank_returns_deepest_nonzero_classified_rank():
    root = KrakenNode(100.0, 100, 0, "R", 1, "root")
    bacteria = KrakenNode(100.0, 100, 0, "D", 2, "Bacteria", depth=1)
    escherichia = KrakenNode(100.0, 100, 4, "G", 561, "Escherichia", depth=2)
    ecoli = KrakenNode(20.0, 20, 20, "S", 562, "Escherichia coli", depth=3)

    root.add_child(bacteria)
    bacteria.add_child(escherichia)
    escherichia.add_child(ecoli)

    assert lowest_rank(make_report(root, bacteria, escherichia, ecoli)) == "S"


def test_lowest_rank_handles_genus_only_unclassified_only_and_root_only_reports():
    genus_root = KrakenNode(100.0, 100, 0, "R", 1, "root")
    genus = KrakenNode(100.0, 100, 100, "G", 561, "Escherichia", depth=1)
    genus_root.add_child(genus)

    unclassified = KrakenNode(100.0, 100, 100, "U", 0, "unclassified")
    root_only = KrakenNode(100.0, 100, 0, "R", 1, "root")

    assert lowest_rank(make_report(genus_root, genus)) == "G"
    assert lowest_rank(make_report(unclassified)) is None
    assert lowest_rank(make_report(root_only)) is None


def test_rank_read_counts_rejects_unknown_modes():
    report = KrakenReport()

    try:
        rank_read_counts(report, mode="all-the-reads")
    except ValueError as exc:
        assert str(exc) == "mode must be one of: cumulative, exact"
    else:
        raise AssertionError("rank_read_counts should reject unknown modes")
