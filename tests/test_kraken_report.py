from kraut.models.kraken_data import KrakenReport


def test_parse_preserves_unclassified_and_tree_relationships(synthetic_report_path):
    report = KrakenReport.from_file(str(synthetic_report_path))

    assert report.unclassified is not None
    assert report.unclassified.tax_id == 0
    assert report.unclassified.name == "unclassified"

    assert report.root is not None
    assert report.root.name == "root"
    assert [child.name for child in report.root.children] == ["Bacteria", "Archaea"]

    bacteria = report.nodes[2]
    pseudomonadota = report.nodes[1224]
    escherichia = report.nodes[561]
    ecoli = report.nodes[562]
    ealbertii = report.nodes[61645]

    assert bacteria.parent is report.root
    assert pseudomonadota.parent is bacteria
    assert escherichia.parent is pseudomonadota
    assert [child.name for child in escherichia.children] == [
        "Escherichia coli",
        "Escherichia albertii",
    ]
    assert ecoli.depth == 4
    assert ealbertii.depth == 4


def test_round_trip_preserves_kraken_indentation(synthetic_report_path):
    report = KrakenReport.from_file(str(synthetic_report_path))

    assert report.to_string() == synthetic_report_path.read_text()


def test_min_count_filter_prunes_low_count_nodes(synthetic_report_path):
    report = KrakenReport.from_file(str(synthetic_report_path))

    filtered = report.to_string(min_count=20)

    assert "unclassified" not in filtered
    assert "root" in filtered
    assert "Bacteria" in filtered
    assert "Pseudomonadota" in filtered
    assert "Escherichia coli" in filtered
    assert "Escherichia albertii" not in filtered
    assert "Archaea" not in filtered


def test_min_fraction_filter_uses_kraken_percent_column(synthetic_report_path):
    report = KrakenReport.from_file(str(synthetic_report_path))

    filtered = report.to_string(min_fract=50.0)

    assert "root" in filtered
    assert "Bacteria" in filtered
    assert "Pseudomonadota" in filtered
    assert "Escherichia" not in filtered
