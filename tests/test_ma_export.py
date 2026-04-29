from collections import Counter
from pathlib import Path

import pandas as pd
import pytest
import typer

from kraut.commands import ma_export
from kraut.models.ma_export import export_microbiome_analyst, generate_pseudo_labels


def write_ma_report(
    path: Path,
    multiplier: int = 1,
    genus_taxid: int = 561,
    genus_name: str = "Escherichia",
) -> None:
    taxa = [
        ("D", 2, "Bacteria"),
        ("K", 1000, "Pseudomonadati"),
        ("P", 1224, "Proteobacteria"),
        ("C", 1236, "Gammaproteobacteria"),
        ("O", 91347, "Enterobacterales"),
        ("F", 543, "Enterobacteriaceae"),
        ("G", genus_taxid, genus_name),
        ("S", 562, "Escherichia coli"),
    ]
    lines = [
        "  0.00\t0\t0\tU\t0\tunclassified\n",
        f"100.00\t{multiplier * 1000}\t0\tR\t1\troot\n",
    ]
    for idx, (rank, taxid, name) in enumerate(taxa, start=1):
        indent = "  " * idx
        lines.append(
            f" 99.00\t{multiplier * (100 + idx)}\t{multiplier * idx}"
            f"\t{rank}\t{taxid}\t{indent}{name}\n"
        )
    lines.append(
        f"  1.00\t{multiplier * 999}\t{multiplier * 999}"
        "\tS1\t83333\t                  Escherichia coli K-12\n"
    )
    path.write_text("".join(lines))


def test_ma_export_writes_microbiomeanalyst_files(tmp_path):
    alpha = tmp_path / "alpha.krep.tsv"
    beta = tmp_path / "beta.tsv"
    metadata = tmp_path / "metadata.tsv"
    outdir = tmp_path / "ma"
    write_ma_report(alpha, multiplier=1)
    write_ma_report(beta, multiplier=2)
    metadata.write_text("Sample\tCol1\nbeta\tB\nextra\tX\nalpha\tA\n")

    ma_export.run(
        input_files=[alpha, beta],
        outdir=outdir,
        metadata=metadata,
        metadata_sample_col="Sample",
        pseudo_col="Group",
        metric="LVL",
    )

    counts = pd.read_csv(outdir / "counts.csv")
    assert list(counts.columns) == ["#NAME", "alpha", "beta"]
    assert counts.to_dict("records")[0] == {"#NAME": "Feat_1", "alpha": 1, "beta": 2}
    assert counts.to_dict("records")[-1] == {
        "#NAME": "Feat_8",
        "alpha": 8,
        "beta": 16,
    }
    assert "Feat_9" not in counts["#NAME"].values

    taxonomy = pd.read_csv(outdir / "taxonomy.csv")
    kingdom = taxonomy[taxonomy["#TAXONOMY"] == "Feat_2"].iloc[0]
    assert kingdom["Domain"] == "d__Pseudomonadati"
    assert kingdom["Phylum"] == "p__"

    species = taxonomy[taxonomy["#TAXONOMY"] == "Feat_8"].iloc[0]
    assert species.to_dict() == {
        "#TAXONOMY": "Feat_8",
        "Domain": "d__Bacteria",
        "Phylum": "p__Proteobacteria",
        "Class": "c__Gammaproteobacteria",
        "Order": "o__Enterobacterales",
        "Family": "f__Enterobacteriaceae",
        "Genus": "g__Escherichia",
        "Species": "s__Escherichia coli",
    }

    reformatted_metadata = pd.read_csv(outdir / "metadata.csv")
    assert list(reformatted_metadata.columns) == ["#NAME", "Col1", "Group"]
    assert reformatted_metadata.to_dict("records") == [
        {"#NAME": "alpha", "Col1": "A", "Group": "A"},
        {"#NAME": "beta", "Col1": "B", "Group": "A"},
    ]

    tree = (outdir / "tree.nwk").read_text()
    for idx in range(1, 9):
        assert f"Feat_{idx}:1" in tree
    assert "Feat_9" not in tree
    assert tree.endswith(";\n")


def test_ma_export_can_use_tot_counts_and_generated_metadata(tmp_path):
    alpha = tmp_path / "alpha.tsv"
    beta = tmp_path / "beta.tsv"
    outdir = tmp_path / "ma"
    write_ma_report(alpha, multiplier=1)
    write_ma_report(beta, multiplier=2)

    ma_export.run(
        input_files=[alpha, beta],
        outdir=outdir,
        metadata=None,
        metadata_sample_col=None,
        pseudo_col="Random_label",
        metric="TOT",
    )

    counts = pd.read_csv(outdir / "counts.csv")
    assert counts[counts["#NAME"] == "Feat_8"].to_dict("records") == [
        {"#NAME": "Feat_8", "alpha": 108, "beta": 216}
    ]

    metadata = pd.read_csv(outdir / "metadata.csv")
    assert metadata.to_dict("records") == [
        {"#NAME": "alpha", "Random_label": "A"},
        {"#NAME": "beta", "Random_label": "A"},
    ]


def test_ma_export_metadata_sample_column_must_exist(tmp_path):
    alpha = tmp_path / "alpha.tsv"
    beta = tmp_path / "beta.tsv"
    metadata = tmp_path / "metadata.csv"
    outdir = tmp_path / "ma"
    write_ma_report(alpha, multiplier=1)
    write_ma_report(beta, multiplier=2)
    metadata.write_text("Wrong,Col1\nalpha,A\nbeta,B\n")

    with pytest.raises(typer.Exit) as excinfo:
        ma_export.run(
            input_files=[alpha, beta],
            outdir=outdir,
            metadata=metadata,
            metadata_sample_col="Sample",
            pseudo_col="Random_label",
            metric="LVL",
        )

    assert excinfo.value.exit_code == 1
    assert not outdir.exists()


def test_generate_pseudo_labels_avoids_singletons_when_possible():
    labels, warnings = generate_pseudo_labels(1)
    assert labels == ["A"]
    assert warnings

    for sample_count, expected in [
        (2, ["A", "A"]),
        (3, ["A", "A", "A"]),
        (5, ["A", "B", "A", "B", "A"]),
    ]:
        labels, warnings = generate_pseudo_labels(sample_count)
        assert labels == expected
        assert not warnings
        assert all(count > 1 for count in Counter(labels).values())


def test_ma_export_detects_conflicting_parentage(tmp_path):
    alpha = tmp_path / "alpha.tsv"
    beta = tmp_path / "beta.tsv"
    write_ma_report(alpha, multiplier=1, genus_taxid=561, genus_name="Escherichia")
    write_ma_report(beta, multiplier=2, genus_taxid=620, genus_name="Shigella")

    with pytest.raises(ValueError, match="conflicting taxonomy"):
        export_microbiome_analyst(
            input_files=[alpha, beta],
            outdir=tmp_path / "ma",
        )
