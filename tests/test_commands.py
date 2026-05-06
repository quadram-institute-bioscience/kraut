import gzip
from pathlib import Path

import pandas as pd
import pytest
import typer
from typer.testing import CliRunner

from kraut import __version__
import kraut.alpha_diversity as alpha_diversity
from kraut.commands import (
    alpha_cmd,
    beta as beta_cmd,
    dendrogram as dendrogram_cmd,
    list_reads,
    make_table_cmd,
    plot_multi,
    plot_single,
    ranks_cmd,
    split_table_cmd,
)
from kraut.commands.cli import app

runner = CliRunner()


def write_report(
    path: Path,
    species_clade_count: int,
    species_taxon_count: int,
) -> None:
    path.write_text(
        " 10.00\t10\t10\tU\t0\tunclassified\n"
        " 90.00\t90\t0\tR\t1\troot\n"
        f" 90.00\t90\t{90 - species_taxon_count}\tD\t2\t  Bacteria\n"
        f" 80.00\t{species_clade_count}\t0\tG\t561\t    Escherichia\n"
        f" 70.00\t{species_clade_count}\t{species_taxon_count}"
        "\tS\t562\t      Escherichia coli\n"
    )


def write_raw_output(path: Path) -> None:
    path.write_text(
        "C\tread_parent\t561\t100\t561:100\n"
        "C\tread_child_one\t562\t100\t562:100\n"
        "C\tread_child_two\t61645\t100\t61645:100\n"
        "C\tread_domain\t2\t100\t2:100\n"
        "U\tread_unclassified\t0\t100\t0:100\n"
    )


def run_list_reads(
    raw_output: Path,
    taxon: list[int] | None = None,
    report_file: Path | None = None,
    children: bool = False,
    unclassified: bool = False,
    invert: bool = False,
    output_file: Path | None = None,
) -> None:
    list_reads.run(
        raw_output=raw_output,
        taxon=taxon,
        report_file=report_file,
        children=children,
        unclassified=unclassified,
        invert=invert,
        output_file=output_file,
    )


def test_cli_help_groups_commands_by_category():
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "--version" in result.output
    for panel_name in [
        "Report Processing",
        "Diversity Analysis",
        "Visualization",
        "Export and Reads",
    ]:
        assert panel_name in result.output


def test_cli_version_option_prints_package_version():
    result = runner.invoke(app, ["--version"])

    assert result.exit_code == 0
    assert result.output.strip() == f"kraut {__version__}"


def test_list_reads_writes_exact_taxon_matches(tmp_path):
    raw_output = tmp_path / "kraken.out"
    output = tmp_path / "reads.txt"
    write_raw_output(raw_output)

    run_list_reads(raw_output, taxon=[562], output_file=output)

    assert output.read_text().splitlines() == ["read_child_one"]


def test_list_reads_accepts_gzipped_raw_output(tmp_path):
    raw_output = tmp_path / "kraken.out.gz"
    output = tmp_path / "reads.txt"
    with gzip.open(raw_output, "wt") as handle:
        handle.write("C\tread_child_one\t562\t100\t562:100\n")

    run_list_reads(raw_output, taxon=[562], output_file=output)

    assert output.read_text().splitlines() == ["read_child_one"]


def test_list_reads_can_list_only_unclassified_reads(tmp_path, capsys):
    raw_output = tmp_path / "kraken.out"
    write_raw_output(raw_output)

    run_list_reads(raw_output, unclassified=True)

    assert capsys.readouterr().out.splitlines() == ["read_unclassified"]


def test_list_reads_children_include_parent_and_descendants(
    tmp_path,
    synthetic_report_path,
):
    raw_output = tmp_path / "kraken.out"
    output = tmp_path / "reads.txt"
    write_raw_output(raw_output)

    run_list_reads(
        raw_output,
        taxon=[561],
        report_file=synthetic_report_path,
        children=True,
        output_file=output,
    )

    assert output.read_text().splitlines() == [
        "read_parent",
        "read_child_one",
        "read_child_two",
    ]


def test_list_reads_invert_applies_after_combined_match_set(tmp_path):
    raw_output = tmp_path / "kraken.out"
    output = tmp_path / "reads.txt"
    write_raw_output(raw_output)

    run_list_reads(
        raw_output,
        taxon=[562],
        unclassified=True,
        invert=True,
        output_file=output,
    )

    assert output.read_text().splitlines() == [
        "read_parent",
        "read_child_two",
        "read_domain",
    ]


def test_list_reads_requires_a_taxon_or_unclassified(tmp_path):
    raw_output = tmp_path / "kraken.out"
    write_raw_output(raw_output)

    with pytest.raises(typer.Exit) as excinfo:
        run_list_reads(raw_output)

    assert excinfo.value.exit_code == 1


def test_make_table_defaults_cover_simple_merged_report_output(tmp_path):
    alpha = tmp_path / "alpha.krep.tsv"
    beta = tmp_path / "beta.tsv"
    output = tmp_path / "merged.tsv"
    write_report(alpha, species_clade_count=70, species_taxon_count=30)
    write_report(beta, species_clade_count=8, species_taxon_count=5)

    make_table_cmd.run(
        input_files=[alpha, beta],
        output_file=output,
        metric="TOT",
        level="S",
        rank_prefix=False,
        use_taxid=False,
        no_unclassified=False,
        min_perc=0.0,
    )

    df = pd.read_csv(output, sep="\t")

    assert list(df.columns) == ["#Taxon", "alpha", "beta"]
    assert df.to_dict("records") == [
        {"#Taxon": "unclassified", "alpha": 10, "beta": 10},
        {"#Taxon": "Escherichia coli", "alpha": 70, "beta": 8},
    ]


def test_make_table_command_writes_named_samples_from_input_stems(tmp_path):
    alpha = tmp_path / "alpha.krep.tsv"
    beta = tmp_path / "beta.tsv"
    output = tmp_path / "table.tsv"
    write_report(alpha, species_clade_count=70, species_taxon_count=30)
    write_report(beta, species_clade_count=8, species_taxon_count=5)

    make_table_cmd.run(
        input_files=[alpha, beta],
        output_file=output,
        metric="LVL",
        level="S",
        rank_prefix=True,
        use_taxid=False,
        no_unclassified=True,
        min_perc=0.0,
    )

    df = pd.read_csv(output, sep="\t")

    assert list(df.columns) == ["#Taxon", "alpha", "beta"]
    assert df.to_dict("records") == [
        {"#Taxon": "s__Escherichia coli", "alpha": 30, "beta": 5}
    ]


def write_combined_krakentools_table(path: Path) -> None:
    path.write_text(
        "#Number of Samples: 2\n"
        "#S1\treports/alpha.krep.tsv\n"
        "#S2\treports/beta.tsv\n"
        "#perc\ttot_all\ttot_lvl\tS1_all\tS1_lvl\tS2_all\tS2_lvl"
        "\tlvl_type\ttaxid\tname\n"
        "70.0\t100\t10\t70\t7\t30\t3\tS\t562\t  Escherichia coli\n"
        "30.0\t900\t0\t400\t0\t500\t0\tR\t1\troot\n"
    )


def test_split_combine_table_creates_all_and_lvl_outputs(tmp_path):
    input_table = tmp_path / "combined.tsv"
    output_base = tmp_path / "split"
    write_combined_krakentools_table(input_table)

    split_table_cmd.run(
        input_file=input_table,
        output_basename=str(output_base),
        use_taxid=False,
        rank_filter="S",
    )

    all_df = pd.read_csv(f"{output_base}_all.tsv", sep="\t")
    lvl_df = pd.read_csv(f"{output_base}_lvl.tsv", sep="\t")

    assert all_df.to_dict("records") == [
        {"#Taxon": "Escherichia coli", "alpha": 70, "beta": 30}
    ]
    assert lvl_df.to_dict("records") == [
        {"#Taxon": "Escherichia coli", "alpha": 7, "beta": 3}
    ]


def test_split_combine_table_can_use_taxids_as_row_labels(tmp_path):
    input_table = tmp_path / "combined.tsv"
    output_base = tmp_path / "taxid"
    write_combined_krakentools_table(input_table)

    split_table_cmd.run(
        input_file=input_table,
        output_basename=str(output_base),
        use_taxid=True,
        rank_filter=None,
    )

    all_df = pd.read_csv(f"{output_base}_all.tsv", sep="\t")

    assert list(all_df["#Taxon"]) == [562, 1]
    assert list(all_df["alpha"]) == [70, 400]
    assert list(all_df["beta"]) == [30, 500]


def write_plot_report(
    path: Path,
    unclassified_count: int,
    species_counts: dict,
) -> None:
    classified_total = sum(species_counts.values())
    total = unclassified_count + classified_total
    species_taxids = {
        "Escherichia coli": 562,
        "Salmonella enterica": 28901,
        "Lactobacillus crispatus": 47770,
    }
    lines = [
        (
            f"{unclassified_count / total * 100:6.2f}\t{unclassified_count}"
            f"\t{unclassified_count}\tU\t0\tunclassified\n"
        ),
        f"{100.0:6.2f}\t{classified_total}\t0\tR\t1\troot\n",
        (
            f"{classified_total / total * 100:6.2f}\t{classified_total}"
            "\t0\tD\t2\t  Bacteria\n"
        ),
    ]
    for name, count in species_counts.items():
        lines.append(
            f"{count / total * 100:6.2f}\t{count}\t{count}\tS"
            f"\t{species_taxids[name]}\t    {name}\n"
        )
    path.write_text("".join(lines))


def write_bracken_report(path: Path, species_counts: dict) -> None:
    total = sum(species_counts.values())
    lines = [
        f"{100.0:6.2f}\t{total}\t0\tR\t1\troot\n",
        f"{100.0:6.2f}\t{total}\t0\tD\t2\t  Bacteria\n",
    ]
    for idx, (name, count) in enumerate(species_counts.items(), start=1000):
        lines.append(
            f"{count / total * 100:6.2f}\t{count}\t{count}\tS" f"\t{idx}\t    {name}\n"
        )
    path.write_text("".join(lines))


def fake_alpha_diversity(metric, counts, ids=None, validate=True):
    data = counts.to_numpy()
    values = []
    for sample_counts in data:
        observed = int((sample_counts > 0).sum())
        total = int(sample_counts.sum())
        if metric == "observed_features":
            values.append(observed)
        elif metric == "chao1":
            values.append(observed + 1)
        else:
            values.append(total / max(observed, 1))
    return pd.Series(values, index=ids)


def fake_get_alpha_diversity_metrics():
    return [
        "observed_features",
        "shannon",
        "simpson",
        "inv_simpson",
        "pielou_e",
        "dominance",
        "goods_coverage",
        "chao1",
    ]


def test_alpha_command_writes_table_and_plot_for_kraken_and_bracken(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        alpha_diversity,
        "_load_skbio_diversity",
        lambda: (fake_alpha_diversity, fake_get_alpha_diversity_metrics),
    )
    kraken = tmp_path / "alpha.krep.tsv"
    bracken = tmp_path / "beta.brep"
    output_table = tmp_path / "alpha.tsv"
    output_plot = tmp_path / "alpha.svg"
    write_plot_report(
        kraken,
        unclassified_count=10,
        species_counts={
            "Escherichia coli": 70,
            "Salmonella enterica": 20,
        },
    )
    write_bracken_report(
        bracken,
        species_counts={
            "Escherichia coli": 15,
            "Salmonella enterica": 35,
        },
    )

    alpha_cmd.run(
        input_files=[kraken, bracken],
        output_table=output_table,
        plot_file=output_plot,
        rank="S",
        metric="TOT",
        metrics="core",
        add_metrics="chao1",
        width=4.0,
        height=1.5,
        dpi=72,
    )

    df = pd.read_csv(output_table, sep="\t")

    assert list(df.columns) == ["#Metric", "alpha", "beta"]
    assert list(df["#Metric"]) == [
        "observed_features",
        "shannon",
        "simpson",
        "inv_simpson",
        "pielou_e",
        "dominance",
        "goods_coverage",
        "chao1",
    ]
    assert df.loc[df["#Metric"] == "observed_features", ["alpha", "beta"]].to_dict(
        "records"
    ) == [{"alpha": 2.0, "beta": 2.0}]
    assert output_plot.exists()
    assert output_plot.stat().st_size > 0


def test_beta_command_writes_matrix_heatmap_and_pca_from_reports(tmp_path):
    kraken = tmp_path / "alpha.krep.tsv"
    bracken = tmp_path / "beta.brep"
    output_table = tmp_path / "beta.tsv"
    heatmap = tmp_path / "beta.svg"
    pca = tmp_path / "beta_pca.svg"
    write_plot_report(
        kraken,
        unclassified_count=10,
        species_counts={
            "Escherichia coli": 70,
            "Salmonella enterica": 20,
        },
    )
    write_bracken_report(
        bracken,
        species_counts={
            "Escherichia coli": 15,
            "Salmonella enterica": 35,
        },
    )

    beta_cmd.run(
        input_files=[kraken, bracken],
        output_file=output_table,
        plot_file=heatmap,
        pca_file=pca,
        metric="braycurtis",
        width=4.0,
        height=3.0,
        dpi=72,
    )

    df = pd.read_csv(output_table, sep="\t")

    assert list(df.columns) == ["#Sample", "alpha", "beta"]
    assert list(df["#Sample"]) == ["alpha", "beta"]
    assert df.loc[df["#Sample"] == "alpha", "alpha"].iloc[0] == pytest.approx(0.0)
    assert heatmap.exists()
    assert heatmap.stat().st_size > 0
    assert pca.exists()
    assert pca.stat().st_size > 0


def test_beta_command_reads_wide_table_and_writes_stdout(tmp_path, capsys):
    table = tmp_path / "table.tsv"
    table.write_text("#Taxon\talpha\tbeta\tgamma\n" "t1\t10\t0\t10\n" "t2\t0\t10\t10\n")

    beta_cmd.run(input_files=[table], metric="jaccard")
    captured = capsys.readouterr()

    assert "Warning: Jaccard" in captured.err
    assert captured.out.splitlines()[0] == "#Sample\talpha\tbeta\tgamma"
    assert "alpha\t0.0\t1.0\t0.5" in captured.out


def test_beta_command_rejects_unsupported_plot_suffix(tmp_path):
    table = tmp_path / "table.tsv"
    output = tmp_path / "beta.tsv"
    plot = tmp_path / "beta.txt"
    table.write_text("#Taxon\talpha\tbeta\n" "t1\t10\t0\n" "t2\t0\t10\n")

    with pytest.raises(typer.Exit) as exc:
        beta_cmd.run(input_files=[table], output_file=output, plot_file=plot)

    assert exc.value.exit_code == 1
    assert not output.exists()
    assert not plot.exists()


def test_dendrogram_command_writes_plot_with_metadata(tmp_path):
    kraken = tmp_path / "alpha.krep.tsv"
    bracken = tmp_path / "beta.brep"
    output_plot = tmp_path / "dendrogram.svg"
    metadata = tmp_path / "metadata.tsv"
    write_plot_report(
        kraken,
        unclassified_count=10,
        species_counts={
            "Escherichia coli": 70,
            "Salmonella enterica": 20,
        },
    )
    write_bracken_report(
        bracken,
        species_counts={
            "Escherichia coli": 15,
            "Salmonella enterica": 35,
        },
    )
    metadata.write_text("sample\tgroup\nalpha\tA\nbeta\tB\n")

    dendrogram_cmd.run(
        input_files=[kraken, bracken],
        output_file=output_plot,
        distance="braycurtis",
        clustering="average",
        metadata_file=metadata,
        color_by="group",
        width=4.0,
        height=3.0,
        dpi=72,
    )

    assert output_plot.exists()
    assert output_plot.stat().st_size > 0


def test_ranks_command_writes_percentages_that_sum_to_100(tmp_path):
    kraken = tmp_path / "alpha.krep.tsv"
    bracken = tmp_path / "beta.brep"
    output = tmp_path / "ranks.tsv"
    write_report(kraken, species_clade_count=70, species_taxon_count=30)
    write_bracken_report(
        bracken,
        species_counts={
            "Escherichia coli": 15,
            "Salmonella enterica": 35,
        },
    )

    ranks_cmd.run(input_files=[kraken, bracken], output_file=output)

    df = pd.read_csv(output, sep="\t")

    assert list(df["#Rank"]) == ["U", "R", "D", "K", "P", "C", "O", "F", "G", "S"]
    assert df["alpha"].sum() == pytest.approx(100.0)
    assert df["beta"].sum() == pytest.approx(100.0)
    assert df.loc[df["#Rank"] == "U", "alpha"].iloc[0] == pytest.approx(10.0)
    assert df.loc[df["#Rank"] == "D", "alpha"].iloc[0] == pytest.approx(60.0)
    assert df.loc[df["#Rank"] == "S", "alpha"].iloc[0] == pytest.approx(30.0)
    assert df.loc[df["#Rank"] == "U", "beta"].iloc[0] == pytest.approx(0.0)
    assert df.loc[df["#Rank"] == "S", "beta"].iloc[0] == pytest.approx(100.0)


def test_ranks_command_can_write_counts(tmp_path):
    kraken = tmp_path / "alpha.krep.tsv"
    output = tmp_path / "ranks.tsv"
    write_report(kraken, species_clade_count=70, species_taxon_count=30)

    ranks_cmd.run(input_files=[kraken], output_file=output, counts=True)

    df = pd.read_csv(output, sep="\t")

    assert list(df["#Rank"]) == ["U", "R", "D", "K", "P", "C", "O", "F", "G", "S"]
    assert df.loc[df["#Rank"] == "U", "alpha"].iloc[0] == 10
    assert df.loc[df["#Rank"] == "R", "alpha"].iloc[0] == 0
    assert df.loc[df["#Rank"] == "D", "alpha"].iloc[0] == 60
    assert df.loc[df["#Rank"] == "G", "alpha"].iloc[0] == 0
    assert df.loc[df["#Rank"] == "S", "alpha"].iloc[0] == 30


@pytest.mark.parametrize("suffix", [".html", ".png"])
def test_ranks_command_can_write_stacked_rank_plot(tmp_path, suffix):
    kraken = tmp_path / "alpha.krep.tsv"
    bracken = tmp_path / "beta.brep"
    output = tmp_path / "ranks.tsv"
    plot = tmp_path / f"ranks{suffix}"
    write_report(kraken, species_clade_count=70, species_taxon_count=30)
    write_bracken_report(
        bracken,
        species_counts={
            "Escherichia coli": 15,
            "Salmonella enterica": 35,
        },
    )

    ranks_cmd.run(input_files=[kraken, bracken], output_file=output, plot_file=plot)

    assert output.exists()
    assert plot.exists()
    assert plot.stat().st_size > 0
    if suffix == ".html":
        html = plot.read_text().lower()
        assert "<html" in html
        assert "rank composition" in html


def test_ranks_command_rejects_unsupported_plot_suffix(tmp_path):
    kraken = tmp_path / "alpha.krep.tsv"
    output = tmp_path / "ranks.tsv"
    plot = tmp_path / "ranks.txt"
    write_report(kraken, species_clade_count=70, species_taxon_count=30)

    with pytest.raises(typer.Exit) as exc:
        ranks_cmd.run(input_files=[kraken], output_file=output, plot_file=plot)

    assert exc.value.exit_code == 1
    assert not output.exists()
    assert not plot.exists()


@pytest.mark.parametrize("suffix", [".html", ".png", ".pdf", ".svg"])
def test_plot_single_writes_supported_outputs(tmp_path, suffix):
    input_file = tmp_path / "alpha.tsv"
    output_file = tmp_path / f"single{suffix}"
    write_plot_report(
        input_file,
        unclassified_count=10,
        species_counts={
            "Escherichia coli": 70,
            "Salmonella enterica": 20,
        },
    )

    plot_single.run(
        input_file=input_file,
        output_file=output_file,
        min_perc=0.0,
        width=4.0,
        height=3.0,
        dpi=72,
    )

    assert output_file.exists()
    assert output_file.stat().st_size > 0
    if suffix == ".html":
        assert "<html" in output_file.read_text().lower()


@pytest.mark.parametrize("suffix", [".html", ".png"])
def test_plot_multi_writes_supported_outputs(tmp_path, suffix):
    alpha = tmp_path / "alpha.krep.tsv"
    beta = tmp_path / "beta.tsv"
    output_file = tmp_path / f"multi{suffix}"
    write_plot_report(
        alpha,
        unclassified_count=10,
        species_counts={
            "Escherichia coli": 70,
            "Salmonella enterica": 20,
        },
    )
    write_plot_report(
        beta,
        unclassified_count=5,
        species_counts={
            "Escherichia coli": 10,
            "Salmonella enterica": 85,
        },
    )

    plot_multi.run(
        input_files=[alpha, beta],
        output_file=output_file,
        min_perc=0.0,
        top_taxa=1,
        width=4.0,
        height=3.0,
        dpi=72,
    )

    assert output_file.exists()
    assert output_file.stat().st_size > 0


@pytest.mark.parametrize("suffix", [".html", ".png"])
def test_plot_multi_writes_bubble_outputs(tmp_path, suffix):
    alpha = tmp_path / "alpha.tsv"
    beta = tmp_path / "beta.tsv"
    output_file = tmp_path / f"bubble{suffix}"
    write_plot_report(
        alpha,
        unclassified_count=10,
        species_counts={
            "Escherichia coli": 70,
            "Salmonella enterica": 20,
        },
    )
    write_plot_report(
        beta,
        unclassified_count=5,
        species_counts={
            "Escherichia coli": 10,
            "Salmonella enterica": 85,
        },
    )

    plot_multi.run(
        input_files=[alpha, beta],
        output_file=output_file,
        min_perc=0.0,
        kind="bubble",
        width=4.0,
        height=3.0,
        dpi=72,
    )

    assert output_file.exists()
    assert output_file.stat().st_size > 0
    if suffix == ".html":
        html = output_file.read_text().lower()
        assert "<html" in html
        assert "scatter" in html


def test_plot_multi_rejects_unsupported_kind(tmp_path):
    alpha = tmp_path / "alpha.tsv"
    output_file = tmp_path / "multi.png"
    write_plot_report(
        alpha,
        unclassified_count=10,
        species_counts={"Escherichia coli": 90},
    )

    with pytest.raises(typer.Exit) as exc:
        plot_multi.run(input_files=[alpha], output_file=output_file, kind="cloud")

    assert exc.value.exit_code == 1
    assert not output_file.exists()


def test_plot_single_rejects_unsupported_output_suffix(tmp_path):
    input_file = tmp_path / "alpha.tsv"
    output_file = tmp_path / "single.txt"
    write_plot_report(
        input_file,
        unclassified_count=10,
        species_counts={"Escherichia coli": 90},
    )

    with pytest.raises(typer.Exit) as exc:
        plot_single.run(input_file=input_file, output_file=output_file)

    assert exc.value.exit_code == 1
    assert not output_file.exists()
