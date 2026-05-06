import gzip
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, List, Optional, TextIO

import typer

from kraut.models.kraken_data import KrakenNode, KrakenReport


@contextmanager
def _open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            yield handle
    else:
        with path.open() as handle:
            yield handle


def run(
    raw_output: Path = typer.Argument(..., help="Kraken raw output TSV"),
    taxon: Optional[List[int]] = typer.Option(
        None,
        "--taxon",
        "-t",
        help="TaxID to match. Can be repeated.",
    ),
    report_file: Optional[Path] = typer.Option(
        None,
        "--report",
        "-r",
        help="Kraken report TSV. Required with --children.",
    ),
    children: bool = typer.Option(
        False,
        "--children",
        help="Include each selected taxon and its descendants from the report.",
    ),
    unclassified: bool = typer.Option(
        False,
        "--unclassified",
        help="Include unclassified reads.",
    ),
    invert: bool = typer.Option(
        False,
        "--invert",
        "-v",
        help="Output reads that do not match the selection.",
    ),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file (default: stdout).",
    ),
):
    """
    List read names matching taxon selections from Kraken raw output.
    """
    if not raw_output.exists():
        typer.echo(f"Error: Raw output file {raw_output} does not exist.", err=True)
        raise typer.Exit(code=1)

    selected_taxa = set(taxon or [])
    if not selected_taxa and not unclassified:
        typer.echo("Error: Provide at least one --taxon or --unclassified.", err=True)
        raise typer.Exit(code=1)

    if children:
        if report_file is None:
            typer.echo("Error: --children requires --report.", err=True)
            raise typer.Exit(code=1)
        if not report_file.exists():
            typer.echo(f"Error: Report file {report_file} does not exist.", err=True)
            raise typer.Exit(code=1)
        selected_taxa = _expand_with_children(selected_taxa, report_file)

    if output_file:
        with output_file.open("w") as output_handle:
            _write_matching_reads(
                raw_output,
                selected_taxa,
                unclassified,
                invert,
                output_handle,
            )
    else:
        _write_matching_reads(
            raw_output,
            selected_taxa,
            unclassified,
            invert,
            None,
        )


def _expand_with_children(taxa: set[int], report_file: Path) -> set[int]:
    report = KrakenReport.from_file(str(report_file))
    expanded = set(taxa)
    for tax_id in taxa:
        node = report.nodes.get(tax_id)
        if node:
            _collect_descendant_taxa(node, expanded)
    return expanded


def _collect_descendant_taxa(node: KrakenNode, taxa: set[int]) -> None:
    taxa.add(node.tax_id)
    for child in node.children:
        _collect_descendant_taxa(child, taxa)


def _write_matching_reads(
    raw_output: Path,
    selected_taxa: set[int],
    include_unclassified: bool,
    invert: bool,
    output_handle: Optional[TextIO],
) -> None:
    with _open_text(raw_output) as input_handle:
        for line in input_handle:
            read = _parse_raw_output_line(line)
            if read is None:
                continue
            status, read_name, tax_id = read
            matched = tax_id in selected_taxa or (
                include_unclassified and (status == "U" or tax_id == 0)
            )
            if matched == invert:
                continue

            if output_handle:
                output_handle.write(f"{read_name}\n")
            else:
                typer.echo(read_name)


def _parse_raw_output_line(line: str) -> Optional[tuple[str, str, int]]:
    line = line.rstrip("\n")
    if not line:
        return None

    fields = line.split("\t")
    if len(fields) < 3:
        return None

    try:
        tax_id = int(fields[2])
    except ValueError:
        return None

    return fields[0], fields[1], tax_id
