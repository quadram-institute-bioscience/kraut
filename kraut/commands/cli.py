import typer

from kraut import __version__

from . import alpha_cmd
from . import beta
from . import dendrogram
from . import list_reads
from . import ma_export
from . import make_mpa_table_cmd
from . import make_table_cmd
from . import plot_multi
from . import plot_single
from . import ranks_cmd
from . import single_report
from . import split_table_cmd
from . import table_summary_cmd

HELP_CONTEXT = {"help_option_names": ["-h", "--help"]}

REPORTS_PANEL = "Report Processing"
DIVERSITY_PANEL = "Diversity Analysis"
VISUALIZATION_PANEL = "Visualization"
EXTRAS_PANEL = "Extras"

app = typer.Typer(
    name="kraut",
    help="A tool for parsing and manipulating Kraken2 reports",
    add_completion=False,
    context_settings=HELP_CONTEXT,
    invoke_without_command=True,
)


def version_callback(value: bool) -> None:
    if value:
        typer.echo(f"kraut {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    ctx: typer.Context,
    version: bool = typer.Option(
        False,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Show the version and exit.",
    ),
) -> None:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()


app.command(
    name="single-report",
    context_settings=HELP_CONTEXT,
    rich_help_panel=REPORTS_PANEL,
)(single_report.run)
app.command(
    name="make-table",
    context_settings=HELP_CONTEXT,
    rich_help_panel=REPORTS_PANEL,
)(make_table_cmd.run)
app.command(
    name="split-combine-table",
    context_settings=HELP_CONTEXT,
    rich_help_panel=REPORTS_PANEL,
)(split_table_cmd.run)
app.command(
    name="table-summary",
    context_settings=HELP_CONTEXT,
    rich_help_panel=REPORTS_PANEL,
)(table_summary_cmd.run)

app.command(
    name="alpha",
    context_settings=HELP_CONTEXT,
    rich_help_panel=DIVERSITY_PANEL,
)(alpha_cmd.run)
app.command(
    name="beta",
    context_settings=HELP_CONTEXT,
    rich_help_panel=DIVERSITY_PANEL,
)(beta.run)
app.command(
    name="ranks",
    context_settings=HELP_CONTEXT,
    rich_help_panel=DIVERSITY_PANEL,
)(ranks_cmd.run)

app.command(
    name="plot-single",
    context_settings=HELP_CONTEXT,
    rich_help_panel=VISUALIZATION_PANEL,
)(plot_single.run)
app.command(
    name="plot-multi",
    context_settings=HELP_CONTEXT,
    rich_help_panel=VISUALIZATION_PANEL,
)(plot_multi.run)
app.command(
    name="dendrogram",
    context_settings=HELP_CONTEXT,
    rich_help_panel=VISUALIZATION_PANEL,
)(dendrogram.run)

app.command(
    name="ma-export",
    context_settings=HELP_CONTEXT,
    rich_help_panel=EXTRAS_PANEL,
)(ma_export.run)
app.command(
    name="list-reads",
    context_settings=HELP_CONTEXT,
    rich_help_panel=EXTRAS_PANEL,
)(list_reads.run)
app.command(
    name="make-mpa-table",
    context_settings=HELP_CONTEXT,
    rich_help_panel=EXTRAS_PANEL,
)(make_mpa_table_cmd.run)

if __name__ == "__main__":
    app()
