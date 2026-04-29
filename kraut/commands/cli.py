import typer
from .import single_report
from .import merge_reports
from .import make_table_cmd
from .import plot_multi
from .import plot_single
from .import split_table_cmd
from .import alpha_cmd
from .import ranks_cmd

app = typer.Typer(
    name="kraut",
    help="A tool for parsing and manipulating Kraken2 reports",
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)

app.command(name="single-report", context_settings={"help_option_names": ["-h", "--help"]})(single_report.run)
app.command(name="merge-reports", context_settings={"help_option_names": ["-h", "--help"]})(merge_reports.run)
app.command(name="make-table", context_settings={"help_option_names": ["-h", "--help"]})(make_table_cmd.run)
app.command(name="alpha", context_settings={"help_option_names": ["-h", "--help"]})(alpha_cmd.run)
app.command(name="ranks", context_settings={"help_option_names": ["-h", "--help"]})(ranks_cmd.run)
app.command(name="plot-single", context_settings={"help_option_names": ["-h", "--help"]})(plot_single.run)
app.command(name="plot-multi", context_settings={"help_option_names": ["-h", "--help"]})(plot_multi.run)
app.command(name="split-combine-table", context_settings={"help_option_names": ["-h", "--help"]})(split_table_cmd.run)

if __name__ == "__main__":
    app()
