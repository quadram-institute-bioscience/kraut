import typer
from typing import List, Optional
from pathlib import Path
from kraut.models.kraken_data import KrakenReport
from kraut.models.multi_report import MultiKrakenReport

# No app instantiation here, just the function
def run(
    input_files: List[Path] = typer.Argument(..., help="Input Kraken report files"),
    output_file: Optional[Path] = typer.Option(None, "--output", "-o", help="Output file (default: stdout)"),
    metric: str = typer.Option("TOT", "--metric", "-m", help="Metric: LVL (specific), TOT (cumulative), PERCENTAGE"),
    level: str = typer.Option("S", "--rank", "-r", help="Taxonomic rank (K,P,C,O,F,G,S,ALL)"),
    rank_prefix: bool = typer.Option(False, "--rank-prefix", "-p", help="Add rank prefix (e.g. s__Species)"),
    use_taxid: bool = typer.Option(False, "--taxid", help="Use TaxID instead of name"),
    no_unclassified: bool = typer.Option(False, "--no-unclassified", help="Exclude unclassified reads") # Implicitly requested as useful?
):
    """
    Generate a table from multiple reports with customizable formatting.
    """
    multi_report = MultiKrakenReport()
    
    for p in input_files:
        if not p.exists():
             typer.echo(f"Warning: File {p} does not exist, skipping.", err=True)
             continue
             
        sample_name = p.stem
        if sample_name.endswith(".krep"):
             sample_name = Path(sample_name).stem
             
        report = KrakenReport.from_file(str(p))
        multi_report.add_report(report, sample_name)
        
    result = multi_report.to_tsv(
        metric=metric, 
        level=level, 
        no_unclassified=no_unclassified,
        use_taxid=use_taxid,
        rank_prefix=rank_prefix
    )
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(result)
    else:
        typer.echo(result)
