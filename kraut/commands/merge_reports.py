import typer
from typing import List, Optional
from pathlib import Path
from kraut.models.kraken_data import KrakenReport
from kraut.models.multi_report import MultiKrakenReport

# app = typer.Typer(help="Merge multiple Kraken reports")
# Removed app to allow direct registration in cli.py

def run(
    input_files: List[Path] = typer.Argument(..., help="Input Kraken report files (KREP)"),
    output_file: Optional[Path] = typer.Option(None, "--output", "-o", help="Output file (default: stdout)"),
    metric: str = typer.Option("COUNTS", "--metric", "-m", help="Metric to report: COUNTS or PERCENTAGES"),
    level: str = typer.Option("S", "--level", "-l", help="Taxonomic level to report (default: S)"),
    no_unclassified: bool = typer.Option(False, "--no-unclassified", help="Exclude unclassified reads"),
):
    """
    Merges multiple Kraken reports into a single table.
    """
    multi_report = MultiKrakenReport()
    
    for p in input_files:
        if not p.exists():
             typer.echo(f"Warning: File {p} does not exist, skipping.", err=True)
             continue
             
        # Infer sample name from filename (basename without extension)
        # Assuming extension is .tsv based on prompt examples, but let's be safe
        sample_name = p.stem
        if sample_name.endswith(".krep"): # Handle .krep extension if present
             sample_name = Path(sample_name).stem
             
        report = KrakenReport.from_file(str(p))
        multi_report.add_report(report, sample_name)
        
    result = multi_report.to_tsv(metric=metric, level=level, no_unclassified=no_unclassified)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(result)
    else:
        typer.echo(result)
