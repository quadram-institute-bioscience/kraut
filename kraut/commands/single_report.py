import typer
from typing import Optional
from pathlib import Path
from kraut.models.kraken_data import KrakenReport

# app = typer.Typer(help="Parse and print a single Kraken report")
# Removed app to allow direct registration in cli.py

def run(
    input_file: Path = typer.Option(..., "--input", "-i", help="Input Kraken report file (KREP)"),
    output_file: Optional[Path] = typer.Option(None, "--output", "-o", help="Output file (default: stdout)"),
    min_fract: float = typer.Option(0.0, "--min-fract", "-m", help="Minimum fraction of reads to keep a taxon"),
    min_count: int = typer.Option(0, "--min-count", "-c", help="Minimum count of reads to keep a taxon"),
    min_level: Optional[str] = typer.Option(None, "--min-level", "-l", help="Minimum level to keep (K, P, C, O, F, G, S)"),
    max_level: Optional[str] = typer.Option(None, "--max-level", "-L", help="Maximum level to keep (e.g. do not print below Species)"),
):
    """
    Parses a Kraken report and prints it in text format.
    """
    if not input_file.exists():
        typer.echo(f"Error: Input file {input_file} does not exist.", err=True)
        raise typer.Exit(code=1)

    report = KrakenReport.from_file(str(input_file))
    
    # Generate output string with strictly formatted columns
    result = report.to_string(
        min_fract=min_fract,
        min_count=min_count,
        min_level=min_level,
        max_level=max_level
    )
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(result)
    else:
        typer.echo(result, nl=False)
