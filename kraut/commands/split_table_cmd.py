import typer
import pandas as pd
from typing import Optional
from pathlib import Path
import sys

def run(
    input_file: Path = typer.Option(..., "--input", "-i", help="Input table (KrakenTools combined format)"),
    output_basename: str = typer.Option(..., "--output", "-o", help="Output basename"),
    use_taxid: bool = typer.Option(False, "--taxid", help="Use TaxID as first column"),
    rank_filter: Optional[str] = typer.Option(None, "--rank", "-r", help="Filter by rank code (e.g. S)"),
):
    """
    Split a KrakenTools combined table into two tables (ALL and LVL) with clean headers.
    """
    if not input_file.exists():
        typer.echo(f"Error: Input file {input_file} does not exist.", err=True)
        raise typer.Exit(code=1)

    # 1. Parse Header
    sample_map = {}
    header_line_index = 0
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith("#S"):
                # Format: #S1\tpath/to/file.tsv
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    s_id = parts[0][1:] # remove #
                    filepath = parts[1]
                    # Use stem
                    sample_name = Path(filepath).stem
                    if sample_name.endswith(".krep"):
                        sample_name = Path(sample_name).stem
                    sample_map[s_id] = sample_name
            
            if line.startswith("#perc") or line.startswith("perc"):
                header_line_index = i
                break
    
    # 2. Read Table
    try:
        df = pd.read_csv(input_file, sep='\t', header=header_line_index)
    except Exception as e:
        typer.echo(f"Error reading table: {e}", err=True)
        raise typer.Exit(code=1)
        
    # Clean header: remove leading # if present in column names
    df.columns = [c.replace('#', '') for c in df.columns]
    
    # Identify relevant columns
    # We expect columns like S1_all, S1_lvl
    
    # Filter by rank if needed
    if rank_filter:
        if 'lvl_type' in df.columns:
            df = df[df['lvl_type'] == rank_filter]
    
    # Prepare ID column
    id_col = 'taxid' if use_taxid else 'name'
    if id_col not in df.columns:
        # Maybe "name" is indented, stripped handled by pandas if sep=\t?
        # Actually name usually has spaces.
        pass

    # Function to create subtable
    def create_subtable(suffix):
        cols_to_keep = [id_col]
        rename_map = {id_col: '#Taxon'}
        
        for s_id, s_name in sample_map.items():
            col_name = f"{s_id}_{suffix}"
            if col_name in df.columns:
                cols_to_keep.append(col_name)
                rename_map[col_name] = s_name
        
        sub_df = df[cols_to_keep].copy()
        sub_df.rename(columns=rename_map, inplace=True)
        
        # Clean name column if it's name (strip indentation)
        if id_col == 'name' and '#Taxon' in sub_df.columns:
             sub_df['#Taxon'] = sub_df['#Taxon'].str.strip()
             
        return sub_df

    # 3. Split and Save
    df_all = create_subtable('all')
    df_lvl = create_subtable('lvl')
    
    out_all = f"{output_basename}_all.tsv"
    out_lvl = f"{output_basename}_lvl.tsv"
    
    df_all.to_csv(out_all, sep='\t', index=False)
    df_lvl.to_csv(out_lvl, sep='\t', index=False)
    
    typer.echo(f"Created {out_all} and {out_lvl}")
