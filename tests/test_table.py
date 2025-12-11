import pytest
from pathlib import Path
from kraut.models.kraken_data import KrakenReport, KrakenNode
from kraut.models.multi_report import MultiKrakenReport
import pandas as pd
import io

# Helper to create synthetic report
def create_report(counts):
    # counts: dict {taxid: (clade_count, taxon_count, rank, name)}
    r = KrakenReport()
    r.nodes = {}
    for tid, (cc, tc, rank, name) in counts.items():
        node = KrakenNode(0, cc, tc, rank, tid, name)
        r.nodes[tid] = node
        if tid == 1: r.root = node
        if tid == 0: r.unclassified = node
    
    # Manually link children for traversal if needed, but MultiKrakenReport primarily iterates nodes if we add manually? 
    # Actually add_report traverses from root. So we must link.
    # Simple hierarchy: 1 -> others
    if 1 in r.nodes:
         for tid, node in r.nodes.items():
             if tid not in [0, 1]:
                 r.nodes[1].add_child(node)
                 
    return r

def test_make_table_logic():
    r1 = create_report({
        1: (100, 0, 'R', 'root'), 
        2: (50, 50, 'S', 'E. coli'),
        0: (10, 10, 'U', 'unclassified')
    })
    
    mr = MultiKrakenReport()
    mr.add_report(r1, "S1")
    
    # Test TOT metric (clade counts)
    df = mr.to_dataframe(metric='TOT', level='S') 
    assert not df.empty
    row = df[df['#Taxon'] == 'E. coli'].iloc[0]
    assert row['S1'] == 50
    
    # Test prefixes
    df = mr.to_dataframe(metric='TOT', level='S', rank_prefix=True)
    assert 's__E. coli' in df['#Taxon'].values
    
    # Test TaxID
    df = mr.to_dataframe(metric='TOT', level='S', use_taxid=True)
    assert '2' in df['#Taxon'].values
