import pytest
from pathlib import Path
from kraut.models.kraken_data import KrakenReport, KrakenNode
from kraut.models.multi_report import MultiKrakenReport

# Sample data path
SAMPLE_DATA = Path(__file__).parent / "input/kraken-reports/Segatella_copri.tsv"

def test_parsing_and_reproduction():
    """
    Test that we can parse a file and reproduce it exactly.
    """
    if not SAMPLE_DATA.exists():
        pytest.skip(f"Sample data {SAMPLE_DATA} not found")
        
    report = KrakenReport.from_file(str(SAMPLE_DATA))
    
    assert report.root is not None
    assert report.root.name == "root"
    assert report.root.tax_id == 1
    
    # Check reproduction
    original_content = SAMPLE_DATA.read_text()
    generated_content = report.to_string()
    
    # Debugging: print differences if assertion fails
    if original_content != generated_content:
        print("Original:\n", original_content[:500])
        print("Generated:\n", generated_content[:500])
        
    assert generated_content == original_content

def test_node_structure():
    node = KrakenNode(100.0, 10, 0, "R", 1, "root")
    child = KrakenNode(50.0, 5, 5, "D", 2, "Bacteria", depth=1)
    node.add_child(child)
    
    assert len(node.children) == 1
    assert node.children[0].parent == node
    assert node.children[0].depth == 1
    
    # Check string formatting indentation
    # Root has depth 0 -> 0 spaces
    # Child has depth 1 -> 2 spaces
    assert "\troot" in node.to_string()
    assert "\t  Bacteria" in child.to_string()

def test_merging():
    # Create synthetic reports
    r1 = KrakenReport()
    n1_root = KrakenNode(100, 100, 0, "R", 1, "root")
    n1_bac = KrakenNode(50, 50, 50, "D", 2, "Bacteria", depth=1, parent=n1_root)
    n1_root.add_child(n1_bac)
    r1.root = n1_root
    r1.nodes = {1: n1_root, 2: n1_bac}
    
    r2 = KrakenReport()
    n2_root = KrakenNode(100, 200, 0, "R", 1, "root")
    n2_bac = KrakenNode(50, 100, 100, "D", 2, "Bacteria", depth=1, parent=n2_root)
    n2_root.add_child(n2_bac)
    r2.root = n2_root
    r2.nodes = {1: n2_root, 2: n2_bac}
    
    mr = MultiKrakenReport()
    mr.add_report(r1, "sample1")
    mr.add_report(r2, "sample2")
    
    df = mr.to_dataframe(metric="COUNTS", level="D")
    
    print(df)
    
    assert not df.empty
    # Updated to reflect new MultiKrakenReport behavior (returns #Taxon column)
    # Check for Bacteria by name
    row = df[df['#Taxon'] == "Bacteria"].iloc[0]
    
    # Or check we have correct counts
    assert row['sample1'] == 50
    assert row['sample2'] == 100
    assert row['sample1'] == 50
    assert row['sample2'] == 100
