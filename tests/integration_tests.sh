#!/bin/bash

# Integration tests for kraut
# Assumes 'kraut' is installed and available in the path

set -e

echo "Starting integration tests..."

# Ensure we have input data
INPUT_DIR="${1:-kraken}"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory \"$INPUT_DIR\" not found. Please run this script from the project root or provide a valid path."
    exit 1
fi

TEST_OUT_DIR="tests_output"
mkdir -p "$TEST_OUT_DIR"

echo "[1/4] Testing single-report..."
kraut single-report -i "$INPUT_DIR/genome/Segatella_copri.tsv" -o "$TEST_OUT_DIR/single_report.txt"
if [ -s "$TEST_OUT_DIR/single_report.txt" ]; then
    echo "  PASS: Output created"
else
    echo "  FAIL: Output empty or missing"
    exit 1
fi
# Check if output is equal to input
diff "$TEST_OUT_DIR/single_report.txt" "$INPUT_DIR/genome/Segatella_copri.tsv"
if [ $? -ne 0 ]; then
    echo "  FAIL: Output does not match input"
    exit 1
fi

echo "[2/4] Testing make-table default output..."
kraut make-table "$INPUT_DIR/shredded-mixes-filtered"/*.tsv -o "$TEST_OUT_DIR/merged.tsv"
if [ -s "$TEST_OUT_DIR/merged.tsv" ]; then
    echo "  PASS: Output created"
else
    echo "  FAIL: Output empty or missing"
    exit 1
fi

echo "[3/4] Testing make-table options..."
kraut make-table "$INPUT_DIR/shredded-mixes-filtered"/*.tsv -o "$TEST_OUT_DIR/table_tot.tsv" --metric TOT --rank S --rank-prefix
if [ -s "$TEST_OUT_DIR/table_tot.tsv" ]; then
    echo "  PASS: Output created"
else
    echo "  FAIL: Output empty or missing"
    exit 1
fi

echo "[4/4] Testing split-combine-table..."
# Ensure input file exists (created in previous steps or assumed)
if [ -f "$INPUT_DIR/shredded-mix-merge.tsv" ]; then
    kraut split-combine-table -i "$INPUT_DIR/shredded-mix-merge.tsv" -o "$TEST_OUT_DIR/split" --rank S --taxid
    if [ -s "$TEST_OUT_DIR/split_all.tsv" ] && [ -s "$TEST_OUT_DIR/split_lvl.tsv" ]; then
        echo "  PASS: Outputs created"
    else
        echo "  FAIL: Split outputs missing"
        exit 1
    fi
else
    echo "  SKIP: Input for split-combine-table not found"
fi

echo "All tests passed successfully!"
rm -rf "$TEST_OUT_DIR"
