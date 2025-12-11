
```bash
# Combine Kraken Reports
combine_kreports.py \
	-o krakentools-combined/combine_kreports_kraken.tsv -r kraken-reports/*tsv

# Combine Bracken Report
combine_kreports.py \
	-o krakentools-combined/combine_kreports_bracken.tsv -r bracken
-report/*brep

# Combine Bracken Outputs
combine_bracken_outputs.py \
	-o krakentools-combined/combine_bracken_outputs.tsv --files bracken-output/* 
```
