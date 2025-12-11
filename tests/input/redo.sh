#!/bin/bash
set -euxo pipefail

kraut make-table -o kraut/make-table-kraken.tsv kraken-reports/*tsv

kraut make-table -o kraut/make-table-bracken.tsv bracken-report/*brep

kraut split-combine-table -i krakentools-combined/combine_kreports_kraken.tsv -o kraut/split_kraken_combined

kraut split-combine-table -i krakentools-combined/combine_kreports_bracken.tsv -o kraut/split_bracken_combined
