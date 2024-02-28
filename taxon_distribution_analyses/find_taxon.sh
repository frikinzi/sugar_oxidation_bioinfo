#!/bin/bash
gene_file=$1
output=$2
# this script outputs a tab separated table with the MGnify accession of the genome and its corresponding taxon
while IFS= read -r name; do
    echo finding "$name"
    grep "$name" MGnify_accession_Lineage.txt >> $gene_file
done < $output
