#!/bin/bash
gene=$1
faa=$2

awk -v seq_name=$gene 'BEGIN {print_seq = 0} /^>/ {if (print_seq) {exit} print_seq = ($0 ~ seq_name)} {if (print_seq) print}' $faa >> gard_aa.fa
