#!/bin/bash
module load FastTree
module load clustalw

set -e
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

file=$1
input_file=$file

# Output files
aligned_file=${file%.*}_aligned.fasta
trimmed_file=${aligned_file%.*}_aligned_cleaned.fasta

# Perform alignment with clustalw
clustalw $input_file -OUTPUT=FASTA -OUTFILE=$aligned_file

# Trim the alignment with Goalign
goalign clean sites --char=GAP -c 0.97 -i $aligned_file -o $trimmed_file

# Run FastTree on the trimmed alignment
FastTreeMP "$trimmed_file" > "$trimmed_file".tree
