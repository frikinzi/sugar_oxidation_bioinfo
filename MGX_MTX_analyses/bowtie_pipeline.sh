#!/bin/bash
module load bowtie
#perform for each batch for parallel processing
batch=$1

#get a list of all prefixes
prefixes=($(ls -1 2018-05-04-unprocessed/"${batch}"/*_R1.fastq.gz | sed 's/_R1.fastq.gz//'))

# perform mapping for every file
for prefix in "${prefixes[@]}"; do

# map reads
bowtie2 -x ibd_genes -1 ${prefix}_R1.fastq.gz -2 ${prefix}_R2.fastq.gz -S sam/${prefix}.sam

# convert sam to bam
samtools view -bS sam/${prefix}.sam > bam/${prefix}.bam

# create tables with mapped reads count
samtools flagstat bam/${prefix}.bam -O tsv > tables_mgx/${prefix}.tsv

done

# creates big table of reads per sample for R. assume folder is named tables_mgx
python3 create_reads_table.py
