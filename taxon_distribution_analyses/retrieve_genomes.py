'''
Script to grab the nucleotide sequences of the significant genes

All of the nucleotide sequences will be output in separate files
'''

import pandas as pd
import os
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq

gut_genomes = pd.read_csv("genomes-all_metadata.tsv", sep="\t")

list_files = open("ibd_genes.txt", "r")

lines = list_files.readlines()

list_files.close()

for i in lines:
	genome_name = i.split("_0")[0].strip()
	result_row = gut_genomes[gut_genomes["Genome"] == genome_name]
	mgx = result_row
	mgnify = result_row.at[result_row.index[0], 'MGnify_accession']
	genome_dir = f'/data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/{mgnify[:-2]}/{mgnify}/genome/{mgnify}.faa'
	gff_dir = f'/data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/{mgnify[:-2]}/{mgnify}/genome/{mgnify}.gff'
	fna_dir = f'/data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/{mgnify[:-2]}/{mgnify}/genome/{mgnify}.fna'
	os.system(f'grep -A 1 {i.strip()} {genome_dir} >> results.fa') # amino acids
	grep_command = f"grep -i {i.strip()} {gff_dir}"
	result = subprocess.check_output(grep_command, shell=True, text=True)
	print(result)
	result_row = result.split("\t")
	start_ind = int(result_row[3])-1
	end_ind = int(result_row[4])
	contig = result_row[0]
	rev = result_row[6]
	os.system(f"")
	ntd = open("tmp.bed", "w")
	ntd.write(f"{contig}\t{str(start_ind)}\t{end_ind}\t1\t{rev}\t{rev}\n")
	ntd.close()
	os.system(f"fastaFromBed -fi {fna_dir} -fo {i.strip()}.fa -bed tmp.bed -s")


