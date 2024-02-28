import pandas as pd
import os
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
import subprocess

# builds amino acid file given a list of genes

gut_genomes = pd.read_csv("genomes-all_metadata.tsv", sep="\t")

list_files = open("gard_gudd_genes.txt", "r")

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
	#os.system(f'grep -A 1 {i.strip()} {genome_dir} >> aa_garl.fa') # amino acids
	os.system(f'./get_aa_genes.sh {i.strip()} {genome_dir.strip()}')


