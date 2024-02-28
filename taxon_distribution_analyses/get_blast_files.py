import os

#simple script to extract the blast tables from the raw ProkFunFind output given an output folder with the genomes with the function
for i in os.listdir("gudgar_significant"):
        name = (i.split(".")[1] + "\n").strip()
        os.system("cp gudd_gard_garr/ggg." + name + ".tsv significant_hits_garr/")
