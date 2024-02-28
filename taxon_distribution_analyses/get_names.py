import os

#outputs a txt file list of genomes that have function hits given a folder
newfile = open("genomes_with_significant_hits_garr.txt", "w")
for i in os.listdir("gudgar_significant"):
	newfile.write(i.split(".")[1] + "\n")

newfile.close()
