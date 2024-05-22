library(ggplot2)
library(gggenes)

genes <- read.table("operons/GCA_000312365.2_neighborhood.txt", header=FALSE, sep="\t",
                    col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))

genes2 <- read.table("operons/GCA_002858805.1_neighborhood.txt", header=FALSE, sep="\t",
                    col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))
gene3 <- read.table("operons/GCA_000479045.1_neighborhood.txt", header=FALSE, sep="\t",
                    col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))

gene4 <- read.table("operons/C_clos.neighborhood.txt",header=FALSE,col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))
gene5 <- read.table("operons/GCA_003697165.2_neighborhood.txt",sep="\t",header=FALSE,col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))
ecoli_operon2 <-read.table("operons/e_coli_operon2.txt",sep="\t",header=FALSE,col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))

gene6 <- read.table("operons/salmonella_neighborhood.txt",sep="\t",header=FALSE,col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))

gene7 <- read.table("operons/salmonella_operon2.txt",sep="\t",header=FALSE,col.names=c("Organism", "Query_Length", "Subject_ID", "Subject_Operon", "Start", "End", "Strand", "Annotation"))


combined_df <- rbind(genes, genes2, gene3, gene4, gene5,ecoli_operon2,gene6, gene7)

genes_plot <- ggplot(combined_df, aes(xmin = Start, xmax = End, y = Subject_Operon, forward = Strand == "+", fill = Annotation)) +
  geom_gene_arrow() +  # Removed the fill aesthetic
  geom_gene_label(aes(label = Annotation), size = 3) + 
  facet_wrap(~ Organism, scales = "free_y", ncol = 1) +  # Separate plots for each organism
  theme_minimal() +
  scale_fill_manual(values = c("garD" = "#FF7F00", "garK" = "#B15928", "garR" = "#6A3D9A", "gudD" = "#33A02C", "gudL" = "#E31A1C", "garL" = "#1F78B4", "garABC" = "#F781BF", "garP" = "darkturquoise", "gudX" = "goldenrod", "gudP" = "darkolivegreen")) +
  labs(title = "Gene Operons", x = "Genomic Position", y = "") +
  theme(axis.text.y = element_text(angle = 0), strip.background = element_blank(), strip.text.x = element_text(face = "bold"))

ggsave(file="neighborhood.svg", plot=genes_plot, width=6, height=10)
