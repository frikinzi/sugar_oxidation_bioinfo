library(ggtree)
library(ComplexHeatmap)


tree <- read.tree("/Users/angela/Downloads/polymorphum.newick")

pa_matrix <- read.csv("/Users/angela/Downloads/polymorphum_presence_absence.csv", sep="\t",row.names=1)

heatmap_data <- as.matrix(pa_matrix)

tree_plot <- ggtree(tree)

heatmap_data <- ifelse(heatmap_data == 1, "present", "absent")

transposed <- t(heatmap_data)

heat <- gheatmap(tree_plot, transposed, width=4, font.size=10)+
  scale_fill_manual(breaks=c("absent", "present"), 
                    values=c("lightgray", "steelblue"), name="gene")+
  theme(
    legend.title = element_text(size = 12),  # Adjust the legend title size
    legend.text = element_text(size = 12),   # Adjust the legend text/label size
    legend.key.size = unit(1.5, "cm")        # Adjust the size of the legend keys
  )

ggsave("test.svg", plot = heat, width = 10, height = 25, dpi = 300)

