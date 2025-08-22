library(phytools)
library(treeio)
library(phylotree)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(microViz)
library(RColorBrewer)
library(phangorn)
library(tidytree)
library(dplyr)

# use pacman package to install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, data.table, tidyverse, ggpubr, vegan, forcats, 
               gridGraphics, reshape, phytools, ggtree, ggtreeExtra, tidytree,
               RColorBrewer, phangorn, treeio, ggnewscale)

# import tree and metadata 
tree_bacteria <- read.raxml("results/raxml/pseudoalteromonas/pangenome/RAxML_bipartitionsBranchLabels.pseudoalteromonas-pangenome")
metadata <- as.data.frame(read.table("metadata/pseudoalteromonas-large-dataset.txt", sep = "\t", header = T, fill = T, ))

# get list of colors 
list_color_palettes = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, list_color_palettes$maxcolors, rownames(list_color_palettes)))

# assign one color for each phylogroup and keep the same theme throughout the analysis
colors.variable <- data.frame(
  phylogroup = unique(metadata$Phylogroup),
  mycolors = color_vector[1:57]
)


tree_plot <- ggtree(tree_bacteria) %<+% metadata +
  ggplot2::ylim(-1, 250) +
  ggplot2::xlim(0, 3) +
  geom_tiplab(aes(color = Clade), size = 1.25, offset = 0.25, align = T, linetype = 0) +
  scale_color_manual(values = c("Pigmented" = "darkorange", "Nonpigmented" = "purple")) +
  new_scale_color() +
  geom_tippoint(aes(shape = Lifestyle, color = Lifestyle)) +
  scale_color_manual(values = c("Host-associated" = "red", "Free-living" = "blue")) +
  guides(color = guide_legend(nrow = 1)) +
  geom_text2(aes(subset=(isTip==F & bootstrap < 80), label=bootstrap), size = 2, color = "black", hjust = 1, vjust = -.5) 

tree_plot_phylogroups <- tree_plot + new_scale_fill() +
  geom_fruit(geom=geom_tile, mapping=aes(fill = Phylogroup), width=.1, offset=0.1) +
  scale_fill_manual(values = color_vector) +
  guides(fill = guide_legend(ncol = 2), shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) +
  theme(legend.position=c(0.85, 0.5), legend.text = element_text(size = 8), legend.title = element_text(size = 10, face = "bold", hjust = 0.5)) 

tree_plot_phylogroups

png("figures/raxml/2024-10-10-tree-bootstraphs-not-rooted-2.png", res = 300, width = 8.5, height = 11, units = "in")
print(tree_plot_phylogroups)
dev.off()

