# use pacman package to install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, data.table, tidyverse, ggpubr, vegan, forcats, 
               gridGraphics, reshape, svglite, dplyr, tidyr, magrittr, ggtree,
               ggnewscale)

tree <- ape::read.tree("results/raxml/pseudoalteromonas/pangenome/RAxML_bipartitionsBranchLabels.pseudoalteromonas-pangenome")
metadata <- as.data.frame(read.table("metadata/pseudoalteromonas-large-dataset.txt", sep = "\t", header = T, fill = T, ))

scoary_host_pigmented <- read.table("results/scoary2-pigmented-pirate/traits/host-associated/pirate-scoary-pigmented-host-associated.txt", header = T, sep = "\t", row.names = 1)
scoary_host_nonpigmented <- read.table("results/scoary2-nonpigmented-pirate/traits/host-associated/pirate-scoary-nonpigmented-host-associated.txt", header = T, sep = "\t", row.names = 1)

# get list of colors 
list_color_palettes = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, list_color_palettes$maxcolors, rownames(list_color_palettes)))

# assign one color for each phylogroup and keep the same theme throughout the analysis
colors.variable <- data.frame(
  phylogroup = unique(metadata$Phylogroup),
  mycolors = color_vector[1:57]
)

tree.unroot <- unroot(tree)
tree.midpoint <- midpoint(tree.unroot)

p <- ggtree(tree.midpoint, ladderize = T, right = F, aes(color = Phylogroup)) %<+% metadata +
  scale_color_manual(values = color_vector, guide = "none") +
  new_scale_color() +
  # geom_tiplab(size = 1, offset = 3, align = T, linesize = 0) +
  ggplot2::xlim(0, 3.75) +
  ggplot2::ylim(-75, 240) +
  geom_tippoint(aes(color = Lifestyle, x = 0.915), size = 0.5) + 
  scale_color_manual(values = c("Host-associated" = "red", "Free-living" = "blue"), guide = "none") +
  theme(legend.position = 'right',
        plot.margin = unit(c(2,1,2,1), "mm")) 


# plot vitamin biosynthesis
p2 <- gheatmap(p, scoary_host_pigmented, width=0.1, colnames = T, font.size = 1.5, colnames_angle = 90, hjust = 1, color = "grey90", 
               offset = 0.02) +
  scale_fill_manual(values=c("white", "orange"), guide = "none") + new_scale_fill() 

p3 <- gheatmap(p2, scoary_host_nonpigmented, width=2, colnames = T, font.size = 1.5, colnames_angle = 90, hjust = 1, color = "grey90", 
               offset = 0.20) +
  scale_fill_manual(values=c("white", "purple"), guide = "none") + new_scale_fill() 
p3

png("figures/metabolic/2025-06-14-scoary-summary.png", res = 300, width = 14, height = 14, units = "in")
print(p3)
dev.off()
