library(PhyloOrchard)

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(adegenet, data.table, ape, ggtree, phangorn, ggnewscale, ggtreeExtra, ggimage, cowplot, ggpubr, gridExtra,
               ggplot2, picante, aplot, tidyr, phytools, RColorBrewer, microViz, PhyloOrchard)


# get list of colors 
list_color_palettes = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, list_color_palettes$maxcolors, rownames(list_color_palettes)))

# assign one color for each phylogroup and keep the same theme throughout the analysis
colors.variable <- data.frame(
  phylogroup = unique(metadata$Phylogroup),
  mycolors = color_vector[1:57]
)

tree <- ape::read.tree("results/raxml/pseudoalteromonas/pangenome/RAxML_bipartitionsBranchLabels.pseudoalteromonas-pangenome")
metadata <- as.data.frame(read.table("metadata/pseudoalteromonas-large-dataset.txt", sep = "\t", header = T, fill = T, ))

vitamins <- read.table("results/metabolic/METABOLIC-vitamins.txt", header = T, sep = "\t", row.names = 1)
secretion_systems <- read.table("results/metabolic/METABOLIC-secretion-system.txt", header = T, sep = "\t", row.names = 1)
nutrient_cycles <- read.table("results/metabolic/METABOLIC-nutrient-cycles.txt", header = T, sep = "\t", row.names = 1)
fermentation <- read.table("results/metabolic/METABOLIC-fermentation.txt", header = T, sep = "\t", row.names = 1)
cazymes <- read.table("results/metabolic/METABOLIC-cazymes.txt", sep = "\t", header = T, row.names = 1)
#drug_resistance <- read.table("results/05-metabolic/METABOLIC-drug-resistance.txt", header = T, sep = "\t", row.names = 1)

tree.unroot <- unroot(tree)
tree.midpoint <- midpoint(tree.unroot)

p <- ggtree(tree.midpoint, ladderize = T, right = F, aes(color = Phylogroup), branch.length = "none") %<+% metadata +
  scale_color_manual(values = color_vector, guide = "none") +
  new_scale_color() +
#  geom_tiplab(size = 2.5, offset = 0.8, align = T, linesize = 0) +
  ggplot2::xlim(0, 50) +
  ggplot2::ylim(-30, 240) +
  geom_tippoint(aes(color = Lifestyle), size = 0.5) + 
  scale_color_manual(values = c("Host-associated" = "red", "Free-living" = "blue"), guide = "none") +
  theme(legend.position = 'right',
        plot.margin = unit(c(0,0,0,0), "mm")) 

# plot cazyme/peptidases summary
#p2 <- p + geom_fruit(data = cazymes, geom = geom_col(), mapping=aes(y=Samples, x=Abundance, fill = Type),
#           position = position_dodgex()) + scale_fill_manual(values = c("black", "grey50") )+ new_scale_fill() 

# plot vitamin biosynthesis
p2 <- gheatmap(p, vitamins, width=0.09, colnames = T, font.size = 1, colnames_angle = 90, hjust = 1, color = "grey70", 
               offset = 0.02) +
  scale_fill_manual(values=c("white", "salmon4", "salmon3", "salmon2", "salmon1", "salmon"), guide = "none") + new_scale_fill() 
#  theme(legend.position = "none") +
#  theme(legend.text = element_text(size=6), legend.key.size = unit(3, 'mm'), legend.title = element_text(size = 8))
  
# secrection systems
p3 <- gheatmap(p2, secretion_systems, width=0.12, colnames = T, font.size = 1, colnames_angle = 90, hjust = 1, color = "grey70", 
               offset = 0.11) +
  scale_fill_manual(values=c("white", "seagreen4", "seagreen3", "seagreen2", "seagreen1", "seagreen"), guide = "none") + new_scale_fill()
#  theme(legend.position = "none") +
#  theme(legend.text = element_text(size=6), legend.key.size = unit(3, 'mm'), legend.title = element_text(size = 8))
  
# nutrient cycles
p4 <- gheatmap(p3, nutrient_cycles, width=0.15, colnames = T, font.size = 1, colnames_angle = 90, hjust = 1, color = "grey70", 
               offset = 0.23) +
  scale_fill_manual(values=c("white", "darkblue", "darkcyan", "cyan3", "cyan2", "cyan1", "lightcyan"), guide = "none") + new_scale_fill( )
 # theme(legend.position = "none") +
#  theme(legend.text = element_text(size=6), legend.key.size = unit(3, 'mm'), legend.title = element_text(size = 8))
  
# fermentation
p5 <- gheatmap(p4, fermentation, width=0.20, colnames = T, font.size = 1, colnames_angle = 90, hjust = 1, color = "grey70", 
               offset = 0.38) +
  scale_fill_manual(values=c("white", "antiquewhite4", "antiquewhite3", "antiquewhite2", "antiquewhite1", "burlywood3"), guide = "none") + new_scale_fill() 
  #theme(legend.position = "none") +
  #theme(legend.text = element_text(size=6), legend.key.size = unit(3, 'mm'), legend.title = element_text(size = 8))
  
#cazymes
p6 <- gheatmap(p5, cazymes, width=2, colnames = T, font.size = 1, colnames_angle = 90, hjust = 1, color = "white", 
               offset = 0.58) +
  scale_fill_gradient(low = "white", high = "darkblue") +
#  theme(legend.position = "none") +
  theme(legend.text = element_text(size=6), legend.key.size = unit(12, 'mm'), legend.title = element_text(size = 8))

# drug resistance
#p7 <- gheatmap(p6, drug_resistance, width=0.75, colnames = T, font.size = 3, colnames_angle = 90, hjust = 1, color = "grey70", 
#               offset = 33.75) +
#  scale_fill_manual(values=c("white", "grey90", "grey80", "grey70", "grey60", "grey50", "grey40", "grey30", "grey20", "grey10", "grey0", "grey85" , "grey65","grey55")) + new_scale_fill() +
#  theme(legend.position = "none") +
#  theme(legend.text = element_text(size=6), legend.key.size = unit(3, 'mm'), legend.title = element_text(size = 8))

png("figures/metabolic/circular-poster.png", res = 300, width = 11.5, height = 11.5, units = "in")
print(p6)
dev.off()


#cazymes
p6 <- gheatmap(p5, cazymes, width=2, colnames = T, font.size = 1.5, colnames_angle = 90, hjust = 1, color = "white", 
               offset = 0.75) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  #theme(legend.position = "none") +
  theme(legend.text = element_text(size=6), legend.key.size = unit(2, 'mm'), legend.title = element_text(size = 8))


#cazymes <- read.table("results/metabolic/METABOLIC-cazymes.txt", sep = "\t", header = T, row.names = 1)
#cazymes.long <- gather(cazymes, cazyme, abundance, factor_key = T, GH001:PL042)
#cazymes_heatplot <- ggplot(cazymes.long, aes(x=reorder(cazyme, -abundance), y=SampleID)) + 
#  geom_tile(aes(fill=abundance)) + scale_fill_viridis_c() + 
#  scale_fill_gradient(low = "white", high = "darkblue", limits=c(0,15)) +
#  theme_minimal() + xlab(NULL) + ylab(NULL) +
#  theme(axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face ="bold", size = 15),
#        legend.key.height= unit(2, 'cm'))

#cazyme_final_figure <- cazymes_heatplot %>% insert_left(p)
#cazyme_final_figure

#png("figures/metabolic/figure-cazymes.png", res = 300, width = 15, height = 20, units = "in")
#print(cazyme_final_figure)
#dev.off()

#peptidases
#peptidases <- read.table("results/metabolic/METABOLIC-", sep = "\t", header = T)
#peptidases.long <- gather(peptidases, peptidase, abundance, factor_key = T, PL001:PL042)
#peptidase_heatplot <- ggplot(peptidases.long, aes(x=reorder(peptidase, -abundance), y=Sample)) + 
#  geom_tile(aes(fill=abundance)) + scale_fill_viridis_c() + 
#  scale_fill_gradient(low = "white", high = "darkblue", limits=c(0,15)) +
#  theme_minimal() + xlab(NULL) + ylab(NULL) +
#  theme(axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face ="bold", size = 10),
#        legend.key.height= unit(2, 'cm'))

#peptidase_final_figure <- peptidase_heatplot %>% insert_left(p7, width = 20)
#peptidase_final_figure

png("figures/metabolic/figure-peptidases.png", res = 300, width = 25, height = 10, units = "in")
print(p2)
dev.off()


#RPKM
#rpmk_final_figure <- rpkm_heatplot %>% insert_left(p7, width = 1)
#rpmk_final_figure

#png("figures-final/metabolic-pathways/figure-rpkm_metagenomes_read_mapping.png", res = 300, width = 12, height = 15, units = "in")
#print(rpmk_final_figure)
#dev.off()
#png("figures-final/metabolic-pathways/figure-peptidases.png", res = 300, width = 25, height = 10, units = "in")
#print(peptidase_final_figure)
#dev.off()
