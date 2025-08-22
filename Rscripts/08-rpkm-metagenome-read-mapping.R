library(adegenet)
library(ape)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(ggtreeExtra)
library(ggimage)
library(microViz)
library(tidyr)
library(aplot)
library(ggdendroplot)
library(RColorBrewer)
library(sets)
library(dplyr)
library(purrr)

# use papurrr# use pacman package to install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(adegenet, ape, ggtree, phangorn, ggnewscale, ggtreeExtra, 
               ggimage, microViz, tidyr, aplot, ggdendro, RColorBrewer, sets)


# get list of colors 
list_color_palettes = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, list_color_palettes$maxcolors, rownames(list_color_palettes)))

# assign one color for each phylogroup and keep the same theme throughout the analysis
colors.variable <- data.frame(
  phylogroup = unique(metadata$Phylogroup),
  mycolors = color_vector[1:57]
)

tree <- ape::read.tree("results/raxml/pseudoalteromonas/pangenome/RAxML_bestTree.pseudoalteromonas-pangenome")
rpkm <- read.table("results/rrap/nematode_metagenomes_rpkm_noLog-replicates-only.csv", sep = ",", header = T, check.names = F)
metagenome_metadata <- read.table("metadata/metagenomes_metadata.tsv", sep = "\t", header = T)

#midroot the tree
tree.unroot <- unroot(tree)
tree.midpoint <- midpoint(tree.unroot)

#convert rpkm values wide to long form
rpkm.long <- gather(rpkm, metagenome, RPKM, factor_key = T, "metoncholaimus.1":"adoncholaimus.11")
rpkm.metadata <- merge.data.frame(rpkm.long, metagenome_metadata, by.x = "metagenome", by.y = "nematodeID")

phylo_tree <- ggtree(tree.midpoint, ladderize = T, right = F, aes(color = Phylogroup)) %<+% metadata +
  scale_color_manual(values = color_vector, guide = "none") +
  new_scale_color() +
  geom_tiplab(size = 1.5, offset = 0, align = T, linesize = 0) +
  ggplot2::xlim(0, 1) +
  ggplot2::ylim(-30, 240) +
  geom_tippoint(aes(color = Lifestyle, x = 0.915), size = 0.5) + 
  scale_color_manual(values = c("Host-associated" = "red", "Free-living" = "blue"), guide = "none") +
  theme(legend.position = 'right',
        plot.margin = unit(c(2,1,2,1), "mm")) 

rpkm_heatplot <- ggplot(rpkm.metadata, aes(x=reorder(metagenome, -RPKM), y=ACC)) + 
  geom_tile(aes(fill=RPKM)) +
  scale_fill_gradient(low = "white", high = "blue", space = "Lab" , limits=c(0.1, 7), breaks=seq(1, 7, by=1), na.value = 'white') +
  theme_void() + xlab(NULL) + ylab(NULL) +
  facet_grid( .~Order + Family, scales = "free", space = "free") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        legend.key.height = unit(2, 'cm')) 

rpkm_heatplot

rpmk_final_figure <- rpkm_heatplot %>% insert_left(phylo_tree, width = 3) 
rpmk_final_figure

png("figures/rrap/2024-11-22-rpkm_metagenomes_read_mapping-5.png", res = 300, width = 8, height = 11.5, units = "in")
print(rpmk_final_figure)
dev.off()


rpkm.onch <- rpkm.metadata %>%
  filter(Family == "Oncholaimidae") %>%
  filter(RPKM > 0.01) %>%
  group_by(metagenome) %>%
  summarize(sum(RPKM))
mean(rpkm.onch$`sum(RPKM)`)

rpkm.lepto <- rpkm.metadata %>%
  filter(Family == "Leptosomatidae") %>%
  filter(RPKM > 0.01) %>%
  group_by(metagenome) %>%
  summarize(sum(RPKM))
mean(rpkm.lepto$`sum(RPKM)`)

rpkm.thora <- rpkm.metadata %>%
  filter(Family == "Thoracostomopsidae") %>%
  filter(RPKM > 0.01) %>%
  group_by(metagenome) %>%
  summarize(sum(RPKM))
mean(rpkm.thora$`sum(RPKM)`)

rpkm.sphar <- rpkm.metadata %>%
  filter(Family == "Sphaerolaimidae") %>%
  filter(RPKM > 0.01) %>%
  group_by(metagenome) %>%
  summarize(sum(RPKM)) 
mean(rpkm.sphar$`sum(RPKM)`)
