# use pacman package to install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, ggpubr, data.table)

checkm <- fread("results/checkm/checkm-statistics.txt", sep = "\t", header = T)
metadadata <- fread("metadata/pseudoalteromonas-large-dataset.tsv", sep = "\t", header = T)

#merge checkm with metadata
checkm_metadata <- merge(checkm, metadadata, by.x = "Bin Id", by.y = "SampleID")

# GC Content
codingDensityGC <- ggplot(checkm_metadata, aes(x=`Lifestyle`, y=`GC`, fill = Lifestyle)) + 
  geom_boxplot(outlier.colour="black", outlier.size=2) +
  scale_fill_manual(values = c("Host-associated" = "darkred", "Free-living" = "dodgerblue4")) +
  ylim(30,55) +
  stat_compare_means(label.x = 1.1, label.y = 55) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 
  

# Coding Density
codingDensityLifestyle <- ggplot(checkm_metadata, aes(x=`Lifestyle`, y=`Coding density`, fill = Lifestyle)) + 
  geom_boxplot(outlier.colour="black", outlier.size=2) +
  scale_fill_manual(values = c("Host-associated" = "darkred", "Free-living" = "dodgerblue4")) +
  ylim(85,95) +
  stat_compare_means(label.x = 1.1, label.y = 95) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 

# Number of Genes
numberGenesLifestyle <- ggplot(checkm_metadata, aes(x=`Lifestyle`, y=`# predicted genes`, fill = Lifestyle)) + 
  geom_boxplot(outlier.colour="black", outlier.size=2) +
  scale_fill_manual(values = c("Host-associated" = "darkred", "Free-living" = "dodgerblue4")) +
  ylim(3250,6000) + ylab("Number of Genes") +
  stat_compare_means(label.x = 1.1, label.y = 6000) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

# Genome Size 
genomeSizeLifestyle <- ggplot(checkm_metadata, aes(x=`Lifestyle`, y=`Genome size (bp)`, fill = Lifestyle)) + 
  geom_boxplot(outlier.colour="black", outlier.size=2) +
  scale_fill_manual(values = c("Host-associated" = "darkred", "Free-living" = "dodgerblue4")) +
  #ylim(3500000,6800000) +
  scale_y_continuous(labels=function(x)x/1000000) + ylab("Genome Size (Mbp)") +
  stat_compare_means(label.x = 1.1, label.y = 6800000) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 


# correlation genome size and CDS
correlation <- ggplot(checkm_metadata, aes(x=`Genome size (bp)`, y=`# predicted genes`, color=Lifestyle)) + geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + theme_bw() + ylab("Number of Genes") + xlab("Genome Size (Mbp)") +
  scale_color_manual(values = c("Host-associated" = "darkred", "Free-living" = "dodgerblue4")) +
  scale_x_continuous(labels=function(x)x/1000000) + theme(legend.position = "bottom") +
  theme_bw()

genomeInfo <- ggarrange(genomeSizeLifestyle, numberGenesLifestyle, codingDensityLifestyle, codingDensityGC, 
          ncol = 4, legend = "none")

genomeInfoCorrelation <- ggarrange(genomeInfo, correlation, nrow = 2, common.legend = T, align = "h", legend = "bottom")



png("figures/checkm/2024-10-10-genome-size-density-gc.png", res = 300, width = 11.8, height = 8, units = "in")
print(genomeInfoCorrelation)
dev.off()

