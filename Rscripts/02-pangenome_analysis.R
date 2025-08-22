# use pacman package to install and load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, data.table, tidyverse, ggpubr, vegan, forcats, 
               gridGraphics, reshape, svglite)

# load functions
source("scripts/Rscripts/functions/functions.R")

# pangenomic analysis of pseudoalteromonas genomes
pangenome <- as.data.frame(fread("results/pirate/PIRATE.gene_families.tsv", sep = "\t", header = T))
metadata <- as.data.frame(fread("metadata/pseudoalteromonas-large-dataset.txt", sep = "\t", header = T))


# CREATE PRESENCE/ANBSENCE MATRIX WITH METADATA
rownames(metadata) <- metadata$SampleID 

df <- create_matrix(pangenome, metadata)
df_host <- subset(df, Lifestyle == "Host-associated")
df_free <- subset(df, Lifestyle == "Free-living")

n_isolates <- nrow(df)
n_host <- nrow(df_host)
n_free <- nrow(df_free)

# GET PANGENOME DISTRIBUTION PLOTS BY FAMILY
# 45598 is the number of gene families according to the PIRATE output and indicate number of columns in dataframe
df_gene_distribution <- distribution_gene_families(df, 45598) # number of gene families
df_gene_distribution_plot <- plot_distribution_gene_families(df_gene_distribution) +
  theme(legend.position = "bottom") +
  ylim(0,20000) + xlim(0,235) + labs(y = "n gene families \n", x = "\n n genomes in family")

df_host_gene_distribution <- distribution_gene_families(df_host, 45598)
df_host_gene_distribution_plot <- plot_distribution_gene_families(df_host_gene_distribution) + 
  theme(legend.position = "none") +
  ylim(0,20000) + xlim(0,150)

df_free_gene_distribution <- distribution_gene_families(df_free, 45598)
df_free_gene_distribution_plot <- plot_distribution_gene_families(df_free_gene_distribution) +
  theme(legend.position = "none") +
  ylim(0, 20000) + xlim(0, 95)


# GET PANGENOME DISTRIBUTION PIE PLOT
n_pangenome_all <- df_gene_families(df_gene_distribution, n_isolates)
n_pangenome_host <- df_gene_families(df_host_gene_distribution, n_host)
n_pangenome_free <- df_gene_families(df_free_gene_distribution, n_free)

pie_chart_all <- plot_piechart(n_pangenome_all) 
pie_chart_host <- plot_piechart(n_pangenome_host)
pie_chart_free <- plot_piechart(n_pangenome_free)

# PLOT ACCUMULATION CURVES 
df_sans_meta <- df[1:45598]
df_sans_meta_host <- df_host[1:45598]
df_sans_meta_free <- df_free[1:45598]

df_sans_meta_accy <- df_sans_meta[, !colSums(df_sans_meta) >= (234 * .95)]
df_sans_meta_accy_sans_uniq <- df_sans_meta_accy[, colSums(df_sans_meta_accy) > 1 ]
df_sans_meta_accy_host <- df_sans_meta_host[, !colSums(df_sans_meta_host) >= (143 * .95)]
df_sans_meta_accy_sans_uniq_host <- df_sans_meta_accy_host[, colSums(df_sans_meta_accy_host) > 1 ]
df_sans_meta_accy_free <- df_sans_meta_free[, !colSums(df_sans_meta_free) >= (91 * .95)]
df_sans_meta_accy_sans_uniq_free <- df_sans_meta_accy_free[, colSums(df_sans_meta_accy_free) > 1 ]

# calculate accumulation curves
df_sans_meta_accy_curve <- specaccum(df_sans_meta_accy, method = "random", permutations = 100)
df_sans_meta_accy_sans_uniq_curve <- specaccum(df_sans_meta_accy_sans_uniq, method = "random", permutations = 100)
df_sans_meta_accy_curve_host <- specaccum(df_sans_meta_accy_host, method = "random", permutations = 100)
df_sans_meta_accy_sans_uniq_curve_host <- specaccum(df_sans_meta_accy_sans_uniq_host, method = "random", permutations = 100)
df_sans_meta_accy_curve_free <- specaccum(df_sans_meta_accy_free, method = "random", permutations = 100)
df_sans_meta_accy_sans_uniq_curve_free <- specaccum(df_sans_meta_accy_sans_uniq_free, method = "random", permutations = 100)

# accumulation curves of accessory + singelton genes
par(mar=c(6,6,1,1))
plot(df_sans_meta_accy_curve, col="#00000E", lwd=1, ci.lty=0, ci.type ="poly", ci.col="#00000E50", xlab="n Genomes in Gene Families", ylab="n Gene Families \n", ylim = c(0,50000), las = 1)
plot(df_sans_meta_accy_curve_host,  col="#800080", lwd=1, ci.lty=0,ci.type ="poly", ci.col="#80008E50", add = T)
plot(df_sans_meta_accy_curve_free,  col="#008B8B", lwd=1, ci.lty=0, ci.type ="poly", ci.col="#008B8B50", add = T)
legend(1, 50000, legend=c("Genus-wide pangenome", "Host-associated pangenome", "Free-living pangenome"), col=c("#00000E", "#800080", "#008B8B"), lwd=5, lty=1, cex=0.6, bty = "n")
accumulation_curve_pan_pot <- recordPlot()
dev.off()

# accumulation curves of accessory genes
par(mar=c(6,6,1,1))
plot(df_sans_meta_accy_sans_uniq_curve, col="#00000E", lwd=1, ci.lty=0, ci.type ="poly", ci.col="#00000E50", xlab="n Genomes in Gene Families", ylab="n Gene Families \n", ylim = c(0,35000), las = 1)
plot(df_sans_meta_accy_sans_uniq_curve_host,  col="#800080", lwd=1, ci.lty=0,ci.type ="poly", ci.col="#80008E70", add = T)
plot(df_sans_meta_accy_sans_uniq_curve_free,  col="#008B8B", lwd=1, ci.lty=0, ci.type ="poly", ci.col="#008B8B70", add = T)
legend(1, 35000, legend=c("Genus-wide pangenome", "Host-associated pangenome", "Free-living pangenome"), col=c("#00000E", "#800080", "#008B8B"), lwd=5, lty=1, cex=0.6, bty = "n")
accumulation_curve_acc_plot <- recordPlot()
dev.off()

# plot figure aa
final_pie_plot <- ggarrange(pie_chart_all, pie_chart_host, pie_chart_free, nrow = 1, ncol = 3)
ggsave("figures/pirate/2024-10-10-pie-charts-core-accessory-unique-genes.svg", width = 10, height =5, units = "in")
print(final_pie_plot)
dev.off()

ggsave("figures/pirate/2024-10-10-gene-distribution-plot.svg", width = 10, height =5, units = "in")
print(df_gene_distribution_plot)
dev.off()

final_curves <- ggarrange(accumulation_curve_pan_pot, NULL, accumulation_curve_acc_plot, nrow = 3, heights  = c(1,0.2,1)) 
ggsave("figures/pirate/2024-10-10-gene-accumulation-curves.svg", width = 5, height = 10, units = "in")
print(final_curves)
dev.off()

# Chi-squared test 
# Number of Core, Accessory, and Unique Genes in Host and Free Living Isolates
M <- as.table(rbind(c(1723, 20992, 13002), c(1758, 15394, 11749)))
dimnames(M) <- list(lifestyle = c("Host", "Free"),
                    freq = c("Core","Accessory", "Unique"))
Xsq <- chisq.test(M)  # Prints test summary


