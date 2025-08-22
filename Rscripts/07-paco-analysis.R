library(TreeDist)
library(paco)
library(ape)
library(treeio)

#####################################
##### PACo and GRF tests for Animalia
#####################################
tree.H <- unroot(read.tree("results/raxml/host/Animalia_phylum.nwk"))
tree.B <- unroot(read.tree("results/raxml/pseudoalteromonas/pangenome/RAxML_bestTree.pseudoalteromonas-pangenome"))
am_matrix <- read.delim("results/raxml/associations/animalia-ass.txt", header = F)
am_matrix_links <- as.matrix(read.delim("results/raxml/associations/animalia-ass-list.txt", header = T, row.names = 1))

tree_trim <- keep.tip(tree.H, c("Annelida", "Arthropoda", "Chordata", "Cnidaria", "Ctenophora", "Echinodermata",
                                      "Mollusca", "Nematoda-1", "Porifera"))

rename_host <- data.frame(old = c("Annelida", "Arthropoda", "Chordata", "Cnidaria", 
                                  "Ctenophora", "Echinodermata", "Mollusca", "Nematoda-1", "Porifera"),
                          new= c("Annelida", "Arthropoda", "Chordata", "Cnidaria", 
                                 "Ctenophora", "Echinodermata", "Mollusca", "Nematoda", "Porifera"))

# Generalized Robinson-Foulds 
tree.BR <- rename_taxa(tree.B, am_matrix, key = 1, value = 2)
tree.HR <- rename_taxa(tree_trim, rename_host, key = 1, value = 2)

tree.BR <- drop.tip(tree.BR, c("Seawater"))

distance <- TreeDistance(tree.HR, tree.BR)
distance ## GRF = 0.88
VisualizeMatching(PhylogeneticInfoDistance, tree.HR, tree.BR)

# PACo analysis
gtree <- cophenetic(tree.HR)
ltree <- cophenetic(tree.B)

am_matrix_links <- subset(am_matrix_links, select = -c(Seawater))

D <- prepare_paco_data(H = gtree, P = ltree, HP = am_matrix_links)
D <- add_pcoord(D, correction = 'cailliez')

D <- PACo(D, nperm = 1000, seed = 10, method = "quasiswap", symmetric = T)
D <- paco_links(D)
res <- residuals_paco(D$proc)

# high cophylogenetic signal despite weak taxonomic relationship
# aka high host-switching
D$gof # p-value = 0; m2 = 0.93; r2 = 0.07 high cophylogenetic signal 


#####################################
##### PACo and GRF tests for Mollusca
#####################################
tree.H <- unroot(read.tree("results/raxml/host/mollusca-host-genus.nwk"))
tree.B <- unroot(read.tree("results/raxml/pseudoalteromonas/mollusca/RAxML_bestTree.pseudoalteromonas-mollusca"))
am_matrix <- read.delim("results/raxml/associations/mollusca-ass.txt", header = F)
am_matrix_links <- as.matrix(read.delim("results/raxml/associations/mollusca-ass_link.txt", header = T, row.names = 1))

# Generalized Robinson-Foulds 

tree.BR <- rename_taxa(tree.B, am_matrix, key = 1, value = 2)
distance <- TreeDistance(tree.H, tree.BR)
distance ## GRF = 0.61
VisualizeMatching(PhylogeneticInfoDistance, tree.H, tree.BR)

# PACo analysis
gtree <- cophenetic(tree.H)
ltree <- cophenetic(tree.B)

D <- prepare_paco_data(H = gtree, P = ltree, HP = am_matrix_links)
D <- add_pcoord(D, correction = 'cailliez')

D <- PACo(D, nperm = 1000, seed = 10, method = "quasiswap", symmetric = T)
D <- paco_links(D)
res <- residuals_paco(D$proc)

# high cophylogenetic signal despite weak taxonomic relationship?? 
# aka high host-switching
D$gof # p-value = 0.013; m2 = 0.60; r2 = 0.40 high cophylogenetic signal 


#####################################
##### PACo and GRF tests for Mollusca - Bivalvia
#####################################
tree.H <- unroot(read.tree("results/raxml/host/bivalvia-list.nwk"))
tree.B <- unroot(read.tree("results/raxml/pseudoalteromonas/mollusca/RAxML_bestTree.pseudoalteromonas-bivalvia"))
am_matrix <- read.delim("results/raxml/associations/bivalvia-ass.txt", header = F)
am_matrix_links <- as.matrix(read.delim("results/raxml/associations/bivalvia-ass_link.txt", sep = "\t", header = T, row.names = 1))

# Generalized Robinson-Foulds 
tree.BR <- rename_taxa(tree.B, am_matrix, key = 1, value = 2)
distance <- TreeDistance(tree.H, tree.B)
distance ## GRF = 0.75
VisualizeMatching(PhylogeneticInfoDistance, tree.H, tree.BR)

# PACo analysis
gtree <- cophenetic(tree.H)
ltree <- cophenetic(tree.B)

D <- prepare_paco_data(H = gtree, P = ltree, HP = am_matrix_links)
D <- add_pcoord(D, correction = 'cailliez')

D <- PACo(D, nperm = 100, seed = 10, method = "swap", symmetric = T)
D <- paco_links(D)
res <- residuals_paco(D$proc)

# high cophylogenetic signal despite weak taxonomic relationship?? 
# aka high host-switching
D$gof # p-value = 0.03; m2 = 0.66; r2 = 0.34 high cophylogenetic signal 


#####################################
##### PACo and GRF tests for Porifera
#####################################

tree.H <- unroot(read.tree("results/raxml/host/Porifera_class.nwk"))
tree.B <- unroot(read.tree("results/raxml/pseudoalteromonas/porifera/RAxML_bestTree.pseudoalteromonas-porifera"))
am_matrix <- read.delim("results/raxml/associations/porifera-ass.txt", header = F)
am_matrix_links <- as.matrix(read.delim("results/raxml/associations/porifera-ass_link.txt", header = T, row.names = 1))

# Generalized Robinson-Foulds 
tree.BR <- rename_taxa(tree.B, am_matrix, key = 1, value = 2)
distance <- TreeDistance(tree.H, tree.BR)
distance ## GRF = 0.93
VisualizeMatching(MutualClusteringInformation, tree.H, tree.BR)

# PACo
gtree <- cophenetic(tree.H)
ltree <- cophenetic(tree.B)

D <- prepare_paco_data(H = gtree, P = ltree, HP = am_matrix_links)
D <- add_pcoord(D, correction = 'cailliez')

D <- PACo(D, nperm = 1000, seed = 10, method = "quasiswap", symmetric = T)
D <- paco_links(D)
res <- residuals_paco(D$proc)

# high cophylogenetic signal despite weak taxonomic relationship?? 
# aka high host-switching
D$gof # p-value = 0.4702; m2 = 0.86; r2 = 0.14 no cophylogenetic signal 


#####################################
##### PACo and GRF tests for Cnidarians
#####################################

tree.H <- unroot(read.tree("results/raxml/host/scleractinia-list.nwk"))
tree.B <- unroot(read.tree("results/raxml/pseudoalteromonas/cnidaria/RAxML_bestTree.pseudoalteromonas-scleractinia"))
am_matrix <- read.delim("results/raxml/associations/scleractinia-ass.txt", header = F)
am_matrix_links <- as.matrix(read.delim("results/raxml/associations/scleractinia-ass_list.txt", header = T, row.names = 1))

# Generalized Robinson-Foulds 
tree.BR <- rename_taxa(tree.B, am_matrix, key = 1, value = 2)
distance <- TreeDistance(tree.H, tree.BR)
distance ## GRF = 0.85
VisualizeMatching(MutualClusteringInformation, tree.H, tree.BR)

# PACo
gtree <- cophenetic(tree.H)
ltree <- cophenetic(tree.B)

D <- prepare_paco_data(H = gtree, P = ltree, HP = am_matrix_links)
D <- add_pcoord(D, correction = 'cailliez')

D <- PACo(D, nperm = 1000, seed = 10, method = "quasiswap", symmetric = T)
D <- paco_links(D)
res <- residuals_paco(D$proc)

# high cophylogenetic signal despite weak taxonomic relationship?? 
# aka high host-switching
D$gof # p-value = 0.011; m2 = 0.89; r2 = 0.11 low cophylogenetic signal

#####################################
##### PACo and GRF tests for Cnidarians (Faviidae)
#####################################

tree.H <- unroot(read.tree("results/raxml/host/faviidae-list.nwk"))
tree.B <- unroot(read.tree("results/raxml/pseudoalteromonas/cnidaria/RAxML_bestTree.pseudoalteromonas-faviidae"))
am_matrix <- read.delim("results/raxml/associations/faviidae-ass.txt", header = F)
am_matrix_links <- as.matrix(read.delim("results/raxml/associations/faviidae-ass_link.txt", header = T, row.names = 1))

# Generalized Robinson-Foulds 
tree.BR <- rename_taxa(tree.B, am_matrix, key = 1, value = 2)
distance <- TreeDistance(tree.H, tree.BR)
distance ## GRF = 0.26
VisualizeMatching(MutualClusteringInformation, tree.H, tree.BR)

# PACo
gtree <- cophenetic(tree.H)
ltree <- cophenetic(tree.B)

D <- prepare_paco_data(H = gtree, P = ltree, HP = am_matrix_links)
D <- add_pcoord(D, correction = 'cailliez')

D <- PACo(D, nperm = 1000, seed = 10, method = "quasiswap", symmetric = T)
D <- paco_links(D)
res <- residuals_paco(D$proc)

# high cophylogenetic signal despite weak taxonomic relationship?? 
# aka high host-switching
D$gof # p-value = 0.025; m2 = 0.59; r2 = .41 # high cophylogenetic signal

#####################################
##### PACo and GRF tests for Nematoda
#####################################
library(vegan)

tree.H <- unroot(read.tree("results/raxml/host/RAxML_bestTree.nematoda-metagenomes"))
rpkm.B <- read.table("results/rrap/nematode_metagenomes_rpkm_noLog-replicates-only.csv", sep = ",", header = T, check.names = F)
rpkm.B[rpkm.B < 0.1] <- 0
tree.B <- as.matrix(vegdist(t(rpkm.B[,-1]), method = "bray"))

#tree.H <- as.matrix(dist(t(rpkm[,-1]), method = "bray-curtis", diag = T, upper = T))
am_matrix_links <- as.matrix(read.delim("results/raxml/associations/nematoda-ass_link.txt", header = T, row.names = 1))

# PACo
#gtree <- cophenetic(tree.H)
ltree <- cophenetic(tree.H)

D <- prepare_paco_data(H = tree.B, P = ltree, HP = am_matrix_links)
D <- add_pcoord(D)

D <- PACo(D, nperm = 1000, seed = 10, method = "quasiswap", symmetric = T)
D <- paco_links(D)
res <- residuals_paco(D$proc)

# high cophylogenetic signal despite weak taxonomic relationship?? 
# aka high host-switching
D$gof # p-value = 0.00; m2 = 0.80; r2 = .20 #  cophylogenetic signal


# https://academic.oup.com/sysbio/article/73/3/613/7628010
# https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347(21)00179-8#f0020
