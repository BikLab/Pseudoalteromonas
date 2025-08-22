library("phytools")
library("phangorn")

## for cnidaria
# load trees
tree.H <- read.tree("results/raxml/host/scleractinia-host-genus.nwk")
tree.S <- read.tree("results/raxml/pseudoalteromonas/mollusca/raxml-10000boot/RAxML_bestTree.pseudoalteromonas-bivalvia")
am_matrix <- read.delim("results/raxml/associations/bivalvia-ass_link.txt", header = F)

tree.SR <- midpoint.root(unroot(tree.S))

# rotate nodes to optimize matching
obj <- cophylo(tree.SR, tree.H, assoc = am_matrix, rotate = T, tangle = "both")

# get list of colors 
list_color_palettes = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, list_color_palettes$maxcolors, rownames(list_color_palettes)))

# assign one color for each phylogroup and keep the same theme throughout the analysis
colors.variable <- data.frame(
  phylogroup = unique(obj$assoc$V2),
  mycolors = color_vector[1:16]
)

png(filename = "figures/cophylogeny/scleractinia-cophylogeny-2.png",  width = 11, 
    height = 8.5, units = "in", res = 300)

plot(obj,type="phylogram",fsize=0.75 ,part=0.4, tip.len=0,
     link.type="curved",link.lty="solid",link.lwd=3,
     link.col=make.transparent("black",0.5), pts=1, ftype="off")

tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[1]])),"1"),
                  cex=0.3,piecol="darkred")

tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[2]])),"1"),
                  cex=0.3,piecol="black",which="right")

dev.off()

## for cnidaria
tree.H <- read.tree("results/raxml/host/cnidaria-faviidae-host-genus.nwk")
tree.S <- read.tree("results/raxml/pseudoalteromonas/cnidaria/RAxML_bestTree.pseudoalteromonas-faviidae")
am_matrix <- read.delim("results/raxml/associations/faviidae-ass.txt", header = F)

tree.SR <- midpoint.root(unroot(tree.S))

obj <- cophylo(tree.SR, tree.H, assoc = am_matrix, rotate = T, tangle = "left")

png(filename = "figures/cophylogeny/faviidae-cophylogeny-2.png",  width = 11, 
    height = 8.5, units = "in", res = 300)

plot(obj,type="phylogram",fsize=0.75 ,part=0.4, tip.len=0,
     link.type="curved",link.lty="solid",link.lwd=3,
     link.col=make.transparent("black",0.5), pts=1)
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[1]])),"1"),
                  cex=0.3,piecol="darkred")
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[2]])),"1"),
                  cex=0.3,piecol="black",which="right")

dev.off()

## for bivalvia
tree.H <- read.tree("results/raxml/hosts/bivalvia-host-genus.nwk")
tree.S <- read.tree("results/raxml/pseudoalteromonas/mollusca/RAxML_bestTree.pseudoalteromonas-bivalvia")
am_matrix <- read.delim("results/raxml/associations/bivalvia-ass.txt", header = F)

tree.SR <- midpoint.root(unroot(tree.S))

obj <- cophylo(tree.SR, tree.H, assoc = am_matrix, rotate = T, tangle = "left")

png(filename = "figures/cophylogeny/bivalvia-cophylogeny.png",  width = 11, 
    height = 8.5, units = "in", res = 300)

plot(obj,type="phylogram",fsize=0.75 ,part=0.4, tip.len=0,
     link.type="curved",link.lty="solid",link.lwd=3,
     link.col=make.transparent("black",0.5), pts=1, ftype="off")
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[1]])),"1"),
                  cex=0.3,piecol="darkred")
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[2]])),"1"),
                  cex=0.3,piecol="black",which="right")

dev.off()

## for mollusca
tree.H <- read.tree("results/raxml/hosts/mollusca-host-genus.nwk")
tree.S <- read.tree("results/raxml/pseudoalteromonas/mollusca/RAxML_bestTree.pseudoalteromonas-mollusca")
am_matrix <- read.delim("results/raxml/associations/mollusca-ass.txt", header = F)

tree.SR <- midpoint.root(unroot(tree.S))

obj <- cophylo(tree.SR, tree.H, assoc = am_matrix, rotate = T, tangle = "both")

png(filename = "figures/cophylogeny/mollusca-cophylogeny-2.png",  width = 11, 
    height = 8.5, units = "in", res = 300)

plot(obj,type="phylogram",fsize=0.75 ,part=0.4, tip.len=0,
     link.type="curved",link.lty="solid",link.lwd=3,
     link.col=make.transparent("black",0.5), pts=1)
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[1]])),"1"),
                  cex=0.3,piecol="darkred")
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[2]])),"1"),
                  cex=0.3,piecol="black",which="right")

dev.off()

## for porifera
tree.H <- read.tree("results/raxml/hosts/Porifera_class.nwk")
tree.S <- read.tree("/results/raxml/pseudoalteromonas/")
am_matrix <- read.delim("results/raxml/hosts/porifera-host-bacterial-association.txt", header = F)
tree.SR <- midpoint.root(unroot(tree.S))

obj <- cophylo(tree.SR, tree.H, assoc = am_matrix, rotate = F, tangle = "both")
cophyloplot(tree.SR, tree.H, assoc = am_matrix,
            length.line = 4, space = 28, gap = 3)
png(filename = "figures/cophylogeny/porifera-cophylogeny-2.png",  width = 11, 
    height = 8.5, units = "in", res = 300)

plot(obj,type="phylogram",fsize=0.75 ,part=0.4, tip.len=0,
     link.type="curved",link.lty="solid",link.lwd=3,
     link.col=make.transparent("black",0.5), pts=1, ftype = "off")
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[1]])),"1"),
                  cex=0.3,piecol="darkred")
tiplabels.cophylo(pie=to.matrix(rep(1,Ntip(obj$trees[[2]])),"1"),
                  cex=0.3,piecol="black",which="right")

dev.off()
