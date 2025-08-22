# import fastani table 
fastani <- read.csv("results/", sep = ",", header = T)
fastani <- fastani[, 1:3]

# convert from long to wide format
fastani_wide <- cast(fastani, reference~querry)
rownames(fastani_wide) <- fastani_wide[,1]
fastani_wide[,1] <- NULL

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

fastani_reorder <- reorder_cormat(as.matrix(fastani_wide))
write.table(fastani_reorder, file='tables/2024-10-10-fastani-pairwise-comparisons.txt',
            quote=FALSE, sep='\t', col.names = NA)

fastani_long <- melt(fastani_reorder, na.rm = T)

# Get upper triangle of the correlation matrix
fastani_heatmap <- ggplot(data = fastani_long, aes(x=X1, y=X2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", high = "darkblue", mid = "#8080C5", midpoint = 0.90, 
                       limits=c(0.75,1), name = "fastANI \n") +
  theme_minimal() + 
  ylab("Reference \n") + xlab("\n Query") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 4, hjust = 1),
        axis.text.y = element_text(size = 4), 
        legend.key.height = unit(1, 'in'), legend.text = element_text(size = 12), legend.title = element_text(size = 20, face = "bold")) +
  coord_fixed()

fastani_heatmap
png("figures/drep/2024-10-10-ani-heatmap.png", res = 300, width = 15, height = 15, units = "in")
print(fastani_heatmap)
dev.off()

