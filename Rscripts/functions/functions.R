create_matrix <- function(data, meta) {
  rownames(data) <- data$allele_name 
  df_matrix <- ifelse(data[21:256]=="", 0, 1)
  df_matrix_t <- as.data.frame(t(df_matrix))
  
  df_matrix_t_meta <- merge(x = df_matrix_t, y = meta, by = 0)
  
  rownames(df_matrix_t_meta) <- df_matrix_t_meta$Row.names
  df_matrix_t_meta <- subset(df_matrix_t_meta, select = -c(Row.names))
  
}

df_gene_families <- function(data, n_isolates) {
  n_core_genes <- sum(data$count > (n_isolates * .95))
  n_accessory_genes <- sum(data$count < (n_isolates * .95) & data$count > 1)
  n_unique_genes <- sum(data$count == 1)
  n_core_accessory_unique <- data.frame(c(n_core_genes, n_accessory_genes, n_unique_genes), c("Core Genes", "Accessory Genes", "Unique Genes"))
  names(n_core_accessory_unique) <- c('count', 'group')
  return(n_core_accessory_unique)
}

distribution_gene_families <- function(df, n_genes) {
  gene_distribution <- data.frame(matrix(ncol = 2, nrow = n_genes))
  colnames(gene_distribution) <- c('gene_type', 'count')
  
  for (i in 1:n_genes){
    gene_count <- sum(df[, i])
    if (gene_count == 1) {
      gene_distribution$gene_type[i] = "Singleton"
      gene_distribution$count[i] = gene_count
    } else if (gene_count >= nrow(df) * .95) {
      gene_distribution$gene_type[i] = "Core"
      gene_distribution$count[i] = gene_count
    } else {
      gene_distribution$gene_type[i] = "Accessory"
      gene_distribution$count[i] = gene_count
    }
  }
  
  gene_distribution <- subset(gene_distribution, count > 0)
  
  return(gene_distribution)
}

plot_distribution_gene_families <- function(df) {
  p <- df %>%
    ggplot(aes(x = count)) +
    geom_bar(aes(color = gene_type, fill = gene_type),
             position = "identity") +
#    ggtitle("Pangenome Distribution by Gene Family") +
#    ylab("n Gene Families") + xlab("n Genomes in Family") +
    theme_bw() +
    theme(
      plot.title = element_text(size=15), legend.title = element_blank()) +
    scale_color_manual(values = c("#BA9141", "#806633", "#316A6E")) +
    scale_fill_manual(values = c("#BA9141", "#806633", "#316A6E")) 
  p
  
}

plot_piechart <- function(data) {
  
  df <- data
  df$fraction <- df$count / sum(df$count)
  df$ymax <- cumsum(df$fraction) 
  df$ymin <- c(0, head(df$ymax, n=-1))
  df$labelPosition <- (df$ymax + df$ymin) / 2
  df$label <- paste0(df$group, "\n", df$count)
  
  plot <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
    geom_rect() +
    geom_text(x=1.25, aes(y=labelPosition, label=label, color=group), size=3) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#BA9141", "#806633", "#316A6E")) +
    scale_color_manual(values = c("#BA9141", "#806633", "#316A6E"))
  
  return(plot)
  
}
