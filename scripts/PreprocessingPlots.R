library(ggplot2)
library(gridExtra)

plotKNNCutoffGenes <- function(genes_nn_dist, threshold) {
  toPlot <- as.data.frame(genes_nn_dist)
  toPlot$idx <- seq_len(nrow(toPlot))
  ggplot(toPlot, aes(x = idx, y = genes_nn_dist)) + geom_line() + geom_hline(yintercept = threshold, colour = "red")
}

plotTopMADGenes <- function(metadata_, top_genes) {

  toPlot <- data.frame(log2_mad = log2(metadata_$genes_mad + 1),
                       gene = names(metadata_$genes_mad))
  toPlot$idx <- 1:nrow(toPlot)
  toPlot$pass_filter <- F
  toPlot[toPlot$gene %in% top_genes, "pass_filter"] <- T
  ggplot(toPlot, aes(x = log2_mad, color = pass_filter)) +
    geom_histogram(bins = 100, fill = 'white') +
    theme_minimal()
}

plotMADExpMedianScatter <- function(metadata_, top_genes) {
  toPlot <- data.frame(log2_mad = log2(metadata_$genes_mad + 1),
                       log2_median = log2(metadata_$genes_median + 1),
                       gene = names(metadata_$genes_mad))
  toPlot$idx <- 1:nrow(toPlot)
  toPlot$pass_filter <- F
  toPlot[toPlot$gene %in% top_genes, "pass_filter"] <- T
  ggplot(toPlot, aes(x = log2_mad, y = log2_median, color = pass_filter)) +
    geom_point() +
    theme_minimal()
}

plotDistances <- function(metadata_, filter_genes = 0, filter_samples = 0) {
  toPlot_Genes <- data.frame(Distance = metadata_$distance_genes)
  rownames(toPlot_Genes) <- names(metadata_$distance_genes)
  toPlot_Genes$idx <- 1:nrow(toPlot_Genes)

  toPlot_Samples <- data.frame(Distance = metadata_$distance_samples)
  rownames(toPlot_Samples) <- names(metadata_$distance_samples)
  toPlot_Samples$idx <- 1:nrow(toPlot_Samples)

  pltGenes <- ggplot(toPlot_Genes, aes(x = idx, y = Distance)) +
    geom_point(size = 0.1) +
    geom_line() +
    theme_minimal()

  if (filter_genes > 0) {
    pltGenes <- pltGenes + geom_vline(xintercept = filter_genes,
                                      color = "red", linetype = "dashed")
  }

  pltSamples <- ggplot(toPlot_Samples, aes(x = idx, y = Distance)) +
    geom_point(size = 0.1) +
    geom_line() +
    theme_minimal()

  if (filter_samples > 0) {
    pltSamples <- pltSamples + geom_vline(xintercept = filter_samples,
                                          color = "red", linetype = "dashed")
  }

  grid.arrange(pltGenes, pltSamples)
}

plotMergedDistances <- function(metadata_, filter_genes = 0, filter_samples = 0) {
  toPlot_Genes <- data.frame(Distance = metadata_$merged_distance_genes)
  rownames(toPlot_Genes) <- names(metadata_$merged_distance_genes)
  toPlot_Genes$idx <- 1:nrow(toPlot_Genes)

  toPlot_Samples <- data.frame(Distance = metadata_$merged_distance_samples)
  rownames(toPlot_Samples) <- names(metadata_$merged_distance_samples)
  toPlot_Samples$idx <- 1:nrow(toPlot_Samples)

  pltGenes <- ggplot(toPlot_Genes, aes(x = idx, y = Distance)) +
    geom_point(size = 0.1) +
    geom_line() +
    theme_minimal()

  if (filter_genes > 0) {
    pltGenes <- pltGenes + geom_vline(xintercept = filter_genes,
                                      color = "red", linetype = "dashed")
  }

  pltSamples <- ggplot(toPlot_Samples, aes(x = idx, y = Distance)) +
    geom_point(size = 0.1) +
    geom_line() +
    theme_minimal()

  if (filter_samples > 0) {
    pltSamples <- pltSamples + geom_vline(xintercept = filter_samples,
                                          color = "red", linetype = "dashed")
  }

  grid.arrange(pltGenes, pltSamples)
}

plotSVD <- function(metadata_, save_rds, components = 50) {
  svd_ <- svd(metadata_$V_row)
  vars <- svd_$d^2
  vars <- cumsum(vars / sum(vars))
  df <- data.frame(dimension = 1:length(vars), variance = vars, filtered = FALSE)
  saveRDS(df, save_rds)
  colors <- 1
  ggplot(data = df, aes(x = dimension, y = variance)) +
    geom_point(aes(y = variance, color = filtered), size = 0.5, alpha = 1) +
    geom_line(aes(y = variance, color = filtered, group = filtered),
              size = 0.5, alpha = 1) +
    scale_color_manual(values = c("#999999", "#E41A1C")[1:colors]) +
    theme_minimal(base_size = 8) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
          legend.position = "none") +
    scale_x_continuous(minor_breaks = 1:components,
                       limits = c(0, components))
}

plotSVDMerged <- function(metadata_, svd_before_f, components = 50) {
  svd_ <- svd(metadata_$V_row)
  vars <- svd_$d^2
  vars <- cumsum(vars / sum(vars))
  df <- data.frame(dimension = 1:length(vars), variance = vars, filtered = FALSE)
  df$processing <- "after"
  svd_before <- readRDS(svd_before_f)
  svd_before$processing <- "before"
  toPlot <- rbind(svd_before, df)
  toPlot$processing <- factor(toPlot$processing, levels = c("before", "after"))
  colors <- 2
  ggplot(data = toPlot, aes(x = dimension, y = variance)) +
    geom_point(aes(y = variance, color = processing), size = 0.5, alpha = 1) +
    geom_line(aes(y = variance, color = processing, group = processing),
              size = 0.5, alpha = 1) +
    scale_color_manual(values = c("#999999", "#E41A1C")[1:colors]) +
    theme_minimal(base_size = 8) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid')) +
    scale_x_continuous(minor_breaks = 1:components,
                       limits = c(0, components))
}
