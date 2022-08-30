library(Rcpp)
source("scripts/SinkhornNNLSLinseedC.R")
source("scripts/PreprocessingPlots.R")
sourceCpp("scripts/pipeline.cpp")

filterZeroMAD <- function(dataset) {
  mad_ <- apply(dataset, 1, mad)
  dataset[mad_ > 0,]
}

get_normalized_svd_projections <- function(
    self,
    dims,
    keep_genes = NULL,
    keep_samples = NULL
) {
    if (is.null(keep_genes)) {
        keep_genes <- rownames(self$V_row)
    }
    if (is.null(keep_samples)) {
        keep_samples <- colnames(self$V_row)
    }
    V_row_flt <- self$V_row[
        rownames(self$V_row) %in% keep_genes,
        colnames(self$V_row) %in% keep_samples
    ]
    V_column_flt <- self$V_column[
        rownames(self$V_column) %in% keep_genes,
        colnames(self$V_column) %in% keep_samples
    ]
    svd_ <- svd(V_row_flt)
    S <- t(svd_$u[, dims])
    R <- t(svd_$v[, dims])
    S[1,] <- -S[1,]
    R[1,] <- -R[1,]
    projOmega <- as.data.frame(t(S %*% V_column_flt))[, seq_len(length(dims))]
    colnames(projOmega) <- paste0("dim_", dims)
    projX <- as.data.frame(V_row_flt %*% t(R))[, seq_len(length(dims))]
    colnames(projX) <- paste0("dim_", dims)
    list(projX, projOmega)
}

calculate_knn_cutoff <- function(self, k_neighbours = 20, svd_dims = NULL) {
  if (is.null(svd_dims)) {
      genes <- self$V_row
  } else {
      genes <- get_normalized_svd_projections(self, seq_len(svd_dims))[[1]]
  }
  sort(dbscan::kNNdist(genes, k = k_neighbours))
}

data_ <- readRDS(snakemake@params[["dataset"]])
data_ <- filterZeroMAD(data_)
tmp_snk <- SinkhornNNLSLinseed$new(dataset = snakemake@config[["dataset"]],
                                   path = snakemake@output[[1]],
                                   data = data_,
                                   analysis_name = snakemake@config[["analysis_name"]],
                                   cell_types = snakemake@config[["cell_types"]])
print(tmp_snk$cell_types)
if (!is.null(snakemake@config[["top_mad"]])) {
  mad_limit <- min(tmp_snk$M, snakemake@config[["top_mad"]])
  top_genes <- names(sort(tmp_snk$genes_mad, decreasing = T)[1:mad_limit])
  png(snakemake@output[["top_mad"]])
  print(plotTopMADGenes(tmp_snk, top_genes))
  dev.off()
  tmp_snk$selectTopGenes(mad_limit)
} else {
  min_mad <- 0
  min_median <- 0
  if (!is.null(snakemake@config[["min_mad"]])) {
    min_mad <- snakemake@config[["min_mad"]]
  }
  if (!is.null(snakemake@config[["min_median"]])) {
    min_median <- snakemake@config[["min_median"]]
  }
  print(min_mad)
  print(min_median)
  top_genes1 <- names(tmp_snk$genes_mad[tmp_snk$genes_mad > min_mad])
  top_genes2 <- names(tmp_snk$genes_median[tmp_snk$genes_median >= min_median])
  top_genes <- intersect(top_genes1, top_genes2)
  png(snakemake@output[["top_mad"]])
  print(plotTopMADGenes(tmp_snk, top_genes1))
  dev.off()

  png(snakemake@output[["top_mad_joined"]])
  print(plotTopMADGenes(tmp_snk, top_genes))
  dev.off()

  png(snakemake@output[["mad_med_scatter"]])
  print(plotMADExpMedianScatter(tmp_snk, top_genes))
  dev.off()

  tmp_snk$filterByMADExpMedian(min_mad, min_median)
}
cat(paste0("Genes: ", tmp_snk$M), paste0("Samples: ", tmp_snk$N),
    file = snakemake@output[["metadata"]], sep = "\n")
tmp_snk$scaleDataset(snakemake@config[["scale_iterations"]])
png(snakemake@output[["svd_before"]])
plotSVD(tmp_snk, snakemake@output[["svd_before_plot"]])
dev.off()

if (snakemake@config[["knn_filter"]]) {
  genes_nn_dist <- calculate_knn_cutoff(
    tmp_snk,
    20,
    NULL
  )
  keep_n_genes <- tmp_snk$M - snakemake@config[["filter_genes"]]
  print("KEEP_N_GENES:")
  print(keep_n_genes)
  stopifnot(keep_n_genes > 0)
  threshold <- genes_nn_dist[keep_n_genes]
  png(snakemake@output[["distance_before"]])
  print(plotKNNCutoffGenes(genes_nn_dist, threshold))
  dev.off()
  keep_genes <- names(genes_nn_dist[genes_nn_dist < threshold])
  tmp_snk$filterGenes(keep_genes)
  tmp_snk$calculateDistances()
  genes_nn_dist_after <- calculate_knn_cutoff(
    tmp_snk,
    20,
    NULL
  )
  png(snakemake@output[["distance_after"]])
  print(plotKNNCutoffGenes(genes_nn_dist_after, 0))
  dev.off()
} else {
  tmp_snk$getSvdProjectionsNew()
  tmp_snk$calculateDistances()
  png(snakemake@output[["distance_before"]])
  plotDistances(tmp_snk, snakemake@config[["filter_genes"]], snakemake@config[["filter_samples"]])
  dev.off()
  tmp_snk$filterByDistance(filter_genes = snakemake@config[["filter_genes"]],
                           filter_samples = snakemake@config[["filter_samples"]],
                           iterations = snakemake@config[["scale_iterations"]])
  tmp_snk$calculateDistances()
  png(snakemake@output[["distance_after"]])
  plotDistances(tmp_snk, 0, 0)
  dev.off()
}

png(snakemake@output[["svd_after"]])
plotSVDMerged(tmp_snk, snakemake@output[["svd_before_plot"]])
dev.off()
saveRDS(tmp_snk, file = snakemake@output[[1]])