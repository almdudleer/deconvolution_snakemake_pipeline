library(patchwork)

filterZeroMAD <- function(dataset) {
  mad_ <- apply(dataset, 1, mad)
  dataset[mad_ > 0,]
}

calculate_knn_distances <- function(self, k_neighbours = 20, svd_dims = NULL, samples = FALSE) {
  if (is.null(svd_dims)) {
    if (!samples) {
      data <- self$V_row
    } else {
      data <- t(self$V_row)
    }
  } else {
    proj <- get_normalized_svd_projections(self, seq_len(svd_dims))
    if (!samples) {
      data <- proj[[1]]
    } else {
      data <- proj[[2]]
    }
  }
  sort(dbscan::kNNdist(data, k = k_neighbours))
}

calculate_zero_distances <- function(data) {
  sqrt(rowSums(data^2))
}

plotDistancesToZero <- function(
  annotations,
  col = NULL,
  thresholds = list(
    zero_distance = NULL,
    knn_distance = NULL,
    plane_distance = NULL
  ),
  bins = 200,
  highlight_items = NULL
) {
  columns <- c("knn_distance", "zero_distance", "plane_distance")
  if (!is.null(highlight_items)) {
    annotations$highlight <- highlight_items
    col <- "highlight"
  }

  distr_plots <- lapply(columns, function(colname) {
    ggplot(
      annotations,
      aes_string(x = colname, fill = col)
    ) +
      geom_histogram(bins = bins, aes(y = (..count..)/sum(..count..))) +
      geom_vline(xintercept=thresholds[[colname]], colour = "red") +
      ylab("frequency")
  })

  sorted_plots <- lapply(columns, function(colname) {
    ggplot(
      annotations[order(annotations[, colname]),],
      aes_string(x = 1:nrow(annotations), y = colname, col = col)
    ) +
      geom_hline(yintercept=thresholds[[colname]], colour = "red") +
      geom_point(size = 0.1) +
      xlab("order")
  })

  if (!is.null(col)) {
    arranged_annotations <- annotations[order(annotations[, col]),]
  } else {
    arranged_annotations <- annotations
  }

  vs_plots <- apply(combn(columns, 2), 2, function(colnames) {
    ggplot(
      arranged_annotations,
      aes_string(x = colnames[[1]], y = colnames[[2]], col = col)
    ) +
      geom_vline(xintercept=thresholds[[colnames[[1]]]], colour = "red") +
      geom_hline(yintercept=thresholds[[colnames[[2]]]], colour = "red") +
      geom_point(size = 0.1)
  })

  wrap_plots(c(distr_plots, sorted_plots, vs_plots), nrow = 3, ncol = 3, guides = "collect")
}

createGeneAnnotations <- function(data, sim = F) {
  annotations <- data.frame(gene_name = rownames(data))
  annotations$RPLS <- grepl('^RPL|^RPS', annotations$gene_name)
  genes_annotation_lists <- list(
    HK = "datasets/gene_annotation_lists/HK.txt",
    CC = "datasets/gene_annotation_lists/CC.txt",
    CODING = "datasets/gene_annotation_lists/CODING_ENSEMBL.txt"
  )
  for (name in names(genes_annotation_lists)) {
    path <- genes_annotation_lists[[name]]
    genes_string <- readChar(path, file.info(path)$size)
    genes_list <- stringr::str_split(genes_string, " ")[[1]]
    annotations[name] <- annotations$gene_name %in% genes_list
  }
  if (sim) {
    annotations$CODING <- TRUE
  }
  annotations$mad <- apply(data, 1, mad)
  annotations$mad_gt_2 <- annotations$mad > 2
  annotations$mad_gt_1 <- annotations$mad > 1
  annotations$mad_gt_0 <- annotations$mad > 0
  annotations
}

init_tmp_lo2 <- function(data, cell_types, dataset_name) {
  lo2 <- SinkhornNNLSLinseed$new(dataset = dataset_name,
                                 path = paste0("datasets/", dataset_name, ".rds"),
                                 data = data,
                                 analysis_name = "TEST",
                                 cell_types = cell_types)
  lo2$scaleDataset(30)
  lo2
}

# TODO: replace with list, then, with ExpressionSet or our own object
calculate_distances <- function(data, gene_annotations, sample_annotations, dataset_name, cell_types, svd_dims = 20) {
  new_gene_annotations <- gene_annotations
  new_sample_annotations <- sample_annotations

  lo2 <- init_tmp_lo2(data, cell_types, dataset_name)
  lo2$getSvdProjectionsNew()
  lo2$calculateDistances()

  new_gene_annotations$plane_distance <- lo2$distance_genes[new_gene_annotations$gene_name]
  new_sample_annotations$plane_distance <- lo2$distance_samples[new_sample_annotations$sample_name]

  proj <- get_normalized_svd_projections(lo2, 1:cell_types)
  projected_genes <- proj$projX[, 2:cell_types]
  projected_samples <- proj$projOmega[, 2:cell_types]

  new_gene_annotations$knn_distance <- calculate_knn_distances(lo2, 20, svd_dims)[new_gene_annotations$gene_name]
  new_gene_annotations$zero_distance <- calculate_zero_distances(projected_genes)[new_gene_annotations$gene_name]

  new_sample_annotations$knn_distance <- calculate_knn_distances(
    lo2, 20, svd_dims, samples = TRUE
  )[new_sample_annotations$sample_name]
  new_sample_annotations$zero_distance <- calculate_zero_distances(
    projected_samples
  )[new_sample_annotations$sample_name]

  list(
    gene_annotations = new_gene_annotations,
    sample_annotations = new_sample_annotations
  )
}

knn_cutoff <- function(data, annotations, threshold, samples = F) {
  new_annotations <- annotations[annotations$knn_distance < threshold,]
  if (samples) {
    new_data <- data[,new_annotations$sample_name]
  } else {
    new_data <- data[new_annotations$gene_name,]
  }
  list(
    data = new_data,
    annotations = new_annotations
  )
}

plotMadHistograms <- function(gene_annotations) {
  plt1 <- ggplot(gene_annotations, aes(x = log2(mad + 1), color = mad_gt_0)) +
    geom_histogram(bins = 100, fill = 'white') +
    theme_minimal()
  
  plt2 <- gene_annotations %>%
    filter_str(flt) %>%
    ggplot(aes(x = log2(mad + 1), color = as.factor(mad_gt_1 + mad_gt_2))) +
    geom_histogram(bins = 100, fill = 'white') +
    theme_minimal()
  
  grid.arrange(plt1, plt2)
}

filter_str <- function(data, expr, ...) {
  filter(data, eval(rlang::parse_expr(expr)), ...)
}

fit_lo2 <- function(data, cell_types, dataset_name) {
  lo2 <- init_tmp_lo2(data, cell_types, dataset_name)
  lo2$getSvdProjectionsNew()
  lo2$selectInitX()
  run_block(lo2, iterations = 15000)
  show(plotErrors(lo2))
  show(plotPoints2D(lo2))
  show(plotSumToOneChange(lo2))
  show(plotNegProportionsChange(lo2))
  show(plotNegBasisChange(lo2))
  show(plotAbundance(lo2))
  lo2
}

get_basis_column <- function(lo2) {
  ct_names <- paste0("Cell_type_", 1:lo2$cell_types)
  colnames(lo2$full_basis) <- ct_names
  rownames(lo2$full_proportions) <- ct_names
  toSave <- t(t(lo2$full_basis) / rowSums(t(lo2$full_basis)))
  toSave <- lo2$getFoldChange(toSave)
  toSave <- rbind(c(rep(NA,lo2$cell_types), round(apply(lo2$full_proportions,1,mean),4)),toSave)
  rownames(toSave) <- c("avg_proportions",rownames(lo2$filtered_dataset))
  toSave
}

get_cell_type_markers <- function(lo2, top_genes) {
  cell_types <- lo2$cell_types
  basis_column <- get_basis_column(lo2)
  basis_ <- data.frame(basis_column)
  basis_$gene <- rownames(basis_)
  basis_ <- basis_[-1,]
  genesList <- data.frame(genes=rep(" ", cell_types))
  rownames(genesList) <- paste0("Cell_type_", 1:cell_types)

  for (ct in 1:cell_types) {
    row_name <- paste0("Cell_type_", ct)
    genes <- (basis_ %>% arrange(desc(.data[[row_name]])))[1:top_genes, "gene"]
    genesList[row_name, "genes"] <- paste0(genes, collapse = " ")
  }
  genesList
}