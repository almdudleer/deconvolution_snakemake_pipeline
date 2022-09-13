library(dplyr)

init_lo2 <- function(
    data,
    top_genes,
    cell_types = 3,
    seed = 42,
    plot = TRUE
) {
    lo2 <- SinkhornNNLSLinseed$new(
        dataset = data,
        path = "../data/linseed2_results",
        analysis_name = "deconv_comparison_2",
        cell_types = cell_types,
    )
    lo2$selectTopGenes(top_genes)
    lo2$scaleDataset(iterations = 20)
    lo2$getSvdProjectionsNew(k = cell_types)
    lo2$calculateDistances()
    lo2$selectInitRandom(seed = seed)
    if (plot) {
        plotPointsColored(lo2, "init")
        # do.call("grid.arrange", c(lo2$plotDistances(arrange = F), ncol = 2))
    }
    lo2
}

guess_borders_unit <- function(borders) {
    unit <- "num"
    if (borders[1] < 0) {
        unit <- "mad"
    } else if (borders[2] < 1) {
        unit <- "quantile"
    }
    unit
}

cast_borders_unit <- function(data, borders, unit = "mad") {
    if (unit == "mad") {
        lo <- max(median(data) + borders[1] * mad(data), 0)
        hi <- median(data) + borders[2] * mad(data)
    } else if (unit == "num") {
        sorted <- sort(data)
        lo <- sorted[borders[1]]
        hi <- sorted[borders[2]]
    } else if (unit == "quantile") {
        lohi <- quantile(data,  borders)
        lo <- lohi[[1]]
        hi <- lohi[[2]]
    }
    c(lo, hi)
}

cast_borders <- function(data, borders) {
    cast_borders_unit(data, borders, guess_borders_unit(borders))
}

plotPointsColored <- function(
    self, 
    points="init",
    mid_alpha = 0.5,
    genes_borders = c(-2, 2),
    samples_borders = c(-2, 2),
    highlight_samples = c(),
    highlight_genes = c(),
    dims=3,
    arrange=TRUE
) {
    if (is.null(self$distance_samples)) {
      self$calculateDistances()
    }

    if (points == "init") {
      X <- self$init_X
      Omega <- self$init_Omega
      count_neg_props <- self$init_count_neg_props
      count_neg_basis <- self$init_count_neg_basis
    }
    if (points == "current") {
      X <- self$X
      Omega <- self$Omega
      count_neg_props <- self$count_neg_props
      count_neg_basis <- self$count_neg_basis
    }

    X <- X[,1:dims]
    Omega <- Omega[1:dims,]

    ## plot X
    toPlotX <- as.data.frame(self$V_row %*% t(self$R))[,1:dims]
    colnames(toPlotX) <- c("X","Y","Z")
    colnames(X) <- c("X","Y","Z")
    genes_borders <- cast_borders(self$distance_genes, genes_borders)
    pltX <- toPlotX %>% 
      mutate(
          distance = self$distance_genes[rownames(toPlotX)],
          highlight = rownames(toPlotX) %in% highlight_genes
      ) %>%
      mutate(
          col = case_when(
            highlight ~ "green",
            distance < genes_borders[1] ~ "blue",
            distance > genes_borders[2] ~ "red",
            T ~ alpha("grey", mid_alpha)
          )
      ) %>%
      ggplot(aes(x=Y, y=Z, col=factor(col))) +
        geom_point() + 
        geom_polygon(data=as.data.frame(X), fill=NA, color = "green") + 
        theme_minimal() + 
        scale_color_identity() + 
        ggtitle(paste(nrow(toPlotX), "genes"))
    if (!is.null(count_neg_props)) {
      pltX <- pltX + annotate(
          "text",
          x = Inf,
          y = Inf,
          label = paste0(round(count_neg_props / (self$cell_types*self$N),4)*100,"%"),
          vjust=1,
          hjust=1
      )
    }

    ## plot Omega  
    toPlotOmega <- as.data.frame(t(self$S %*% self$V_column))[,1:dims]
    colnames(toPlotOmega) <- c("X","Y","Z")
    rownames(Omega) <- c("X","Y","Z")
    samples_borders <- cast_borders(self$distance_samples, samples_borders)
    pltOmega <- toPlotOmega %>% 
      mutate(
          distance = self$distance_samples[rownames(toPlotOmega)],
          highlight = rownames(toPlotOmega) %in% highlight_samples
      ) %>%
      mutate(
          col = case_when(
            highlight ~ "green",
            distance < samples_borders[1] ~ "blue",
            distance > samples_borders[2] ~ "red",
            T ~ alpha("grey", mid_alpha)
          )
      ) %>% 
      ggplot(aes(x=Y, y=Z, col=col)) +
        geom_point() + 
        geom_polygon(data=as.data.frame(t(Omega)), fill=NA, color = "green") +
        theme_minimal() +
        scale_color_identity() +
        ggtitle(paste(nrow(toPlotOmega), "samples"))
    if (!is.null(count_neg_basis)) {
      pltOmega <- pltOmega + annotate(
          "text",
          x = Inf,
          y = Inf,
          label = paste0(round(count_neg_basis / (self$cell_types*self$M),4)*100,"%"),
          vjust=1,
          hjust=1
      )
    }
    
    if (arrange) {
        grid.arrange(pltX, pltOmega, nrow=1)
    } else {
        list(pltX, pltOmega)
    }
}

# sim_3$proportions[, colnames(sim_3_data_flt1)])
coercePredTrueProps <- function(pred_props, true_props) {
    pred_props <- pred_props[guessOrder(pred_props, true_props), ]
    colnames(pred_props) <- colnames(true_props)
    rownames(pred_props) <- rownames(true_props)
    list(pred_props, true_props)
}

visualize_lo2 <- function(lo2, true_props) {
    show(lo2$plotErrors())
    plotPointsColored(lo2, "current")
    if (!is.null(true_props)) {
        pred_true_props <- coercePredTrueProps(lo2$full_proportions, true_props)
        show(do.call("plotProportions", pred_true_props) +
          labs(caption = paste("RMSE =", do.call("rmse", pred_true_props))))
    }
}

run_block <- function(
    lo2,
    true_props = NULL,
    coef_der_X = 0.001,
    coef_der_Omega = 0.001,
    coef_hinge_H = 10,
    coef_hinge_W = 1,
    coef_pos_D_h = 0,
    coef_pos_D_w = 0,
    iterations = 1000,
    startWithInit = F
) {
    lo2$runGradientBlock(
        coef_der_X = coef_der_X,
        coef_der_Omega = coef_der_Omega,
        coef_hinge_H = coef_hinge_H,
        coef_hinge_W = coef_hinge_W,
        coef_pos_D_h = coef_pos_D_h,
        coef_pos_D_w = coef_pos_D_w,
        iterations = iterations,
        startWithInit = startWithInit
    )
}

get_normalized_svd_projections <- function(
    self,
    dims = seq_len(3),
    keep_genes = NULL,
    keep_samples = NULL,
    plot = T
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
    if (plot) {
        plot(svd_$d)
    }
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

calc_mad_cutoff <- function(lo2, samples = T, k = 2) {
    if (is.null(lo2$distance_samples)) {
        lo2$calculateDistances()
    }
    if (samples) {
        ds <- lo2$distance_samples
    } else {
        ds <- lo2$distance_genes
    }
    lo <- max(median(ds) - k*mad(ds), 0)
    hi <- median(ds) + k*mad(ds)
    c(lo, hi)
}

do_mad_cutoff <- function(lo2, samples = T, k = 2){
    borders <- calc_mad_cutoff(lo2, samples, k)
    if (samples) {
        plt <- lo2$plotDistances(arrange = FALSE)[[2]]
    } else {
        plt <- lo2$plotDistances(arrange = FALSE)[[1]]
    }
    show(
        plt +
        geom_hline(yintercept=borders[1], col="red") +
        geom_hline(yintercept=borders[2], col="red")
    )
    if (samples) {
        ds <- lo2$distance_samples
    } else {
        ds <- lo2$distance_genes
    }
    names(ds[ds > borders[1] & ds < borders[2]])
}

do_knn_cutoff <- function(lo2, keep_n_genes = 4001, k_neighbours = 20, svd_dims = NULL, cell_types = 3, plot_simplices = TRUE) {
    if (is.null(svd_dims)) {
        genes <- lo2$V_row
    } else {
        genes <- get_normalized_svd_projections(lo2, seq_len(svd_dims))[[1]]
    }
    genes_nn_dist <- sort(dbscan::kNNdist(genes, k = k_neighbours))
    threshold <- genes_nn_dist[keep_n_genes]
    plot(genes_nn_dist, type="l")
    abline(h=threshold, col="red")
    keep_genes <- names(genes_nn_dist[genes_nn_dist < threshold])
    init_lo2(lo2$dataset[keep_genes,], top_genes = length(keep_genes), cell_types = cell_types, plot = plot_simplices)
}