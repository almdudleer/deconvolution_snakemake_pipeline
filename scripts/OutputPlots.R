library(ggplot2)
library(reshape2)
library(gridExtra)
library(uwot)
library(ggpubr)
library(ggrepel)

loadRData <- function(fileName){
  #loads an R file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

plotErrors <- function(metadata,variables = c("deconv_error","lamdba_error","beta_error",
    "D_h_error","D_w_error","total_error")) {
      toPlot <- data.frame(metadata$errors_statistics[,variables])
      toPlot$iteration <- 0:(nrow(metadata$errors_statistics)-1)
      toPlot <- melt(toPlot,id.vars="iteration",measure.vars = variables)
      plt <- ggplot(toPlot,aes(x=iteration,y=log10(value),color=variable)) +
      geom_point(size=0.2) +
      geom_line() + theme_minimal()
      plt
}

plotNegBasisChange <- function(metadata_) {
  total_ <- nrow(metadata_$V_row)*metadata_$cell_types
  toPlot <- as.data.frame(metadata_$errors_statistics[,"neg_basis_count",drop=F])
  last_basis_count <- round(toPlot[nrow(toPlot),"neg_basis_count"] / total_,6) * 100
  toPlot$idx <- 1:nrow(toPlot)
  plt <- ggplot(toPlot,aes(y=neg_basis_count,x=idx)) + geom_line() + 
    theme_minimal() + xlab("Iteration") + ylab("Negative basis") +
    annotate("text",  x=Inf, y = Inf, label = paste0(last_basis_count,"%"), vjust=1, hjust=1)
  plt
}

plotNegProportionsChange <- function(metadata_) {
  total_ <- ncol(metadata_$V_row)*metadata_$cell_types
  toPlot <- as.data.frame(metadata_$errors_statistics[,"neg_props_count",drop=F])
  last_prop_count <- round(toPlot[nrow(toPlot),"neg_props_count"] / total_,6) * 100
  toPlot$idx <- 1:nrow(toPlot)
  plt <- ggplot(toPlot,aes(y=neg_props_count,x=idx)) + geom_line() + theme_minimal() + 
    xlab("Iteration") + ylab("Negative proportions")+
    annotate("text",  x=Inf, y = Inf, label = paste0(last_prop_count,"%"), vjust=1, hjust=1)
  plt
}

plotSumToOneChange <- function(metadata_) {
  toPlot <- as.data.frame(metadata_$errors_statistics[,"sum_d_w",drop=F])
  last_sum_dw <- toPlot[nrow(toPlot),"sum_d_w"]
  toPlot$idx <- 1:nrow(toPlot)
  ggplot(toPlot,aes(y=sum_d_w,x=idx)) + geom_line() + theme_minimal() + 
    xlab("Iteration") + ylab("Sum-to-one") +
    annotate("text",  x=Inf, y = Inf, label = sprintf("%.7e",last_sum_dw), vjust=1, hjust=1)
}

plotPoints2D <- function(metadata_,points="init",dims=3) {
      if (!points %in% c("init","current")) {
        stop("Allowed values for points are 'init', 'current'")
      }

      if (points == "init") {
        X <- metadata_$init_X
        Omega <- metadata_$init_Omega
        count_neg_props <- metadata_$init_count_neg_props
        count_neg_basis <- metadata_$init_count_neg_basis
      }
      if (points == "current") {
        X <- metadata_$final_X
        Omega <- metadata_$final_Omega
        count_neg_props <- metadata_$count_neg_props
        count_neg_basis <- metadata_$count_neg_basis
      }

      X <- X[,1:dims]
      Omega <- Omega[1:dims,]

      ## plot X
        toPlot <- as.data.frame(metadata_$V_row %*% t(metadata_$R))[,1:dims]
        colnames(toPlot) <- c("X","Y","Z")
        colnames(X) <- c("X","Y","Z")
        pltX <- ggplot(toPlot, aes(x=Y, y=Z)) +
          geom_point() + 
          geom_polygon(data=as.data.frame(X), fill=NA, color = "green") +
          theme_minimal()
      if (!is.null(count_neg_props)) {
        pltX <- pltX + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_props / (metadata_$cell_types*metadata_$N),4)*100,"%"), vjust=1, hjust=1)
      }

      ## plot Omega  
      toPlot <- as.data.frame(t(metadata_$S %*% metadata_$V_column))[,1:dims]
      colnames(toPlot) <- c("X","Y","Z")
      rownames(Omega) <- c("X","Y","Z")
      pltOmega <- ggplot(toPlot, aes(x=Y, y=Z)) +
        geom_point() + 
        geom_polygon(data=as.data.frame(t(Omega)), fill=NA, color = "green") +
        theme_minimal()
      
      if (!is.null(count_neg_basis)) {
        pltOmega <- pltOmega + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_basis / (metadata_$cell_types*metadata_$M),4)*100,"%"), vjust=1, hjust=1)
      }

      grid.arrange(pltX,pltOmega,nrow=1)
}

plotUMAP <- function(data_, best_run){
  print("DIM")
  print(dim(data_))
  toPlot <- as.data.frame(umap(t(data_), n_neighbors=5))
  rownames(toPlot) <- colnames(data_)
  toPlot$best <- grepl(paste0("^",best_run), rownames(toPlot))
  toPlot$best_id <- NULL
  init_number <- strsplit(best_run,"_")[[1]]
  init_number <- init_number[length(init_number)]
  toPlot[toPlot$best,"best_id"] <- rownames(toPlot[toPlot$best,])
  toPlot$best_id <- gsub(best_run,init_number,toPlot$best_id)
  ggplot() +  geom_point(data=toPlot[!toPlot$best,],aes(x=V1,y=V2)) + 
    geom_point(data=toPlot[toPlot$best,],aes(x=V1,y=V2,color=best)) +
    geom_text_repel(data=toPlot[toPlot$best,],aes(x=V1,y=V2,label=best_id),
                    min.segment.length = unit(0, 'lines')) +
    theme_minimal() + theme(legend.position = "none") +
    xlab('UMAP1') + ylab('UMAP2')
}

plotAbundance <- function(metadata_){
  colnames(metadata_$full_proportions) <- colnames(metadata_$filtered_dataset)
  rownames(metadata_$full_proportions) <- paste("Cell type",1:metadata_$cell_types)
  toPlot <- melt(metadata_$full_proportions)
  colnames(toPlot) <- c("cell_type","sample","cluster_percent")
  toPlot$cluster_percent <- round(toPlot$cluster_percent*100,2)

  ggbarplot(toPlot, x = "cell_type", y = "cluster_percent",
                fill = "cell_type",  color  ="black",
                position = position_dodge(0.7), add = c("mean_sd", "point"),
                add.params = list(color = "black"), alpha = 0.95) +
  theme(axis.text.x=element_text(angle = 45, hjust=1, size = 12),
        legend.position="none") + xlab("")
}