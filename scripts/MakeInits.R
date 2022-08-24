source("scripts/SinkhornNNLSLinseedC.R")

prepareInit <- function(obj_, ids_X, ids_Omega) {

    V__ <- obj_$S %*% obj_$V_row %*% t(obj_$R)
    
    ## Omega
    Ae <- obj_$V_column[, ids_X]
    init_Omega <- obj_$S %*% Ae
    
    ## X
    Ae <- obj_$V_row[ids_Omega, ]
    init_X <- Ae %*% t(obj_$R)
    
    ## calculate D_w and D_h
    ## vectorizing deconvolution
    vec_mtx <- matrix(0,obj_$cell_types*obj_$cell_types,obj_$cell_types)
    for (col_ in 1:obj_$cell_types) {
      vec_mtx[,col_] <- cbind(c(t(t(init_Omega[,col_])) %*% init_X[col_,]))
    }
    ## adding sum-to-one constraint
    init_D_w <- matrix(nnls(rbind(vec_mtx,init_Omega),rbind(cbind(c(V__)),obj_$B))$x,nrow=obj_$cell_types,ncol=1)
    init_D_h <- init_D_w * (obj_$N/obj_$M)
    
    list(init_Omega=init_Omega,
                 init_X=init_X,
                 init_D_w= init_D_w,
                 init_D_h = init_D_h)
}

tmp_snk <- readRDS(snakemake@input[[1]])
init_strategy <- snakemake@config[["init_strategy"]]
if (is.null(init_strategy)){
  init_strategy <- "SelectRandom"
}
print(init_strategy)
if (! init_strategy %in% c("SelectX","SelectOmega","SelectRandom","SelectCentered"))
  stop("Selected initialization is not allowed. Available values 'SelectX','SelectOmega','SelectRandom','SelectCentered'")
if (init_strategy == "SelectRandom") {
  tmp_1 <- tmp_snk$initWithSubset(10000,snakemake@config[["num_inits"]])
}

if (init_strategy == "SelectX") {
  for (idx in 1:snakemake@config[["num_inits"]]) {
    tmp_snk$selectInitX()
    saveRDS(list(init_Omega=tmp_snk$init_Omega,
                 init_X=tmp_snk$init_X,
                 init_D_w= tmp_snk$init_D_w,
                 init_D_h = tmp_snk$init_D_h),
                 snakemake@output[[idx]])
  }
}

if (init_strategy == "SelectOmega") {
  for (idx in 1:snakemake@config[["num_inits"]]) {
    tmp_snk$selectInitOmega()
    saveRDS(list(init_Omega=tmp_snk$init_Omega,
                 init_X=tmp_snk$init_X,
                 init_D_w= tmp_snk$init_D_w,
                 init_D_h = tmp_snk$init_D_h),
                 snakemake@output[[idx]])
  }
}

if (init_strategy == "SelectCentered") {
  for (idx in 1:snakemake@config[["num_inits"]]) {
    tmp_snk$selectInitRandomCentered()
    saveRDS(list(init_Omega=tmp_snk$init_Omega,
                 init_X=tmp_snk$init_X,
                 init_D_w= tmp_snk$init_D_w,
                 init_D_h = tmp_snk$init_D_h),
                 snakemake@output[[idx]])
  }
}

if (init_strategy == "SelectRandom") {
  tmp_1 <- tmp_snk$initWithSubset(10000,snakemake@config[["num_inits"]])
  ctn <- ncol(tmp_1$idsTableOmega)
  for (idx in 1:snakemake@config[["num_inits"]]) {
    saveRDS(prepareInit(tmp_snk,tmp_1$idsTableX[idx,-ctn],tmp_1$idsTableOmega[idx,-ctn]),
            snakemake@output[[idx]])
  }
}


