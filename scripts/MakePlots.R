source("scripts/OutputPlots.R")

metadata_ <- loadRData(snakemake@input[[1]])

png(snakemake@output[[1]])
plotErrors(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()

png(snakemake@output[[2]])
plotNegProportionsChange(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()

png(snakemake@output[[3]])
plotNegBasisChange(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()

png(snakemake@output[[4]])
plotSumToOneChange(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()
