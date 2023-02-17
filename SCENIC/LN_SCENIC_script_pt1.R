library(SCENIC)

setwd("/gpfs/home/acs9950/singlecell/2022-12-29/SCENIC/")

# Run commented code only if RDS does not exist.
# scenicOptions <- initializeScenic(org="mgi",
#                                   dbDir="resources/cisTarget_databases")
# motifAnnotations_mgi <- motifAnnotations
# scenicOptions <- initializeScenic(org="mgi",
#                                   dbDir="resources/cisTarget_databases")
# saveRDS(scenicOptions, file="resources/scenicOptions.Rds")

scenicOptions <- readRDS("resources/scenicOptions.Rds")

library(dplyr)
library(Seurat)
library(loomR)
library(SeuratDisk)
seurat_obj <- readRDS("../LN/output/2022-12-29_LN_sub1_rmbc_figs.rds")
# note that this causes an error if loom file already exists
# it will not overwrite the existing file and will halt execution
loom <- as.loom(seurat_obj)

library(SCopeLoomR)

# loom <- open_loom("2022-12-29.loom", mode="r+") # does not work even if file exists
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)