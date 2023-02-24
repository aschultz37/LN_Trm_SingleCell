motifAnnotations_mgi <- motifAnnotations
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

export2loom(scenicOptions, exprMat)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
saveRDS(exprMat_log, file="int/exprMat_log.Rds")