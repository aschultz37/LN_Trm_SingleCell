# NOTE: This section is interactive, run in RStudio

library(SCENIC)
library(loomR)
library(SCopeLoomR)
library(dplyr)
library(Seurat)
library(AUCell)

setwd("/gpfs/home/acs9950/singlecell/2022-12-29/SCENIC")

scenicOptions <- readRDS("int/scenicOptions.Rds")
exprMat_log <- readRDS("int/exprMat_log.Rds")

seurat_obj <- readRDS("../LN/output/2022-12-29_LN_sub1_rmbc_figs.rds")

aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

cellInfo <- data.frame(seuratCluster=Idents(seurat_obj))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo$seuratCluster)
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
saveRDS(exprMat_log, file="int/exprMat_log.Rds")

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

nPcs <- c(15)
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), 
                     onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), 
                                             value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), 
                                             value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=TRUE, 
                         varName="seuratCluster", cex=.5)
# Using only "high-confidence" regulons (normally similar)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), 
                                             value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=TRUE, 
                         varName="seuratCluster", cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 15
scenicOptions@settings$defaultTsne$perpl <- 50
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

exprMat <- exprMat_log
scenicOptions@fileNames$output["loomFile",] <- "output/LN_SCENIC_manual.loom"
export2loom(scenicOptions, exprMat)

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))],
                        plots="Expression")

# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", 
                width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

# show several simultaneously (this is broken, documentation out of date)
# par(mfrow=c(1,2))
# 
# regulonNames <- c("Klf2", "Nfkb1")
# cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
# text(-5,-23, attr(cellCol,"red"), col="red", cex=.7)
# text(-10,-18, attr(cellCol,"green"), col="green", cex=.7)
# 
# regulonNames <- list(red=c("Klf2","Nfkb1"),
#                      green=c("Irf7"),
#                      blue=c("Bhlhe40"))
# cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="Binary")

regulons <- loadInt(scenicOptions, "regulons")
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, 
                                             "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Dlx5"]
viewMotifs(tableSubset) 

# regulators of known clusters heatmap
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")


topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

# binarized version of heatmap
minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))

topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)

# cell-type specific regulators
# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), 
               cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


# plot regulons on Seurat
grep("Bhlhe40_extended (12g)", names(regulons), value=T)
Bhlhe40_targets <- list(regulons[["Bhlhe40_extended (12g)"]])
seurat.integrated <- AddModuleScore(seurat_obj, features = Bhlhe40_targets, 
                                    ctrl = 50, name = "Bhlhe40_signature.list")
FeaturePlot(seurat.integrated, features = "Bhlhe40_signature.list1", 
            order = T, pt.size = 2.5) + 
            ggtitle("Bhlhe40 regulon") + 
            scale_colour_gradientn(colors = (brewer.pal(n = 9, name = "YlGn")))



saveRDS(scenicOptions, file="int/scenicOptions.Rds")



