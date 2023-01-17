library(iCellR)
library(cowplot)

# Read Data into Object
LN.data <- load10x("../LN/filtered/")

# Make iCellR Object
LN.obj <- make.obj(LN.data)

# QC Data on Object
LN.obj <- qc.stats(LN.obj)
# plot UMIs, genes and percent mito all at once and in one plot. 
# you can make them individually as well, see the arguments ?stats.plot.
stats.plot(LN.obj,
           plot.type="three.in.one",
           out.name="UMI-plot",
           interactive=FALSE,
           cell.color="slategray3", 
           cell.size=1, 
           cell.transparency=0.5,
           box.color="red",
           box.line.col="green")

# Filter Cells
LN.obj <- cell.filter(LN.obj,
                      min.mito=0,
                      max.mito=0.04,
                      min.genes=100,
                      max.genes=5500,
                      min.umis=0,
                      max.umis=Inf)
# check to see how many cells are left.  
dim(LN.obj@main.data)

# Normalize Data
LN.obj <- norm.data(LN.obj, 
                    norm.method="ranked.glsf", # dif method than Seurat
                    top.rank=500) # best for scRNA-Seq

# Note: scaling data is not necessary because iCellR scales at runtime
#LN.obj <- data.scale(LN.obj)

# Gene Stats
LN.obj <- gene.stats(LN.obj, which.data="main.data")
head(LN.obj@gene.data[order(LN.obj@gene.data$numberOfCells, decreasing=T),])

# Gene Model for Clustering
#see model
make.gene.model(LN.obj, my.out.put="plot",
                dispersion.limit=1.5, 
                base.mean.rank=500,
                gene.num.max=5500,
                no.mito.model=T, 
                mark.mito=T, 
                interactive=F,
                out.name="gene.model")
#write the model data into object
LN.obj <- make.gene.model(LN.obj, my.out.put="data",
                          dispersion.limit=1.5, 
                          base.mean.rank=500,
                          gene.num.max=5500,
                          no.mito.model=T, 
                          mark.mito=T, 
                          interactive=F,
                          out.name="gene.model")
head(LN.obj@gene.model)

# Principal Component Analysis (PCA)
# Note: skip this step if doing batch correction (use iba function)
# using list of all genes instead of LN.obj@gene.model
LN.obj <- run.pca(LN.obj, method="gene.model", 
                  gene.list=LN.obj@gene.data$genes, data.type="main")
opt.pcs.plot(LN.obj)

# Clustering & Dimensionality Reduction (UMAP)
LN.obj <- run.umap(LN.obj, dims = 1:15)
LN.obj <- iclust(LN.obj, sensitivity=150, data.type="pca", dims=1:15)
cluster.plot(LN.obj, plot.type="umap", interactive=F,
             cell.size=0.5, cell.transparency=1, anno.clust=T)
# reorder clusters based on distance (useful for merging clusters)
LN.obj <- clust.ord(LN.obj, top.rank=500, how.to.order="distance")
cluster.plot(LN.obj,plot.type="umap", interactive=F,
             cell.size=0.5, cell.transparency=1, anno.clust=T)

# PseudoTime analysis & Cell Cycle Prediction
#see tutorial document, not relevant right now
#uses KNetL

# Many more QC stats available for clusters
#see tutorial document, not relevant right now
#gene-gene correlation

# Find Marker Genes
marker.genes <- findMarkers(LN.obj, fold.change=2, padjval=0.1)
dim(marker.genes)
head(marker.genes)

# Heatmap
# find top genes
MyGenes <- top.markers(marker.genes, topde=10, min.base.mean=0.2, filt.ambig=F)
MyGenes <- unique(MyGenes)
# main data 
heatmap.gg.plot(LN.obj, gene=MyGenes, interactive=F, 
                cluster.by="clusters", conds.to.plot=NULL)

# Feature Plots
genelist = c("Itgae", "Ccr7", "Klf2", "Cxcr6", "S1pr1")
rm(list = ls(pattern="PL_"))
for(i in genelist){
  MyPlot <- gene.plot(LN.obj, gene = i,
                      interactive = F,
                      cell.size = 0.1,
                      plot.data.type = "umap",
                      data.type = "main",
                      scaleValue = T,
                      min.scale = 0,max.scale = 2.0,
                      cell.transparency = 1)
  NameCol=paste("PL",i,sep="_")
  eval(call("<-", as.name(NameCol), MyPlot))
}

filenames <- ls(pattern="PL_")

B <- cluster.plot(LN.obj, plot.type="umap", interactive=F, 
                  cell.size=0.1, cell.transparency=1, anno.clust=T)
filenames <- c("B", filenames)

png('genes_umap.png', width=15, height=12, units='in', res=300)
plot_grid(plotlist=mget(filenames))
dev.off()

# Comparing Clusters/Differential Expression
diff.res <- run.diff.exp(LN.obj, de.by="clusters", cond.1=c(1), cond.2=c(2))
diff.res1 <- as.data.frame(diff.res)
diff.res1 <- subset(diff.res1, padj<0.05)
head(diff.res1)

# Volcano Plot
volcano.ma.plot(diff.res,
                sig.value = "pval",
                sig.line = 0.05,
                plot.type = "volcano",
                interactive = F)

# Merging Clusters
#LN.obj <- change.clust(LN.obj, change.clust=3, to.clust=2)

# Unmerge/Reset Clusters
#LN.obj <- change.clust(LN.obj, clust.reset=T)

# Renaming Clusters
# reset after this so you can run further analysis
#LN.obj <- change.clust(LN.obj, change.clust=2, to.clust="Trm")

# Removing Clusters
#LN.obj <- clust.rm(LN.obj, clust.to.rm=1)

# Cell Gating Tool
# by gene
#LN.plot <- gene.plot(LN.obj, gene="S1pr1", 
#                     plot.type="scatterplot",
#                     clust.dim=2,
#                     interactive=F)
# by cluster
#LN.plot <- cluster.plot(LN.obj, 
#                        cell.size=1, 
#                        cell.transparency=0.5, 
#                        clust.dim=2, 
#                        interactive=F)
#cell.gating(LN.obj, LN.plot = LN.plot, plot.type = "umap")	
# rename cluster after downloading cell IDs
#LN.obj <- gate.to.clust(LN.obj, my.gate="cellGating.txt", to.clust=10)