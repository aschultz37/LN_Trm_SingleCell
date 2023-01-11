library(dplyr)
library(Seurat)
library(patchwork)

# INITIALIZING SEURAT OBJECT
# Set working directory
#setwd('/gpfs/home/acs9950/singlecell/2022-12-29/')

# Load the dataset
scobj.data <- Read10X(data.dir="LN/filtered/")

# Initialize Seurat object with the raw data
scobj <- CreateSeuratObject(counts=scobj.data, project="2022-12-29_LN", 
                         min.cells=3, min.features=200)
scobj

# QC & FILTERING CELLS
# [[ operator adds column to object metadata (good place to store QC info)
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern="^mt-")

# Visualize QC metrics as violin plot
VlnPlot(scobj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

# FeatureScatter used to visualize feature-feature relationships
# Can also be used for anything calculated by the object (col in metadata, 
# PC scores, etc.)
plot1 <- FeatureScatter(scobj, feature1="nCount_RNA", feature2="percent.mt")
plot2 <- FeatureScatter(scobj, feature1="nCount_RNA", feature2="nFeature_RNA")
plot1 + plot2

scobj <- subset(scobj, subset=(nFeature_RNA > 100 & nFeature_RNA < 5500 
                & percent.mt < 4))

# NORMALIZING DATA
scobj <- NormalizeData(scobj, normalization.method="LogNormalize",
                       scale.factor=10000)

# IDENTIFYING HIGHLY VARIABLE FEATURES
scobj <- FindVariableFeatures(scobj, selection.method="vst", nfeatures=2500)

# find top 10 most variable genes
top10 <- head(VariableFeatures(scobj), 10)

# plot variable features (w/o labels)
plot1 <- VariableFeaturePlot(scobj)
plot2 <- LabelPoints(plot=plot1, points=top10, repel=TRUE, xnudge=0, ynudge=0)
plot1 + plot2

# SCALING DATA
# results of scaling are stores in scobj[["RNA"]]@scale.data
all.genes <- rownames(scobj)
scobj <- ScaleData(scobj, features=all.genes)

# LINEAR DIMENSIONAL REDUCTION
scobj <- RunPCA(scobj, features=VariableFeatures(object=scobj))
#different ways to visualize
print(scobj[["pca"]], dims=1:5, nfeatures=5)
VizDimLoadings(scobj, dims=1:2, reduction="pca")
DimPlot(scobj, reduction="pca")
DimHeatmap(scobj, dims=1:2, cells=500, balanced=TRUE)

# DETERMINE DIMENSIONALITY OF DATASET
# determines # of PC to include (low p val); removes noise
scobj <- JackStraw(scobj, dims=30, num.replicate=100)
scobj <- ScoreJackStraw(scobj, dims=1:30)

# can visually determine dropoff in p-value after certain # PC
JackStrawPlot(scobj, dims=1:30)

# alternatively, can use elbow plot to determine which # PC capture most signal
ElbowPlot(scobj, ndims=30)

# CLUSTERING CELLS
scobj <- FindNeighbors(scobj, dims=1:15)     # choose dim based on PCA
scobj <- FindClusters(scobj, resolution=0.3) # resolution determines # clusters

# View cluster ID of cells
head(Idents(scobj), 5)

# RUN UMAP/tSNE
scobj <- RunUMAP(scobj, dims=1:15) #choose dim based on PCA & FindNeighbors

# note that you can set 'label = TRUE' or use LabelClusters function to
# label the individual clusters
DimPlot(scobj, reduction="umap")
# print number of cells per cluster
table(Idents(scobj))

# save the object so don't have to redo the above pre-processing steps
saveRDS(scobj, file="LN/output/2022-12-29_LN.rds")

# FINDING DIFF. EXPRESSED FEATURES
# FindMarkers can be used to compare clusters/groups vs. e/o or vs. all
# (specified by ident.1, ident.2)
# min.pct requires features to be expressed by amt b/w groups (higher=faster)
# max.cells.per.ident downsamples each identity class (lower=faster)

# find all markers of cluster 2
#cluster2.markers <- FindMarkers(scobj, ident.1=2, min.pct=0.25)
#head(cluster2.markers, n=10)

# find all markers distinguishing cluster 0 and 1
cluster0v1.markers <- FindMarkers(scobj, ident.1=0, ident.2=c(1), min.pct=0.25)
cluster0v1.markers.sig <- subset(cluster0v1.markers, p_val_adj<=0.05)
head(cluster0v1.markers.sig, n=50)
write.csv(cluster0v1.markers.sig, "LN/output/cluster0v1.csv", row.names=TRUE)

# find all markers distinguishing cluster 1 and 2
cluster1v2.markers <- FindMarkers(scobj, ident.1=1, ident.2=c(2), min.pct=0.25)
cluster1v2.markers.sig <- subset(cluster1v2.markers, p_val_adj<=0.05)
head(cluster1v2.markers.sig, n=50)
write.csv(cluster1v2.markers.sig, "LN/output/cluster1v2.csv", row.names=TRUE)

# find all markers distinguishing cluster 2 and 0
cluster2v0.markers <- FindMarkers(scobj, ident.1=2, ident.2=c(0), min.pct=0.25)
cluster2v0.markers.sig <- subset(cluster2v0.markers, p_val_adj<=0.05)
head(cluster2v0.markers.sig, n=50)
write.csv(cluster2v0.markers.sig, "LN/output/cluster2v0.csv", row.names=TRUE)

# find all markers distinguishing cluster 2 from 0 and 1
cluster2v01.markers <- FindMarkers(scobj, ident.1=2, ident.2=c(0, 1), 
                                  min.pct=0.25)
cluster2v01.markers.sig <- subset(cluster2v01.markers, p_val_adj<=0.05)
head(cluster2v01.markers.sig, n=50)
write.csv(cluster2v01.markers.sig, "LN/output/cluster2v01.csv", row.names=TRUE)

# find markers for every cluster compared to all remaining cells
# N.B. only.pos=TRUE reports only the positive markers
scobj.markers <- FindAllMarkers(scobj, only.pos=TRUE, min.pct=0.25,
                                logfc.threshold=0.25)
scobj.markers %>%
    group_by(cluster) %>%
    slice_max(n=4, order_by=avg_log2FC)

cluster0.markers.sig = subset(subset(scobj.markers, cluster=0), p_val_adj<=0.05)
cluster1.markers.sig = subset(subset(scobj.markers, cluster=1), p_val_adj<=0.05)
cluster2.markers.sig = subset(subset(scobj.markers, cluster=2), p_val_adj<=0.05)
cluster3.markers.sig = subset(subset(scobj.markers, cluster=3), p_val_adj<=0.05)
cluster4.markers.sig = subset(subset(scobj.markers, cluster=4), p_val_adj<=0.05)
cluster5.markers.sig = subset(subset(scobj.markers, cluster=5), p_val_adj<=0.05)
cluster6.markers.sig = subset(subset(scobj.markers, cluster=6), p_val_adj<=0.05)
cluster7.markers.sig = subset(subset(scobj.markers, cluster=7), p_val_adj<=0.05)

write.csv(cluster0.markers.sig, "LN/output/cluster0.csv", row.names=TRUE)
write.csv(cluster1.markers.sig, "LN/output/cluster1.csv", row.names=TRUE)
write.csv(cluster2.markers.sig, "LN/output/cluster2.csv", row.names=TRUE)
write.csv(cluster3.markers.sig, "LN/output/cluster3.csv", row.names=TRUE)
write.csv(cluster4.markers.sig, "LN/output/cluster4.csv", row.names=TRUE)
write.csv(cluster5.markers.sig, "LN/output/cluster5.csv", row.names=TRUE)
write.csv(cluster6.markers.sig, "LN/output/cluster6.csv", row.names=TRUE)
write.csv(cluster7.markers.sig, "LN/output/cluster7.csv", row.names=TRUE)

# multiple tests for differential expression are avail, e.g.:
#cluster0.markers <- FindMarkers(scobj, ident.1=0, logfc.threshold=0.25,
#                                test.use="roc", only.pos=TRUE)

# expression probability distributions across clusters
genelist = c("Itgae", "Ccr7", "Klf2", "Cxcr6", "S1pr1")
VlnPlot(scobj, features=genelist)

# plot raw counts
VlnPlot(scobj, features=genelist, slot="counts", log=TRUE)

# visualize feature expression on tSNE or PCA plot
FeaturePlot(scobj, features=genelist)

# N.B. can also try RidgePlot, CellScatter, and DotPlot to view dataset
# ridge plot
RidgePlot(scobj, features=genelist, ncol=2)

# generate expression heatmap for top 10 markers for each cluster
scobj.markers %>%
    group_by(cluster) %>%
    top_n(n=10, wt=avg_log2FC) -> top10
DoHeatmap(scobj, features=top10$gene) + NoLegend()

# ASSIGN CELL TYPE TO CLUSTERS
#new.cluster.ids <- c("Type1", "Type2", "Type3", "Type4", "Type5", "Type6")
#names(new.cluster.ids) <- levels(scobj)
#scobj <- RenameIdents(scobj, new.cluster.ids)
#DimPlot(scobj, reduction="umap", label=TRUE, pt.size=0.5) + NoLegend()

# save again, includes plots
saveRDS(scobj, file="LN/output/2022-12-29_LN_figures.rds")
