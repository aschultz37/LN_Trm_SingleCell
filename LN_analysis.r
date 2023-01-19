library(dplyr)
library(Seurat)
library(patchwork)

# INITIALIZING SEURAT OBJECT
# Set working directory
#setwd('/gpfs/home/acs9950/singlecell/2022-12-29/')

# LOADING THE SEURAT OBJECT
#scobj <- readRDS("LN/output/2022-12-29_LN_figures.rds")

# Load the dataset
scobj.data <- Read10X(data.dir="LN/filtered/")

# Initialize Seurat object with the raw data
scobj <- CreateSeuratObject(counts=scobj.data, project="2022-12-29_LN", 
                         min.cells=3, min.features=200)
scobj

# QC & FILTERING CELLS
# [[ operator adds column to object metadata (good place to store QC info)
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern="^mt-")
# do the same for ribosomal proteins (Rps, Rpl)
scobj[["percent.rps"]] <- PercentageFeatureSet(scobj, pattern="Rps")
scobj[["percent.rpl"]] <- PercentageFeatureSet(scobj, pattern="Rpl")

# Visualize QC metrics as violin plot
VlnPlot(scobj, features=c("nFeature_RNA", "nCount_RNA"), ncol=2)
VlnPlot(scobj, features=c("percent.mt", "percent.rps", "percent.rpl"), ncol=3)

# FeatureScatter used to visualize feature-feature relationships
# Can also be used for anything calculated by the object (col in metadata, 
# PC scores, etc.)
plot1 <- FeatureScatter(scobj, feature1="nCount_RNA", feature2="percent.mt")
plot2 <- FeatureScatter(scobj, feature1="nCount_RNA", feature2="nFeature_RNA")
plot1 + plot2

# nFeature increased to 6000 & percent.mt to < 5
scobj <- subset(scobj, subset=(nFeature_RNA > 100 & nFeature_RNA < 6000 
                & percent.mt < 5 & percent.rps < 20 & percent.rpl < 23))

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
#head(Idents(scobj), 5)

# RUN UMAP/tSNE
scobj <- RunUMAP(scobj, dims=1:15) #choose dim based on PCA & FindNeighbors

# note that you can set 'label = TRUE' or use LabelClusters function to
# label the individual clusters
DimPlot(scobj, reduction="umap")
# print number of cells per cluster
table(Idents(scobj))

# save the object so don't have to redo the above pre-processing steps
saveRDS(scobj, file="LN/output/2022-12-29_LN_rp.rds")

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
#head(cluster0v1.markers.sig, n=50)
write.csv(cluster0v1.markers.sig, "LN/output/cluster0v1.csv", row.names=TRUE)

# find all markers distinguishing cluster 1 and 2
cluster1v2.markers <- FindMarkers(scobj, ident.1=1, ident.2=c(2), min.pct=0.25)
cluster1v2.markers.sig <- subset(cluster1v2.markers, p_val_adj<=0.05)
#head(cluster1v2.markers.sig, n=50)
write.csv(cluster1v2.markers.sig, "LN/output/cluster1v2.csv", row.names=TRUE)

# find all markers distinguishing cluster 2 and 0
cluster2v0.markers <- FindMarkers(scobj, ident.1=2, ident.2=c(0), min.pct=0.25)
cluster2v0.markers.sig <- subset(cluster2v0.markers, p_val_adj<=0.05)
#head(cluster2v0.markers.sig, n=50)
write.csv(cluster2v0.markers.sig, "LN/output/cluster2v0.csv", row.names=TRUE)

# find all markers distinguishing cluster 2 from 0 and 1
cluster2v01.markers <- FindMarkers(scobj, ident.1=2, ident.2=c(0, 1), 
                                  min.pct=0.25)
cluster2v01.markers.sig <- subset(cluster2v01.markers, p_val_adj<=0.05)
#head(cluster2v01.markers.sig, n=50)
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
cluster8.markers.sig = subset(subset(scobj.markers, cluster=8), p_val_adj<=0.05)

write.csv(cluster0.markers.sig, "LN/output/cluster0.csv", row.names=TRUE)
write.csv(cluster1.markers.sig, "LN/output/cluster1.csv", row.names=TRUE)
write.csv(cluster2.markers.sig, "LN/output/cluster2.csv", row.names=TRUE)
write.csv(cluster3.markers.sig, "LN/output/cluster3.csv", row.names=TRUE)
write.csv(cluster4.markers.sig, "LN/output/cluster4.csv", row.names=TRUE)
write.csv(cluster5.markers.sig, "LN/output/cluster5.csv", row.names=TRUE)
write.csv(cluster6.markers.sig, "LN/output/cluster6.csv", row.names=TRUE)
write.csv(cluster7.markers.sig, "LN/output/cluster7.csv", row.names=TRUE)
write.csv(cluster8.markers.sig, "LN/output/cluster8.csv", row.names=TRUE)

# multiple tests for differential expression are avail, e.g.:
#cluster0.markers <- FindMarkers(scobj, ident.1=0, logfc.threshold=0.25,
#                                test.use="roc", only.pos=TRUE)

# expression probability distributions across clusters
genelist = c("Itgae", "Ccr7", "Klf2", "Cxcr6", "S1pr1", 
             "Cd8a", "Thy1", "Ptprc", "Cd3e")
VlnPlot(scobj, features=genelist)

# plot log transformed counts
VlnPlot(scobj, features=genelist, slot="counts", log=TRUE)

# visualize feature expression on tSNE or PCA plot
FeaturePlot(scobj, features=genelist)
# to detect naive T cells
FeaturePlot(scobj, features="Cd44")

# N.B. can also try RidgePlot, CellScatter, and DotPlot to view dataset
# ridge plot
RidgePlot(scobj, features=genelist, ncol=3)

# generate expression heatmap for top 10 markers for each cluster
scobj.markers %>%
    group_by(cluster) %>%
    top_n(n=10, wt=avg_log2FC) -> top10
DoHeatmap(scobj, features=top10$gene) + NoLegend()

# generate a heatmap based on cluster2v01
#choose markers based on avg_log2FC -- clean up Gm and Rp genes first

# ASSIGN CELL TYPE TO CLUSTERS
#new.cluster.ids <- c("Type1", "Type2", "Type3", "Type4", "Type5", "Type6")
#names(new.cluster.ids) <- levels(scobj)
#scobj <- RenameIdents(scobj, new.cluster.ids)
#DimPlot(scobj, reduction="umap", label=TRUE, pt.size=0.5) + NoLegend()

# save again, includes plots
saveRDS(scobj, file="LN/output/2022-12-29_LN_figures_rp.rds")

# REMOVE CLUSTERS AND RE-ANALYZE
#sub1_scobj = readRDS("LN/output/2022-12-29_LN_sub1.rds")
# Remove clusters 4, 7, 8 and keep the remaining for new analysis
sub1_scobj = subset(x=scobj, idents=c(4, 7, 8), invert=TRUE)

# reidentify highly variable features
sub1_scobj <- FindVariableFeatures(sub1_scobj, selection.method="vst", nfeatures=2500)

# find top 10 most variable genes
sub1_top10 <- head(VariableFeatures(sub1_scobj), 10)

# rescale data
# results of scaling are stores in sub1_scobj[["RNA"]]@scale.data
sub1_all.genes <- rownames(sub1_scobj)
sub1_scobj <- ScaleData(sub1_scobj, features=sub1_all.genes)

# re-run pca analysis
sub1_scobj <- RunPCA(sub1_scobj, features=VariableFeatures(object=sub1_scobj))

# determines # of PC to include (low p val); removes noise
sub1_scobj <- JackStraw(sub1_scobj, dims=30, num.replicate=100)
sub1_scobj <- ScoreJackStraw(sub1_scobj, dims=1:30)

# can visually determine dropoff in p-value after certain # PC
JackStrawPlot(sub1_scobj, dims=1:30)

# alternatively, can use elbow plot to determine which # PC capture most signal
ElbowPlot(sub1_scobj, ndims=30)

# recluster cells
sub1_scobj <- FindNeighbors(sub1_scobj, dims=1:25)     # choose dim based on PCA
sub1_scobj <- FindClusters(sub1_scobj, resolution=0.3) # resolution determines # clusters

# run umap
sub1_scobj <- RunUMAP(sub1_scobj, dims=1:25) #choose dim based on PCA/Neighbors
# note that you can set 'label = TRUE' or use LabelClusters function to
# label the individual clusters
custom_color_palette = c("#1B9E77", "#D95F02", "#7570B3", "#CB008A", "#66A61E", 
                         "#E6AB02", "#A6761D", "#666666")
DimPlot(sub1_scobj, reduction="umap", pt.size=2.5, cols=custom_color_palette)
# print number of cells per cluster
table(Idents(sub1_scobj))

# find all markers
sub1_scobj.markers <- FindAllMarkers(sub1_scobj, only.pos=TRUE, min.pct=0.25,
                                logfc.threshold=0.25)

# expression probability distributions across clusters
sub1_genelist = c("Itgae", "Ccr7", "Klf2", "Cxcr6", "S1pr1", 
             "Cd8a", "Thy1", "Ptprc", "Cd3e")
LN_Trm_genes <- c("Acap1", "Actn2", "Amica1", "Arhgef1", "Atxn7l3b", "Aw112010", 
                  "B4galnt1", "Bcl11b", "Cbx3", "Ccnd2", "Ccr10", "Cd27", "Cd7", 
                  "Cd74", "Chd3", "Cirbp", "Clec2d", "Crot", "Csf1", "Cxcr3", 
                  "Cxcr6", "Eif5", "Evl", "Fam189b", "Fubp1", "Fyb", "Gramd1a", 
                  "Sema4a","Gstp1", "Shisa5", "H2-T23", "Sipa1", "Hmgn1", 
                  "Slfn2", "Hmha1", "Sp100", "Hsp90b1", "Spcs2", "Id2", "Srrm2", 
                  "Ifitm10", "Stap1", "Ikzf3", "Tbc1d10c", "Il16", "Tesc", 
                  "Il18r1", "Tnfaip8", "Il7r", "Tnrc6a", "Irf2bpl", "Tsc22d4", 
                  "Itgae", "Uba52", "Itgal", "Ucp2", "Itm2c", "Wbp1", "Lfng", 
                  "Xist", "Lpar6", "Ypel3", "Lrrc58", "Zbtb7a", "Ltb", "Znrf1", 
                  "Ly6a", "Ly6e", "Ly6g5b", "Malat1", "Mbnl1", "Mrpl52", "Mxd4", 
                  "Mycbp2", "N4bp2l2", "N4bp2l2", "Ndufa3", "Ndufa5", "Nktr", 
                  "Nudcd3", "Ogt", "Pdcd4", "Pdia3", "Pdia6", "Ptpn7", "Ptprc", 
                  "Rapgef6", "Rbpj", "Rgs10", "Rpl15", "Rpl35", "Rpl38", 
                  "Rps28", "Rps29", "Sash3")
cytotoxic_gene_list = c("Prf1", "Gzmb", "Gzmk", "Ccl4", "Ccl5", "Csf1")
traffic_gene_list = c("S1pr1", "Ccr7", "Cxcr4", "Cxcr3", "Cxcr6")
clusters_of_interest = c("0", "1", "2", "3")
VlnPlot(sub1_scobj, features=sub1_genelist, cols=custom_color_palette,
        idents=NULL)
VlnPlot(sub1_scobj, features=cytotoxic_gene_list, cols=custom_color_palette,
        idents=clusters_of_interest)
VlnPlot(sub1_scobj, features=traffic_gene_list, cols=custom_color_palette,
        idents=clusters_of_interest)

# visualize feature expression on tSNE or PCA plot
FeaturePlot(sub1_scobj, features=sub1_genelist)

# visualize LN_Trm genes in feature plots
# FeaturePlot for list of genes of arbitrary length (groups of 12)
final_index <- 0
for(i in 1:length(LN_Trm_genes)){
  if((i %% 12) == 0){
    print(FeaturePlot(sub1_scobj, features=LN_Trm_genes[(i-11):i]))
    final_index <- i
  }
}
FeaturePlot(sub1_scobj, 
            features=LN_Trm_genes[(final_index+1):length(LN_Trm_genes)])


# to detect naive T cells
FeaturePlot(sub1_scobj, features="Cd44")

# N.B. can also try RidgePlot, CellScatter, and DotPlot to view dataset
# ridge plot
RidgePlot(sub1_scobj, features=sub1_genelist, ncol=3)

# generate expression heatmap for top 10 markers for each cluster
sub1_scobj.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> sub1_top10
DoHeatmap(sub1_scobj, features=sub1_top10$gene) + NoLegend()

# find all markers distinguishing cluster 0 and 1
cluster0v1_sub1.markers <- FindMarkers(sub1_scobj, ident.1=0, ident.2=c(1), 
                                       min.pct=0.25)
cluster0v1_sub1.markers.sig <- subset(cluster0v1_sub1.markers, p_val_adj<=0.05)
write.csv(cluster0v1_sub1.markers.sig, "LN/output/cluster0v1_sub1.csv", 
          row.names=TRUE)

# find all markers distinguishing cluster 2 from 0 and 1
cluster2v01_sub1.markers <- FindMarkers(sub1_scobj, ident.1=2, ident.2=c(0, 1), 
                                   min.pct=0.25)
cluster2v01_sub1.markers.sig <- subset(cluster2v01_sub1.markers, p_val_adj<=0.05)
write.csv(cluster2v01_sub1.markers.sig, "LN/output/cluster2v01_sub1.csv", 
          row.names=TRUE)

saveRDS(sub1_scobj, "LN/output/2022-12-29_LN_sub1.rds")
