library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(EnhancedVolcano)

# load copy of desired object
vp_scobj <- readRDS("LN/output/2022-12-29_LN_sub1.rds")

# calculate markers/DEG for desired clusters
cluster2v01_vp.markers <- FindMarkers(vp_scobj, ident.1=2, ident.2=c(0, 1), 
                                        min.pct=0.25)
# filtering steps if wanted
# cluster2v01_vp.markers.pval <- subset(cluster2v01_vp.markers, 
#                                         p_val_adj<=0.01) # filter p_val_adj
# cluster2v01_vp.markers.filter <- subset(cluster2v01_vp.markers.pval, 
#                                          abs(avg_log2FC)>=0.5) # filter log2FC

genes_of_interest = c("Csf1", "Cxcr6") # list of genes to label (selectLab=)
EnhancedVolcano(cluster2v01_vp.markers,
                lab=rownames(cluster2v01_vp.markers),
                x='avg_log2FC',
                y='p_val_adj',
                title="Cluster 2v(0,1)",
                pCutoff=0.01,
                FCcutoff=0.5,
                pointSize=3.0,
                labSize=5.0)
