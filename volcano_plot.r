# library(dplyr)
# library(Seurat)
# library(patchwork)
# library(RColorBrewer)
library(EnhancedVolcano)

# load copy of desired object
vp_scobj <- readRDS("LN/output/2022-12-29_LN_sub1_cc_figs.rds")

# calculate markers/DEG for desired clusters
cluster2vT_vp.markers <- FindMarkers(vp_scobj, 
                                      ident.1=2, 
                                      ident.2=c(0, 1, 3, 5), 
                                      min.pct=0.25)
# filtering steps if wanted
# cluster2vT_vp.markers.pval <- subset(cluster2vT_vp.markers, 
#                                         p_val_adj<=0.01) # filter p_val_adj
# cluster2vT_vp.markers.filter <- subset(cluster2vT_vp.markers.pval, 
#                                          abs(avg_log2FC)>=0.5) # filter log2FC

genes_of_interest = c("Csf1", "Cxcr6") # list of genes to label (selectLab=)
EnhancedVolcano(cluster2vT_vp.markers,
                lab=rownames(cluster2vT_vp.markers),
                x='avg_log2FC',
                y='p_val_adj',
                title="Trm vs Other T",
                pCutoff=0.01,
                FCcutoff=0.5,
                pointSize=3.0,
                labSize=5.0,
                selectLab=NULL)

# cxcr6_deg <- read.csv("~/Downloads/Cluster2_vs_Cluster3.csv")
# p1 <- EnhancedVolcano(cxcr6_deg,
#                 lab=cxcr6_deg$X,
#                 x='log2FoldChange',
#                 y='padj',
#                 title="CXCR6 Trm vs Other T",
#                 pCutoff=0.01,
#                 FCcutoff=0.5,
#                 pointSize=3.0,
#                 labSize=5.0,
#                 selectLab=NULL)
# p1 +
#   ggplot2::coord_cartesian(xlim=c(-8, 8)) +
#   ggplot2::scale_x_continuous(breaks=seq(-8, 8, 2))
