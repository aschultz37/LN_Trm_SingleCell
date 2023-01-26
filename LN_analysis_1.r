# This file removes clusters from initial analysis that are not real/artifacts.
# 4, 7, and 8 are removed.

# library(dplyr)
# library(Seurat)
# library(patchwork)

# REMOVE CLUSTERS AND RE-ANALYZE
sub1_scobj = readRDS("LN/output/2022-12-29_LN_sub1.rds")
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
DimPlot(sub1_scobj, reduction="umap", pt.size=2.5, cols=custom_color_palette)
# print number of cells per cluster
table(Idents(sub1_scobj))

# find all markers
sub1_scobj.markers <- FindAllMarkers(sub1_scobj, only.pos=TRUE, min.pct=0.25,
                                     logfc.threshold=0.25)

# expression probability distributions across clusters
clusters_of_interest = c("0", "1", "2", "3", "5") # exclude cycling cells in 4
VlnPlot(sub1_scobj, features=gen_gene_list, cols=custom_color_palette,
        idents=NULL)
VlnPlot(sub1_scobj, features=cytotoxic_gene_list, cols=custom_color_palette,
        idents=NULL)
VlnPlot(sub1_scobj, features=traffic_gene_list, cols=custom_color_palette,
        idents=NULL)

# visualize feature expression on tSNE or PCA plot
FeaturePlot(sub1_scobj, features=gen_gene_list)
FeaturePlot(sub1_scobj, features=cytotoxic_gene_list)
FeaturePlot(sub1_scobj, features=traffic_gene_list)

# visualize LN_Trm genes in feature plots
# FeaturePlot for list of genes of arbitrary length (groups of 12)
# final_index <- 0
# for(i in 1:length(LN_Trm_genes)){
#   if((i %% 12) == 0){
#     print(FeaturePlot(sub1_scobj, features=LN_Trm_genes[(i-11):i]))
#     final_index <- i
#   }
# }
# FeaturePlot(sub1_scobj, 
#             features=LN_Trm_genes[(final_index+1):length(LN_Trm_genes)])


# to detect naive T cells
# FeaturePlot(sub1_scobj, pt.size=2.5, features="Cd44")

# N.B. can also try RidgePlot, CellScatter, and DotPlot to view dataset
# ridge plot
RidgePlot(sub1_scobj, features=gen_gene_list, cols=custom_color_palette, ncol=3)

# generate expression heatmap for top 10 markers for each cluster
sub1_scobj.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> sub1_top10
DoHeatmap(sub1_scobj, features=sub1_top10$gene, 
          group.colors=custom_color_palette) + NoLegend()

# find all markers distinguishing cluster 0 and 1
cluster0v1_sub1.markers <- FindMarkers(sub1_scobj, ident.1=0, ident.2=c(1), 
                                       min.pct=0.25)
# cluster0v1_sub1.markers.sig <- subset(cluster0v1_sub1.markers, p_val_adj<=0.05)
write.csv(cluster0v1_sub1.markers, "LN/output/cluster0v1_sub1.csv", 
          row.names=TRUE)

# find all markers distinguishing cluster 2 from other T cells (0, 1, 3, 5)
cluster2vTcirc_sub1.markers <- FindMarkers(sub1_scobj, ident.1=2, 
                                        ident.2=c(0, 1, 3, 5), 
                                        min.pct=0.25)
# cluster2vTcirc_sub1.markers.sig <- subset(cluster2vTcirc_sub1.markers, 
#                                           p_val_adj<=0.05)
write.csv(cluster2vTcirc_sub1.markers, "LN/output/cluster2vTcirc_sub1.csv", 
          row.names=TRUE)

# calculate module score for feature expression
sub1_scobj <- AddModuleScore(object=sub1_scobj, features=list(LN_Trm_genes),
                             ctrl=100, name="LN_Trm")
FeaturePlot(sub1_scobj, features="LN_Trm1", pt.size=2.5, 
            cols=rev(brewer.pal(n=11, name="RdBu")))

saveRDS(sub1_scobj, "LN/output/2022-12-29_LN_sub1.rds")
