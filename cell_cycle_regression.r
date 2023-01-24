# library(dplyr)
# library(Seurat)
# library(patchwork)

# Note: Assumes an existing Seurat object
# Object should already have normalized and scaled data, PCA, etc.
scobj_cc <- readRDS("LN/output/2022-12-29_LN_sub1.rds")

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- tolower(cc.genes$s.genes)
g2m.genes <- tolower(cc.genes$g2m.genes)
# correct the case to match the formatting of our genes (Abcd1)
for(i in 1:length(s.genes)){
  substr(s.genes[i],1,1) <- toupper(substr(s.genes[i], 1, 1))
}
for(i in 1:length(g2m.genes)){
  substr(g2m.genes[i],1,1) <- toupper(substr(g2m.genes[i], 1, 1))
}

# set.ident will switch to cell cycles instead of clusters
scobj_cc <- CellCycleScoring(scobj_cc, s.features=s.genes, 
                               g2m.features=g2m.genes, set.ident=FALSE)
RidgePlot(scobj_cc, features=c("Pcna", "Top2a", "Mcm6", "Mki67"), 
          cols=custom_color_palette, ncol=2)

# do regression
scobj_cc <- ScaleData(scobj_cc, vars.to.regress=c("S.Score", "G2M.Score"), 
                      features=rownames(scobj_cc))

saveRDS(scobj_cc, "LN/output/2022-12-29_LN_sub1_cc.rds")

# rerun PCA
scobj_cc <- RunPCA(scobj_cc, features=VariableFeatures(object=scobj_cc))

# determines # of PC to include (low p val); removes noise
scobj_cc <- JackStraw(scobj_cc, dims=30, num.replicate=100)
scobj_cc <- ScoreJackStraw(scobj_cc, dims=1:30)
# can visually determine dropoff in p-value after certain # PC
JackStrawPlot(scobj_cc, dims=1:30)
# alternatively, can use elbow plot to determine which # PC capture most signal
ElbowPlot(scobj_cc, ndims=30)

# recluster cells
scobj_cc <- FindNeighbors(scobj_cc, dims=1:25)     # choose dim based on PCA
scobj_cc <- FindClusters(scobj_cc, resolution=0.3) # resolution determines # clusters

# run umap
scobj_cc <- RunUMAP(scobj_cc, dims=1:25) #choose dim based on PCA/Neighbors
# note that you can set 'label = TRUE' or use LabelClusters function to
# label the individual clusters
DimPlot(scobj_cc, reduction="umap", pt.size=2.5, cols=custom_color_palette)
# print number of cells per cluster
table(Idents(scobj_cc))

# find all markers
scobj_cc.markers <- FindAllMarkers(scobj_cc, only.pos=TRUE, min.pct=0.25,
                                     logfc.threshold=0.25)

# violin plots
VlnPlot(scobj_cc, features=gen_gene_list, cols=custom_color_palette,
        idents=NULL)
VlnPlot(scobj_cc, features=cytotoxic_gene_list, cols=custom_color_palette,
        idents=NULL)
VlnPlot(scobj_cc, features=traffic_gene_list, cols=custom_color_palette,
        idents=NULL)

# feature plots
FeaturePlot(scobj_cc, features=gen_gene_list)
FeaturePlot(scobj_cc, features=cytotoxic_gene_list)
FeaturePlot(scobj_cc, features=traffic_gene_list)

# ridge plot
RidgePlot(scobj_cc, features=gen_gene_list, cols=custom_color_palette, ncol=3)

# generate expression heatmap for top 10 markers for each cluster
scobj_cc_top10 <- head(VariableFeatures(scobj_cc), 10)
scobj_cc.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> scobj_cc_top10
DoHeatmap(scobj_cc, features=scobj_cc_top10$gene, 
          group.colors=custom_color_palette) + NoLegend()

# find all markers distinguishing cluster 0 and 1
cluster0v1_cc.markers <- FindMarkers(scobj_cc, ident.1=0, ident.2=c(1), 
                                       min.pct=0.25)
cluster0v1_cc.markers.sig <- subset(cluster0v1_cc.markers, p_val_adj<=0.05)
write.csv(cluster0v1_cc.markers.sig, "LN/output/cluster0v1_cc.csv", 
          row.names=TRUE)

# find all markers distinguishing cluster 2 from 0 and 1
cluster2v01_cc.markers <- FindMarkers(scobj_cc, ident.1=2, ident.2=c(0, 1), 
                                        min.pct=0.25)
cluster2v01_cc.markers.sig <- subset(cluster2v01_cc.markers, p_val_adj<=0.05)
write.csv(cluster2v01_cc.markers.sig, "LN/output/cluster2v01_cc.csv", 
          row.names=TRUE)

# calculate module score for feature expression
scobj_cc <- AddModuleScore(object=scobj_cc, features=list(LN_Trm_genes),
                             ctrl=100, name="LN_Trm")
FeaturePlot(scobj_cc, features="LN_Trm1", pt.size=2.5, cols=brewer.pal(n=11, name="RdBu"))

saveRDS(scobj_cc, "LN/output/2022-12-29_LN_sub1_cc_figs.rds")