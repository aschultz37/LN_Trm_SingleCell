# Run LN_analysis_init.r first!

# Read in object with finalized data
proc_scobj <- readRDS("LN/output/2022-12-29_LN_sub1_rmbc_figs.rds")

# Copy of object that combines 0, 1, 3, and 5 into cluster 6
comb_scobj <- proc_scobj
comb_scobj_ids <- c("6", "6", "2", "6", "4", "6")
names(comb_scobj_ids) <- levels(comb_scobj)
comb_scobj <- RenameIdents(comb_scobj, comb_scobj_ids)

# UMAP
DimPlot(proc_scobj, reduction="umap", pt.size=2.5, cols=custom_color_palette)
table(Idents(proc_scobj))

# LN Trm ModuleScore
proc_scobj <- AddModuleScore(object=proc_scobj, features=list(LN_Trm_genes),
                             ctrl=100, name="LN_Trm")
FeaturePlot(proc_scobj, features="LN_Trm1", pt.size=2.5, 
            cols=rev(brewer.pal(n=11, name="RdBu")))

# TGFb ModuleScore
proc_scobj <- AddModuleScore(object=proc_scobj, features=list(TGFb_genes),
                             ctrl=100, name="TGFb")
FeaturePlot(proc_scobj, features="TGFb1", pt.size=2.5, 
            cols=rev(brewer.pal(n=11, name="RdBu")))

# Violin Plots
# Trm vs Other T, excl. cycling
clusters_of_interest <- c("2", "6") # exclude cycling cells in 4
violin_colors <- c("#FFFFFF", "#7570B3")
VlnPlot(comb_scobj, features=c("Gzmb", "Gzmk", "Prf1"), 
        cols=violin_colors, idents=clusters_of_interest)
VlnPlot(comb_scobj, features=c("Csf1", "Ccl4", "Ccl5"), 
        cols=violin_colors, idents=clusters_of_interest)
VlnPlot(comb_scobj, features=c("Klf2", "S1pr1", "Tcf7"), 
        cols=violin_colors, idents=clusters_of_interest)

# Feature Plots
# export individually
feature_gene_list <- c("Ccr7", "Sell", "S1pr1", "Gzmb",
                       "Cxcr6", "Itgae", "Csf1")
for(i in 1:length(feature_gene_list)){
  print(FeaturePlot(proc_scobj, features=feature_gene_list[i], pt.size=2.5))
}

# Heatmap
proc_scobj.markers <- FindAllMarkers(proc_scobj, only.pos=TRUE, min.pct=0.25,
                                     logfc.threshold=0.25)
proc_scobj.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> proc_top10
DoHeatmap(proc_scobj, features=proc_top10$gene, 
          group.colors=custom_color_palette) + NoLegend()
