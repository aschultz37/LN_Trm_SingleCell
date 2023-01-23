library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)

sub2_custom_color_palette = c("#1B9E77", "#7570B3", "#66A61E", 
                         "#E6AB02", "#A6761D", "#666666")
sub2_genelist = c("Itgae", "Ccr7", "Klf2", "Cxcr6", "S1pr1", 
                  "Cd8a", "Thy1", "Ptprc", "Cd3e")
LN_Trm_genes <- c("Acap1", "Actn2", "Arhgef1", "Atxn7l3b",
                  "B4galnt1", "Bcl11b", "Cbx3", "Ccnd2", "Ccr10", "Cd27", "Cd7",
                  "Cd74", "Chd3", "Cirbp", "Clec2d", "Crot", "Csf1", "Cxcr3",
                  "Cxcr6", "Eif5", "Evl", "Fam189b", "Fubp1", "Fyb", "Gramd1a",
                  "Sema4a","Gstp1", "Shisa5", "H2-T23", "Sipa1", "Hmgn1",
                  "Slfn2", "Sp100", "Hsp90b1", "Spcs2", "Id2", "Srrm2",
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
traffic_gene_list = c("S1pr1", "Sell", "Ccr7", "Cxcr4", "Cxcr3", "Cxcr6")
sub2_clusters_of_interest = c("0", "2")

# in a new object that is a copy of sub1_scobj, combine clusters 0,1,3
sub2_scobj <- readRDS("LN/output/2022-12-29_LN_sub1.rds")
sub2_cluster_ids <- c("0", "0", "2", "0", "4", "5")
names(sub2_cluster_ids) <- levels(sub2_scobj)
sub2_scobj <- RenameIdents(sub2_scobj, sub2_cluster_ids)

# label the individual clusters
DimPlot(sub2_scobj, reduction="umap", pt.size=2.5, cols=sub2_custom_color_palette)
# print number of cells per cluster
table(Idents(sub2_scobj))

# find top 10 most variable genes
sub2_top10 <- head(VariableFeatures(sub2_scobj), 10)

# find all markers
sub2_scobj.markers <- FindAllMarkers(sub2_scobj, only.pos=TRUE, min.pct=0.25,
                                     logfc.threshold=0.25)

# violin plots
VlnPlot(sub2_scobj, features=sub2_genelist, cols=sub2_custom_color_palette,
        idents=NULL)
VlnPlot(sub2_scobj, features=cytotoxic_gene_list, cols=sub2_custom_color_palette,
        idents=sub2_clusters_of_interest)
VlnPlot(sub2_scobj, features=traffic_gene_list, cols=sub2_custom_color_palette,
        idents=sub2_clusters_of_interest)

# generate expression heatmap for top 10 markers for each cluster
sub2_scobj.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> sub2_top10
DoHeatmap(sub2_scobj, features=sub2_top10$gene, cols=sub2_custom_color_palette) + NoLegend()

# find all markers distinguishing cluster Trm and 0
cluster2v0_sub2.markers <- FindMarkers(sub2_scobj, ident.1=2, ident.2=0, 
                                       min.pct=0.25)
cluster2v0_sub2.markers.sig <- subset(cluster2v0_sub2.markers, p_val_adj<=0.05)
write.csv(cluster2v0_sub2.markers.sig, "LN/output/cluster2v0_sub2.csv", 
          row.names=TRUE)

saveRDS(sub2_scobj, "LN/output/2022-12-29_LN_sub2.rds")