# This file has includes and static data used across LN analysis
# Run this script prior to any LN_analysis*.r files, cell_cycle_regression.r,
# and volcano_plot.r.

library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)

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