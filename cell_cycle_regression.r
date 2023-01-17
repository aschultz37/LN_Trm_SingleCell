library(Seurat)

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
RidgePlot(scobj_cc, features=c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol=2)

# do regression
scobj_cc <- ScaleData(scobj_cc, vars.to.regress=c("S.Score", "G2M.Score"), 
                      features=rownames(scobj_cc))

saveRDS(scobj_cc, "LN/output/2022-12-29_LN_cc.rds")