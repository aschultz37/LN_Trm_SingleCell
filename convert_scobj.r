library(scater)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(patchwork)

# Convert Seurat object to SingleCellExperiment
scobj.sce <- as.SingleCellExperiment(scobj)

# Convert SingleCellExperiment to Seurat object
sceobj.seurat <- as.Seurat(sceobj, counts="counts", data="logcounts")

# Convert Seurat object to loom
scobj.loom <- as.loom(scobj, filename="name.loom", verbose=FALSE)
# call when done with loom data
scobj.loom$close_all()