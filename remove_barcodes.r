# This script removes cells with given barcodes from a Seurat object

#replace sub1_scobj with whatever object you want to remove from, or readRDS()
rmbc_scobj <- sub1_scobj

barcodes_to_remove = as.character(read.csv("clust7cells.csv")$X7)

rmbc_scobj <- rmbc_scobj[, !(rmbc_scobj@assays$RNA@data@Dimnames[[2]] 
                             %in% barcodes_to_remove)]

# This is old code that might be a useful reference later but for here is
# slower/less efficient.
# # read .csv containing barcodes (or create it from another source)
# barcodes_to_remove_df = read.csv("clust7cells.csv")
# 
# # at this point, make sure the only column is the list of barcodes (no index)
# # this line may or may not be necessary depending on input format
# # ensure the column name after $ is correct in below line
# barcodes_to_remove_df$X <- NULL
# 
# # convert dataframe to class "character"
# # ensure the column name after $ is correct in all below lines
# barcodes_to_remove = character(length(barcodes_to_remove_df$X7))
# for(i in 1:length(barcodes_to_remove_df$X7)){
#   barcodes_to_remove[i] = as.character(barcodes_to_remove_df$X7[i])
# }

# # locate the barcodes in the data (find index)
# for(i in 1:length(rmbc_scobj@assays$RNA@data@Dimnames[[2]])){
#   if(rmbc_scobj@assays$RNA@data@Dimnames[[2]][i] %in% barcodes_to_remove){
#     print(i)
#   }
# }

# # remove all the cells with barcodes in barcodes_to_remove
# rmbc_scobj <- rmbc_scobj[, !(rmbc_scobj@assays$RNA@data@Dimnames[[2]] %in% barcodes_to_remove)]

# now, if wanted, put this data back into old object to re-analyze
# sub1_scobj <- rmbc_scobj
