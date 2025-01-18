# Recluster Seurat Object [CCBR] [scRNA-seq] (576fe688-c445-48c6-a40a-033da171c149): v22
ReclusterSO_CD8_GzmbposKi67neg <- function(FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg, MetadataTable_CD8T_Gzmb_Mki67_copied) {
    #image: png

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    
    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Basic Parameters:
    seurat_object = FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg
    meta <- MetadataTable_CD8T_Gzmb_Mki67_copied
    reduction <- "umap"

    #Old Clustering Column Parameters:
    old_columns_to_save <- c()
    prepend_text = "old"

    #Reclustering Parameters:
    number_of_pcs = 30
    cluster_resolution_low_range <- 0.2
    cluster_resolution_high_range <- 6
    cluster_resolution_range_bins <- 1
    
    #Advanced Parameters:
    seurat_object_filename <- "seurat_object.rds"

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##
    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    ## Input SO.
    cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:     path <- nidapGetPath(seurat_object,seurat_object_filename)
    SO <- seurat_object

    ## Recluster the SO.
    reclustered_SO <- reclusterSeuratObject(object = SO,
                                            prepend.txt = prepend_text,
                                            old.columns.to.save = old_columns_to_save,
                                            number.of.pcs = number_of_pcs,
                                            cluster.resolution.low.range = cluster_resolution_low_range,
                                            cluster.resolution.high.range = cluster_resolution_high_range,
                                            cluster.resolution.range.bins = cluster_resolution_range_bins,
                                            reduction.type = reduction
                                            )

    ## Print reclustered dimensionality reduction plot.
    print(reclustered_SO$plot)

    ## Output reclustered SO.
# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(reclustered_SO$object)
    return(output_fs)
}

######### Node Execution Steps ##########
print("template_ReclusterSO_CD8_GzmbposKi67neg.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg
var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg<-readRDS(paste0(rds_output,"/var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg.rds"))

if (!('Seurat' %in% class(var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg))) { if (!(class(var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg <- as.data.frame(var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg)}}
#############################


# Processing input variable: var_MetadataTable_CD8T_Gzmb_Mki67_copied
var_MetadataTable_CD8T_Gzmb_Mki67_copied<-readRDS(paste0(rds_output,"/var_MetadataTable_CD8T_Gzmb_Mki67_copied.rds"))

if (!('Seurat' %in% class(var_MetadataTable_CD8T_Gzmb_Mki67_copied))) { if (!(class(var_MetadataTable_CD8T_Gzmb_Mki67_copied) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_CD8T_Gzmb_Mki67_copied <- as.data.frame(var_MetadataTable_CD8T_Gzmb_Mki67_copied)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ReclusterSO_CD8_GzmbposKi67neg<-ReclusterSO_CD8_GzmbposKi67neg(var_FilteredSO_Duallabel_CD8T_Gzmb_Mki67neg,var_MetadataTable_CD8T_Gzmb_Mki67_copied)
saveRDS(var_ReclusterSO_CD8_GzmbposKi67neg, paste0(rds_output,"/var_ReclusterSO_CD8_GzmbposKi67neg.rds"))
