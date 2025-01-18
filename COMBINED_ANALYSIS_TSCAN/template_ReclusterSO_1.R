# Recluster Seurat Object [CCBR] [scRNA-seq] (576fe688-c445-48c6-a40a-033da171c149): v22
ReclusterSO_1 <- function(SubsettedSO, MetadataTable_7) {
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
    seurat_object = SubsettedSO
    meta <- MetadataTable_7
    reduction <- "umap"

    #Old Clustering Column Parameters:
    old_columns_to_save <- c()
    prepend_text = "old"

    #Reclustering Parameters:
    number_of_pcs = 30
    cluster_resolution_low_range <- 0.2
    cluster_resolution_high_range <- 1.0
    cluster_resolution_range_bins <- 0.2
    
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
print("template_ReclusterSO_1.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_SubsettedSO
var_SubsettedSO<-readRDS(paste0(rds_output,"/var_SubsettedSO.rds"))

if (!('Seurat' %in% class(var_SubsettedSO))) { if (!(class(var_SubsettedSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_SubsettedSO <- as.data.frame(var_SubsettedSO)}}
#############################


# Processing input variable: var_MetadataTable_7
var_MetadataTable_7<-readRDS(paste0(rds_output,"/var_MetadataTable_7.rds"))

if (!('Seurat' %in% class(var_MetadataTable_7))) { if (!(class(var_MetadataTable_7) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_7 <- as.data.frame(var_MetadataTable_7)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ReclusterSO_1<-ReclusterSO_1(var_SubsettedSO,var_MetadataTable_7)
saveRDS(var_ReclusterSO_1, paste0(rds_output,"/var_ReclusterSO_1.rds"))
