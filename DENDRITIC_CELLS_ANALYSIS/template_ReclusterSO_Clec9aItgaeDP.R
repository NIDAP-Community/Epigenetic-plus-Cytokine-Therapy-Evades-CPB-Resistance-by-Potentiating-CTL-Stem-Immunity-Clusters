# Recluster Seurat Object [scRNA-seq][CCBR] (576fe688-c445-48c6-a40a-033da171c149): v19
ReclusterSO_Clec9aItgaeDP <- function(FilteredSO_Clec9aItgaeDP, Filteredmetadata_Clec9aItgaeDP) {
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
    seurat_object = FilteredSO_Clec9aItgaeDP
    meta <- Filteredmetadata_Clec9aItgaeDP
    reduction <- "umap"

    #Old Clustering Column Parameters:
    old_columns_to_save <- c()
    prepend_text = "old"

    #Reclustering Parameters:
    number_of_pcs = 30
    cluster_resolution_low_range <- 0.2
    cluster_resolution_high_range <- 1.2
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
print("template_ReclusterSO_Clec9aItgaeDP.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_FilteredSO_Clec9aItgaeDP
var_FilteredSO_Clec9aItgaeDP<-readRDS(paste0(rds_output,"/var_FilteredSO_Clec9aItgaeDP.rds"))

if (!('Seurat' %in% class(var_FilteredSO_Clec9aItgaeDP))) { if (!(class(var_FilteredSO_Clec9aItgaeDP) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_FilteredSO_Clec9aItgaeDP <- as.data.frame(var_FilteredSO_Clec9aItgaeDP)}}
#############################


# Processing input variable: var_Filteredmetadata_Clec9aItgaeDP
var_Filteredmetadata_Clec9aItgaeDP<-readRDS(paste0(rds_output,"/var_Filteredmetadata_Clec9aItgaeDP.rds"))

if (!('Seurat' %in% class(var_Filteredmetadata_Clec9aItgaeDP))) { if (!(class(var_Filteredmetadata_Clec9aItgaeDP) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_Filteredmetadata_Clec9aItgaeDP <- as.data.frame(var_Filteredmetadata_Clec9aItgaeDP)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ReclusterSO_Clec9aItgaeDP<-ReclusterSO_Clec9aItgaeDP(var_FilteredSO_Clec9aItgaeDP,var_Filteredmetadata_Clec9aItgaeDP)
saveRDS(var_ReclusterSO_Clec9aItgaeDP, paste0(rds_output,"/var_ReclusterSO_Clec9aItgaeDP.rds"))
