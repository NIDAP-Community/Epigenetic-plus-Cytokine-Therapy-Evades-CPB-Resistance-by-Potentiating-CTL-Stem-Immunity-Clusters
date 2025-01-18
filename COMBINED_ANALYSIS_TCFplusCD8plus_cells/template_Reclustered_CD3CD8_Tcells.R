# Recluster Seurat Object [CCBR] [scRNA-seq] (576fe688-c445-48c6-a40a-033da171c149): v22
Reclustered_CD3CD8_Tcells <- function(FSO_CD3CD8_Tcells, MDT_FSO_CD3CD8_Tcells) {
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
    seurat_object = FSO_CD3CD8_Tcells
    meta <- MDT_FSO_CD3CD8_Tcells
    reduction <- "umap"

    #Old Clustering Column Parameters:
    old_columns_to_save <- c("SCT_snn_res_0_2","SCT_snn_res_0_4","SCT_snn_res_0_6","SCT_snn_res_0_8","SCT_snn_res_1","SCT_snn_res_1_2")
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
print("template_Reclustered_CD3CD8_Tcells.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_FSO_CD3CD8_Tcells
var_FSO_CD3CD8_Tcells<-readRDS(paste0(rds_output,"/var_FSO_CD3CD8_Tcells.rds"))

if (!('Seurat' %in% class(var_FSO_CD3CD8_Tcells))) { if (!(class(var_FSO_CD3CD8_Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_FSO_CD3CD8_Tcells <- as.data.frame(var_FSO_CD3CD8_Tcells)}}
#############################


# Processing input variable: var_MDT_FSO_CD3CD8_Tcells
var_MDT_FSO_CD3CD8_Tcells<-readRDS(paste0(rds_output,"/var_MDT_FSO_CD3CD8_Tcells.rds"))

if (!('Seurat' %in% class(var_MDT_FSO_CD3CD8_Tcells))) { if (!(class(var_MDT_FSO_CD3CD8_Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MDT_FSO_CD3CD8_Tcells <- as.data.frame(var_MDT_FSO_CD3CD8_Tcells)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_Reclustered_CD3CD8_Tcells<-Reclustered_CD3CD8_Tcells(var_FSO_CD3CD8_Tcells,var_MDT_FSO_CD3CD8_Tcells)
saveRDS(var_Reclustered_CD3CD8_Tcells, paste0(rds_output,"/var_Reclustered_CD3CD8_Tcells.rds"))
