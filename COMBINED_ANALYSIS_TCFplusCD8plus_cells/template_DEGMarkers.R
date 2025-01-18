# DEG with Find Markers [scRNA-seq][CCBR] (54b6dd44-e233-4fb8-9bae-6ab0cb46e399): v93
DEGMarkers <- function(ReclusteredSO,SampleNames_Reclustered,MetadataTable_Reclustered) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    seurat_object <- ReclusteredSO
    metadata_table <- MetadataTable_Reclustered
    samples <- 'c("LN_aPD1","LN_ENT","LN_ENT+aPD1","LN_ENT+N803","LN_ENT+N803+aPD1","LN_N803","LN_N803+aPD1","LN_PBS","T_aPD1","T_ENT","T_ENT+aPD1","T_ENT+N803","T_ENT+N803+aPD1","T_N803","T_N803+aPD1","T_PBS")'
    parameter_to_test <- "SCT_snn_res_0_4"
    contrasts <- c("0-1","0-2","0-3")
    test_to_use <- "MAST"
    log_fc_threshold <- 0
    use_spark <- FALSE
    assay_to_use <- "SCT"
    use_log_2 <- TRUE

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##
    
    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

   ## Load SO 
        cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:         fs <- seurat_object$fileSystem()
# auto removed:         path <- fs$get_path("seurat_object.rds", 'r')
        SO <- seurat_object
        print(SO)

results <- degGeneExpressionMarkers(object = SO,
                                     samples = samples,
                                     contrasts = contrasts,
                                     parameter.to.test = parameter_to_test,
                                     test.to.use = test_to_use,
                                     log.fc.threshold = log_fc_threshold,
                                     use.spark = use_spark,
                                     assay.to.use = assay_to_use)

    return(results$df)

}

######### Node Execution Steps ##########
print("template_DEGMarkers.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ReclusteredSO
var_ReclusteredSO<-readRDS(paste0(rds_output,"/var_ReclusteredSO.rds"))

if (!('Seurat' %in% class(var_ReclusteredSO))) { if (!(class(var_ReclusteredSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ReclusteredSO <- as.data.frame(var_ReclusteredSO)}}
#############################


# Processing input variable: var_SampleNames_Reclustered
var_SampleNames_Reclustered<-readRDS(paste0(rds_output,"/var_SampleNames_Reclustered.rds"))

if (!('Seurat' %in% class(var_SampleNames_Reclustered))) { if (!(class(var_SampleNames_Reclustered) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_SampleNames_Reclustered <- as.data.frame(var_SampleNames_Reclustered)}}
#############################


# Processing input variable: var_MetadataTable_Reclustered
var_MetadataTable_Reclustered<-readRDS(paste0(rds_output,"/var_MetadataTable_Reclustered.rds"))

if (!('Seurat' %in% class(var_MetadataTable_Reclustered))) { if (!(class(var_MetadataTable_Reclustered) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_Reclustered <- as.data.frame(var_MetadataTable_Reclustered)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_DEGMarkers<-DEGMarkers(var_ReclusteredSO,var_SampleNames_Reclustered,var_MetadataTable_Reclustered)
saveRDS(var_DEGMarkers, paste0(rds_output,"/var_DEGMarkers.rds"))
