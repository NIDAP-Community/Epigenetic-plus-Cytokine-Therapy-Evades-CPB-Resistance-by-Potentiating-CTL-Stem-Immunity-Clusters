# Dot Plot of Genes by Metadata [CCBR] [scRNA-seq] (79573b27-8a93-4f22-9863-993be1a44fc1): v37
DotGenesClusters <- function(NameClustersSO,My_InputForDotplot,MetadataTable_NameClusters) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    
    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Primary Inputs:
    input_dataset <- NameClustersSO
    metadata_table <- MetadataTable_NameClusters
    category_labels_and_genes_table <-  My_InputForDotplot

    #Category Parameters:
    metadata_category_to_plot <- "Clusternames"
    category_labels_to_plot <- "My_ClusterNames"
    genes_column <- "My_GeneNames"

    #Output Parameters:
    return_percentage_cells_expressing <- TRUE
    seurat_object_filename <- "seurat_object.rds"

    #Visualization Parameters:
    reverse_plot <- FALSE
    reverse_categories <- FALSE
    dot_color <- "darkblue"
    
    ## -------------------------------- ##
    ## Parameter Misspecifation Errors  ##
    ## -------------------------------- ##

    if(class(input_dataset) != "FoundryTransformInput"){
# auto removed:         stop("Input should be Seurat object in rds file format. Rerun previous step with latest released version to produce correct input format for Seurat Object")
    }

    if(sum(category_labels_and_genes_table[[category_labels_to_plot]] %in% metadata_table[[metadata_category_to_plot]]) == 0){
        stop(paste0("At least some category element in input table column: ", category_labels_to_plot, " should match Seurat Object metadata category column: ",metadata_category_to_plot))
    }

    ## -------------------------------- ##
    ## Functions                        ##
    ## -------------------------------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

# auto removed:     path <- nidapGetPath(input_dataset,seurat_object_filename)
    so <- input_dataset

    #Change metadata tables to match NIDAP:
    metadata.df <- so@meta.data
    colnames(metadata.df) <- gsub("\\.","_",colnames(metadata.df))
    so@meta.data <- metadata.df

    cells <- category_labels_and_genes_table[[category_labels_to_plot]]
    cells <- as.factor(cells[!is.na(cells)])

### Temporary fix - ignore extra labels:
library(Seurat)
ExtraValue <- sum(!cells %in% unique(metadata.df[[metadata_category_to_plot]]))
if (ExtraValue > 0) {
         missinglab2 <- cells[!cells %in% unique(metadata.df[[metadata_category_to_plot]])]
         warning(sprintf("There are %s additional elements in your input categories\n that are missing from your metadata table: ", ExtraValue))
         missinglab2 <- cat(paste(as.character(missinglab2), collapse = "\n"))
#
    cells <- cells[cells %in% unique(metadata.df[[metadata_category_to_plot]])]
#
    }
### End of the fix  

    genes <- category_labels_and_genes_table[[genes_column]]
    genes <- genes[!is.na(genes)]

    results <- dotPlotMet(object = so,
                       metadata = metadata_category_to_plot,
                       cells = cells,
                       markers = genes,
                       plot.reverse = reverse_plot,
                       cell.reverse.sort = reverse_categories,
                       dot.color = dot_color)

    print(results$plot)

    if(return_percentage_cells_expressing == TRUE){
        return(results$pct)
    } else {
        return(results$exp)
    }
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

######### Node Execution Steps ##########
print("template_DotGenesClusters.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_NameClustersSO
var_NameClustersSO<-readRDS(paste0(rds_output,"/var_NameClustersSO.rds"))

if (!('Seurat' %in% class(var_NameClustersSO))) { if (!(class(var_NameClustersSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_NameClustersSO <- as.data.frame(var_NameClustersSO)}}
#############################


# Processing input variable: var_My_InputForDotplot
var_My_InputForDotplot<-readRDS(paste0(rds_output,"/var_My_InputForDotplot.rds"))

if (!('Seurat' %in% class(var_My_InputForDotplot))) { if (!(class(var_My_InputForDotplot) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_My_InputForDotplot <- as.data.frame(var_My_InputForDotplot)}}
#############################


# Processing input variable: var_MetadataTable_NameClusters
var_MetadataTable_NameClusters<-readRDS(paste0(rds_output,"/var_MetadataTable_NameClusters.rds"))

if (!('Seurat' %in% class(var_MetadataTable_NameClusters))) { if (!(class(var_MetadataTable_NameClusters) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_NameClusters <- as.data.frame(var_MetadataTable_NameClusters)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_DotGenesClusters<-DotGenesClusters(var_NameClustersSO,var_My_InputForDotplot,var_MetadataTable_NameClusters)
saveRDS(var_DotGenesClusters, paste0(rds_output,"/var_DotGenesClusters.rds"))
