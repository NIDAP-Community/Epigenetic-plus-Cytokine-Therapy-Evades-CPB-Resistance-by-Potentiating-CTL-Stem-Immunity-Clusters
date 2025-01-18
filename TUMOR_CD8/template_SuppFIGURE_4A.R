# Dot Plot of Genes by Metadata [CCBR] [scRNA-seq] (79573b27-8a93-4f22-9863-993be1a44fc1): v37
SuppFIGURE_4A <- function(FilteredSO_CD8CD3,NewDataset_7,MetadataTable_CD8cd3) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    
    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Primary Inputs:
    input_dataset <- FilteredSO_CD8CD3
    metadata_table <- MetadataTable_CD8cd3
    category_labels_and_genes_table <-  NewDataset_7

    #Category Parameters:
    metadata_category_to_plot <- "orig_ident"
    category_labels_to_plot <- "My_ClusterNames"
    genes_column <- "My_GeneNames"

    #Output Parameters:
    return_percentage_cells_expressing <- TRUE
    seurat_object_filename <- "seurat_object.rds"

    #Visualization Parameters:
    reverse_plot <- TRUE
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
print("template_SuppFIGURE_4A.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_FilteredSO_CD8CD3
var_FilteredSO_CD8CD3<-readRDS(paste0(rds_output,"/var_FilteredSO_CD8CD3.rds"))

if (!('Seurat' %in% class(var_FilteredSO_CD8CD3))) { if (!(class(var_FilteredSO_CD8CD3) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_FilteredSO_CD8CD3 <- as.data.frame(var_FilteredSO_CD8CD3)}}
#############################


# Processing input variable: var_NewDataset_7
var_NewDataset_7<-readRDS(paste0(rds_output,"/var_NewDataset_7.rds"))

if (!('Seurat' %in% class(var_NewDataset_7))) { if (!(class(var_NewDataset_7) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_NewDataset_7 <- as.data.frame(var_NewDataset_7)}}
#############################


# Processing input variable: var_MetadataTable_CD8cd3
var_MetadataTable_CD8cd3<-readRDS(paste0(rds_output,"/var_MetadataTable_CD8cd3.rds"))

if (!('Seurat' %in% class(var_MetadataTable_CD8cd3))) { if (!(class(var_MetadataTable_CD8cd3) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_CD8cd3 <- as.data.frame(var_MetadataTable_CD8cd3)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_SuppFIGURE_4A<-SuppFIGURE_4A(var_FilteredSO_CD8CD3,var_NewDataset_7,var_MetadataTable_CD8cd3)
saveRDS(var_SuppFIGURE_4A, paste0(rds_output,"/var_SuppFIGURE_4A.rds"))
