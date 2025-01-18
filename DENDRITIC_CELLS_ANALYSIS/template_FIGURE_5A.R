# Dotplot of Gene Expression by Metadata [scRNA-seq][CCBR] (79573b27-8a93-4f22-9863-993be1a44fc1): v30
FIGURE_5A <- function(FilteredSO_Clec9aItgaeDP,NewDataset_2,Filteredmetadata_Clec9aItgaeDP) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    
    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Primary Inputs:
    input_dataset <- FilteredSO_Clec9aItgaeDP
    metadata_table <- Filteredmetadata_Clec9aItgaeDP
    category_labels_and_genes_table <-  NewDataset_2

    #Basic Parameters:
    metadata_category_to_plot <- "orig_ident"
    category_labels_to_plot <- "My_ClusterNames"
    genes_column <- "My_GeneNames"
    return_percentage_cells_expressing <- TRUE
    
    #Visualization Parameters:
    reverse_plot <- TRUE
    reverse_categories <- FALSE
    dot_color <- "darkred"
    seurat_object_filename <- "seurat_object.rds"
    
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
print("template_FIGURE_5A.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_FilteredSO_Clec9aItgaeDP
var_FilteredSO_Clec9aItgaeDP<-readRDS(paste0(rds_output,"/var_FilteredSO_Clec9aItgaeDP.rds"))

if (!('Seurat' %in% class(var_FilteredSO_Clec9aItgaeDP))) { if (!(class(var_FilteredSO_Clec9aItgaeDP) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_FilteredSO_Clec9aItgaeDP <- as.data.frame(var_FilteredSO_Clec9aItgaeDP)}}
#############################


# Processing input variable: var_NewDataset_2
var_NewDataset_2<-readRDS(paste0(rds_output,"/var_NewDataset_2.rds"))

if (!('Seurat' %in% class(var_NewDataset_2))) { if (!(class(var_NewDataset_2) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_NewDataset_2 <- as.data.frame(var_NewDataset_2)}}
#############################


# Processing input variable: var_Filteredmetadata_Clec9aItgaeDP
var_Filteredmetadata_Clec9aItgaeDP<-readRDS(paste0(rds_output,"/var_Filteredmetadata_Clec9aItgaeDP.rds"))

if (!('Seurat' %in% class(var_Filteredmetadata_Clec9aItgaeDP))) { if (!(class(var_Filteredmetadata_Clec9aItgaeDP) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_Filteredmetadata_Clec9aItgaeDP <- as.data.frame(var_Filteredmetadata_Clec9aItgaeDP)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_FIGURE_5A<-FIGURE_5A(var_FilteredSO_Clec9aItgaeDP,var_NewDataset_2,var_Filteredmetadata_Clec9aItgaeDP)
saveRDS(var_FIGURE_5A, paste0(rds_output,"/var_FIGURE_5A.rds"))
