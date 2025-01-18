# Rename Clusters by Cell Type [CCBR] [scRNA-seq] (7a47246e-1fa7-4400-9393-f22849758f92): v51
NameClustersSO <- function(ModScoreSO, MetadataTable_ModScore, Cluster_Identities_table) {

    ## This function maps seurat cluster ids to annotated cell types

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
        
    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Primary Inputs:
    seurat_object <- ModScoreSO  
    cluster_identities_table <- Cluster_Identities_table

    #Basic Parameters:
    cluster_number_column <- "My_ClusterNumbers"
    cluster_name_column <- "My_ClusterNames"
    cluster_column_from_metadata <- "orig_ident"
    celltype_column_from_metadata <- "MS_Celltype"

    #Advanced Parameters:
    order_clusters_by <- c()
    order_celltypes_by <- c()
	interactive <- FALSE
    seurat_object_filename <- "seurat_object.rds"
		
	## -------------------------------- ##
	## Parameter Misspecifation Errors  ##
	## -------------------------------- ##
	## -------------------------------- ##
	## Functions                        ##
	## -------------------------------- ##
	## --------------- ##
	## Main Code Block ##
	## --------------- ##

	# load data
# auto removed: 	path <- nidapGetPath(seurat_object,seurat_object_filename)
    SO <- seurat_object
	
    #Set up the parameters 

    cluster_numbers <- cluster_identities_table[[cluster_number_column]]
    cluster_names <- cluster_identities_table[[cluster_name_column]]
    
    #Change cluster column to match SO (NIDAPism)
    metadata_df <- SO@meta.data
	colnames(metadata_df) <- gsub("\\.", "_", colnames(metadata_df))
	SO.renamed <- SO
    SO.renamed@meta.data <- metadata_df 

############## A new ways of handling labels
library(Seurat)
library(dplyr)
# Prefix
Prefix <- "Cluster_"

# Check for entries in cluster_number_column that are not in the Seurat object and print a warning if there are missing clusters
missing_clusters <- setdiff(cluster_numbers, unique(metadata_df[[cluster_column_from_metadata]]))
if (length(missing_clusters) > 0) { warning("The following cluster numbers are not found in the Seurat object: ", paste(missing_clusters, collapse = ", ")) }

# Create a named vector for easy lookup of new cluster names
new_cluster_names <- setNames(cluster_names, cluster_numbers)
existing_clusters <- metadata_df[[cluster_column_from_metadata]]

# Adding new 'Clusternames' column to metadata.df
metadata_new <- metadata_df %>%
  mutate(Clusternames = ifelse(existing_clusters %in% cluster_numbers,
                               new_cluster_names[as.character(existing_clusters)],
                               paste0(Prefix, existing_clusters)))

cluster_numbers <- unique(metadata_new[[cluster_column_from_metadata]])
cluster_names <- unique(metadata_new$Clusternames)
##############

    so_result <- nameClusters(object = SO.renamed,
                        cluster.numbers = cluster_numbers,
                        cluster.names = cluster_names,
                        cluster.column = cluster_column_from_metadata,
                        labels.column = celltype_column_from_metadata,
                        order.clusters.by = order_clusters_by,
                        order.celltypes.by = order_celltypes_by,
                        interactive = interactive)

	
    print(so_result$plot)

    #Add new metadata column to original so:
	SO@meta.data$Clusternames <- so_result$object@meta.data$Clusternames

    #Print figure and save output file
# auto removed: 	output <- new.output()
# auto removed: 	output_fs <- output$fileSystem()
return(SO)
	return(SO)
}

######### Node Execution Steps ##########
print("template_NameClustersSO.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ModScoreSO
var_ModScoreSO<-readRDS(paste0(rds_output,"/var_ModScoreSO.rds"))

if (!('Seurat' %in% class(var_ModScoreSO))) { if (!(class(var_ModScoreSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ModScoreSO <- as.data.frame(var_ModScoreSO)}}
#############################


# Processing input variable: var_MetadataTable_ModScore
var_MetadataTable_ModScore<-readRDS(paste0(rds_output,"/var_MetadataTable_ModScore.rds"))

if (!('Seurat' %in% class(var_MetadataTable_ModScore))) { if (!(class(var_MetadataTable_ModScore) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_ModScore <- as.data.frame(var_MetadataTable_ModScore)}}
#############################


# Processing input variable: var_Cluster_Identities_table
var_Cluster_Identities_table<-readRDS(paste0(rds_output,"/var_Cluster_Identities_table.rds"))

if (!('Seurat' %in% class(var_Cluster_Identities_table))) { if (!(class(var_Cluster_Identities_table) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_Cluster_Identities_table <- as.data.frame(var_Cluster_Identities_table)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_NameClustersSO<-NameClustersSO(var_ModScoreSO,var_MetadataTable_ModScore,var_Cluster_Identities_table)
saveRDS(var_NameClustersSO, paste0(rds_output,"/var_NameClustersSO.rds"))
