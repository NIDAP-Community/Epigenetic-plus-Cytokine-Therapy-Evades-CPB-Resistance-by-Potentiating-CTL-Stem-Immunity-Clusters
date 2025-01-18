# Cell Type Annotation with SingleR [CCBR] [scRNA-seq] (ba316b6a-e424-49ab-a8af-32af8c85529f): v132
CellTypesSO <- function(CombNormSO,CellDex_snapshotDate_2021_10_19, MetadataTable_CombNorm) {

## --------- ##
## Libraries ##
## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    library(grid)
    library(gridExtra)
    library(gridBase)
    library(cowplot)
    

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

# Primary Inputs:
CombinedSO <-CombNormSO
Metadata <- MetadataTable_CombNorm
RDS <- CellDex_snapshotDate_2021_10_19

# Basic Parameters:
species <- "Mouse"
useClusters <- FALSE
ClusterColumn <- c()[1]

# Visualization Parameters:
Reduction <- "umap"
Legend_Dot_Size <- 2

# Advanced Parameters:
doFineTuning <- FALSE
Number_of_cells_per_partition <- 400

## -------------------------------- ##
## Errors                           ##
## -------------------------------- ##

## -------------------------------- ##
## Functions                        ##
## -------------------------------- ##

 
## --------------- ##
## Main Code Block ##
## --------------- ##

# Loading Seurat object
    cat("1. Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:     fs <- CombinedSO$fileSystem()
# auto removed:     path <- fs$get_path("seurat_object.rds", 'r')
    so <- CombinedSO

# Loading CellDex annotation
#cell.dex <- list(HPCA, BP, mousernaseq, immgen)
    cat("2. Reading pre-saved CellDex Data from dataset: CellDexDatabase.rds\n\n")
# auto removed:     fs2 <- RDS$fileSystem()
# auto removed:     path2 <- fs2$get_path("CellDexDatabase.rds", 'r')
    CellDexData <- RDS

Reduction = "umap"
if (useClusters)
  {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData,
                            use.clusters = ClusterColumn)
  } else {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData)

  }
    gc()

## Print figures
    print(anno$p1)
    print(anno$p2)

Reduction = "tsne"
if (useClusters)
  {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData,
                            use.clusters = ClusterColumn)
  } else {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData)

  }
    gc()

## Print figures
    print(anno$p1)
    print(anno$p2)

## Save the annotated Seurat object
# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(anno$object)

# auto removed:     return(NULL)
}

######### Node Execution Steps ##########
print("template_CellTypesSO.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_CombNormSO
var_CombNormSO<-readRDS(paste0(rds_output,"/var_CombNormSO.rds"))

if (!('Seurat' %in% class(var_CombNormSO))) { if (!(class(var_CombNormSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_CombNormSO <- as.data.frame(var_CombNormSO)}}
#############################


# Processing input variable: var_CellDex_snapshotDate_2021_10_19
var_CellDex_snapshotDate_2021_10_19<-readRDS(paste0(rds_output,"/var_CellDex_snapshotDate_2021_10_19.rds"))

if (!('Seurat' %in% class(var_CellDex_snapshotDate_2021_10_19))) { if (!(class(var_CellDex_snapshotDate_2021_10_19) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_CellDex_snapshotDate_2021_10_19 <- as.data.frame(var_CellDex_snapshotDate_2021_10_19)}}
#############################


# Processing input variable: var_MetadataTable_CombNorm
var_MetadataTable_CombNorm<-readRDS(paste0(rds_output,"/var_MetadataTable_CombNorm.rds"))

if (!('Seurat' %in% class(var_MetadataTable_CombNorm))) { if (!(class(var_MetadataTable_CombNorm) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_CombNorm <- as.data.frame(var_MetadataTable_CombNorm)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_CellTypesSO<-CellTypesSO(var_CombNormSO,var_CellDex_snapshotDate_2021_10_19,var_MetadataTable_CombNorm)
saveRDS(var_CellTypesSO, paste0(rds_output,"/var_CellTypesSO.rds"))
