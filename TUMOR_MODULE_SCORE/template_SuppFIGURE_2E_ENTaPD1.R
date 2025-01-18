# Color by Metadata [CCBR] [scRNA-seq] (0fac593c-01e9-4c02-96d0-cc4876c4dfef): v163
SuppFIGURE_2E_ENTaPD1 <- function(ModScoreSO,SampleNames_CellTypes,MetadataTable_ModScore) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    library(grid)
    library(gridExtra)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    # Input Parameters:
    Seurat_Object <- ModScoreSO
    Sample_Table <- SampleNames_CellTypes
    Metadata_Table <- MetadataTable_ModScore
    
    # Samples and Metadata Parameters
    Samples_to_Include <- 'c("ENT+aPD1")'
    Metadata_to_Plot <- 'c("MS_Celltype")'

    # Visualization Parameters:
    Number_of_Columns_for_Final_Image <- 0
    Show_Labels <- FALSE
    Dot_Size <- 1
   
    # Legend Parameters
    Legend_Text_Size <- 1
    Legend_Position <- "right"

    # Summarization Parameters
    Columns_to_Summarize <- c()
    Summarization_Cut_Off <- 5

    # Advanced Parameters:
    Save_the_Entire_Dataset <- FALSE
    Use_CITE_seq <- FALSE

   

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
# auto removed:     fs <- Seurat_Object$fileSystem()
# auto removed:     path <- fs$get_path("seurat_object.rds", 'r')
    SO <- Seurat_Object
    print(SO)

All_Reduction_Types <- c("umap", "tsne", "pca")
for(Reduction_Type in All_Reduction_Types){

results <- plotMetadata(object = SO,
  samples.to.include = Samples_to_Include,
  metadata.to.plot = Metadata_to_Plot,
  columns.to.summarize = Columns_to_Summarize,
  summarization.cut.off = Summarization_Cut_Off,
  reduction.type = Reduction_Type,
  use.cite.seq = Use_CITE_seq,
  show.labels = Show_Labels,
  legend.text.size = Legend_Text_Size,
  legend.position = Legend_Position,
  dot.size = Dot_Size
  ) 

## Print Graphic output

matched_commas <- gregexpr(",", Metadata_to_Plot, fixed = TRUE)
n_commas <- length(matched_commas[[1]])
    if (Number_of_Columns_for_Final_Image == 0) {
        n = ceiling((n_commas+1)^0.5)
    } else {
        n = Number_of_Columns_for_Final_Image
    }
 do.call("grid.arrange", c(results$plot, ncol=n))

}

## Save dataset if requested

if (Save_the_Entire_Dataset){
# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(results$object)
# auto removed:     return(NULL)
    } else {
    return(results$object@meta.data)
    }
}

######### Node Execution Steps ##########
print("template_SuppFIGURE_2E_ENTaPD1.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ModScoreSO
var_ModScoreSO<-readRDS(paste0(rds_output,"/var_ModScoreSO.rds"))

if (!('Seurat' %in% class(var_ModScoreSO))) { if (!(class(var_ModScoreSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ModScoreSO <- as.data.frame(var_ModScoreSO)}}
#############################


# Processing input variable: var_SampleNames_CellTypes
var_SampleNames_CellTypes<-readRDS(paste0(rds_output,"/var_SampleNames_CellTypes.rds"))

if (!('Seurat' %in% class(var_SampleNames_CellTypes))) { if (!(class(var_SampleNames_CellTypes) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_SampleNames_CellTypes <- as.data.frame(var_SampleNames_CellTypes)}}
#############################


# Processing input variable: var_MetadataTable_ModScore
var_MetadataTable_ModScore<-readRDS(paste0(rds_output,"/var_MetadataTable_ModScore.rds"))

if (!('Seurat' %in% class(var_MetadataTable_ModScore))) { if (!(class(var_MetadataTable_ModScore) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_ModScore <- as.data.frame(var_MetadataTable_ModScore)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_SuppFIGURE_2E_ENTaPD1<-SuppFIGURE_2E_ENTaPD1(var_ModScoreSO,var_SampleNames_CellTypes,var_MetadataTable_ModScore)
saveRDS(var_SuppFIGURE_2E_ENTaPD1, paste0(rds_output,"/var_SuppFIGURE_2E_ENTaPD1.rds"))
