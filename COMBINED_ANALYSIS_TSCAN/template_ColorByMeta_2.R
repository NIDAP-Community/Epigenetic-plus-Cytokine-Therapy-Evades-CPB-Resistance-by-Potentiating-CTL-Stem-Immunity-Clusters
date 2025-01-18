# Color by Metadata [CCBR] [scRNA-seq] (0fac593c-01e9-4c02-96d0-cc4876c4dfef): v163
ColorByMeta_2 <- function(RSO_ExcludeNK_CD8Tcells,RSN_ExcludeNK_CD8Tcells,RMDT_ExcludeNK_CD8Tcells) {

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
    Seurat_Object <- RSO_ExcludeNK_CD8Tcells
    Sample_Table <- RSN_ExcludeNK_CD8Tcells
    Metadata_Table <- RMDT_ExcludeNK_CD8Tcells
    
    # Samples and Metadata Parameters
    Samples_to_Include <- 'c("LN_aPD1","LN_ENT","LN_ENT+aPD1","LN_ENT+N803","LN_ENT+N803+aPD1","LN_N803","LN_N803+aPD1","LN_PBS","T_aPD1","T_ENT","T_ENT+aPD1","T_ENT+N803","T_ENT+N803+aPD1","T_N803","T_N803+aPD1","T_PBS")'
    Metadata_to_Plot <- 'c("My_Variable_1")'

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
print("template_ColorByMeta_2.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_RSO_ExcludeNK_CD8Tcells
var_RSO_ExcludeNK_CD8Tcells<-readRDS(paste0(rds_output,"/var_RSO_ExcludeNK_CD8Tcells.rds"))

if (!('Seurat' %in% class(var_RSO_ExcludeNK_CD8Tcells))) { if (!(class(var_RSO_ExcludeNK_CD8Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RSO_ExcludeNK_CD8Tcells <- as.data.frame(var_RSO_ExcludeNK_CD8Tcells)}}
#############################


# Processing input variable: var_RSN_ExcludeNK_CD8Tcells
var_RSN_ExcludeNK_CD8Tcells<-readRDS(paste0(rds_output,"/var_RSN_ExcludeNK_CD8Tcells.rds"))

if (!('Seurat' %in% class(var_RSN_ExcludeNK_CD8Tcells))) { if (!(class(var_RSN_ExcludeNK_CD8Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RSN_ExcludeNK_CD8Tcells <- as.data.frame(var_RSN_ExcludeNK_CD8Tcells)}}
#############################


# Processing input variable: var_RMDT_ExcludeNK_CD8Tcells
var_RMDT_ExcludeNK_CD8Tcells<-readRDS(paste0(rds_output,"/var_RMDT_ExcludeNK_CD8Tcells.rds"))

if (!('Seurat' %in% class(var_RMDT_ExcludeNK_CD8Tcells))) { if (!(class(var_RMDT_ExcludeNK_CD8Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RMDT_ExcludeNK_CD8Tcells <- as.data.frame(var_RMDT_ExcludeNK_CD8Tcells)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ColorByMeta_2<-ColorByMeta_2(var_RSO_ExcludeNK_CD8Tcells,var_RSN_ExcludeNK_CD8Tcells,var_RMDT_ExcludeNK_CD8Tcells)
saveRDS(var_ColorByMeta_2, paste0(rds_output,"/var_ColorByMeta_2.rds"))
