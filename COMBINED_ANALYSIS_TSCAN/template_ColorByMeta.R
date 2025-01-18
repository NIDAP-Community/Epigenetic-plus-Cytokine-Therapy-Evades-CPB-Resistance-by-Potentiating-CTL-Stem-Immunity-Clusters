# Color by Metadata [CCBR] [scRNA-seq] (0fac593c-01e9-4c02-96d0-cc4876c4dfef): v163
ColorByMeta <- function(FSO_ExcludeNK_CD8Tcells,SampleNames_CellTypes,FMDT_ExcludeNK_CD8Tcells) {

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
    Seurat_Object <- FSO_ExcludeNK_CD8Tcells
    Sample_Table <- SampleNames_CellTypes
    Metadata_Table <- FMDT_ExcludeNK_CD8Tcells
    
    # Samples and Metadata Parameters
    Samples_to_Include <- 'c("LN_aPD1","LN_ENT","LN_ENT+aPD1","LN_ENT+N803","LN_ENT+N803+aPD1","LN_N803","LN_N803+aPD1","LN_PBS","T_aPD1","T_ENT","T_ENT+aPD1","T_ENT+N803","T_ENT+N803+aPD1","T_N803","T_N803+aPD1","T_PBS")'
    Metadata_to_Plot <- 'c("mouseRNAseq_main","immgen_main")'

    # Visualization Parameters:
    Number_of_Columns_for_Final_Image <- 0
    Show_Labels <- FALSE
    Dot_Size <- 0.01
   
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
print("template_ColorByMeta.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_FSO_ExcludeNK_CD8Tcells
var_FSO_ExcludeNK_CD8Tcells<-readRDS(paste0(rds_output,"/var_FSO_ExcludeNK_CD8Tcells.rds"))

if (!('Seurat' %in% class(var_FSO_ExcludeNK_CD8Tcells))) { if (!(class(var_FSO_ExcludeNK_CD8Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_FSO_ExcludeNK_CD8Tcells <- as.data.frame(var_FSO_ExcludeNK_CD8Tcells)}}
#############################


# Processing input variable: var_SampleNames_CellTypes
var_SampleNames_CellTypes<-readRDS(paste0(rds_output,"/var_SampleNames_CellTypes.rds"))

if (!('Seurat' %in% class(var_SampleNames_CellTypes))) { if (!(class(var_SampleNames_CellTypes) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_SampleNames_CellTypes <- as.data.frame(var_SampleNames_CellTypes)}}
#############################


# Processing input variable: var_FMDT_ExcludeNK_CD8Tcells
var_FMDT_ExcludeNK_CD8Tcells<-readRDS(paste0(rds_output,"/var_FMDT_ExcludeNK_CD8Tcells.rds"))

if (!('Seurat' %in% class(var_FMDT_ExcludeNK_CD8Tcells))) { if (!(class(var_FMDT_ExcludeNK_CD8Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_FMDT_ExcludeNK_CD8Tcells <- as.data.frame(var_FMDT_ExcludeNK_CD8Tcells)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ColorByMeta<-ColorByMeta(var_FSO_ExcludeNK_CD8Tcells,var_SampleNames_CellTypes,var_FMDT_ExcludeNK_CD8Tcells)
saveRDS(var_ColorByMeta, paste0(rds_output,"/var_ColorByMeta.rds"))
