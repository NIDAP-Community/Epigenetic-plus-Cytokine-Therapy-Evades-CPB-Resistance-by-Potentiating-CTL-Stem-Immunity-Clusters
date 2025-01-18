# Color by Genes [CCBR] [scRNA-seq] (82c9ec27-5cd8-4fd2-87a9-2a0f83e2ac3a): v108
SuppFIGURE_3A_PBS <- function(ReclusterSO_CD8cd3,SampleNames_Reclustered_CD8CD3) {
    
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##
    
    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    library(grid)
    library(gridExtra)
    library(scales)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    # Primary Inputs:
    seurat_object <- ReclusterSO_CD8cd3

    # Basic Parameters:
    samples_to_include <- 'c("PBS")'
    gene <- c("Gzmb","Tcf7","Il7r","Cd69","Itgae","Lef1","Sell","Cd44")

    # Visualization Parameters:
    reduction_type <- "umap"
    number_of_rows <- 2
    color <- "red"
    point_size <- 2
    point_shape <- 16
    point_transparency <- 2
    image_type <- "png"

    # Advanced Parameters:
    save_the_entire_dataset <- FALSE
    use_cite_seq_data <- FALSE
    assay <- "SCT"

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##
    
    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    # Loading Seurat Object
    cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:     fs <- seurat_object$fileSystem()
# auto removed:     path <- fs$get_path("seurat_object.rds", 'r')
    SO <- seurat_object
    print(SO)

ColorByGene.result <- colorByGene(object = SO,
                                    samples.to.include = samples_to_include,
                                    gene = gene,
                                    reduction.type = reduction_type,
                                    number.of.rows = number_of_rows,
                                    return.seurat.object = save_the_entire_dataset,
                                    color = color,
                                    point.size = point_size,
                                    point.shape = point_shape,
                                    point.transparency = point_transparency,
                                    use.cite.seq.data = use_cite_seq_data,
                                    assay = assay)

# Preparing Graphic Output
    if (number_of_rows == 0) {
        n = ceiling(length(ColorByGene.result$plot)^0.5)
    } else {
        n = number_of_rows
    }

do.call("grid.arrange", c(ColorByGene.result$plot, nrow=n))

# Saving dataset (if requested)
if(save_the_entire_dataset){
    return(ColorByGene.result$object)
  }
  else{
    gene = as.data.frame(gene) 
    return(gene)
  }
}

######### Node Execution Steps ##########
print("template_SuppFIGURE_3A_PBS.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ReclusterSO_CD8cd3
var_ReclusterSO_CD8cd3<-readRDS(paste0(rds_output,"/var_ReclusterSO_CD8cd3.rds"))

if (!('Seurat' %in% class(var_ReclusterSO_CD8cd3))) { if (!(class(var_ReclusterSO_CD8cd3) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ReclusterSO_CD8cd3 <- as.data.frame(var_ReclusterSO_CD8cd3)}}
#############################


# Processing input variable: var_SampleNames_Reclustered_CD8CD3
var_SampleNames_Reclustered_CD8CD3<-readRDS(paste0(rds_output,"/var_SampleNames_Reclustered_CD8CD3.rds"))

if (!('Seurat' %in% class(var_SampleNames_Reclustered_CD8CD3))) { if (!(class(var_SampleNames_Reclustered_CD8CD3) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_SampleNames_Reclustered_CD8CD3 <- as.data.frame(var_SampleNames_Reclustered_CD8CD3)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_SuppFIGURE_3A_PBS<-SuppFIGURE_3A_PBS(var_ReclusterSO_CD8cd3,var_SampleNames_Reclustered_CD8CD3)
saveRDS(var_SuppFIGURE_3A_PBS, paste0(rds_output,"/var_SuppFIGURE_3A_PBS.rds"))
