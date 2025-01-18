# Color by Genes [CCBR] [scRNA-seq] (82c9ec27-5cd8-4fd2-87a9-2a0f83e2ac3a): v108
ColorByGene_1 <- function(RSO_ExcludeNK_CD8Tcells,RSN_ExcludeNK_CD8Tcells) {
    
    
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
    seurat_object <- RSO_ExcludeNK_CD8Tcells

    # Basic Parameters:
    samples_to_include <- 'c("LN_aPD1","LN_ENT","LN_ENT+aPD1","LN_ENT+N803","LN_ENT+N803+aPD1","LN_N803","LN_N803+aPD1","LN_PBS","T_aPD1","T_ENT","T_ENT+aPD1","T_ENT+N803","T_ENT+N803+aPD1","T_N803","T_N803+aPD1","T_PBS")'
    gene <- c("Gzmm","Tcf7","Il7r","Ly6a","Ccl5","Lef1","Sell","Cd44")

    # Visualization Parameters:
    reduction_type <- "umap"
    number_of_rows <- 2
    color <- "red"
    point_size <- 1
    point_shape <- 16
    point_transparency <- 4
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
print("template_ColorByGene_1.R #########################################################################")
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

# Saving Output Dataset
invisible(graphics.off())
var_ColorByGene_1<-ColorByGene_1(var_RSO_ExcludeNK_CD8Tcells,var_RSN_ExcludeNK_CD8Tcells)
saveRDS(var_ColorByGene_1, paste0(rds_output,"/var_ColorByGene_1.rds"))
