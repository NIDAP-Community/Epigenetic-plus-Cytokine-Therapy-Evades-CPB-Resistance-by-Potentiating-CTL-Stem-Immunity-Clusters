# Get Metadata Table [CCBR] [scRNA-seq] (6a8139b7-45b4-4c6a-8648-c8ec34e6fc60): v71
RMDT_ExcludeNK_CD8Tcells <- function(RSO_ExcludeNK_CD8Tcells) {
    
## This function extracts Metadata Table from a Seurat Object. 

## --------- ##
## Libraries ##
## --------- ##
    
library(nidapFunctions)
nidapLoadPackages(c("SCWorkflow", "magrittr", "tibble"))
library(ggplot2)

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

seurat_object <- RSO_ExcludeNK_CD8Tcells
return_cell_embeddings <- FALSE
seurat_object_filename <- "seurat_object.rds"

## -------------------------------- ##
## Generating PostIt label: start   ##
## -------------------------------- ##
# text params
    text = "Seurat\nMetadata"
    text_color = 'Black'
    text_size = 5
    text_fontface = "italic"
    text_fontfamily = 'Arial'
    text_angle = 0
    custom_text_color = c()

# background params
    fill_color = 'Yellow1'
    fill_alpha = 0.3
    border_alpha = 1
    border_width = 0
    custom_fill = c()

# image params
    image_width = 750
    image_height = 250
    image_resolution = 300  

# output value
    tag = "CCBR"

# paints    
    paints = c(        
        White="#FFFFFF", 
        Gray1 = "#EDEFF2", Gray2="#D3D3D3", Gray3="#999999", Black="#000000",
        Red1="#F44E3B", Red2="#D33115",Red3="#9F0500",
        Orange1="#FE9200", Orange2="#E27300", Orange3="#C45100",
        Yellow1="#FCDC00", Yellow2="#FCC400", Yellow3="#FB9E00",
        YellowGreen1="#DBDF00", YellowGreen2="#B0BC00", Yellowgreen3="#808900",
        Green1="#A4DD00", Green2="#68BC00", Green3="#194D33",
        Teal1="#68CCCA",Teal2="#16A5A5", Teal3="#0C797D",
        Blue1="#73D8FF",Blue2="#009CE0",Blue3="#0062B1",
        Purple1="#AEA1FF", Purple2="#7B64FF", Purple3="#653294",
        Magenta1="#FDA1FF", Magenta2="#FA28FF", Magenta3="#AB149E")    

#image: png
    png(filename="RMDT_ExcludeNK_CD8Tcells.png", width=image_width, height=image_height,      units="px", pointsize=4, bg="white", res=image_resolution, type="cairo") 
    
    fill_color = paints[fill_color]
    if(!is.null(custom_fill)) { fill_color = custom_fill }
    border_color = fill_color

    text_color = paints[text_color]
    if(!is.null(custom_text_color)) { text_color = custom_text_color }

    print(ggplot() + annotate("text", x = 1, y = 1, size = text_size, col=text_color, label = text, angle=text_angle, family=text_fontfamily, fontface=text_fontface) + theme_void() + theme(plot.background = element_rect(fill=alpha(fill_color, fill_alpha), colour = alpha(border_color, border_alpha), size = border_width, linetype = 1)))

## -------------------------------- ##
## Finished PostIt generation       ##
## -------------------------------- ##

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
cat(sprintf("\nReading Seurat Object from dataset: %s\n\n", paste(seurat_object_filename, collapse = ", ")))
# auto removed: path <- nidapGetPath(seurat_object, seurat_object_filename)
so <- seurat_object
print(so)

# get cell embeddings
if (return_cell_embeddings) {
  reds = names(so@reductions)
  cat(sprintf("\n2. Extracting cell embeddings from Seurat Object: %s\n\n", paste(reds, collapse=", ")))
  for (i in seq_along(reds)) {
    so = Seurat::AddMetaData(so, as.data.frame(so@reductions[[i]]@cell.embeddings))
    }
} else {
    cat("\n2. Skipping extracting cell embeddings from Seurat Object\n\n")
}

# extract meta.data from Seurat Object
cat("3. Extracting Metadata Table from Seurat Object\n\n")
if ("meta.data" %in% slotNames(so)) {
  met.df <- so@meta.data
} else {
  met.df <- so$RNA@meta.data
}

# if no barcode column, get barcodes from the rownames.
if (!("Barcode" %in% colnames(met.df))) {
  met.df %>% rownames_to_column("Barcode") -> met.df
}
    
# Return metadata table.
print(summary(met.df))
return(met.df)
    
}

######### Node Execution Steps ##########
print("template_RMDT_ExcludeNK_CD8Tcells.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_RSO_ExcludeNK_CD8Tcells
var_RSO_ExcludeNK_CD8Tcells<-readRDS(paste0(rds_output,"/var_RSO_ExcludeNK_CD8Tcells.rds"))

if (!('Seurat' %in% class(var_RSO_ExcludeNK_CD8Tcells))) { if (!(class(var_RSO_ExcludeNK_CD8Tcells) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RSO_ExcludeNK_CD8Tcells <- as.data.frame(var_RSO_ExcludeNK_CD8Tcells)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_RMDT_ExcludeNK_CD8Tcells<-RMDT_ExcludeNK_CD8Tcells(var_RSO_ExcludeNK_CD8Tcells)
saveRDS(var_RMDT_ExcludeNK_CD8Tcells, paste0(rds_output,"/var_RMDT_ExcludeNK_CD8Tcells.rds"))
