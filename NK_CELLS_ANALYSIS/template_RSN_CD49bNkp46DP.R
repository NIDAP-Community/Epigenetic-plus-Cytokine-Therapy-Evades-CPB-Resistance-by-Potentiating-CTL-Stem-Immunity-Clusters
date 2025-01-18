# Get Sample Names [CCBR] [scRNA-seq] (51ea15ce-2be0-415f-902b-3c86175eb6cd): v57
RSN_CD49bNkp46DP <- function(RSO_CD49bNkp46DP) {
    
    ## This function extracts Sample Names from a Seurat Object. 

## --------- ##
## Libraries ##
## --------- ##

library(Seurat)
library(tidyverse)
library(ggplot2)

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

Seurat_Object <- RSO_CD49bNkp46DP

## -------------------------------- ##
## Generating PostIt label: start   ##
## -------------------------------- ##
# text params
    text = "Seurat\nSample Names"
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
    png(filename="RSN_CD49bNkp46DP.png", width=image_width, height=image_height,      units="px", pointsize=4, bg="white", res=image_resolution, type="cairo") 
    
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

cat("1. Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed: fs <- Seurat_Object$fileSystem()
# auto removed: path <- fs$get_path("seurat_object.rds", 'r')
SO <- Seurat_Object
print(SO)

# extract metadata table and sample names
cat("\n2. Extracting Sample Names from metadata table: ")
if ("meta.data" %in% slotNames(SO)) {
  if ("orig.ident" %in% colnames(SO@meta.data)) {
    dim.df <-
      as.data.frame(t(as.matrix(table(
        SO@meta.data$orig.ident
      ))))
  } else {
    dim.df <-
      as.data.frame(t(as.matrix(table(
        SO@meta.data$orig_ident
      ))))
  }
  colnames(dim.df) <-
    lapply(colnames(dim.df), function(x)
      gsub(".h5", "", x))
} else {
  if ("orig.ident" %in% colnames(SO$RNA@meta.data)) {
    dim.df <-
      as.data.frame(t(as.matrix(table(
        SO$RNA@meta.data$orig.ident
      ))))
  } else {
    dim.df <-
      as.data.frame(t(as.matrix(table(
        SO$RNA@meta.data$orig_ident
      ))))
  }
  dim.df <-
    as.data.frame(t(as.matrix(table(
      SO$RNA@meta.data$orig.ident
    ))))
  colnames(dim.df) <-
    lapply(colnames(dim.df), function(x)
      gsub(".h5", "", x))
}

# output data.frame
cat(sprintf("%g sample names returned", ncol(dim.df)))
return(dim.df)    

}

######### Node Execution Steps ##########
print("template_RSN_CD49bNkp46DP.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_RSO_CD49bNkp46DP
var_RSO_CD49bNkp46DP<-readRDS(paste0(rds_output,"/var_RSO_CD49bNkp46DP.rds"))

if (!('Seurat' %in% class(var_RSO_CD49bNkp46DP))) { if (!(class(var_RSO_CD49bNkp46DP) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RSO_CD49bNkp46DP <- as.data.frame(var_RSO_CD49bNkp46DP)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_RSN_CD49bNkp46DP<-RSN_CD49bNkp46DP(var_RSO_CD49bNkp46DP)
saveRDS(var_RSN_CD49bNkp46DP, paste0(rds_output,"/var_RSN_CD49bNkp46DP.rds"))
