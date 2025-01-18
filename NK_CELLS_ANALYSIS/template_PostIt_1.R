# Post-It [CCBR] (e5b1eab9-72a2-4de9-b8eb-943618524152): v21
PostIt_1 <- function() {

library(ggplot2)

# text params

    text = "Isolate NK cells\n(Cd49b and Nkp46 dp)"
    text_color = 'Black'
    text_size = 40
    text_fontface = "plain"
    text_fontfamily = 'Arial'
    text_angle = 0
    custom_text_color = c()

 # background params
    
    fill_color = 'Purple1'
    fill_alpha = 0.3
    border_alpha = 1
    border_width = 0
    custom_fill = c()

# image params
    image_width = 8500
    image_height = 3000
    image_resolution = 300  

# output value
    tag = "CCBR"

    # RUN ====   

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
    png(filename="PostIt_1.png", width=image_width, height=image_height,      units="px", pointsize=4, bg="white", res=image_resolution, type="cairo") 
    
    fill_color = paints[fill_color]
    if(!is.null(custom_fill)) { fill_color = custom_fill }
    border_color = fill_color

    text_color = paints[text_color]
    if(!is.null(custom_text_color)) { text_color = custom_text_color }

    print(ggplot() + annotate("text", x = 1, y = 1, size = text_size, col=text_color, label = text, angle=text_angle, family=text_fontfamily, fontface=text_fontface) + theme_void() + theme(plot.background = element_rect(fill=alpha(fill_color, fill_alpha), colour = alpha(border_color, border_alpha), size = border_width, linetype = 1)))
    
    return(data.frame("Label"=text, "Tag" = tag, Time=sprintf("Sys.time:%s EDT",Sys.time())))
    
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

install_bioconductor_package <- function(pkg) {
}

######### Node Execution Steps ##########
print("template_PostIt_1.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################
# Saving Output Dataset
invisible(graphics.off())
var_PostIt_1<-PostIt_1()
saveRDS(var_PostIt_1, paste0(rds_output,"/var_PostIt_1.rds"))
