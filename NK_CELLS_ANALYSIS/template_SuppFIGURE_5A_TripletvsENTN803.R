# Volcano Plot - Enhanced [CCBR] [scRNA-seq] [Bulk] (0c91aa57-0f76-4513-a063-5f9263d65727): v77
SuppFIGURE_5A_TripletvsENTN803 <- function(DEGMarkers_CD49bNkp46DP) {
    # image: png

    # Changelog
    # 2022-09-14 Rearranged structure and description
    # 2020-10-29 Add support for pval == 0

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(stringr)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tidyr)

    library(EnhancedVolcano)
    # For interactive plot:
    library(plotly)
    library(grid)
    

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Basic Parameters:
    df.orig <- DEGMarkers_CD49bNkp46DP
    label.col <- "Gene"
    sig.col <- c("C_ENT+N803+aPD1_vs_ENT+N803_pval")
    pCutoff  = 0.05
    lfc.col <- c("C_ENT+N803+aPD1_vs_ENT+N803_logFC")
    FCcutoff = 0.585
   
    
    #Label Parameters
    value_to_sort_the_output_dataset <- "p-value"
    no_genes_to_label <- 30
    use_only_addition_labels <- TRUE
    additional_labels <- "Xcl1"
    labSize <- 10    

    #Title and Axis labels Parameters
    change_sig_name <- "p-value"
    change_lfc_name <- "log2FC"
    title <- "ENTN803aPD1_vs_ENTN803"
    #subtitle <- ""
    use_custom_lab <- FALSE
            
    #Plot Parameters
    ylim <- 0
    custom_xlim <- ""
    xlim_additional <- 0
    ylim_additional <- 0
    axisLabSize <- 24
    pointSize <- 4

    #Image Parameters
    imageWidth = 3000
    imageHeight = 3000
    dpi = 300

  
    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    rank <- list()
    for(i in 1:length(lfc.col)){
        lfccol <- lfc.col[i]
        sigcol <- sig.col[i]
        columns_of_interest <- c(label.col,lfc.col[i],sig.col[i])
        df <- df.orig %>% dplyr::select(one_of(columns_of_interest)) %>% 
                        mutate(!!sym(lfccol) := replace_na(!!sym(lfccol), 0)) %>%
                        mutate(!!sym(sigcol) := replace_na(!!sym(sigcol), 1)) 
                        #mutate(.data[[lfc.col[i]]] = replace_na(.data[[lfc.col[i]]], 0)) %>%
                        #mutate(.data[[sig.col[i]]] = replace_na(.data[[sig.col[i]]], 1)) 
        if (use_custom_lab==TRUE){
          if (nchar(change_lfc_name)==0){lfc_name = lfc.col[i]}
          if (nchar(change_sig_name)==0){sig_name = sig.col[i]}
          colnames(df) <- c(label.col,change_lfc_name, sig_name)
        } else {
          lfc_name = lfc.col[i]
          sig_name = sig.col[i]
        }
    
        group <- gsub("_pval|p_val_","",sig_name)
        rank[[i]] <- -log10(df[[sig_name]]) * sign(df[[lfc_name]]) 
        names(rank)[i] <- paste0("C_",group,"_rank")
        
        cat(paste0("Genes in initial dataset: ", nrow(df),"\n"))

        #Select top genes by logFC or Significance
            
        if (value_to_sort_the_output_dataset=="fold-change") {
            df <- df %>% dplyr::arrange(desc(.data[[lfc_name]]))
        } else if (value_to_sort_the_output_dataset=="p-value") {
            df <- df %>% dplyr::arrange(.data[[sig_name]]) 
        }
      
        genes_to_label <- as.character(df[1:no_genes_to_label,label.col])
#        additional_labels <- unlist(str_split(additional_labels,","))
 ## Modifying Additional Labels List:
 ## Replace commas with spaces and split the string
split_values <- unlist(strsplit(gsub(",", " ", additional_labels), " "))
additional_labels <- split_values[split_values != ""]

        filter <- additional_labels %in% df[,label.col]
        missing_labels <- additional_labels[!filter]
        additional_labels <- additional_labels[filter]

        if(length(missing_labels) > 0){
            cat("Could not find:\n")
            print(missing_labels)
        }

        if(use_only_addition_labels){
            genes_to_label <- additional_labels
        }else{
            genes_to_label <- unique(append(genes_to_label, additional_labels))
        }

        significant = vector(length = nrow(df))
        significant[] = "Not significant"
        significant[which(abs(df[,2]) > FCcutoff)] = "Fold change only"
        significant[which(df[,3] < pCutoff)] = "Significant only"
        significant[which(abs(df[,2]) > FCcutoff & df[,3] < pCutoff)] = "Significant and fold change"
        print(table(significant))
        
        # fix pvalue == 0
        shapeCustom <- rep(19,nrow(df))
        maxy <-  max(-log10(df[[sig_name]]), na.rm=TRUE)
        if(ylim > 0){
            maxy <- ylim
        }
      
        cat(paste0("Maxy: ",maxy,"\n"))
        if(maxy == Inf){
            # Sometimes, pvalues == 0
            keep <- df[[sig_name]] > 0
            df[[sig_name]][!keep] <- min(df[[sig_name]][keep])
            shapeCustom[!keep] <- 17

            maxy <- -log10(min(df[[sig_name]][keep]))
            cat("Some p-values equal zero. Adjusting y-limits.\n")
            cat(paste0("Maxy adjusted: ",maxy,"\n"))

        }

        # By default, nothing will be greater than maxy. User can set this value lower
        keep <- -log10(df[[sig_name]]) <= maxy
        df[[sig_name]][!keep] <- maxy
        shapeCustom[!keep] <- 17

        names(shapeCustom)<- rep("Exact",length(shapeCustom))
        names(shapeCustom)[shapeCustom == 17] <- "Adjusted"
        
        #Remove if nothin' doin'
        if(all(shapeCustom == 19)){
            shapeCustom <- NULL
        }
        
        maxy <- ceiling(maxy)

        if (grepl("log",lfc.col[i]) ){
                xlab <- bquote(~Log[2]~ "fold change")
        } else {
            xlab <- "Fold change"
        }
        if (grepl("adj",sig.col[i])){
            ylab <- bquote(~-Log[10]~ "FDR")
        } else {
            ylab <- bquote (~-Log[10]~ "p-value")
        }
        if(use_custom_lab){
            if(lfc_name != lfc.col[i]){
                xlab <- gsub("_"," ",lfc_name)
            }
            if (sig_name != sig.col[i]){ 
                ylab <- gsub("_"," ",sig_name)
            }
        }
  
## X-axis custom range change:
if (custom_xlim == "") {
        xlim=c(floor(min(df[,lfc_name])) - xlim_additional,ceiling(max(df[,lfc_name]))+ xlim_additional)
} else if (grepl(",", custom_xlim) == FALSE) {
    xlim=c(-1*as.numeric(trimws(custom_xlim)), as.numeric(trimws(custom_xlim)))
} else {
    split_values <- strsplit(custom_xlim, ",")[[1]]

    # Trim whitespace and convert to numeric values
    x_min <- as.numeric(trimws(split_values[1]))
    x_max <- as.numeric(trimws(split_values[2]))

    xlim <- c(x_min, x_max)
}

        p <- EnhancedVolcano( df,x=lfc_name,y=sig_name,
                              lab=df[,label.col],
                              selectLab = genes_to_label,
                              title=title, #CHANGE NW: See line 78
                              subtitle <- group,
                              xlab=xlab,
                              ylab=ylab,
                              xlim=xlim,
                              ylim=c(0, maxy + ylim_additional),
                              pCutoff=pCutoff,
                              FCcutoff=FCcutoff,
                              axisLabSize=axisLabSize,
                              labSize=labSize,
                              pointSize=pointSize,
                              shapeCustom=shapeCustom
                              )
        print(p)   
 ## Adding interactive plot
    px <- p +
        xlab("Fold Change") +  # Simplify x-axis label
        ylab("Significance") +    # Simplify y-axis label
        theme_minimal()    

    interactive_plot <-
      ggplotly(px) %>%
      layout(hovermode = 'closest') %>%
      style(hoverinfo = "text", text = df[[label.col]])
    grid.newpage()
    print(interactive_plot)
    grid.newpage()
    ## The end of addition
   
    }
   
    df.final <- cbind(df.orig, do.call(cbind, rank))
    return(df.final)
}

#################################################
## Global imports and functions included below ##
#################################################

######### Node Execution Steps ##########
print("template_SuppFIGURE_5A_TripletvsENTN803.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_DEGMarkers_CD49bNkp46DP
var_DEGMarkers_CD49bNkp46DP<-readRDS(paste0(rds_output,"/var_DEGMarkers_CD49bNkp46DP.rds"))

if (!('Seurat' %in% class(var_DEGMarkers_CD49bNkp46DP))) { if (!(class(var_DEGMarkers_CD49bNkp46DP) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_DEGMarkers_CD49bNkp46DP <- as.data.frame(var_DEGMarkers_CD49bNkp46DP)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_SuppFIGURE_5A_TripletvsENTN803<-SuppFIGURE_5A_TripletvsENTN803(var_DEGMarkers_CD49bNkp46DP)
saveRDS(var_SuppFIGURE_5A_TripletvsENTN803, paste0(rds_output,"/var_SuppFIGURE_5A_TripletvsENTN803.rds"))
