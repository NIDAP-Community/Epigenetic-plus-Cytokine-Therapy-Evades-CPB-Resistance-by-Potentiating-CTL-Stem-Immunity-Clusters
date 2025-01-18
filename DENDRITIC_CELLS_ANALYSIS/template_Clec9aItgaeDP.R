# Cell Annotation with Co-Expression [CCBR] [scRNA-seq] (4452781a-7015-4dad-923c-61e2a183855f): v136
Clec9aItgaeDP <- function(ReclusterSO_DualLabel_CD11cNcr1,SampleNames_DualLabel_CD11cNcr1) {  # produces 2 color channels and the overlay

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages(c("SCWorkflow","Seurat","tibble","dplyr","ggplot2","ggExtra","grid","gridExtra","scales"))

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Primary inputs:
    seurat_object <- ReclusterSO_DualLabel_CD11cNcr1

    #Basic Parameters:
    samples <- c("aPD1","ENT","ENT+aPD1","ENT+N803","ENT+N803+aPD1","N803","N803+aPD1","PBS")
    marker1 <- "Clec9a"
    marker_1_threshold <- 0.02
    marker_1_type <- "SCT"
    marker2 <- "Itgae"
    marker_2_threshold <- 0.02
    marker_2_type <- "SCT"
    
    #Filter Parameters:
    filter_data <- TRUE
    parameter_name <- "Clec9aItgaeDP"
    marker.1.filter.direction <- "greater than"
    marker.2.filter.direction <- "greater than"
    apply_filter_1 <- TRUE
    apply_filter_2 <- TRUE
    filter_condition <- TRUE
    
    #Visualization Parameters:
    data_reduction <- "umap"
    add_marker_thresholds <- TRUE
    point_size <- 0.5
    point_shape <- 16
    point_transparency <- 0.5

    #Advanced Parameters:
    trim_marker_1 <- FALSE
    trim_marker_2 <- FALSE
    pre_scale_trim <- 0.99    
    display_unscaled_values <- FALSE
    seurat_object_filename <- "seurat_object.rds"

    ## -------------------------------- ##
    ## Functions                        ##
    ## -------------------------------- ##
    
# auto removed:     path <- nidapGetPath(seurat_object, seurat_object_filename)
    so <- seurat_object
    
    #In case of NIDAPism:
    colnames(so@meta.data)[colnames(so@meta.data) == "orig_ident"] <- "orig.ident"

    so.result <- dualLabeling(object = so,
                           samples = samples,
                           marker.1 = marker1,
                           marker.2 = marker2,
                           marker.1.type = marker_1_type,
                           marker.2.type = marker_2_type,
                           data.reduction = data_reduction,
                           point.size = point_size,
                           point.shape = point_shape,
                           point.transparency = point_transparency,
                           add.marker.thresholds = add_marker_thresholds,
                           marker.1.threshold = marker_1_threshold,
                           marker.2.threshold = marker_2_threshold,
                           filter.data = filter_data,
                           marker.1.filter.direction = marker.1.filter.direction,
                           marker.2.filter.direction = marker.2.filter.direction,
                           apply.filter.1 = apply_filter_1,
                           apply.filter.2 = apply_filter_2,
                           filter.condition = filter_condition,
                           parameter.name = parameter_name,
                           trim.marker.1 = trim_marker_1,
                           trim.marker.2 = trim_marker_2,
                           pre.scale.trim = pre_scale_trim,
                           display.unscaled.values = display_unscaled_values) 

    #First plot showing tSNE or UMAP 
if (data_reduction=='tsne'|data_reduction=='umap') {      
      
      grid.draw(so.result$plot)

    } else if (data_reduction=='both'){
      
       grid.draw(so.result$plot_umap)
       grid.newpage()
       grid.draw(so.result$plot_tsne)
    }

    #First plot showing Density Heatmap
    grid.newpage()
    grid.draw(so.result$plot_densityHM)

    #Second plot showing numbers of cells annotated for filtering
    g <- so.result$plot_table
    g$width <- g$width * 2
    g$height <- g$height * 2
    grid.newpage()
    grid.draw(g)

# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(so.result$object)

    return(so.result$object)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

dualLabeling <- function (object, 
                          samples, 
                          marker.1, 
                          marker.2, 
                          marker.1.type = "SCT", 
                          marker.2.type = "SCT", 
                          data.reduction = "both", 
                          point.size = 0.5, 
                          point.shape = 16, 
                          point.transparency = 0.5, 
                          add.marker.thresholds = TRUE, 
                          marker.1.threshold = 0.5, 
                          marker.2.threshold = 0.5, 
                          filter.data = TRUE, 
                          marker.1.filter.direction = "greater than", 
                          marker.2.filter.direction = "greater than", 
                          apply.filter.1 = TRUE, 
                          apply.filter.2 = TRUE, 
                          filter.condition = TRUE, 
                          parameter.name = "My_CoExp", 
                          trim.marker.1 = FALSE, 
                          trim.marker.2 = FALSE, 
                          pre.scale.trim = 0.99, 
                          display.unscaled.values = FALSE) 
{
    
  #### Error Messages ####
  
  #Errors for genes not available in dataset/slot
    if (!(marker.1 %in% rownames(object))) {
        stop(sprintf("%s is not found in dataset", marker.1))
    }
    if (!(marker.2 %in% rownames(object))) {
        stop(sprintf("%s is not found in dataset", marker.2))
    }
    if (!(marker.1.type %in% names(object@assays))) {
        stop(sprintf("%s slot is not found in dataset", marker.1.type))
    }
    if (!(marker.2.type %in% names(object@assays))) {
        stop(sprintf("%s slot is not found in dataset", marker.2.type))
    }
    if (data.reduction=='both') {
      if(sum(c('tsne','umap')%in%names(object@reductions))<2){
        rdctns=names(object@reductions)[names(object@reductions)%in%c('umap','tsne')]
      stop(sprintf("Object does not contain both umap and tsne reductions.
       Change Data Reduction parameter to %s",
                   paste(rdctns,collapse=' or ')
                   )
           )
      }
    }else{ 
      if (data.reduction%in%names(object@reductions)==F){
        stop(sprintf("Object does not contain %s reduction. \n    ",
                     data.reduction),
             sprintf("Change Data Reduction parameter to %s",
                c('tsne','umap')[c('tsne','umap')%in%names(object@reductions)])
        )
      }
    }
  #### Functions ####
  
  #Function for drawing overlay images for umap/tsne:
    .ggOverlay <- function(so.sub, df, marker.1, marker.2,reduction) {
        df <- df %>% arrange(mark1.scale)
        
        xmin <- min(df$dr1) - 0.1 * min(df$dr1)
        xmax <- max(df$dr1) + 0.1 * min(df$dr1)
        
        #ggplot for umap/tsne (marker 1)
        gg.z1 <- ggplot(df, aes(dr1, dr2)) + 
          geom_point(
            color = rgb(
            red = df$mark1.scale, 
            green = 0, 
            blue = 0
            ), 
            shape = point.shape, 
            size = point.size, 
            alpha = point.transparency
          ) + 
          theme_classic() + 
          xlab(paste0(reduction, "-1")) + 
          ylab(paste0(reduction, "-2")) + 
          ggtitle(marker.1) + 
          coord_fixed()
        
        df <- df %>% arrange(mark2.scale)
        
        #ggplot for umap/tsne (marker 2)
        gg.z2 <- ggplot(df, aes(dr1, dr2)) + 
          geom_point(
            color = rgb(
              red = 0, 
              green = df$mark2.scale, 
              blue = 0
              ), 
            shape = point.shape, 
            size = point.size, 
            alpha = point.transparency
          ) + 
          theme_classic() + 
          xlab(paste0(reduction, "-1")) + 
          ylab(paste0(reduction, "-2")) + 
          ggtitle(marker.2) + 
            coord_fixed()
        
        df <- df %>% 
          mutate(avg = mark2.scale + mark1.scale) %>% 
            arrange(avg)
        
        #ggplot for umap/tsne (marker 1 & marker 2)
        gg <- ggplot(df, aes(dr1, dr2)) + 
          geom_point(
            color = rgb(
              red = df$mark1.scale, 
              green = df$mark2.scale, 
              blue = 0
            ), 
            shape = point.shape, 
            size = point.size, 
            alpha = point.transparency
          ) + 
            theme_classic() + 
          xlab(paste0(reduction, "-1")) + 
          ylab(paste0(reduction, "-2")) + 
          ggtitle("Combined") + 
          coord_fixed()
        
        return(list(gg.z1, gg.z2, gg))
    }

    #Function for plotting expression data in xy overlay format:
    .ggOverlay2 <- function(so.sub, df, marker.1, marker.2) {
        df <- df %>% arrange(mark1.scale)
        # Create unscaled axis labels
        
        if (display.unscaled.values == TRUE) {
            label1.min <- 
              paste("unscaled min:", round(min(mark1),digits = 2))
            label1.max <- paste("unscaled max:", round(max(mark1), digits = 2))
            label1 <- 
              paste(as.character(marker.1), label1.min, label1.max, sep = "\n")
            
            label2.min <- 
              paste("unscaled min:", round(min(mark2), digits = 2))
            label2.max <- 
              paste("unscaled max:", round(max(mark2), digits = 2))
            label2 <- 
              paste(as.character(marker.2), label2.min, label2.max, sep = "\n")
        } else {
          label1 <- as.character(marker.1)
          label2 <- as.character(marker.2)
        }
        
        #ggplot for scatter plot (marker 1)
        gg.z1 <- ggplot(df, aes(mark1.scale, mark2.scale)) + 
            geom_point(
              color = rgb(
                red = df$mark1.scale, 
                green = 0, 
                blue = 0
              ), 
              shape = 20, 
              size = point.size
            ) + 
            theme_classic() + 
            xlab(label1) + 
            ylab(label2) + 
            coord_fixed()
        
        df <- df %>% arrange(mark2.scale)
        
        #ggplot for scatter plot (marker 2)
        gg.z2 <- ggplot(df, aes(mark1.scale, mark2.scale)) + 
            geom_point(
              color = rgb(
                red = 0, 
                green = df$mark2.scale, 
                blue = 0
              ), 
              shape = 20, 
              size = point.size
            ) + 
            theme_classic() + 
            xlab(label1) + 
            ylab(label2) + 
            coord_fixed()
        
        df <- df %>% 
          mutate(avg = mark2.scale + mark1.scale) %>% 
          arrange(avg)
        
        #ggplot for scatter plot (marker 1 & marker.2)
        gg <- ggplot(df, aes(mark1.scale, mark2.scale)) + 
          geom_point(
            color = rgb(
              red = df$mark1.scale, 
              green = df$mark2.scale, 
              blue = 0
            ), 
            shape = 20, 
            size = point.size
          ) + 
          theme_classic() + 
          xlab(label1) + 
          ylab(label2) + 
          coord_fixed()
        
        if (add.marker.thresholds == TRUE) {
            gg.z1 <- gg.z1 + 
              geom_vline(xintercept = t1, linetype = "dashed") + 
              geom_hline(yintercept = t2, linetype = "dashed")
            gg.z2 <- gg.z2 + 
              geom_vline(xintercept = t1, linetype = "dashed") + 
              geom_hline(yintercept = t2, linetype = "dashed")
            gg <- gg + 
              geom_vline(xintercept = t1, linetype = "dashed") + 
              geom_hline(yintercept = t2, linetype = "dashed")
        }
        
        return(list(gg.z1, gg.z2, gg))
    }
    
    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##
    
    # Load and subset data using sample names
    
    if ("active.ident" %in% slotNames(object)) {
        sample.name <- as.factor(object@meta.data$orig.ident)
        names(sample.name) <- names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        so.sub <- subset(object, ident = samples)
    } else {
        sample.name <- as.factor(object@meta.data$orig.ident)
        names(sample.name) <- names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        so.sub <- subset(object, ident = samples)
    }
    t1 <- marker.1.threshold
    t2 <- marker.2.threshold
    
    
    #Select marker 1 values and scale:
    # use data slot for Spatial assay
    if(marker.1.type == "Spatial"){
        mark1 <- so.sub@assays[[marker.1.type]]@data[marker.1,]
    } else {
      mark1 <- so.sub@assays[[marker.1.type]]@scale.data[marker.1,]
    }
    if (trim.marker.1 == TRUE) {
        q1 <- quantile(mark1, pre.scale.trim)
        q0 <- quantile(mark1, 1 - pre.scale.trim)
        mark1[mark1 < q0] <- q0
        mark1[mark1 > q1] <- q1
    }
    mark1.scale <- rescale(mark1, to = c(0, 1))
    
    #Select marker 2 values and scale:
    if(marker.2.type == "Spatial"){
        mark2 <- so.sub@assays[[marker.2.type]]@data[marker.2, ]
    } else {
    mark2 <- so.sub@assays[[marker.2.type]]@scale.data[marker.2, ]
    }

    if (trim.marker.2 == TRUE) {
        q1 <- quantile(mark2, pre.scale.trim)
        q0 <- quantile(mark2, 1 - pre.scale.trim)
        mark2[mark2 < q0] <- q0
        mark2[mark2 > q1] <- q1
    }
    mark2.scale <- rescale(mark2, to = c(0, 1))
    
   
    
    if (data.reduction=='tsne'|data.reduction=='umap') {
      
      #Draw Plots:
      df <- data.frame(
        cbind(
          dr1 = so.sub@reductions[[data.reduction]]@cell.embeddings[,1], 
          dr2 = so.sub@reductions[[data.reduction]]@cell.embeddings[,2], 
          mark1.scale, 
          mark2.scale
        )
      )
      
      gg.list <- .ggOverlay(so.sub, df, marker.1, marker.2,data.reduction)
      gg.list2 <- .ggOverlay2(so.sub, df, marker.1, marker.2)
      
      grob <- 
        arrangeGrob(gg.list[[1]], 
                    gg.list[[2]], 
                    gg.list[[3]],
                    gg.list2[[1]], 
                    gg.list2[[2]], 
                    gg.list2[[3]], 
                    ncol = 3)
      
      
    } else if (data.reduction=='both'){
      
      #Draw Plots:
      df.u <- data.frame(
        cbind(
          dr1 = so.sub@reductions[['umap']]@cell.embeddings[,1], 
          dr2 = so.sub@reductions[['umap']]@cell.embeddings[,2], 
          mark1.scale, 
          mark2.scale
        )
      )
      gg.list.u <- .ggOverlay(so.sub, df.u, marker.1, marker.2,'umap')
      
      
      #Draw Plots:
      df.t <- data.frame(
        cbind(
          dr1 = so.sub@reductions[['tsne']]@cell.embeddings[,1], 
          dr2 = so.sub@reductions[['tsne']]@cell.embeddings[,2], 
          mark1.scale, 
          mark2.scale
        )
      )
      gg.list.t <- .ggOverlay(so.sub, df.t, marker.1, marker.2,'tsne')
      
      
      ## df.u being used is arbitrary. only difference between df.u and df.t is 
      ## reduction columns which are not being used by ggOverlay2
      df=df.u[,colnames(df.u)%in%c('dr1','dr2')==F]
      
      gg.list2 <- .ggOverlay2(so.sub, df, marker.1, marker.2)
     
      
      grob.u <- 
        arrangeGrob(gg.list.u[[1]], 
                    gg.list.u[[2]], 
                    gg.list.u[[3]],
                    gg.list2[[1]], 
                    gg.list2[[2]], 
                    gg.list2[[3]], 
                    ncol = 3)
      grob.t <- 
        arrangeGrob(gg.list.t[[1]], 
                    gg.list.t[[2]], 
                    gg.list.t[[3]],
                    gg.list2[[1]], 
                    gg.list2[[2]], 
                    gg.list2[[3]], 
                    ncol = 3)
      
    } else {
      stop("Incorrect selection for data.reduction: use either umap,tsne or both")
    }
    
    x = df$mark1.scale
    y = df$mark2.scale

        df_heatmap <- data.frame(
          x = x, 
          y = y, 
          d <- densCols(
            x, 
            y, 
            nbin = 1000, 
            bandwidth = 1, 
            colramp <- 
              colorRampPalette(rev(rainbow(10, end = 4/6)))
          )
        )
        
  p <- ggplot(df_heatmap) + 
    geom_point(aes(x, y, col = d), size = 1) + 
    scale_color_identity() + xlab(marker.1) + ylab(marker.2) + 
    theme_bw() + 
    geom_vline(xintercept = t1, linetype = "dashed") + 
    geom_hline(yintercept = t2, linetype = "dashed") 
            
  p2 <- ggMarginal(p, df_heatmap, x = marker.1, y = marker.2, type = "density")

    grobHM <- 
      arrangeGrob(p2,ncol=1,nrow=1)
  
    
    #Applying Filters to Data using Thresholds:
    if (filter.data == TRUE && 
        (apply.filter.1 == TRUE | apply.filter.2 == TRUE)
        ){
        df <- df %>% mutate(sample = so.sub@meta.data$orig.ident) %>% 
            mutate(cellbarcode = rownames(so.sub@meta.data))
        
        if (marker.1.filter.direction == "greater than") {
            ind1 <- df$mark1.scale > t1
        } else {
            ind1 <- df$mark1.scale < t1
        }
        
        cat("\n")
        print("Marker 1 filter:")
        print(sum(ind1))
        
        if (marker.2.filter.direction == "greater than") {
            ind2 <- df$mark2.scale > t2
        } else {
            ind2 <- df$mark2.scale < t2
        }
        
        cat("\n")
        print("Marker 2 filter:")
        print(sum(ind2))
        
        if (apply.filter.1 == TRUE) {
            if (apply.filter.2 == TRUE) {
                if (filter.condition == TRUE) {
                  df <- df[c(ind1 & ind2), ]
                } else {
                  df <- df[c(ind1 | ind2), ]
                }
            } else {
                df <- df[ind1, ]
            }
        } else {
            if (apply.filter.2) {
                df <- df[ind2, ]
            }
        }
        
        
        # Print out numbers of cells that meet threshold cutoffs for marker 1,
        # marker 2 and for either intersection or union of 2 thresholds:
        
        colnames(df)[3:4] <- c(marker.1, marker.2)
        so.sub.df <- so.sub@meta.data %>% 
          mutate(x = case_when(
            rownames(so.sub@meta.data) %in% df$cellbarcode ~ TRUE, 
            TRUE ~ FALSE
          ))
        
        colnames(so.sub.df) <- sub("x", parameter.name, 
                                   colnames(so.sub.df))
        
        data.filt <- 
          as.data.frame.matrix(table(so.sub.df[[parameter.name]], 
                                     so.sub.df$orig.ident))
        data.filt$Total <- rowSums(data.filt)
        data.filt <- data.filt %>% rownames_to_column("Passed Filter")
        
        # Add a title:
        if (filter.condition == TRUE) {
            cond = "and"
        } else {
            cond = "or"
        }
        
        if (apply.filter.1 == TRUE) {
          if (apply.filter.2 == TRUE) {
                titlename <- paste(
                  "Number of cells that pass filters:\n", 
                  marker.1, 
                  marker.1.filter.direction, 
                  marker.1.threshold, 
                  cond, 
                  marker.2, 
                  marker.2.filter.direction, 
                  marker.2.threshold
                )
            } else {
                titlename <- paste(
                  "Number of cells that pass filter:\n", 
                  marker.1, 
                  marker.1.filter.direction, 
                  marker.1.threshold
                )
            }
        } else {
            titlename <- paste(
              "Number of cells that pass filter:\n",
              marker.2, 
              marker.2.filter.direction, 
              marker.2.threshold
            )
        }
        
        title <- 
          textGrob(
            titlename,
            y = 1,
            vjust = 1, 
            gp = gpar(fontsize = 15)
          )
        grid.table <- tableGrob(data.filt, rows = NULL)
        g <- arrangeGrob(grid.table, top = title)
        g$heights[[2]] <- unit(0.5, "npc") - max(g$grobs[[1]]$heights)
        
        rownames(so.sub.df) <- rownames(so.sub@meta.data)
        so.sub@meta.data <- so.sub.df
    } else {
        g <- textGrob("No filtering thresholds applied")
    }
    
    
    if (data.reduction=='tsne'|data.reduction=='umap') {
      
      result.list <- list(object = so.sub, 
                          plot = grob,
                          plot_densityHM = grobHM,
                          plot_table = g)
      
      
    } else if (data.reduction=='both'){
      
      result.list <- list(object = so.sub, 
                          plot_tsne = grob.t,
                          plot_umap = grob.u,
                          plot_densityHM = grobHM,
                          plot_table = g)
    }

    
    return(result.list)
}

######### Node Execution Steps ##########
print("template_Clec9aItgaeDP.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ReclusterSO_DualLabel_CD11cNcr1
var_ReclusterSO_DualLabel_CD11cNcr1<-readRDS(paste0(rds_output,"/var_ReclusterSO_DualLabel_CD11cNcr1.rds"))

if (!('Seurat' %in% class(var_ReclusterSO_DualLabel_CD11cNcr1))) { if (!(class(var_ReclusterSO_DualLabel_CD11cNcr1) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ReclusterSO_DualLabel_CD11cNcr1 <- as.data.frame(var_ReclusterSO_DualLabel_CD11cNcr1)}}
#############################


# Processing input variable: var_SampleNames_DualLabel_CD11cNcr1
var_SampleNames_DualLabel_CD11cNcr1<-readRDS(paste0(rds_output,"/var_SampleNames_DualLabel_CD11cNcr1.rds"))

if (!('Seurat' %in% class(var_SampleNames_DualLabel_CD11cNcr1))) { if (!(class(var_SampleNames_DualLabel_CD11cNcr1) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_SampleNames_DualLabel_CD11cNcr1 <- as.data.frame(var_SampleNames_DualLabel_CD11cNcr1)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_Clec9aItgaeDP<-Clec9aItgaeDP(var_ReclusterSO_DualLabel_CD11cNcr1,var_SampleNames_DualLabel_CD11cNcr1)
saveRDS(var_Clec9aItgaeDP, paste0(rds_output,"/var_Clec9aItgaeDP.rds"))
