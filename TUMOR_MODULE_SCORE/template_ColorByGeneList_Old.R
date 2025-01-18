# Color by Gene Lists [CCBR] [scRNA-seq] (d71ed4e6-a25d-4f66-a186-27c00a50a703): v125
ColorByGeneList_Old <- function(CellTypesSO,SampleNames_CellTypes, My_GeneListsForModScores) {
    
    #image: png
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Input Parameters:
    seurat_object <- CellTypesSO
    marker_list = My_GeneListsForModScores

    #Basic Parameters:
    manual_genes <- c()
    cells_of_interest <- c("T_cell","CD8_T","CD4_T","Treg","Macrophage","M1","M2","Neutrophil","Monocyte","cDC","B_cell","NK")
    protein_presence <- FALSE

    #Plot Parameters:
    assay_to_plot <- "SCT"
    slot_to_plot <- "scale.data"
    reduction_type = "umap"
    point_transparency <- 0.5
    point_shape <- 16
    cite_seq <- FALSE
    consolidated_marker_fig <- TRUE

    #Filesave Parameters:
    save_the_entire_dataset <- FALSE
    seurat_object_filename <- "seurat_object.rds"

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

# auto removed:     path <- nidapGetPath(seurat_object,seurat_object_filename)
    so <- seurat_object

    colorByMarkerTable(object = so,
                        manual.genes = manual_genes,
                        marker.table = marker_list,
                        cells.of.interest = cells_of_interest,
                        protein.presence = protein_presence,
                        assay = assay_to_plot,
                        slot = slot_to_plot,
                        reduction.type = reduction_type,
                        point.transparency = point_transparency,
                        point.shape = point_shape,
                        cite.seq = FALSE)

# auto removed: return(NULL)
}

#################################################
## Global imports and functions included below ##
#################################################

colorByMarkerTable <- function (object, samples.subset, samples.to.display, 
    manual.genes = c(), marker.table, 
    cells.of.interest, protein.presence = FALSE, assay = "SCT", slot = "scale.data",
    reduction.type = "umap", point.transparency = 0.5, point.shape = 16, 
    cite.seq = FALSE){ 

        library(ggplot2)
        library(Seurat)
        library(stringr)
        library(grid)
        library(gridExtra)
        
        .plotMarkers <- function(markers) {
            if (is.na(markers) == TRUE) {
                g <- ggplot() + theme_void()
                return(g)
            } else {
                markers.mat = GetAssayData(object, assay = assay, slot = slot)[markers, ]
                markers.quant = quantile(markers.mat[markers.mat > 
                    1], probs = c(0.1, 0.5, 0.9))
                markers.mat[markers.mat > markers.quant[3]] = markers.quant[3]
                markers.mat[markers.mat < markers.quant[1]] = 0
                if (!(cite.seq)) {
                    if (reduction.type == "tsne") {
                    p1 <- DimPlot(object, reduction = "tsne", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$tSNE_1, 
                        umap2 = p1$data$tSNE_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else if (reduction.type == "umap") {
                    p1 <- DimPlot(object, reduction = "umap", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$UMAP_1, 
                        umap2 = p1$data$UMAP_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else {
                    p1 <- DimPlot(object, reduction = "pca", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$PC_1, 
                        umap2 = p1$data$PC_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    }
                } else {
                    if (reduction.type == "tsne") {
                    p1 <- DimPlot(object, reduction = "protein_tsne", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$protein_tsne_1, 
                        umap2 = p1$data$protein_tsne_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else if (reduction.type == "umap") {
                    p1 <- DimPlot(object, reduction = "protein_umap", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$protein_umap_1, 
                        umap2 = p1$data$protein_umap_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else {
                    p1 <- DimPlot(object, reduction = "protein_pca", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$protein_pca_1, 
                        umap2 = p1$data$protein_pca_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    }
                }
                
                g <- ggplot(clusmat, aes(x = umap1, y = umap2, 
                group = ident)) + theme_bw() + theme(legend.title = element_blank()) + 
                ggtitle(markers) + geom_point(aes(color = markers, 
                shape = ident), alpha = point.transparency, 
                shape = point.shape, size = 1) + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank(), 
                legend.text = element_text(size = rel(0.5)), plot.background=element_rect(fill="white", colour="white"), rect = element_rect(fill = 'white')) +
                scale_color_gradient(limits = c(0, markers.quant[3]), 
                    low = "lightgrey", high = "red") + xlab("umap-1") + 
                ylab("umap-2")
                return(g)
                
            }
        }
        
        marker.table <- marker.table[cells.of.interest]
        present.marker.ls <- list()
        for (celltype in colnames(marker.table)) {
            print(names(marker.table[celltype]))
            present = lapply(marker.table[[celltype]], function(x) x %in% 
                rownames(GetAssayData(object, assay = assay, slot = slot)))
            absent.genes = unlist(marker.table[[celltype]])[present == 
                FALSE]
            present.genes = unlist(marker.table[[celltype]])[present == 
                TRUE]
            print(paste0("Genes not present: ", paste0(absent.genes, 
                collapse = ",")))
            print(paste0("Genes present: ", paste0(present.genes, 
                collapse = ",")))
            if (length(present.genes) == 0) {
                print(paste0(names(marker.table[celltype]), " genes were not found in object and will not be analyzed"))
            } else {
                present.marker.ls[[celltype]] <- present.genes
            }
        }
        
        padded.ls <- lapply(present.marker.ls, `length<-`, max(lengths(present.marker.ls)))
        markers.from.list <- do.call(cbind, padded.ls)
        markers.present = unlist(markers.from.list)

        if (!length(markers.present) > 0) {
            print("No markers found in dataset")
# auto removed:             return(NULL)
        }

        # Store images for making consolidated figure
        cons.gg.storage <- list()

         # Plot arranged figures with padding
        for (cell in colnames(markers.from.list)) {
            title <- cell
            markers.to.analyze <- as.character(markers.from.list[, 
                cell])
            grob <- lapply(markers.to.analyze, function(x) .plotMarkers(x))

            arranged <- gridExtra::arrangeGrob(grobs = grob, ncol = 1, newpage = F, as.table = F, top = ggpubr::text_grob(title, size = 15, face = "bold"))
            
            cons.gg.storage[[cell]] <- arranged
        }
        # Manually specify background 
        background <- rectGrob(gp = gpar(fill = "white", col = NA)) 
        
        cons.fig <- do.call(grid.arrange, c(cons.gg.storage, ncol = ncol(markers.from.list)))
        cons.fig.bkgrd <- grobTree(background, cons.fig)
        grid.draw(cons.fig.bkgrd)
        grid.newpage()

        # Plot individual arranged figures with padding
        for (cell in colnames(markers.from.list)) {
            title <- cell
            markers.to.analyze <- as.character(markers.from.list[, 
                cell])
            grob <- lapply(markers.to.analyze, function(x) .plotMarkers(x))
            arranged <- gridExtra::arrangeGrob(grobs = grob, newpage = F, as.table = F, top = ggpubr::text_grob(title, size = 15, face = "bold"))

            # Combine the background with the arranged plots
            combined <- grobTree(background, arranged)
            # Draw the combined grob with background
            grid.draw(combined)
            grid.newpage()
        }

        # Plot manual genes if not empty
        if(!is.null(manual.genes)){

            # Str-spit and use only present genes
            manual.genes.processed <- str_split(manual.genes, pattern = "[^a-zA-Z0-9]+")[[1]]
            manual.genes.processed <- manual.genes.processed[manual.genes.processed %in% rownames(GetAssayData(object, assay = assay, slot = slot))]

            tryCatch({
            manual.grob <- lapply(manual.genes.processed, function(x) .plotMarkers(x))
            manual.arranged <- gridExtra::arrangeGrob(grobs = manual.grob, newpage = F, as.table = F, top = ggpubr::text_grob("Manual Genes", size = 15, face = "bold"))
            plot(manual.arranged)
            grid.newpage()}, error = function(e) {
            cat(e$message, "\n", "Possible Reason: No manual genes were found in expression matrix")
})
        }
        
# auto removed:         return(NULL)
    }

######### Node Execution Steps ##########
print("template_ColorByGeneList_Old.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_CellTypesSO
var_CellTypesSO<-readRDS(paste0(rds_output,"/var_CellTypesSO.rds"))

if (!('Seurat' %in% class(var_CellTypesSO))) { if (!(class(var_CellTypesSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_CellTypesSO <- as.data.frame(var_CellTypesSO)}}
#############################


# Processing input variable: var_SampleNames_CellTypes
var_SampleNames_CellTypes<-readRDS(paste0(rds_output,"/var_SampleNames_CellTypes.rds"))

if (!('Seurat' %in% class(var_SampleNames_CellTypes))) { if (!(class(var_SampleNames_CellTypes) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_SampleNames_CellTypes <- as.data.frame(var_SampleNames_CellTypes)}}
#############################


# Processing input variable: var_My_GeneListsForModScores
var_My_GeneListsForModScores<-readRDS(paste0(rds_output,"/var_My_GeneListsForModScores.rds"))

if (!('Seurat' %in% class(var_My_GeneListsForModScores))) { if (!(class(var_My_GeneListsForModScores) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_My_GeneListsForModScores <- as.data.frame(var_My_GeneListsForModScores)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ColorByGeneList_Old<-ColorByGeneList_Old(var_CellTypesSO,var_SampleNames_CellTypes,var_My_GeneListsForModScores)
saveRDS(var_ColorByGeneList_Old, paste0(rds_output,"/var_ColorByGeneList_Old.rds"))
