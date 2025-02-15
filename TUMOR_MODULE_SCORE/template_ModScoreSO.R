# Cell Annotation with Module Scores [CCBR] [scRNA-seq] (10cf059e-0bd0-4a1c-9a3a-65b8688dab23): v209
ModScoreSO <- function(CellTypesSO,My_GeneListsForModScores,My_MultiLevelClassesForModScores) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Primary Inputs:
    seurat_object = CellTypesSO
    marker_table <- My_GeneListsForModScores
    levels_data_frame <- My_MultiLevelClassesForModScores

    #Basic Parameters:
    celltype_to_analyze <- c("T_cell","CD8_T","CD4_T","Treg","Macrophage","M1","M2","Neutrophil","Monocyte","cDC","NK","B_cell")
    signature_threshold <- c("T_cell 0.35","CD8_T 0.2","CD4_T 0.12","Treg 0.1","Macrophage 0.4","M1 0.25","M2 0.4","Neutrophil 0.3","Monocyte 0.4","cDC 0.15","NK 0.125","B_cell 0.05")
    general_class <- c("T_cell","Macrophage","Neutrophil","Monocyte","cDC","NK","B_cell")
    multi_level_class <- TRUE

    #Plot Parameters:
    reduction = "umap"
    nbins <- 24
    gradient_density_font_size <- 6
    violinplot_font_size <- 6
    step_size <- 0.1
    
    #Advanced Parameters:
    seurat_object_filename <- "seurat_object.rds"

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

# auto removed:     path <- nidapGetPath(seurat_object,seurat_object_filename)
    so <- seurat_object

    MS_res <- modScore(object = so,
             marker.table = marker_table,
             use_columns = celltype_to_analyze,
             ms_threshold = signature_threshold,
             general.class = general_class,
             multi.lvl = multi_level_class,
             lvl.df = levels_data_frame,
             reduction = reduction,
             nbins = nbins,
             gradient.ft.size = gradient_density_font_size,
             violin.ft.size = violinplot_font_size,
             step.size = step_size)

# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(MS_res)

# auto removed:     return(NULL)
}
#################################################
## Global imports and functions included below ##
#################################################

modScore <- function(object, marker.table, ms_threshold, use_columns,
    general.class, multi.lvl = FALSE, lvl.df, 
    reduction = "tsne", nbins = 10, gradient.ft.size = 6, 
    violin.ft.size = 6, step.size = 0.1) 
{
    library(Seurat)
    library(gridExtra)
    library(grid)
    library(dplyr)
    library(stringr)
    library(ggplot2)

    # Function for separating and calling cells by bimodal thresholds
    .modScoreCall <- function(ms.meta, numeric_threshold, reject) {
        thres.ls <- list()
        for (i in 1:ncol(ms.meta)) {
            thres.ls[[i]] <- rep(numeric_threshold[i], nrow(ms.meta))
        }
        thres.df <- data.frame(matrix(unlist(thres.ls), nrow = nrow(ms.meta)))
        thres.filter <- ms.meta > thres.df
        ms.meta.filt <- ms.meta * thres.filter
        max.col.vec <- max.col(ms.meta.filt)
        zero.filt <- as.integer(!apply(ms.meta.filt, 1, function(find_zero_rows) all(find_zero_rows == 0)))
        final.filt <- (max.col.vec * zero.filt) + 1
        append.name <- c(reject, names(ms.meta))
        dupl.data <- ms.meta
        dupl.data[, "MS_Celltype"] <- append.name[final.filt]
        return(dupl.data)
    }

    # Upstream processing
    # String split celltype_thresholds - numeric portion
    numeric_threshold <- sapply(stringr::str_split(ms_threshold, " "), function(x) as.numeric(x[2]))

    if (!"Barcode" %in% colnames(object@meta.data)) {
        object@meta.data$Barcode <- rownames(object@meta.data)
    }
    colnames(object@meta.data) <- gsub("orig_ident", "orig.ident", 
        colnames(object@meta.data))
    
    # Marker table processing
    marker.table <- marker.table[,use_columns]
    marker.tab <- unlist(marker.table)
    celltypes <- sapply(str_split(ms_threshold, " "), function(x) as.character(x[1]))

    if (any(!celltypes %in% use_columns)){
        unmatched_celltypes <- celltypes[!celltypes %in% use_columns]
        celltype_mismatch_message <- paste0("Labels from thresholds does not match columns from marker table: ",paste(unmatched_celltypes, collapse = ", "))
        stop(celltype_mismatch_message) 
    }

    marker = select(marker.table, celltypes)
    marker.list = as.list(marker)
    if (sum(unlist(marker.list) %in% rownames(object@assays$SCT@data)) == 
        0) {
        stop("No genes from list was found in data")
    }
    if (length(numeric_threshold) != length(celltypes)) {
        if (sum(numeric_threshold) == 0) {
            numeric_threshold <- rep(0, length(celltypes))
            print("Module Score threshold set to zero - outputing preliminary data")
        } else {
            stop("Threshold length does not match # celltypes to analyze")
        }
    }

    # For each celltype, print out present / nonpresent genes, calculate MS and generate plots
    names(numeric_threshold) <- celltypes
    figures <- list()
    exclude_cells <- c()
    h = 0
    j = 1
    for (h in seq_along(marker.list)) {
        print(names(marker.list[h]))
        present = lapply(marker.list[[h]], function(x) x %in% 
            rownames(object@assays$SCT@data))
        absentgenes = unlist(marker.list[[h]])[present == FALSE]
        absentgenes = absentgenes[is.na(absentgenes) == F]
        presentgenes = unlist(marker.list[[h]])[present == TRUE]
        presentgenes = presentgenes[is.na(presentgenes) == F]
        print(paste0("Genes not present: ", paste0(absentgenes, 
            collapse = ",")))
        print(paste0("Genes present: ", paste0(presentgenes, 
            collapse = ",")))
        if (length(presentgenes) == 0) {
            print(paste0(names(marker.list[h]), " genes were not found in object and will not be analyzed"))
            exclude_cells[j] <- h
            j = j + 1
        }
    }
    # End of check present / absent genes

    if (length(exclude_cells) > 0) {
        marker.list <- marker.list[-exclude_cells]
    } else {
        marker.list <- marker.list
    }

    # Calculate MS, make density plots
    for (celltype_name in names(marker.list)) {
        object = AddModuleScore(object, marker.list[celltype_name], 
            name = celltype_name, nbin = nbins, assay = "SCT")
        m = paste0(celltype_name, "1")
        object@meta.data[[m]] <- scales::rescale(object@meta.data[[m]], 
            to = c(0, 1))

        # Do plots for just general (parent) celltypes
        if (celltype_name %in% general.class){
        clusid = object@meta.data[[m]]
        d <- density(clusid)

        # Make dimension reduction tables and plots
        if (reduction == "tsne") {
            p1 <- DimPlot(object, reduction = "tsne", group.by = "ident")
        } else if (reduction == "umap") {
            p1 <- DimPlot(object, reduction = "umap", group.by = "ident")
        } else {
            p1 <- DimPlot(object, reduction = "pca", group.by = "ident")
        }

        if (reduction == "tsne") {
            clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$tSNE_1, umap2 = p1$data$tSNE_2, clusid = as.numeric(object@meta.data[[m]]))
            } else if (reduction == "umap") {
            clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$UMAP_1, umap2 = p1$data$UMAP_2, clusid = as.numeric(object@meta.data[[m]]))
            } else {
            clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$PC_1, 
                umap2 = p1$data$PC_2, clusid = as.numeric(object@meta.data[[m]]))}

        clusmat <- mutate(clusmat, sample_clusid = clusmat$clusid)
        umap.pos <- clusmat %>% group_by(clusid) %>% dplyr::summarise(umap1.mean = mean(umap1), umap2.mean = mean(umap2))
        title = as.character(m)
        clusmat <- clusmat %>% dplyr::arrange(clusid)
        clusid.df <- data.frame(id = object@meta.data$orig.ident, 
            ModuleScore = object@meta.data[[m]])

        g <- ggplot(clusmat, aes(x = umap1, y = umap2)) + theme_bw() + 
            theme(legend.title = element_blank()) + geom_point(aes(colour = sample_clusid), alpha = 0.5, shape = 20, size = 1) + scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = scales::rescale(c(0, 
            numeric_threshold[celltype_name]/2, numeric_threshold[celltype_name], (numeric_threshold[celltype_name] + 1)/2, 
            1), limits = c(0, 1))) + guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) + theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank()) + xlab("tsne-1") + ylab("tsne-2")

        g1 <- RidgePlot(object, features = m, group.by = "orig.ident") + 
            theme(legend.position = "none", title = element_blank(), 
                axis.text.x = element_text(size = gradient.ft.size)) + 
            geom_vline(xintercept = numeric_threshold[celltype_name], linetype = "dashed", 
                color = "red3") + scale_x_continuous(breaks = seq(0, 
            1, step.size))

        g2 <- ggplot(clusid.df, aes(x = id, y = ModuleScore)) + 
            geom_violin(aes(fill = id)) + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = rel(0.8)), legend.position = "top", axis.text.y = element_text(size = violin.ft.size)) + 
            guides(colour = guide_legend(override.aes = list(size = 5, 
                alpha = 1))) + geom_hline(yintercept = numeric_threshold[celltype_name], 
            linetype = "dashed", color = "red3") + scale_y_continuous(breaks = seq(0, 1, step.size))

        g3 <- ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + 
            xlab("ModuleScore") + ylab("Density") + geom_line() + 
            geom_segment(aes(xend = d$x, yend = 0, colour = x)) + 
            scale_y_log10() + scale_color_gradientn(colours = c("blue4", 
            "lightgrey", "red"), values = scales::rescale(c(0, 
            numeric_threshold[celltype_name]/2, numeric_threshold[celltype_name], (numeric_threshold[celltype_name] + 1)/2, 
            1), limits = c(0, 1))) + geom_vline(xintercept = numeric_threshold[celltype_name], 
            linetype = "dashed", color = "red3") + geom_vline(xintercept = numeric_threshold[celltype_name], linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0, 1, step.size)) + theme(legend.title = element_blank(), 
            axis.text.x = element_text(size = 6))

        figures[[celltype_name]] = arrangeGrob(g, g1, g2, g3, ncol = 2, top = textGrob(paste0(celltype_name," (General Class)"), gp = gpar(fontsize = 14, fontface = "bold")))
        }
    }

    # Rename MS columns - get rid of "1" at the end
    colnames(object@meta.data)[colnames(object@meta.data) %in% 
        paste0(names(marker.list), 1)] <- names(marker.list)

    # Heirarchical classification: general.class > subtypes
    general.class <- general.class[general.class %in% colnames(object@meta.data)]
    trunc.meta.gen <- object@meta.data[general.class]
    gen.thrs.vec <- numeric_threshold[general.class]
    call.res <- .modScoreCall(trunc.meta.gen, gen.thrs.vec, reject = "unknown")
    call.res$Barcode <- rownames(call.res)

    if (multi.lvl) {
        for (k in 1:ncol(lvl.df)) {

            sub.class.call <- list()
            store.sub.class <- lvl.df[[k]][!is.na(lvl.df[[k]])]
            parent.class <- unique(gsub("(.*)-(.*)", "\\1", store.sub.class))

            for (parent in parent.class) {
                sub.class <- store.sub.class[grepl(parent, store.sub.class)]
                children_class <- gsub("(.*)-(.*)", "\\2", sub.class)
                parents <- call.res$Barcode[call.res$MS_Celltype == 
                  parent]
                trunc.meta.parent <- object@meta.data[parents, 
                  ] %>% select(children_class)

                  gap_ind <- which(names(figures) == parent)

                  # Stop hierarchical classification in case no parent cell can be called
                  if (nrow(trunc.meta.parent) == 0){
                      stop(paste0("No ",parent," can be called in ","level ",k-1," classification, try setting more lenient thresholds"))}

                for (child in children_class) {
                  plot.title <- paste("Density plot for", child, 
                    "Module Scores within", parent, "population", 
                    sep = " ")

                    figures <- append(figures, list(NA), after = gap_ind)

                  figures[[gap_ind+1]] <- ggplot(trunc.meta.parent, aes_string(x = child)) + geom_density() + ggtitle(plot.title) + geom_vline(xintercept = numeric_threshold[child], linetype = "dashed", color = "red3") + theme_classic()
                  names(figures)[gap_ind+1] <- child
                  }

                trunc.meta.no.parent <- call.res[!call.res$MS_Celltype == 
                  parent, ]
                non.parent <- rownames(trunc.meta.no.parent)
                child.thres.vec <- numeric_threshold[children_class]

            sub.class.call[[match(parent, parent.class)]] <- .modScoreCall(trunc.meta.parent, child.thres.vec, reject = parent) %>% select(MS_Celltype)}

            sub.class.call <- do.call(rbind, sub.class.call)
            sub.class.call$Barcode <- rownames(sub.class.call)
            call.res$temp.call <- sub.class.call$MS_Celltype[match(call.res$Barcode, sub.class.call$Barcode)]
            call.res <- call.res %>% mutate(MS_Celltype = case_when(is.na(temp.call) ~ MS_Celltype, TRUE ~ temp.call))
            call.res$temp.call <- NULL
        }
    }

    object@meta.data$MS_Celltype <- call.res$MS_Celltype[match(object@meta.data$Barcode, call.res$Barcode)]

    lapply(figures, plot)

    return(object)
}

######### Node Execution Steps ##########
print("template_ModScoreSO.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_CellTypesSO
var_CellTypesSO<-readRDS(paste0(rds_output,"/var_CellTypesSO.rds"))

if (!('Seurat' %in% class(var_CellTypesSO))) { if (!(class(var_CellTypesSO) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_CellTypesSO <- as.data.frame(var_CellTypesSO)}}
#############################


# Processing input variable: var_My_GeneListsForModScores
var_My_GeneListsForModScores<-readRDS(paste0(rds_output,"/var_My_GeneListsForModScores.rds"))

if (!('Seurat' %in% class(var_My_GeneListsForModScores))) { if (!(class(var_My_GeneListsForModScores) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_My_GeneListsForModScores <- as.data.frame(var_My_GeneListsForModScores)}}
#############################


# Processing input variable: var_My_MultiLevelClassesForModScores
var_My_MultiLevelClassesForModScores<-readRDS(paste0(rds_output,"/var_My_MultiLevelClassesForModScores.rds"))

if (!('Seurat' %in% class(var_My_MultiLevelClassesForModScores))) { if (!(class(var_My_MultiLevelClassesForModScores) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_My_MultiLevelClassesForModScores <- as.data.frame(var_My_MultiLevelClassesForModScores)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ModScoreSO<-ModScoreSO(var_CellTypesSO,var_My_GeneListsForModScores,var_My_MultiLevelClassesForModScores)
saveRDS(var_ModScoreSO, paste0(rds_output,"/var_ModScoreSO.rds"))
