# Trajectory Analysis with TSCAN [CCBR] [scRNA-seq] (5a69e946-b73a-4cfe-96e0-54922455583b): v11
Traj_2 <- function(ReclusterSO_1,MetadataTable_8) {
 ## Trajectory analysis
   
 ## Basic Parameters:
    seurat_object = ReclusterSO_1
    UserFacing_MetaData = MetadataTable_8
    SelectedColumnName <- "SCT_snn_res_0_2"
 ## Advanced Parameters:
    seurat_object_filename <- "seurat_object.rds"
    custom_gene_list <- "Tcf7 Sell Lef1 Ccl5 Gzmb Nkg7 S100a4 S100a6 Prf1 Il7r Mki67"

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")

    library(Seurat)
    library(ggplot2)
    library(grid)
    library(gridExtra)

    library(scran)
    library(scater)
    library(TSCAN)
    library(mclust)
    #library(igraph)

 ## Modifying Custom Gene List:
 ## Replace commas with spaces and split the string
split_values <- unlist(strsplit(gsub(",", " ", custom_gene_list), " "))
custom_gene_list <- split_values[split_values != ""]

 ## Input SO.
    cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:     path <- nidapGetPath(seurat_object,seurat_object_filename)
    so <- seurat_object

## Workaround Dots and Underscores
RealMetadata <- so@meta.data
LABEL <- colnames(RealMetadata[match(SelectedColumnName,names(UserFacing_MetaData))-1])

## Cleaning a gene List
GeneList <- rownames(so@assays$RNA@counts)
print("The following genes were used:")
print(intersect(custom_gene_list, GeneList)) 

print("The following genes were not found:")
print(setdiff(custom_gene_list, GeneList))

custom_gene_list <- intersect(custom_gene_list, GeneList)

## Part I: Trajectory

Idents(so) <- LABEL
dimplot <- DimPlot(so)
sce <- as.SingleCellExperiment(so)
colLabels(sce) <- colData(sce)[[LABEL]] 

by.cluster <- aggregateAcrossCells(sce, ids=colLabels(sce))
centroids <- reducedDim(by.cluster, "PCA")

# Set clusters=NULL as we have already aggregated above.
mst <- TSCAN::createClusterMST(centroids, clusters=NULL)

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")
plot_tsne1 <- plotTSNE(sce, colour_by="label") + geom_line(data=line.data, mapping=aes(x=tSNE_1, y=tSNE_2, group=edge))

# Lets also try UMAP:
line.data.umap <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP")
plot_umap1 <- plotUMAP(sce, colour_by="label") + geom_line(data=line.data.umap, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="PCA")
tscan.pseudo <- orderCells(map.tscan, mst)
##head(tscan.pseudo)

common.pseudo <- averagePseudotime(tscan.pseudo) 

plot_tsne2 <- plotTSNE(sce, colour_by=I(common.pseudo), 
                       text_by="label", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=tSNE_1, y=tSNE_2, group=edge)) +
  ggtitle("TSCAN-based MST")

plot_umap2 <- plotUMAP(sce, colour_by=I(common.pseudo), 
                       text_by="label", text_colour="red") +
  geom_line(data=line.data.umap, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge)) +
  ggtitle("TSCAN-based MST")

## Part II: Pseudotime

pseudo.all <- TSCAN::quickPseudotime(sce, use.dimred="PCA")
#head(pseudo.all$ordering)

pseudo.mnn <- TSCAN::quickPseudotime(sce, use.dimred="PCA", with.mnn=TRUE)
mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)
plot_tsne3 <- plotTSNE(sce, colour_by=I(mnn.pseudo), text_by="label", text_colour="red") +
  geom_line(data=pseudo.mnn$connected$TSNE, mapping=aes(x=tSNE_1, y=tSNE_2, group=edge)) +
  ggtitle("TSCAN with MNN distances-based MST")
plot_umap3 <- plotUMAP(sce, colour_by=I(mnn.pseudo), text_by="label", text_colour="red") +
  geom_line(data=pseudo.mnn$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge)) +
  ggtitle("TSCAN with MNN distances-based MST")
              

## Part III
### From Chapter 10.3.2, Changes along a trajectory
sce$Pseudotime <- pathStat(tscan.pseudo)[,1]
#sce$TSCAN.first <- pathStat(tscan.pseudo)[,1]
#sce$TSCAN.second <- pathStat(tscan.pseudo)[,2]
pseudo <- TSCAN::testPseudotime(sce, pseudotime=tscan.pseudo[,1])[[1]]
sorted <- pseudo[order(pseudo$p.value),]

up.left <- sorted[sorted$logFC < 0,]
up.left$SYMBOL <- rownames(up.left)
best <- head(up.left$SYMBOL, 10)

TopDecrease <- plotExpression(sce, features=best, x="Pseudotime", colour_by="label", ncol = 5) + ggtitle("Expression of the top 10 genes that decrease in expression with increasing pseudotime")

up.right <- sorted[sorted$logFC > 0,]
up.right$SYMBOL <- rownames(up.right)
best <- head(up.right$SYMBOL, 10)
TopIncrease <- plotExpression(sce, features=best, x="Pseudotime", colour_by="label", ncol = 5) + ggtitle("Expression of the top 10 genes that increase in expression with increasing pseudotime")

on.first.path <- !is.na(sce$Pseudotime)

## Printing Pretty Pics:
print(dimplot)

set.seed(1)
plot(pseudo.all$mst) 

print(plot_tsne1)
grid.arrange(plot_tsne2, plot_tsne3, ncol = 2)

print(plot_umap1)
grid.arrange(plot_umap2, plot_umap3, ncol = 2)

plot(TopDecrease)
plot(TopIncrease)
plotHeatmap(sce[,on.first.path], order_columns_by="Pseudotime", colour_columns_by="label", features=head(up.right$SYMBOL, 50), center=TRUE, main="Expression of the top 50 genes that increase in expression with increasing pseudotime")
plotHeatmap(sce[,on.first.path], order_columns_by="Pseudotime", colour_columns_by="label", features=head(up.left$SYMBOL, 50), center=TRUE, main="Expression of the top 50 genes that decrease in expression with increasing pseudotime")

if(length(custom_gene_list) > 0)
   {
     CustomScatter <- plotExpression(sce, features=custom_gene_list, x="Pseudotime", colour_by="label", ncol = 5) + 
          ggtitle("Expression of genes with increasing pseudotime")
     plot(CustomScatter)

     plotHeatmap(sce[,on.first.path],
          order_columns_by="Pseudotime",
          colour_columns_by="label",
          features=custom_gene_list,
          center=TRUE,
          main="Expression of genes with increasing pseudotime")
   }

# auto removed: return(NULL)

}

#################################################
## Global imports and functions included below ##
#################################################

# R.cache depends on a home folder
Sys.setenv(R_USER_CACHE_DIR = Sys.glob(file.path(R.home())));

######### Node Execution Steps ##########
print("template_Traj_2.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ReclusterSO_1
var_ReclusterSO_1<-readRDS(paste0(rds_output,"/var_ReclusterSO_1.rds"))

if (!('Seurat' %in% class(var_ReclusterSO_1))) { if (!(class(var_ReclusterSO_1) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ReclusterSO_1 <- as.data.frame(var_ReclusterSO_1)}}
#############################


# Processing input variable: var_MetadataTable_8
var_MetadataTable_8<-readRDS(paste0(rds_output,"/var_MetadataTable_8.rds"))

if (!('Seurat' %in% class(var_MetadataTable_8))) { if (!(class(var_MetadataTable_8) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_MetadataTable_8 <- as.data.frame(var_MetadataTable_8)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_Traj_2<-Traj_2(var_ReclusterSO_1,var_MetadataTable_8)
saveRDS(var_Traj_2, paste0(rds_output,"/var_Traj_2.rds"))
