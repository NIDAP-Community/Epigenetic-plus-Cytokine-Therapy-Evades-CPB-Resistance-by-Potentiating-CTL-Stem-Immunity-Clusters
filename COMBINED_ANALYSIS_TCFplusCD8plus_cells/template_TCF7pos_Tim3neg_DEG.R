# DE with Find Markers [CCBR] [scRNA-seq] (54b6dd44-e233-4fb8-9bae-6ab0cb46e399): v100
TCF7pos_Tim3neg_DEG <- function(RSO_TCF7pos_Tim3neg,RSN_TCF7pos_Tim3neg,RMDT_TCF7pos_Tim3neg) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    seurat_object <- RSO_TCF7pos_Tim3neg
    metadata_table <- RMDT_TCF7pos_Tim3neg

# Samples and Assay parameters
    samples <- 'c("LN_aPD1","LN_ENT","LN_ENT+aPD1","LN_ENT+N803","LN_ENT+N803+aPD1","LN_N803","LN_N803+aPD1","LN_PBS","T_aPD1","T_ENT","T_ENT+aPD1","T_ENT+N803","T_ENT+N803+aPD1","T_N803","T_N803+aPD1","T_PBS")'
    assay_to_use <- "SCT"

# Contrasts parameters
    parameter_to_test <- "orig_ident"
    contrasts <- c("LN_N803-LN_PBS","LN_aPD1-LN_PBS","LN_ENT+N803-LN_PBS","LN_ENT+N803+aPD1-LN_PBS","LN_ENT+N803+aPD1-LN_N803","LN_ENT+N803+aPD1-LN_aPD1","LN_ENT+N803+aPD1-LN_ENT+aPD1","LN_ENT+N803+aPD1-LN_ENT+N803","LN_ENT+N803+aPD1-LN_N803+aPD1","LN_N803+aPD1-LN_PBS","T_ENT-T_PBS","T_N803-T_PBS","T_aPD1-T_PBS","T_ENT+aPD1-T_PBS","T_ENT+N803-T_PBS","T_N803+aPD1-T_PBS","T_ENT+N803+aPD1-T_PBS","T_ENT+N803+aPD1-T_ENT","T_ENT+N803+aPD1-T_N803","T_ENT+N803+aPD1-T_aPD1","T_ENT+N803+aPD1-T_ENT+aPD1","T_ENT+N803+aPD1-T_ENT+N803","T_ENT+N803+aPD1-T_N803+aPD1")

# Test parameters
    test_to_use <- "MAST"
    log_fc_threshold <- 0

## Presetting to FALSE
    use_spark <- FALSE

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##
    
    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

   ## Load SO 
        cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:         fs <- seurat_object$fileSystem()
# auto removed:         path <- fs$get_path("seurat_object.rds", 'r')
        SO <- seurat_object
        print(SO)

results <- degGeneExpressionMarkers(object = SO,
                                     samples = samples,
                                     contrasts = contrasts,
                                     parameter.to.test = parameter_to_test,
                                     test.to.use = test_to_use,
                                     log.fc.threshold = log_fc_threshold,
                                     use.spark = use_spark,
                                     assay.to.use = assay_to_use)

### ----------------------------------------------- ###
### Josh Edit 2024-08-20 through 2024-08-30 - Begin ###
### ----------------------------------------------- ###

### Purpose of this edit is to bring output table in-line with 
### column name expectations from GSEA downstream. Column names 
### should be more similar to those from Bulk DEG Analysis output
### when we are done.

# Original column names
original_colnames <- colnames(results$df)

# Function to rename columns based on pattern matching
rename_columns <- function(colnames, pattern, suffix) {
  # Identify columns that match the pattern
  matched <- grepl(pattern, colnames)
  
  # Apply the transformation only to the matching columns
  colnames[matched] <- gsub(
    pattern = paste0("^", pattern),             # Match the specific pattern at the start
    replacement = "C_",                         # Replace pattern with "C_"
    x = colnames[matched]
  )
  
  # Append the suffix to the columns that were matched
  colnames[matched] <- gsub(
    pattern = "C_(.*)_vs_(.*)",                 # Capture the parts after "C_"
    replacement = paste0("C_\\1_vs_\\2", suffix), # Reconstruct the name and append the suffix
    x = colnames[matched]
  )
  
  return(colnames)
}

# Apply the renaming function for various patterns with their corresponding suffixes
new_colnames <- original_colnames
new_colnames <- rename_columns(new_colnames, "p_val_adj_", "_adjpval")
new_colnames <- rename_columns(new_colnames, "avg_log2FC_", "_logFC")
new_colnames <- rename_columns(new_colnames, "pct.1_", "_pct1")
new_colnames <- rename_columns(new_colnames, "pct.2_", "_pct2")
new_colnames <- rename_columns(new_colnames, "p_val_", "_pval")

# Update the column names in the dataframe
colnames(results$df) <- new_colnames

# # Function to calculate t-statistics after removing rows with NA values
# calculate_tstat <- function(df, logfc_col, pval_col) {
#   logFC <- df[[logfc_col]]
#   pval <- df[[pval_col]]
  
#   # Ensure columns are numeric
#   logFC <- as.numeric(logFC)
#   pval <- as.numeric(pval)
  
#   # Remove rows with NA values in either column
#   valid_indices <- !is.na(logFC) & !is.na(pval)
#   logFC <- logFC[valid_indices]
#   pval <- pval[valid_indices]
  
#   # Calculate t-statistic
#   tstat <- sign(logFC) * sqrt(qchisq(1 - pval, df = 1))
  
#   # Initialize the tstat column with NA values and fill in the valid rows
#   tstat_col <- rep(NA, nrow(df))
#   tstat_col[valid_indices] <- tstat
  
#   return(tstat_col)
# }

# # Identify contrasts from column names, excluding non-contrast columns
# all_colnames <- colnames(results$df)
# contrast_patterns <- c("pval_", "logFC_", "pct1_", "pct2_", "adjpval_")
# contrast_cols <- all_colnames[grepl("^C_.*_(pval|logFC|pct1|pct2|adjpval)$", all_colnames)]
# contrasts <- unique(gsub("C_(.*?)_(pval|logFC|pct1|pct2|adjpval)$", "\\1", contrast_cols))

# # Exclude the "Gene" column from contrasts
# contrasts <- contrasts[contrasts != "Gene"]

# # Loop through each contrast and calculate the t-statistic
# for (contrast in contrasts) {
#   logfc_col <- paste0("C_", contrast, "_logFC")
#   pval_col <- paste0("C_", contrast, "_pval")
  
#   # Ensure columns exist in the dataframe
#   if (logfc_col %in% colnames(results$df) && pval_col %in% colnames(results$df)) {
#     cat("Processing contrast:", contrast, "\n")
    
#     # Calculate t-statistic
#     tstat <- calculate_tstat(results$df, logfc_col, pval_col)
    
#     # Add t-statistic column to the dataframe
#     results$df[[paste0("C_", contrast, "_tstat")]] <- tstat
#   } else {
#     warning(paste("Columns", logfc_col, "or", pval_col, "do not exist in the dataframe."))
#   }
# }

# View the updated dataframe with updated column names
head(results$df)

# View the updated column names
colnames(results$df)

### --------------------------------------------- ###
### Josh Edit 2024-08-20 through 2024-08-30 - End ###
### --------------------------------------------- ###

    return(results$df)

}

######### Node Execution Steps ##########
print("template_TCF7pos_Tim3neg_DEG.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_RSO_TCF7pos_Tim3neg
var_RSO_TCF7pos_Tim3neg<-readRDS(paste0(rds_output,"/var_RSO_TCF7pos_Tim3neg.rds"))

if (!('Seurat' %in% class(var_RSO_TCF7pos_Tim3neg))) { if (!(class(var_RSO_TCF7pos_Tim3neg) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RSO_TCF7pos_Tim3neg <- as.data.frame(var_RSO_TCF7pos_Tim3neg)}}
#############################


# Processing input variable: var_RSN_TCF7pos_Tim3neg
var_RSN_TCF7pos_Tim3neg<-readRDS(paste0(rds_output,"/var_RSN_TCF7pos_Tim3neg.rds"))

if (!('Seurat' %in% class(var_RSN_TCF7pos_Tim3neg))) { if (!(class(var_RSN_TCF7pos_Tim3neg) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RSN_TCF7pos_Tim3neg <- as.data.frame(var_RSN_TCF7pos_Tim3neg)}}
#############################


# Processing input variable: var_RMDT_TCF7pos_Tim3neg
var_RMDT_TCF7pos_Tim3neg<-readRDS(paste0(rds_output,"/var_RMDT_TCF7pos_Tim3neg.rds"))

if (!('Seurat' %in% class(var_RMDT_TCF7pos_Tim3neg))) { if (!(class(var_RMDT_TCF7pos_Tim3neg) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_RMDT_TCF7pos_Tim3neg <- as.data.frame(var_RMDT_TCF7pos_Tim3neg)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_TCF7pos_Tim3neg_DEG<-TCF7pos_Tim3neg_DEG(var_RSO_TCF7pos_Tim3neg,var_RSN_TCF7pos_Tim3neg,var_RMDT_TCF7pos_Tim3neg)
saveRDS(var_TCF7pos_Tim3neg_DEG, paste0(rds_output,"/var_TCF7pos_Tim3neg_DEG.rds"))
