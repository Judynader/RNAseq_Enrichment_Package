#' Title: Identify Differentially Expressed Genes (DEGs)
#' Description: This function performs differential expression analysis on RNA-seq count data using the edgeR package.
#' It identifies genes that are differentially expressed based on specified thresholds for FDR and logFC.
#' @param count_data A matrix or data frame containing the RNA-seq count data, where rows represent genes and columns represent samples.
#' @param group_labels A vector indicating the group labels for each sample (e.g., control vs treatment).
#' @param fdr_threshold A numeric value representing the false discovery rate threshold. Default is 0.05.
#' @param logFC_threshold A numeric value representing the log-fold change threshold. Default is 1.
#' @return A data frame containing the differentially expressed genes (DEGs) with associated statistics like FDR and logFC.
#' @examples
#' # Example usage
#' count_data <- read.csv("count_data.csv", row.names = 1)
#' group_labels <- c("control", "control", "treatment", "treatment")
#' degs <- identify_DEGs(count_data, group_labels)
identify_DEGs <- function(count_data, group_labels, fdr_threshold = 0.05, logFC_threshold = 1) {
    library(edgeR)
    # Create a DGEList object
    dge <- DGEList(counts = count_data, group = group_labels)
    # Normalize data
    dge <- calcNormFactors(dge)
    # Estimate dispersion
    dge <- estimateDisp(dge)
    # Fit a generalized linear model
    fit <- glmFit(dge)
    # Perform likelihood ratio test
    results <- glmLRT(fit)
    # Extract top genes based on FDR and logFC
    degs <- topTags(results, n = Inf)$table
    degs <- degs[degs$FDR < fdr_threshold & abs(degs$logFC) > logFC_threshold, ]
    return(degs)
}

de_analysis <- function(count_data, group_labels) {
  library(edgeR)
  dge <- DGEList(counts = count_data, group = group_labels)
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge)
  fit <- glmFit(dge)
  results <- glmLRT(fit)
  degs <- topTags(results)$table
  return(degs)
}
