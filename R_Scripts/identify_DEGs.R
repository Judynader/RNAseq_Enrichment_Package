# Function to identify differentially expressed genes (DEGs)
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
