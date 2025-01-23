#' Filter low-expressed genes
#'
#' This function filters genes with low expression from the dataset.
#'
#' @param data A dataframe containing gene expression data.
#' @param threshold Numeric value for minimum expression threshold. Default is 10.
#'
#' @return A dataframe with filtered genes.
#' @export
filter_genes <- function(data, threshold = 10) {
    return(data[data$expression >= threshold, ])
}
