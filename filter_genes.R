# Function to filter low-expressed genes
filter_genes <- function(count_data, threshold) {
    filtered_data <- count_data[rowSums(count_data > threshold) > 0, ]
    return(filtered_data)
}
