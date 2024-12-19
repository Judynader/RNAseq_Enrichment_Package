#' RNA-seq Analysis and Enrichment
#'
#'@description
#'This function performs RNA-seq data analysis including filtering low-expressed genes,
#' differential expression analysis, and over-representation analysis using GO and KEGG. 
#' It filters DEGs based on a FDR value < 0.05 and absolute logFC value of > 1.3.
#'
#' @param counts_file A string representing the path to the RNA-seq count data text file.
#'                    Default is "E-MTAB-2523.counts.txt".
#' @param sample_table_file A string representing the path to the sample table text file.
#'                          Default is "E-MTAB-2523_sample table.txt".
#' @param output_dir A string representing the directory where output files will be saved.
#' @param fdr_threshold A numeric value for the FDR threshold for filtering DE genes.
#' @param log2fc_threshold A numeric value for the log2 fold change threshold for filtering DE genes.
#'
#' @import DESeq2
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import openxlsx
#' @import ggplot2
#' @import enrichplot
#'
#' @export
rna_seq_analysis <- function(counts_file = "E-MTAB-2523.counts.txt",
                             sample_table_file = "E-MTAB-2523_sample table.txt",
                             output_dir = "data/processed/",
                             fdr_threshold = 0.05,
                             log2fc_threshold = 1.3) {

  # Load necessary libraries
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(openxlsx)
  library(ggplot2)
  library(enrichplot)

  # Create the data/processed/output directory
  dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
  output_dir <- "data/processed/"

  # Step 1: Import RNA-seq count data
  counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")  # Tab-separated values
  sample_table <- read.table(sample_table_file, header = TRUE, row.names = 1, sep = "\t", fill = TRUE)  # Tab-separated values

  # Step 2: Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sample_table,
                                design = ~ disease)  # Assuming 'condition' is the column in sample_table

  # Step 3: Filter low-expressed genes
  dds <- dds[rowSums(counts(dds)) >= 10, ]

  # Step 4: Perform differential expression analysis
  dds <- DESeq(dds)
  res <- results(dds)

  # Step 5: Filter results based on FDR and log2 fold change
  filtered_res <- res[which(res$padj < fdr_threshold & abs(res$log2FoldChange) > log2fc_threshold), ]

  # Step 6: Calculate total logCPM total for each gen
  # Get normalised counts
  norm_counts <- counts(dds, normalized = TRUE)

  # Calculate the total CPM for each gene (average CPM over all samples)
  cpm_total <- rowMeans(norm_counts) / colSums(norm_counts) * 1e6

  # Calculate logCPM (log2 de CPM + 1)
  log_cpm_total <- log2(cpm_total + 1)

  # Step 7: Convert gene symbols to Entrez IDs for KEGG analysis
  gene_ids <- rownames(filtered_res)
  entrez_ids <- bitr(gene_ids, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Step 8: Perform over-representation analysis
  go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05)

  kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                             organism = 'hsa',
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05)

  # Step 9: Create the final results table table with the requested values.
  # The table will contain: Gene names, LogFC, logCPM, PValue y FDR.

  final_table <- data.frame(
    Gene = rownames(filtered_res),
    LogFC = filtered_res$log2FoldChange,
    logCPM = log_cpm_total[rownames(filtered_res)],
    PValue = filtered_res$pvalue,
    FDR = filtered_res$padj
  )

  # Step 10: Export the table with differentially expressed genes to Excel
  write.xlsx(final_table, file.path(output_dir, "differentially_expressed_genes.xlsx"))

  # Step 11: Create plots to visualize the results of the enrichment analysis
  dotplot(go_results) + ggtitle("GO Enrichment Analysis")
  ggsave(file.path(output_dir, "go_enrichment_dotplot.png"))

  dotplot(kegg_results) + ggtitle("KEGG Enrichment Analysis")
  ggsave(file.path(output_dir, "kegg_enrichment_dotplot.png"))


  # Step 12: Create a similarity matrix and treeplot for GO results
  go_sim <- pairwise_termsim(go_results)
  treeplot(go_sim) + ggtitle("GO Term Similarity Tree")
  ggsave(file.path(output_dir, "go_term_similarity_tree.png"))
}
