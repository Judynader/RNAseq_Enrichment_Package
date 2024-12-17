#' RNA-seq Analysis and Enrichment
#'
#' This function performs RNA-seq data analysis including filtering low-expressed genes,
#' differential expression analysis, and over-representation analysis using GO and KEGG.
#'
#' @param counts_file A string representing the path to the RNA-seq count data CSV file.
#' @param sample_table_file A string representing the path to the sample table CSV file.
#' @param output_dir A string representing the directory where output files will be saved.
#' @param fdr_threshold A numeric value for the FDR threshold for filtering DE genes (default is 0.05).
#' @param log2fc_threshold A numeric value for the log2 fold change threshold for filtering DE genes (default is 1).
#'
#' @import DESeq2
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import openxlsx
#' @import ggplot2
#' @import enrichplot
#'
#' @export
rna_seq_analysis <- function(counts_file, sample_table_file, output_dir,
                             fdr_threshold = 0.05, log2fc_threshold = 1) {

  # Load necessary libraries
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(openxlsx)
  library(ggplot2)
  library(enrichplot)

  # Step 1: Import RNA-seq count data
  counts <- read.csv(counts_file, row.names = 1)
  sample_table <- read.csv(sample_table_file)

  # Step 2: Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sample_table,
                                design = ~ condition)  # Assuming 'condition' is the column in sample_table

  # Step 3: Filter low-expressed genes
  dds <- dds[rowSums(counts(dds)) >= 10, ]

  # Step 4: Perform differential expression analysis
  dds <- DESeq(dds)
  res <- results(dds)

  # Step 5: Filter results based on FDR and log2 fold change
  filtered_res <- res[which(res$padj < fdr_threshold & abs(res$log2FoldChange) > log2fc_threshold), ]

  # Step 6: Convert gene symbols to Entrez IDs for KEGG analysis
  gene_ids <- rownames(filtered_res)
  entrez_ids <- bitr(gene_ids, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Step 7: Perform over-representation analysis
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

  # Step 8: Export the table with differentially expressed genes to Excel
  write.xlsx(as.data.frame(filtered_res), file.path(output_dir, "differentially_expressed_genes.xlsx"))

  # Step 9: Create plots to visualize the results of the enrichment analysis
  dotplot(go_results) + ggtitle("GO Enrichment Analysis")
  ggsave(file.path(output_dir, "go_enrichment_dotplot.png"))

  dotplot(kegg_results) + ggtitle("KEGG Enrichment Analysis")
  ggsave(file.path(output_dir, "kegg_enrichment_dotplot.png"))

  # Step 10: Create a similarity matrix and treeplot for GO results
  go_sim <- pairwise_termsim(go_results)
  treeplot(go_sim) + ggtitle("GO Term Similarity Tree")
  ggsave(file.path(output_dir, "go_term_similarity_tree.png"))
}
