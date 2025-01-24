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
#' @import edgeR
#'
#' @export

ImportData <- function(counts_file, sample_file,
                       fdr_threshold = 0.05, log2fc_threshold = 1) {

  # Load necessary libraries
  library(edgeR)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(openxlsx)
  library(DESeq2)
  library(ggplot2)

  #Step 1: Import data and make factors
  counts_table <- read.table(counts_file,
                             header = T,
                             row.names = 1,
                             sep = '\t')

  sample_table <- read.table(sample_file,
                             header = T,
                             row.names = 1,
                             sep = '\t')

  #Creating factor variables
  disease <- factor(sample_table$disease)
  sex <- factor(sample_table$sex)
  individual <- factor(sample_table$individual)

  #Step 2: Filter low values
  meanlog2CPM <- rowMeans(log2(cpm(counts_table + 1)))
  counts_table <- counts_table[meanlog2CPM > 1, ]

  #Step 3: Create DESeq Data Set
  dds <- DESeqDataSetFromMatrix(countData = counts_table,
                                colData = sample_table,
                                design = ~ disease)
  #Step 4: Perform DEG analysis
  dds2 <- DESeq(dds)
  res <- results(dds2)

  # Step 5: Filter results based on FDR and log2 fold change
  filtered_res <- res[which(res$padj < fdr_threshold & abs(res$log2FoldChange) > log2fc_threshold), ]

  # Step 6: Convert gene symbols to Entrez IDs for KEGG analysis
  gene_ids <- rownames(filtered_res)
  entrez_ids <- bitr(gene_ids,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

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
  write.xlsx(data.frame(filtered_res), rowNames = T, file = "DEGs.xlsx")

  # Step 9: Create plots to visualize the results of the enrichment analysis
  # Saves dotplot, one for go and one for kegg
  dotplot(go_results) + ggtitle("GO Enrichment Analysis")
  ggsave("go_enrichment_dotplot.png")

  dotplot(kegg_results) + ggtitle("KEGG Enrichment Analysis")
  ggsave("kegg_enrichment_dotplot.png")

  # Step 10: Create a similarity matrix and treeplot for GO results
  # saves treeplot made by go
  go_sim <- pairwise_termsim(go_results)
  treeplot(go_sim) + ggtitle("GO Term Similarity Tree")
  ggsave("go_term_similarity_tree.png")
}
