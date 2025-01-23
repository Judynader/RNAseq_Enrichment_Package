#Install necessary packages
packageVec <- c(
	"oligo",
	"SummarizedExperiment",
	"DESeq2",
	"limma",
	"genefilter",
	"pheatmap",
	"edgeR",
	"cluster",
	"Biobase",
	"BiocGenerics",
	"bioDist",
	"IRanges",
	"GenomicRanges",
	"openxlsx",
	"data.table",
	"ggplot2",
	"rmarkdown",
	"pd.hugene.1.0.st.v1",
	"hugene10sttranscriptcluster.db",
	"ReactomePA",
	"enrichR",
	"biomaRt",
	"org.Hs.eg.db",
	"org.Ss.eg.db",
	"devtools"
)

if(!require(BiocManager)) install.packages("BiocManager")
for(p in packageVec)
  if(!suppressWarnings(require(p, character.only = TRUE)))
    BiocManager::install(p, update = FALSE)
devtools::install_gitlab(repo = "wolftower/biostudies")
