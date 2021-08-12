# Aims:
# - Analyse the Djurec NPF/CAF data in the same way as the ECM-Targeted RNA-Seq panel data

# Notes:
# - Experiment on the GEO: GSE106901
# - The data provided by the group is in FPKM format which is not compatible with the Qiaseq data.
# - Therefore original counts were obtained from the SRA


## Raw Counts ----
library("Rsubread")
#Raw counts are extracted using a standard Rsubread pipeline
buildindex(basename = "mouse_index", reference = "./Data/DjurecData/GRCm38.primary_assembly.genome.fa.gz")
fastqFiles <- list(
  "./Data/DjurecData/SRR6292394.fastq",
  "./Data/DjurecData/SRR6292395.fastq",
  "./Data/DjurecData/SRR6292396.fastq",
  "./Data/DjurecData/SRR6292397.fastq",
  "./Data/DjurecData/SRR6292399.fastq",
  "./Data/DjurecData/SRR6292401.fastq",
  "./Data/DjurecData/SRR6292402.fastq",
  "./Data/DjurecData/SRR6292398.fastq"
)
align(index = "mouse_index", readfile1 = fastqFiles, nthreads = 10)
bamFiles <- c(
  "./Data/DjurecData/SRR6292394.fastq.subread.BAM",
  "./Data/DjurecData/SRR6292395.fastq.subread.BAM",
  "./Data/DjurecData/SRR6292396.fastq.subread.BAM",
  "./Data/DjurecData/SRR6292397.fastq.subread.BAM",
  "./Data/DjurecData/SRR6292399.fastq.subread.BAM",
  "./Data/DjurecData/SRR6292401.fastq.subread.BAM",
  "./Data/DjurecData/SRR6292402.fastq.subread.BAM",
  "./Data/DjurecData/SRR6292398.fastq.subread.BAM"
)
features <- featureCounts(
  bamFiles,
  isGTFAnnotationFile = T,
  annot.ext = "./Data/DjurecData/Mus_musculus.GRCm38.101.gtf.gz",
  GTF.attrType = "gene_id",
  nthreads = 10
)
djurecRawCounts <- data.frame(features$counts)
write.csv(djurecRawCounts, file = "./Data/Djurec NPFCAF Count Matrix.csv")


## DESeq2 Differential Expression Analysis ----
#The NPF/CAF data is normalised using DESeq2, the same way as the Qiaseq data.
library("DESeq2")
djurecRawCounts <- read.csv("./Data/Djurec NPFCAF Count Matrix.csv", stringsAsFactors = FALSE, row.names = 1)
djurecMetadata <- data.frame(
  row.names = colnames(djurecRawCounts),
  condition = factor(c("NPF", "NPF", "NPF", "CAF", "CAF", "CAF", "CAF", "CAF"))
)
djurecDESeq <- DESeqDataSetFromMatrix(djurecRawCounts, colData = djurecMetadata, design = ~ condition)
djurecDESeq <- DESeq(djurecDESeq)

#The contrast argument provides the correct CAF vs NPF comparison
djurecDERes <- data.frame(results(djurecDESeq, contrast = c("condition", "CAF", "NPF")))

#Gene symbols are pulled from biomaRt and merged with the DE results.
library("biomaRt")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
conversion <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
  filters = "ensembl_gene_id",
  values = rownames(djurecRawCounts),
  mart = mouse, uniqueRows = TRUE
)
conversion <- conversion[!duplicated(conversion$ensembl_gene_id), ]
row.names(conversion) <- conversion$ensembl_gene_id
djurecDERes <- merge(conversion, djurecDERes, by = "row.names")
djurecDERes <- data.frame(djurecDERes, row.names = 1)

#The DESeq2-normalised counts are extracted for the heatmap.
djurecNorm <- data.frame(counts(djurecDESeq, normalized = TRUE))

#Results are saved to a CSV for later use.
write.csv(djurecDERes, "./Data/Djurec NPFCAF DE Results.csv")
write.csv(djurecNorm, "./Data/Djurec NPFCAF Normalised Counts.csv")
