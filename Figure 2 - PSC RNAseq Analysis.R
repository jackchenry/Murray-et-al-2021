# Aims:
# - Analyse ECM-Targeted RNA-Seq panel data.
# - Look for differentially expressed genes.

# Variables
# - Samples are either Wt or KO with respects to PKN2
# - Samples are either untreated (basal) or treated with TGF-b

library("tidyr")
library("DESeq2")
library("biomaRt")

#Read data into the environment
csv <- read.csv("./Data/Raw Qiaseq Counts.csv", stringsAsFactors = FALSE)
#Only keep the interesting data - sample, ID and read columns
rawData <- data.frame(
  sample = as.factor(csv$read.set),
  gene = csv$gene.id,
  reads = as.numeric(gsub(",","",csv$reads))
)
#Remove gDNA controls (all 0 or 1)
rawData <- rawData[which(substr(rawData$gene, start = 1, stop = 3) == "ENS"), ]
#Data transformed into a count matrix
rawCounts <- data.frame(pivot_wider(rawData, names_from = sample, values_from = reads), row.names = 1)

#The reference genes provided by Qiaseq are used to normalise the sequencing results
referenceGenes <- c(
  "ENSMUSG00000034601",
  "ENSMUSG00000043323",
  "ENSMUSG00000024383",
  "ENSMUSG00000028651",
  "ENSMUSG00000022771",
  "ENSMUSG00000031706",
  "ENSMUSG00000039759",
  "ENSMUSG00000051390",
  "ENSMUSG00000033961"
)



## DESeq2 ----
#The sample group information is set up in a metadata dataframe
metadata <- data.frame(
  sampleID = colnames(rawCounts),
  condition = as.factor(c(
   "KO_Basal", "KO_Basal", "KO_Basal",
   "KO_TGF", "KO_TGF", "KO_TGF",
   "Wt_Basal", "Wt_Basal", "Wt_Basal",
   "Wt_TGF", "Wt_TGF", "Wt_TGF")),
  genotype = as.factor(c(
   "KO", "KO", "KO",
   "KO", "KO", "KO",
   "Wt", "Wt", "Wt",
   "Wt", "Wt", "Wt")),
  stimulation = as.factor(c(
   "Basal", "Basal", "Basal",
   "TGF", "TGF", "TGF",
   "Basal", "Basal", "Basal",
   "TGF", "TGF", "TGF")),
  stringsAsFactors = FALSE
)

#Differential expression analysis is performed using DESeq2
#DESeqData object created from raw count matrix and metadata
DESeqData <- DESeqDataSetFromMatrix(rawCounts, colData = metadata, design = ~ condition)

#Size factors are estimated using the reference genes
DESeqData <- estimateSizeFactors(DESeqData, controlGenes = (rownames(DESeqData) %in% referenceGenes))

#Other DESeq things happen
DESeqData <- DESeq(DESeqData)

#DESeq2-Normalised counts matrices are saved for use in heatmaps etc
normalisedCounts <- data.frame(counts(DESeqData, normalized = TRUE))
basalNormalisedCounts <- normalisedCounts[ , metadata[metadata$stimulation == "Basal", "sampleID"]]
tgfNormalisedCounts <- normalisedCounts[ , metadata[metadata$stimulation == "TGF", "sampleID"]]

write.csv(basalNormalisedCounts, file = "./Data/PSCs Basal Normalised Count Matrix.csv")
write.csv(tgfNormalisedCounts, file = "./Data/PSCs TGF Normalised Count Matrix.csv")

#DESeq calculates results. Results are obtained for each TGF treated and untreated condition.
basalDEResults <- results(DESeqData, c("condition", "KO_Basal", "Wt_Basal"))
tgfDEResults <- results(DESeqData, c("condition", "KO_TGF", "Wt_TGF"))



## Gene Information and Saving Results ----
#ID information downloaded using Biomart
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
geneInfo <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
  filters = "ensembl_gene_id",
  values = rownames(rawCounts),
  mart = mouse
)
geneInfo <- data.frame(geneInfo, row.names = "ensembl_gene_id")

#GeneInfo merged with DE results
basalResults <- merge(geneInfo, data.frame(basalDEResults), by = "row.names")
tgfResults <- merge(geneInfo, data.frame(tgfDEResults), by = "row.names")
basalResults <- data.frame(basalResults, row.names = "Row.names")
tgfResults <- data.frame(tgfResults, row.names = "Row.names")

#DE results ordered by stat
basalResults <- basalResults[order(basalResults$stat), ]
tgfResults <- tgfResults[order(tgfResults$stat), ]

#Results saved to a csv for later use
write.csv(basalResults, file = "./Data/PSCs Basal DE Results.csv")
write.csv(tgfResults, file = "./Data/PSCs TGF DE Results.csv")



## Log-Fold-Shrink Data ----
#Some log-fold shrink methods are applied for use in plotting the volcano plots
#DESeq is performed again using Wt-Basal as the reference level. This is done to set up the right coefs for the lfcShrink
metadata$condition <- relevel(metadata$condition, ref = "Wt_Basal")
#Similar DESeq things to before
basalShrinkDESeqData <- DESeqDataSetFromMatrix(rawCounts, colData = metadata, design = ~ condition)
basalShrinkDESeqData <- estimateSizeFactors(basalShrinkDESeqData,
                                            controlGenes = (rownames(basalShrinkDESeqData) %in% referenceGenes))
basalShrinkDESeqData <- DESeq(basalShrinkDESeqData)
#Results are shrunk using apeglm
basalShrinkRes <- lfcShrink(basalShrinkDESeqData, coef = "condition_KO_Basal_vs_Wt_Basal", type = "apeglm")
#Results are merged with gene info
basalShrinkRes <- merge(geneInfo, data.frame(basalShrinkRes), by = "row.names")
basalShrinkRes <- data.frame(basalShrinkRes, row.names = "Row.names")


#The same thing is done with the TGF data
metadata$condition <- relevel(metadata$condition, ref = "Wt_TGF")
#Similar DESeq things to before
tgfShrinkDESeqData <- DESeqDataSetFromMatrix(rawCounts, colData = metadata, design = ~ condition)
tgfShrinkDESeqData <- estimateSizeFactors(tgfShrinkDESeqData,
                                          controlGenes = (rownames(tgfShrinkDESeqData) %in% referenceGenes))
tgfShrinkDESeqData <- DESeq(tgfShrinkDESeqData)
#Results are shrunk using apeglm
tgfShrinkRes <- lfcShrink(tgfShrinkDESeqData, coef = "condition_KO_TGF_vs_Wt_TGF", type = "apeglm")
#Results are merged with gene info
tgfShrinkRes <- merge(geneInfo, data.frame(tgfShrinkRes), by = "row.names")
tgfShrinkRes <- data.frame(tgfShrinkRes, row.names = "Row.names")

#Results saved to a csv for later use
write.csv(basalShrinkRes, file = "./Data/PSCs Basal Shrunk DE Results.csv")
write.csv(tgfShrinkRes, file = "./Data/PSCs TGF Shrunk DE Results.csv")


