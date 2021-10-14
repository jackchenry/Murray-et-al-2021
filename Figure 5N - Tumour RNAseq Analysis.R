# Aims:
# - Perform differential expression analysis on the pancreatic tumour data

# Variables
# - Data is either normal pancreas or pancreatic tumour samples - We are only interested in the pancreatic tumour samples
# - Samples are either Wt or KO with respects to PKN2



## Preparation ----
#Load the required packages
library("DESeq2")
library("biomaRt")

#Read data into the environment
raw <- read.csv("./Data/Raw InVivo Counts.csv", stringsAsFactors = FALSE, row.names = "Gene_ID")



##Differential Expression Analysis ----
metadata <- data.frame(
  sample = colnames(raw),
  genotype = ifelse(substr(colnames(raw), 1, 2) == "WT", "WT", "KO"),
  type = c("Tum", "Panc", "Tum", "Tum", "Panc", "Panc", "Tum", "Panc", "Tum", "Tum", "Tum", "Tum", "Panc", "Tum", "Tum", "Panc", "Tum")
  )

#Creating a count matrix for just the tumour samples
tumIDs <- data.frame(metadata[metadata$type == "Tum", c("sample", "genotype")])
tumIDs$genotype <- factor(tumIDs$genotype, levels = c("WT", "KO"))
tumCounts <- raw[ , tumIDs$sample]

#DESeq2 data set created
tumDESeqData <- DESeqDataSetFromMatrix(tumCounts, colData = tumIDs, design = ~ genotype)

#Differential expression is performed
tumDESeqData <- DESeq(tumDESeqData)
tumDEResults <- results(tumDESeqData)

#Normalised counts are extracted and saved
tumNormalisedCounts <- counts(tumDESeqData, normalized = TRUE)
write.csv(tumNormalisedCounts, file = "./Data/InVivo Tumour Normalised Count Matrix.csv")

#The raw tumour counts are also saved in a separate file
tumRawCounts <- counts(tumDESeqData, normalized = FALSE)
write.csv(tumRawCounts, file = "./Data/InVivo Tumour Raw Count Matrix.csv" )



##Gene Information ----
#Gene names are downloaded using biomaRt
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
geneInfo <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                  filters = "ensembl_gene_id", values = rownames(tumDEResults), mart = mouse, uniqueRows = TRUE)

#Gene information is merged with the differential expression results
geneInfo <- geneInfo[!duplicated(geneInfo$ensembl_gene_id), ]
rownames(geneInfo) <- geneInfo$ensembl_gene_id
geneInfo$ensembl_gene_id <- NULL

tumResults <- merge(geneInfo, data.frame(tumDEResults), by = "row.names", all.y = TRUE)
rownames(tumResults) <- tumResults$Row.names
tumResults$Row.names <- NULL

#Counts are ordered by stat
tumResults <- tumResults[order(tumResults$stat), ]



##Human Orthologues Information
#Human biomart created
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Obtains human ortholougue information using LDS
tumConversionTable <- getLDS(
  attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = rownames(tumResults), mart = mouse,
  attributesL = c("ensembl_gene_id", "hgnc_symbol"), martL = human, uniqueRows = TRUE)

#Column names are renamed
colnames(tumConversionTable) <- c("mouse_ensembl_gene_id", "human_ensembl_gene_id", "hgnc_symbol")

#Duplicated ids are removed
tumConversionTable <- tumConversionTable[!duplicated(tumConversionTable$mouse_ensembl_gene_id), ]

#Human annotated data is merged with the tumour differential expression data
tumResults$mouse_ensembl_gene_id <- rownames(tumResults)
tumResults <- merge(tumResults, tumConversionTable, by = "mouse_ensembl_gene_id", all.y = TRUE)

#Data is saved into a csv for later use
write.csv(tumResults, file = "./Data/InVivo Tumour DE Results.csv")
