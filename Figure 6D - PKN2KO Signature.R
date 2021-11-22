# Aims:
# - Create the PKN2KO signature
# - Convert the PKN2KO signature to human equivalent


## Preparation ----
#Required packages are loaded
library("biomaRt")

#The required data is loaded
pscDE <- read.csv("./Data/PSCs Basal DE Results.csv", stringsAsFactors = FALSE, row.names = 1)
tumDE <- read.csv("./Data/InVivo Tumour DE Results.csv", stringsAsFactors = FALSE, row.names = 1)

#Only the significant genes are selected
pscDE_signif <- pscDE[pscDE$padj < 0.05 & !is.na(pscDE$padj), ]



##Creation of the signature
#Selects only genes that concur between the PSC and tumour data
selection <- ifelse(
  pscDE_signif[, "log2FoldChange"] < 0 &
    tumDE[match(rownames(pscDE_signif), rownames(tumDE)), "log2FoldChange"] < 0, TRUE,
  ifelse(pscDE_signif[, "log2FoldChange"] > 0 &
           tumDE[match(rownames(pscDE_signif), rownames(tumDE)), "log2FoldChange"] > 0, TRUE,
         FALSE)
)

#Selected the genes
currentSig <- pscDE_signif[selection, ]
currentSig <- currentSig[!is.na(currentSig$mgi_symbol), ]
#Selects signature based on combined stat from PSC and tumour data
PKN2KOsig <- currentSig[
  abs(currentSig$stat) > 3.40 &
  abs(tumDE[match(rownames(currentSig), rownames(tumDE)), "stat"]) > 1.23,
  c("mgi_symbol", "description")]

#Adds the ensembl gene id to it's own column
PKN2KOsig$mouse_ensembl_gene_id <- row.names(PKN2KOsig)
row.names(PKN2KOsig) <- NULL



##Converting MGI symbols to HGNC Symbols
#Mouse and human ensemble marts
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Uses getLDS to convert from mouse to human
conversionTable <- getLDS(
  attributes = c("ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = PKN2KOsig$mouse_ensembl_gene_id,
  mart = mouse,
  attributesL = c("ensembl_gene_id", "hgnc_symbol"),
  martL = human,
  uniqueRows=T
)

#Rename columns to match PKN2KOsig table
colnames(conversionTable) <- c("mouse_ensembl_gene_id", "human_ensembl_gene_id", "hgnc_symbol")
#Human names and Ids are merged with the mouse data from PKN2KOsig
PKN2KOsig <- merge(PKN2KOsig, conversionTable, by = "mouse_ensembl_gene_id")

write.csv(PKN2KOsig, "./Data/PKN2KO Signature.csv", row.names = FALSE)
