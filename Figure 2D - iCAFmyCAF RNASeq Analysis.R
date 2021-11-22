# Aims:
# Analyse the Tuveson iCAF/myCAF data in the same way as the ECM-Targeted RNA-Seq panel data

library("DESeq2") #Used to perform differential expression analysis

#Raw count matrix was obtained from the GEO: GSE93313
tuvesonCSV <- read.csv("./Data/TuvesonCAF.csv", stringsAsFactors = FALSE, row.names = 1)

#The iCAF/myCAF data is normalised using DESeq2, the same way as the Qiaseq data.
tuvesonMetadata <- data.frame(
  row.names = colnames(tuvesonCSV),
  condition = factor(c("2D", "2D", "3D", "3D", "Transwell", "Transwell", "Transwell", "Transwell"))
)
tuvesonDESeq <- DESeqDataSetFromMatrix(tuvesonCSV, colData = tuvesonMetadata, design = ~ condition)
tuvesonDESeq <- DESeq(tuvesonDESeq)

#The contrast argument provides the correct iCAF vs myCAF comparison
tuvesonDERes <- data.frame(results(tuvesonDESeq, contrast = c("condition", "Transwell", "2D")))

#The normalised counts are extracted for the heatmap.
tuvesonNorm <- data.frame(counts(tuvesonDESeq, normalized = TRUE))

write.csv(tuvesonDERes, "./Data/Tuveson iCAFmyCAF DE Results.csv")
write.csv(tuvesonNorm, "./Data/Tuveson iCAFmyCAF Normalised Counts.csv")