# Aims:
# - Perform GSEA analysis on Hallmark pathways on tumours
# - Perform GSEA analysis on Hallmark pathways on TCGA-PAAD PKN2KO Signature split data
# - Plot a graph comparing the two



## Preparation ----
#Required packages are loaded
library("DESeq2")
library("biomaRt")
library("fgsea")

#Required data loaded into the environment
PAADRaw <-
PAAD_PKN2KOScore <- read.csv("./Data/TCGA-PAAD PKN2KO Signature Data.csv", stringsAsFactors = FALSE)

tumDE <- read.csv("./Data/InVivo Tumour DE Results.csv", stringsAsFactors = FALSE, row.names = 1)



##TCGA-PAAD PKN2KO Score Differential Expression




##TCGA-PAAD PKN2KO Score Rank
PAADRank <- PAADDE$stat
names(PAADRank) <- PAADDE$Gene.stable.ID.1
PAADRank <- PAADRank[order(-PAADRank)]

print(paste0("Number of genes removed is: ", length(PAADRank[is.na(PAADRank)]), " out of ", length(PAADRank)))
PAADRank <- PAADRank[!is.na(PAADRank)]











##Tumours Rank
tumRank <- tumDE$stat
names(tumRank) <- tumDE$Gene.stable.ID.1
tumRank <- tumRank[order(-tumRank)]

print(paste0("Number of genes removed is: ", length(tumRank[is.na(tumRank)]), " out of ", length(tumRank)))
tumRank <- tumRank[!is.na(tumRank)]







#Saving useful data
#Tumour DE with the hcgn symbols is saved
write.csv(tumDE, "./Data/InVivo Tumour DE Results.csv", row.names = FALSE)


