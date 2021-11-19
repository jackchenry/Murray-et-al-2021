# Aims:
# - Perform GSEA analysis on Hallmark pathways on tumours
# - Perform GSEA analysis on Hallmark pathways on TCGA-PAAD PKN2KO Signature split data
# - Plot a graph comparing the two



## Preparation ----
#Required packages are loaded
library("DESeq2") #To perform differential expression on the TCGA-PAAD patient splits
library("biomaRt") #To get gene labels for the TCGA-PAAD data
library("fgsea") #To perform the GSEA analysis
library("dplyr") #Used for case_when()
library("ggplot2") #Used to plot the final graph

#Required data loaded into the environment
#Non-normalised counts are required for DESeq2 differential expression
rawPAAD <- read.csv(
  "./../common/RawCountMatrices/TCGA-PAAD-Raw-Primary_Tumor.csv",
  stringsAsFactors = FALSE, row.names = 1
)
#These are the genes that make up the PKN2KO score performed previously in Fig6D
PAAD_PKN2KOScore <- read.csv("./Data/TCGA-PAAD PKN2KO Signature Data.csv", stringsAsFactors = FALSE)

#This is the PKN2KO tumour differential expression performed previously in FigS5
tumDE <- read.csv("./Data/InVivo Tumour DE Results.csv", stringsAsFactors = FALSE, row.names = 1)

#This is the Hallmark geneset file downloaded from the MSigDB
hallmarkGeneset <- gmtPathways("./Data/Genesets/h.all.v7.2.symbols.gmt")

#These are a set of manually curated nicely formatted Hallmark pathway names
hallmarkNames <- read.csv("./Data/HallmarkNames.csv", stringsAsFactors = FALSE, header = FALSE)$V1



##TCGA-PAAD PKN2KO Score Differential Expression ----
#Col data for DESeq2 is set up using PKN2KO score splits
colData <- PAAD_PKN2KOScore[, c("patient", "split")]
rownames(colData) <- colData$patient
colData$patient <- NULL
#Required to set up the correct comparision in DESeq2
colData$split <- factor(colData$split, levels = c("Low", "High"))

#Differential expression is formed
DESeqData <- DESeqDataSetFromMatrix(rawPAAD, colData = colData, design = ~ split)
DESeqData <- DESeq(DESeqData)
PAADDEResults <- results(DESeqData)
head(PAADDEResults)



##Collecting gene information for TCGA-PAAD ----
#Genesets are symbols therefore symbols need to be downloaded using biomaRt
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
geneIds <- rownames(rawPAAD)
geneInfo <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id", values = geneIds, mart = human, uniqueRows = TRUE)

#Removes duplicated ensembl Ids and sets up matrix with row names
geneInfo <- geneInfo[!duplicated(geneInfo$ensembl_gene_id), ]
rownames(geneInfo) <- geneInfo$ensembl_gene_id
geneInfo$ensembl_gene_id <- NULL

#Merging DE results with gene information
PAADDE <- merge(data.frame(PAADDEResults), geneInfo, by = "row.names", all.y=TRUE)
rownames(PAADDE) <- PAADDE$Row.names
PAADDE$Row.names <- NULL
PAADDE <- PAADDE[order(PAADDE$stat), ]

#Saving PAAD differential expression results for later use
write.csv(PAADDE, "./Data/TCGA-PAAD PKN2KO Signature Differential Expression Results.csv")



##TCGA-PAAD PKN2KO Score Rank ----
PAADRank <- PAADDE$stat
names(PAADRank) <- PAADDE$hgnc_symbol
PAADRank <- PAADRank[order(-PAADRank)]

print(paste0("Number of genes removed is: ", length(PAADRank[is.na(PAADRank)]), " out of ", length(PAADRank)))
PAADRank <- PAADRank[!is.na(PAADRank)]



##Tumours Rank ----
tumRank <- tumDE$stat
names(tumRank) <- tumDE$hgnc_symbol
tumRank <- tumRank[order(-tumRank)]

print(paste0("Number of genes removed is: ", length(tumRank[is.na(tumRank)]), " out of ", length(tumRank)))
tumRank <- tumRank[!is.na(tumRank)]



##GSEA Analysis ----
#GSEA analysis is performed using fgsea
#An eps of 0 is not required as we are just interested in if p<0.05
hallmarkResultsPAAD <- fgsea(hallmarkGeneset, PAADRank)
hallmarkResultsTum <- fgsea(hallmarkGeneset, tumRank)



##Data preparation for GSEA plot----
##PAAD
#The important information is extracted
GSEAPAAD <- hallmarkResultsPAAD[ , c("pathway", "NES", "padj")]

#The GSEA results are labelled with nicely formatted names which will be used in the graph
GSEAPAAD$pathwayNames <- hallmarkNames

#Case when is used for significant and direction labels
GSEAPAAD$p <- factor(case_when(
    hallmarkResultsPAAD$padj < 0.05 & hallmarkResultsPAAD$NES > 0 ~ "p<0.05 & +NES",
    hallmarkResultsPAAD$padj < 0.05 & hallmarkResultsPAAD$NES < 0 ~ "p<0.05 & -NES",
    hallmarkResultsPAAD$padj >= 0.05 & hallmarkResultsPAAD$NES > 0 ~ "p>0.05 & +NES",
    hallmarkResultsPAAD$padj >= 0.05 & hallmarkResultsPAAD$NES < 0 ~ "p>0.05 & -NES"
  ), levels = c("p<0.05 & +NES", "p>0.05 & +NES", "p>0.05 & -NES", "p<0.05 & -NES")
)
#This will be used for faceting
GSEAPAAD$study <- "TCGA-PAAD"

##Tum
#The important information is extracted
GSEATum <- hallmarkResultsTum[ , c("pathway", "NES", "padj")]

#The GSEA results are labelled with nicely formatted names which will be used in the graph
GSEATum$pathwayNames <- hallmarkNames

#Case when is used for significant and direction labels
GSEATum$p <- factor(case_when(
    hallmarkResultsTum$padj < 0.05 & hallmarkResultsTum$NES > 0 ~ "p<0.05 & +NES",
    hallmarkResultsTum$padj < 0.05 & hallmarkResultsTum$NES < 0 ~ "p<0.05 & -NES",
    hallmarkResultsTum$padj >= 0.05 & hallmarkResultsTum$NES > 0 ~ "p>0.05 & +NES",
    hallmarkResultsTum$padj >= 0.05 & hallmarkResultsTum$NES < 0 ~ "p>0.05 & -NES"
  ), levels = c("p<0.05 & +NES", "p>0.05 & +NES", "p>0.05 & -NES", "p<0.05 & -NES")
)

#This will be used for faceting
GSEATum$study <- "Tumours"

##Merging
GSEAComb <- rbind(GSEATum, GSEAPAAD)



##Final Plot ----
#The pathways are rearranged based on their combined NES
#The different studies are faceted to allow comparisons between the two
hallmarkCombPlot <- ggplot(GSEAComb, aes(reorder(pathwayNames, NES), NES)) +
  geom_col(aes(fill = p)) +
  scale_fill_manual(values=c("#b2182b" ,"#fddbc7", "#d1e5f0" , "#2166ac")) +
  facet_wrap(~ study) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA") +
  theme_bw()
hallmarkCombPlot
