# Aims:
# - Plot individual GSEA plots for some selected pathways

#Selected Hallmark pathways are:
# - EMT
# - Inflammatory response
# - IL6/JAK/STAT3 signalling

#Selected custom pathways are:
# - myCAF up and down - top and bottom 200 genes from the Ohlund et al., 2017 paper
# - iCAF up and down - top and bottom 200 genes from the Ohlund et al., 2017 paper
# - CD105+ vs CD105- - top and bottom 200 genes from the Hutton et al., 2021 paper



##Preparation ----
#The required packages are loaded
library("fgsea") #Used for Gene-set enrichment analysis
library("DESeq2") #Used to perform differential expression on TCGA-PAAD grade data to create genesets.
library("biomaRt")

#The genesets are loaded into the environment. Genesets obtained from the
hallmark <- gmtPathways("./Data/Genesets/h.all.v7.2.symbols.gmt")

#The Ohlund paper DE data is loaded into the environment. Obtained from the GEO: GSE93313.
myCAFDE <- read.csv("./Data/GSE93313_DESeq_2Dvs3D.csv", stringsAsFactors = FALSE) #2dvs3D = myCAF
iCAFDE <- read.csv("./Data/GSE93313_DESeq_3DvsTranswell.csv", stringsAsFactors = FALSE) #3Dvs4D = iCAF

#The Hutton paper DE data is loaded into the environment. Obtained from the GEO: GSE176057.
cd105CAF <- read.csv("./Data/GSE176057_Raw_data_and_in_vitro_PaF_CD105pos_vs_CD105neg_comparison_1_.csv", stringsAsFactors = FALSE)

#The human TCGA-PAAD clinical data is read into the environment. Data was obtained from Liuetal et al., 2016.
paadClinic <- read.csv("./Data/TCGA-PAAD Clinical Data.csv", stringsAsFactors = FALSE)
#
paadRaw <- read.csv("./Data/TCGA-PAAD Raw Primary Tumour Counts.csv", stringsAsFactors = FALSE, row.names = 1)

#
tumDE <- read.csv("./Data/InVivo Tumour DE Results.csv", stringsAsFactors = FALSE, row.names = 1)


##Construction of Ohlund Genesets ----
#Only the significant results are kept
myCAFDESignf <- myCAFDE[myCAFDE$padj < 0.05 &  !is.na(myCAFDE$padj), ]
iCAFDESignf <- iCAFDE[iCAFDE$padj < 0.05  & !is.na(iCAFDE$padj), ]
#Remove any genes that have a fold change of 0 or inf
myCAFDESignf <- myCAFDESignf[myCAFDESignf$foldChange != 0 & is.finite(myCAFDESignf$foldChange), ]
iCAFDESignf <- iCAFDESignf[iCAFDESignf$foldChange != 0 & is.finite(iCAFDESignf$foldChange), ]

#Ordered by fc same as in the paper
myCAFDESignf <- myCAFDESignf[order(myCAFDESignf$foldChange), ] #This one is backwards with its DESeq comparison
iCAFDESignf <- iCAFDESignf[order(-iCAFDESignf$foldChange), ]

#Check that the top genes match those that are in the paper
myCAFDESignf$id[1:25]
iCAFDESignf$id[1:25]



##Construction of Hutton Genesets ----
#Genes are already filtered (significance and errors) and ordered
#Renamed so that only the symbol is kept. Important for creating ranks
cd105CAF$GeneID <- substr(cd105CAF$GeneID, 20, 100)



##Construction of the PAAD grade genesets ----
summary(factor(paadClinic$histological_grade)) #Going to be comparing 129 G1/G2 patients with the 53 G3/G4
levels(factor(paadClinic$histological_grade)) #No complicated grades in the data

paadGrade <- data.frame(
  "row.names" = paadClinic$bcr_patient_barcode,
  "grade" = factor(case_when(
    paadClinic$histological_grade == "G1" ~ "G1G2",
    paadClinic$histological_grade == "G2" ~ "G1G2",
    paadClinic$histological_grade == "G3" ~ "G3G4",
    paadClinic$histological_grade == "G4" ~ "G3G4",
   ))
)
#Removes patients with an NA grade
paadGrade <- paadGrade[!is.na(paadGrade$grade), , drop = FALSE]

#Removes weird formatting for the patient names
colnames(paadRaw) <- gsub(".", "-", substr(colnames(paadRaw), 1, 12), fixed = TRUE)
#Removes clinic data for patients that we do not have sequencing data for
paadGrade <- paadGrade[row.names(paadGrade) %in% colnames(paadRaw), , drop = FALSE]
#Re-orders and removes patients that we do not have clinical data for
paadRaw <- paadRaw[, match(row.names(paadGrade), colnames(paadRaw))]

paadDESeq <- DESeqDataSetFromMatrix(countData = paadRaw, colData = paadGrade, design = ~ grade)
paadDESeq <- DESeq(paadDESeq)
paadDERes <- results(paadDESeq)
head(paadDERes) #G3/G4 vs G1/G2
paadDERes <- data.frame(paadDERes)
#Only significant genes kept
paadDEResSignif <- paadDERes[!is.na(paadDERes$padj) & paadDERes$padj < 0.05, ]
#Data is reordered by stat
paadDEResSignif <- paadDEResSignif[order(-paadDEResSignif$stat), ]



##Construction of the final gene-set lists
#Different lists are made for the different labellings of the ranks.
#hgnc symbol rank list
hNameGenesets <- list(
  "hallmarkEMT" = hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  "hallmarkInflammatory" = hallmark$HALLMARK_INFLAMMATORY_RESPONSE,
  "hallmarkIl6Stat3" = hallmark$HALLMARK_IL6_JAK_STAT3_SIGNALING
)
#Human ensembl id rank list
#The top and bottom 200 differentially expressed genes are selected to form the grade gene-sets.
hIdGenesets <- list(
  "highPAADGrade_Up" = row.names(paadDEResSignif[1:200, ]),
  "highPAADGrade_Down" = row.names(paadDEResSignif[nrow(paadDEResSignif):(nrow(paadDEResSignif) - 199), ])
)

#mgi symbol id rank
#The top 200 and bottom iCAF and myCAF differentially expressed genes are selected to form the gene-sets.
#Some of these may not be expressed in our data, hence why the number of genes shown in the graph is different.
mNameGenesets <- list(
  "myCAF_Up" = myCAFDESignf$id[1:200],
  "myCAF_Down" = rev(myCAFDESignf$id)[1:200],
  "iCAF_Up" = iCAFDESignf$id[1:200],
  "iCAF_Down" = rev(iCAFDESignf$id)[1:200],
  "cd105CAF_Up" = cd105CAF$GeneID[1:200],
  "cd105CAF_Down" = rev(cd105CAF$GeneID)[1:200]
)




##Construction of the different ranks
#All NA genes are removed
ranks <- tumDE[!is.na(tumDE$stat), ]
#Ranks are re-orded by stat
ranks <- ranks[order(-ranks$stat), ]

#hgnc symbol rank list
hNameRanks <- ranks$stat
names(hNameRanks) <- ranks$hgnc_symbol

#Human ensembl id rank list
hIdRanks <- ranks$stat
names(hIdRanks) <- ranks$human_ensembl_gene_id

#mgi symbol id rank
mNameRanks <- ranks$stat
names(mNameRanks) <- ranks$mgi_symbol






##Fgsea analysis
hNameRes <- fgsea(hNameGenesets, hNameRanks, eps = 0)
hIdRes <- fgsea(hIdGenesets, hIdRanks, eps = 0)
mNameRes <- fgsea(mNameGenesets, mNameRanks, eps = 0)




##Enrichment plot function
#This function takes in the results from geneset enrichment analysis,
#the genesets and ranks used in the results and returns a list of enrichment plots for each geneset in the list.
plotEnrichmentGraphs <- function(fgseaResults, genesets, ranks, titleDataType){
  enrichmentPlots <- list()
  #for loop to go through each geneset in the geneset list
  for(currentGeneset in names(genesets)){
    p <- fgseaResults[fgseaResults$pathway == currentGeneset, "padj"] #adjusted p value from the results
    signf <- case_when(       #whether this p value is significant
      p < 0.0001 ~ "****",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p >= 0.05 ~ "ns"
    )
    if(fgseaResults[fgseaResults$pathway == currentGeneset, "NES"] >= 0){ #this adjust the position of the text based on whether it is a negative or positive enrichment
      y <- 0.94
      vjust <- 1
      x = 0.95
      hjust = 1
    }else if(fgseaResults[fgseaResults$pathway == currentGeneset, "NES"] < 0){
      y <- 0.06
      vjust <- 0
      x = 0.05
      hjust = 0
    }
    grob <- grobTree(textGrob(paste("nGenes: ", fgseaResults[fgseaResults$pathway == currentGeneset, "size"], #builds up the text to display over the graph
                                    "\n", "NES: ", signif(fgseaResults[fgseaResults$pathway == currentGeneset, "NES"], 3),
                                    "\n", "padj: ", signif(p, 3),
                                    "\n", "Signif: ", signf, sep=""),
                              x = x,  y = y, hjust = hjust, vjust = vjust))
    #plots the enrichment graph
    enrichmentPlot <- plotEnrichment(genesets[[currentGeneset]], ranks, ticksSize = 0.22) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x ="Rank",
           y = "Enrichment Score",
           title = paste("Enrichment Plot: ", currentGeneset, " in ", titleDataType, sep="")) +
      theme_bw() +
      geom_line(size = 1.5, color = "#228C22") +
      annotation_custom(grob)

    #saves a list of graphs
    enrichmentPlots[[currentGeneset]] <- enrichmentPlot
  }
  return(enrichmentPlots)
}




##Enrichment plot results
plotEnrichmentGraphs(hNameRes, hNameGenesets, hNameRanks, "tumour data")
plotEnrichmentGraphs(hIdRes, hIdGenesets, hIdRanks, "tumour data")
plotEnrichmentGraphs(mNameRes, mNameGenesets, mNameRanks, "tumour data")





