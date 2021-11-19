# Aims:
# - Plot individual GSEA plots for some selected pathways

#Selected custom pathways are:
# - myCAF up - top 200 genes from the Ohlund et al., 2017 paper
# - iCAF down - top 200 genes from the Ohlund et al., 2017 paper



##Preparation ----
#The required packages are loaded
library("fgsea") #Used for Gene-set enrichment analysis
library("ggplot2") #Used to plot the ranked plot and modify enrichment plots
library("dplyr") #Used for case_when()
library("grid") #Used for grobs

#The genesets are loaded into the environment. Genesets obtained from the
hallmark <- gmtPathways("./Data/Genesets/h.all.v7.2.symbols.gmt")

#The Ohlund paper DE data is loaded into the environment. Obtained from the GEO: GSE93313.
myCAFDE <- read.csv("./Data/GSE93313_DESeq_2Dvs3D.csv", stringsAsFactors = FALSE) #2dvs3D = myCAF
iCAFDE <- read.csv("./Data/GSE93313_DESeq_3DvsTranswell.csv", stringsAsFactors = FALSE) #3Dvs4D = iCAF

#The RNAseq data from this study is read into the environment.
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



##Construction of the final gene-set lists ----
#mgi symbol id rank
#The top 200 and bottom iCAF and myCAF differentially expressed genes are selected to form the gene-sets.
#Some of these may not be expressed in our data, hence why the number of genes shown in the graph is different.
mNameGenesets <- list(
  "myCAF_Up" = myCAFDESignf$id[1:200],
  "iCAF_Up" = iCAFDESignf$id[1:200]
)



##Construction of the different ranks ----
#All NA genes are removed
ranks <- tumDE[!is.na(tumDE$stat), ]
#Ranks are re-ordered by stat
ranks <- ranks[order(-ranks$stat), ]

#mgi symbol id rank
mNameRanks <- ranks$stat
names(mNameRanks) <- ranks$mgi_symbol



##Fgsea analysis ----
mNameRes <- fgsea(mNameGenesets, mNameRanks, eps = 0)



##Enrichment plot function ----
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



##Enrichment plot results ----
plotEnrichmentGraphs(mNameRes, mNameGenesets, mNameRanks, "tumour data")



##Plotting the ranked plot ----
ribbonPlotData <- data.frame(
  Rank = 1:length(mNameRanks),
  Stat = mNameRanks
)

tumRankedPlot <- ggplot(ribbonPlotData, aes(x = Rank, y = Stat))+
  geom_area() +
  geom_vline(xintercept = length(mNameRanks[mNameRanks > 0]), linetype = "dashed") +
  labs(x = "Tumour Rank", y = "Stat") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(
    plot.margin=margin(0.7,0.4,0.6,0.2,"cm"),
  )
tumRankedPlot
