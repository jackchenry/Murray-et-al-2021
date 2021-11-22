# Aims:
# - Plot individual GSEA plots for some selected pathways

#Selected Hallmark pathways are:
# - EMT
# - Inflammatory response
# - IL6/JAK/STAT3 signalling



##Preparation ----
#The required packages are loaded
library("fgsea") #Used for Gene-set enrichment analysis
library("ggplot2") #Used to plot the ranked plot and modify enrichment plots
library("dplyr") #Used for case_when()
library("grid") #Used for grobs

#The genesets are loaded into the environment. Genesets obtained from the molecular signatures database.
hallmark <- gmtPathways("./Data/Genesets/h.all.v7.2.symbols.gmt")

#The RNAseq data from this study is read into the environment (FigS5).
tumDE <- read.csv("./Data/InVivo Tumour DE Results.csv", stringsAsFactors = FALSE, row.names = 1)



##Construction of the final gene-set lists ----
#hgnc symbol rank list
hNameGenesets <- list(
  "hallmarkEMT" = hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  "hallmarkInflammatory" = hallmark$HALLMARK_INFLAMMATORY_RESPONSE,
  "hallmarkIl6Stat3" = hallmark$HALLMARK_IL6_JAK_STAT3_SIGNALING
)



##Construction of the different ranks ----
#All NA genes are removed
ranks <- tumDE[!is.na(tumDE$stat), ]
#Ranks are re-orded by stat
ranks <- ranks[order(-ranks$stat), ]

#hgnc symbol rank list
hNameRanks <- ranks$stat
names(hNameRanks) <- ranks$hgnc_symbol



##Fgsea analysis ----
hNameRes <- fgsea(hNameGenesets, hNameRanks, eps = 0)



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
plotEnrichmentGraphs(hNameRes, hNameGenesets, hNameRanks, "tumour data")



##Plotting the ranked plot ----
ribbonPlotData <- data.frame(
  Rank = 1:length(hNameRanks),
  Stat = hNameRanks
)

tumRankedPlot <- ggplot(ribbonPlotData, aes(x = Rank, y = Stat))+
  geom_area() +
  geom_vline(xintercept = length(hNameRanks[hNameRanks > 0]), linetype = "dashed") +
  labs(x = "Tumour Rank", y = "Stat") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(
    plot.margin=margin(0.7,0.4,0.6,0.2,"cm"),
  )
tumRankedPlot
