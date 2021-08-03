# Aims:
# - Plot a heatmap of the ECM-Targeted RNA-Seq panel data.

# Notes:
# - Plot an average of the samples for TGF-treated and untreated (basal) samples.
# - Heatmap should be ordered by log-fold-change
# - If log fold change is greater than 1 then the name should be in bold

#Heatmap ploted in ggplot2
library("ggplot2")



## Basal Heatmap ----
#Data is read into the environment
basalResults <- read.csv("./Data/ECM Targeted RNAseq Basal DE Results.csv")

#Only significant results are kept to plot
basalPlotData <- basalResults[!is.na(basalResults$padj) & basalResults$padj < 0.05,
                              c("mgi_symbol", "log2FoldChange", "stat", "padj")]

basalPlot <- ggplot(basalPlotData, aes(reorder(mgi_symbol, log2FoldChange), y = 1, fill = log2FoldChange)) +
  geom_tile(colour = "white") +
  coord_flip() +
  scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-5, 5), breaks=seq(-5,5,by=1)) +
  labs(x="ECM Genes",
       y="",
       title = "\nBasal",
       fill="Log2 Fold\nChange") +
  guides(fill = guide_colorbar(barheight = 16)) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(
    size = rel(ifelse(abs(basalPlotData$log2FoldChange[order(basalPlotData$log2FoldChange)]) > 1, 1.2, 1)),
    face = ifelse(abs(basalPlotData$log2FoldChange[order(basalPlotData$log2FoldChange)]) > 1, "bold", "plain"))
  )
basalPlot



## TGF-stimulated Heatmap ----
#Data is read into the environment
tgfResults <- read.csv("./Data/ECM Targeted RNAseq TGF DE Results.csv")

#Only significant results are kept to plot
tgfPlotData <- tgfResults[!is.na(tgfResults$padj) & tgfResults$padj < 0.05,
                          c("mgi_symbol", "log2FoldChange", "stat", "padj")]

tgfPlot <- ggplot(tgfPlotData, aes(reorder(mgi_symbol, log2FoldChange), y = 1, fill = log2FoldChange)) +
  geom_tile(colour = "white") +
  coord_flip() +
  scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-5, 5), breaks=seq(-5,5,by=1)) +
  labs(x="ECM Genes",
       y="",
       title = "\nTGF",
       fill="Log2 Fold\nChange") +
  guides(fill = guide_colorbar(barheight = 16)) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(
    size = rel(ifelse(abs(tgfPlotData$log2FoldChange[order(tgfPlotData$log2FoldChange)]) > 1, 1.2, 1)),
    face = ifelse(abs(tgfPlotData$log2FoldChange[order(tgfPlotData$log2FoldChange)]) > 1, "bold", "plain"))
  )
tgfPlot
