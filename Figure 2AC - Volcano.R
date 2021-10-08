# Aims:
# - Plot a volcano plot for the ECM-Targeted RNA-Seq panel data.
# - Volcano plot should show Wt vs KO.
# - Separate volcano plots for the basal and TGF-stimulated conditions.
# - Volcano plot should highlight genes that are p<0.05 and l2fc>1

#Volcano plots will be plotted using Enhanced Volcano
library("EnhancedVolcano")



## Basal Volcano ----
basalShrunkResults <- read.csv("./Data/PSCs Basal Shrunk DE Results.csv")
basalVolcano <- EnhancedVolcano(
  basalShrunkResults,
  lab = basalShrunkResults$mgi_symbol,
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05, FCcutoff = 1.0,
  pointSize = 3, labSize = 4, colAlpha = 1, ylim = c(0, 30), xlim = c(-2.5, 2.5),
  col=c("grey", "lightblue", "lightblue", "red3"), legendPosition = "none", gridlines.minor = FALSE,
  drawConnectors = TRUE, widthConnectors = 0.5, labFace = "italic",
  title = "Basal", subtitle = NULL, caption = NULL
)
basalVolcano



## TGF Volcano ----
tgfShrunkResults <- read.csv("./Data/PSCs TGF Shrunk DE Results.csv")
tgfVolcano <- EnhancedVolcano(
  tgfShrunkResults,
  lab = tgfShrunkResults$mgi_symbol,
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05, FCcutoff = 1.0,
  pointSize = 3, labSize = 4, colAlpha = 1, ylim = c(0, 32), xlim = c(-2.8, 5),
  col=c("grey", "lightblue", "lightblue", "red3"), legendPosition = "none", gridlines.minor = FALSE,
  drawConnectors = TRUE, widthConnectors = 0.5, labFace = "italic",
  title = "TGF stimulated", subtitle = NULL, caption = NULL
)
tgfVolcano
