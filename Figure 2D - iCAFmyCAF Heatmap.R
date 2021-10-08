# Aims:
# - Plot a heatmap comparing the Tuveson iCAF/myCAF data with the Basal ECM-Targeted RNA-Seq panel data.
# - Add annotations to the heatmap highlighting similarities and differences between the experiments.

library("ComplexHeatmap")
library("circlize")



## Loading Data ----
#Differential expression results read into the environment from the previous analysis
basalDE <- read.csv("./Data/PSCs Basal DE Results.csv", stringsAsFactors = FALSE, row.names = 1)
tuvesonDE <- read.csv("./Data/Tuveson iCAFmyCAF DE Results.csv", stringsAsFactors = FALSE, row.names = 1)

#Reduces only to the genes that we are interested from the basal data
basalSignifDE <- basalDE[!is.na(basalDE$padj) & basalDE$padj < 0.05, ]
basalTuvesonDE <- tuvesonDE[which(rownames(tuvesonDE) %in% basalSignifDE$mgi_symbol), ]

#Normalised counts read into the environment from the previous analysis
basalNorm <- read.csv("./Data/PSCs Basal Normalised Count Matrix.csv", stringsAsFactors = FALSE, row.names = 1)
tuvesonNorm <- read.csv("./Data/Tuveson iCAFmyCAF Normalised Counts.csv", stringsAsFactors = FALSE, row.names = 1)

#Normalised counts also reduced only to the genes that we are interested in from the basal data
basalNorm <- basalNorm[rownames(basalSignifDE), ]
tuvesonNorm <- tuvesonNorm[basalSignifDE$mgi_symbol, c(1, 2, 5:8)]

#The basal count matrix is renamed to its mgi_symbol
row.names(basalNorm) <- basalSignifDE[match(row.names(basalSignifDE), row.names(basalNorm)), "mgi_symbol"]



## Scaled Heatmap Matrix ----
#A simple function to return the per-gene z-scores for a count matrix
calculateZScore <- function(reads){
  #Transposes the matrix as scale likes to work on columns
  reads <- t(reads)
  zScores <- reads #Copies the size of the matrix
  #Goes through each column and scales it separately
  zScores <- apply(reads, 2, scale)
  zScores <- t(zScores)
  colnames(zScores) <- row.names(reads)
  return(zScores)
}

basalZScores <- calculateZScore(basalNorm)
tuvesonZScores <- calculateZScore(tuvesonNorm)
combZScores <- merge(basalZScores, tuvesonZScores, by = "row.names")
combZScores <- data.frame(combZScores, row.names = 1)



## Heatmap Metadata ----
#How the heatmap is going to be split
heatmapSplit <- data.frame(
  row.names = colnames(combZScores),
  condition = factor(c(
    "KO_PSC", "KO_PSC", "KO_PSC",
    "Wt_PSC", "Wt_PSC", "Wt_PSC",
    "myCAF", "myCAF", "iCAF", "iCAF", "iCAF", "iCAF"),
    levels = c("Wt_PSC", "KO_PSC", "myCAF", "iCAF")
))

#A dataframe is created containing the annotation data
basalPlotData <- data.frame(
  exp = "KO vs Wt",
  gene = basalSignifDE$mgi_symbol,
  padj = basalSignifDE$padj,
  stat = basalSignifDE$stat
)

tuvesonPlotData <- data.frame(
  exp = "iCAFvsMyCAF",
  gene = rownames(basalTuvesonDE),
  padj = basalTuvesonDE$padj,
  stat = basalTuvesonDE$stat
)

combPlotData <- rbind(basalPlotData, tuvesonPlotData)

#The heatmap is ordered by combined stat, but a x2 weight is given to the PSC data since this is the data we are more interested in.
rowOrder <-  levels(reorder(combPlotData$gene, c(
  (combPlotData[combPlotData$exp == "KO vs Wt", "stat"])*2,
  combPlotData[combPlotData$exp == "iCAFvsMyCAF", "stat"])
))



## Heatmap Annotations ----
#Annotations will include bar plots of summary data and a concur yes/no indication
basalMatch <- match(rownames(combZScores), basalSignifDE$mgi_symbol)
tuvesonMatch <- match(rownames(combZScores), rownames(basalTuvesonDE))

concur <- ifelse(
  basalSignifDE[basalMatch, "log2FoldChange"] < 0 &
    basalTuvesonDE[tuvesonMatch, "log2FoldChange"] < 0,
  "True",
  ifelse(basalSignifDE[basalMatch, "log2FoldChange"] > 0 &
           basalTuvesonDE[tuvesonMatch, "log2FoldChange"] > 0,
         "True",
         "False")
)

rowAnnot = rowAnnotation(
  KOvsWt = anno_barplot(
    basalSignifDE[basalMatch, "log2FoldChange"],
    bar_width = 1,
    gp = gpar(col = "white", fill = ifelse(basalSignifDE[basalMatch, "log2FoldChange"] > 0, "#b2182b", "#2166ac")),
    border = TRUE,
    axis_param = list(at = seq(-4, 4, by = 2)), ylim = c(-4, 4), width = unit(3.5, "cm")
  ),
  iCAFvsmyCAF = anno_barplot(
    basalTuvesonDE[tuvesonMatch, "log2FoldChange"],
    bar_width = 1,
    gp = gpar(col = "white", fill = ifelse(
      basalTuvesonDE[tuvesonMatch, "padj"] < 0.05 &
        basalTuvesonDE[tuvesonMatch, "log2FoldChange"] > 0,
      "#b2182b",
      ifelse(basalTuvesonDE[tuvesonMatch, "padj"] >= 0.05 &
               basalTuvesonDE[tuvesonMatch, "log2FoldChange"] > 0,
             "#fddbc7",
             ifelse(basalTuvesonDE[tuvesonMatch, "padj"] < 0.05 &
                      basalTuvesonDE[tuvesonMatch, "log2FoldChange"] < 0,
                    "#2166ac",
                    "#d1e5f0")
      )
    )),
    border = TRUE,
    axis_param = list(at = seq(-8, 8, by = 2)), ylim = c(-8, 8), width = unit(3.5, "cm")
  ),
  Concur = concur,
  col = list(Concur = c("True" = "#55C54B", "False" = "#C54B4B")), show_legend = FALSE
)



## Custom Legends ----
legends <- list(
  Legend(col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
         title = "Per Gene\nz Scores",
         grid_height  = unit(0.7, "cm"),
         grid_width = unit(0.7, "cm"),
         labels_gp = gpar(fontsize = 11),
         title_gp = gpar(fontsize = 11)
  ),
  Legend(labels  = c("p<0.05 & +l2fc", "p>0.05 & +l2fc", "p>0.05 & -l2fc", "p<0.05 & -l2fc"),
         title = "\nDifferential\nExpression",
         legend_gp = gpar(fill = c("#b2182b" ,"#fddbc7", "#d1e5f0" , "#2166ac")),
         grid_height  = unit(0.7, "cm"),
         grid_width = unit(0.7, "cm"),
         labels_gp = gpar(fontsize = 11),
         title_gp = gpar(fontsize = 11)
  ),
  Legend(labels  = c("True", "False"),
         title = "\nConcur",
         legend_gp = gpar(fill = c("#55C54B" ,"#C54B4B")),
         grid_height  = unit(0.7, "cm"),
         grid_width = unit(0.7, "cm"),
         labels_gp = gpar(fontsize = 11),
         title_gp = gpar(fontsize = 11)
  )
)



## Plotting the Heatmap ----
baseHeatmap <- Heatmap(
  name = "per gene z scores",
  as.matrix(combZScores),

  row_order = rowOrder,
  row_names_side = "left",
  cluster_rows = FALSE,

  column_split = heatmapSplit,
  right_annotation = rowAnnot,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_column_names = FALSE,
  cluster_column_slices = FALSE,

  show_heatmap_legend = FALSE
)
heatmap <- draw(baseHeatmap, annotation_legend_list = legends)
heatmap
