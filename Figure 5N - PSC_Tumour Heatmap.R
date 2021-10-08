# Aims:
# - Plot a heatmap comparing the PSC Qiaseq data with the Pancreatic tumour data.
# - Add annotations to the heatmap highlighting similarities and differences between the experiments.



## Preparation ----
#Required packages are loaded
library("ComplexHeatmap")
library("circlize")

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

#Data is read into the environment
pscNorm <- read.csv("./Data/PSCs Basal Normalised Count Matrix.csv", stringsAsFactors = FALSE, row.names = 1)
tumNorm <- read.csv("./Data/InVivo Tumour Normalised Count Matrix.csv", stringsAsFactors = FALSE, row.names = 1)

pscDE <- read.csv("./Data/PSCs Basal DE Results.csv", stringsAsFactors = FALSE, row.names = 1)
tumDE <- read.csv("./Data/InVivo Tumour DE Results.csv", stringsAsFactors = FALSE, row.names = 1)

pscDE_signif <- pscDE[pscDE$padj < 0.05 & !is.na(pscDE$padj), ]
tumDE_signif <- tumDE[tumDE$padj < 0.05 & !is.na(tumDE$padj), ]

#This is the PKN2-KO signature from the paper
PKN2KOSig <- read.csv("./Data/PKN2KO Signature.csv", stringsAsFactors = FALSE)$mgi_symbol



##Scaled Heatmap Matrix ----
#Important genes are extracted
pscNorm_pscSignif <- pscNorm[rownames(pscDE_signif), ]
tumNorm_pscSignif <- tumNorm[rownames(pscDE_signif), ]

#Z-Scores are calculated using the above function
pscZ <- calculateZScore(pscNorm_pscSignif)
tumZ <- calculateZScore(tumNorm_pscSignif)

#PSC and tumour data are merged together
allZ <- merge(tumZ, pscZ, by = "row.names", all.y=TRUE)
rownames(allZ) <- allZ$Row.names
allZ$Row.names <- NULL



##Annotation data ----
pscData <- data.frame(
  "exp" = "PSCs",
  "gene" = rownames(pscDE_signif),
  "l2fc" = pscDE_signif$log2FoldChange,
  "lfcSE" = pscDE_signif$lfcSE,
  "padj" = pscDE_signif$padj,
  "stat" = pscDE_signif$stat)

tumData <- data.frame(
  "exp" = "Tumour",
  "gene" = rownames(tumDE[rownames(allZ), ]),
  "l2fc" = tumDE[rownames(allZ), "log2FoldChange"],
  "lfcSE" = tumDE[rownames(allZ), "lfcSE"],
  "padj" = tumDE[rownames(allZ), "padj"],
  "stat" = tumDE[rownames(allZ), "stat"])

allData <- rbind(pscData, tumData)



## Heatmap Metadata ----
#The heatmap is split into PSC, Tumour, Wt and KO groups to aid with visualisation of general trends
heatmapSplit <- data.frame(
  row.names = colnames(allZ),
  genotpye = c("WT", "WT", "KO", "KO", "KO", "KO", "KO", "WT", "WT", "KO", "WT", "KO", "KO", "KO", "WT", "WT", "WT"),
  type = c(rep("Tumour", 11), rep("Cells", 6))
  )

#Labels are the MGI symbols that are present in the DE data
labels <- pscDE[rownames(allZ), "mgi_symbol"]

#Rows are ordered based on their stat. A x2 weight is given to the Tumours since this is the one we are comparing to the PSC data
rowOrder <- levels(reorder(
  allData$gene,
  c(allData[allData$exp == "PSCs", "stat"], allData[allData$exp == "Tumour", "stat"]*2)
  ))

#The labels are in bold if they belong to the PKN2KO signature
labelFace <- ifelse(labels %in% PKN2KOSig, "bold", "plain")



##Heatmap Annotations ----
#The annotation that shows whether the PSC and tumour data agree
tumConcur <- ifelse(
  pscDE_signif[match(rownames(allZ), rownames(pscDE_signif)), "log2FoldChange"] < 0 &
    tumDE[match(rownames(allZ), rownames(tumDE)), "log2FoldChange"] < 0, "True",
  ifelse(pscDE_signif[match(rownames(allZ), rownames(pscDE_signif)), "log2FoldChange"] > 0 &
           tumDE[match(rownames(allZ), rownames(tumDE)), "log2FoldChange"] > 0, "True",
         "False"))

#The annotation that shows the Log2-fold change obtained from the DE results data
rowAnnot = rowAnnotation(
  KOvsWt_PSCs = anno_barplot(
    pscDE_signif[match(rownames(allZ), rownames(pscDE_signif)), "log2FoldChange"],
    bar_width = 1,
    gp = gpar(col = "white", fill = ifelse(pscDE_signif[match(rownames(allZ), rownames(pscDE_signif)), "log2FoldChange"] > 0, "#b2182b", "#2166ac")),
    border = TRUE,
    axis_param = list(at = seq(-4, 4, by = 2)), ylim = c(-4, 4), width = unit(3.5, "cm")
  ),
  KOvsWt_Tumour = anno_barplot(
    tumDE[match(rownames(allZ), rownames(tumDE)), "log2FoldChange"],
    bar_width = 1,
    gp = gpar(col = "white", fill = ifelse((tumDE[match(rownames(allZ), rownames(tumDE)), "padj"] < 0.05 & !is.na(tumDE[match(rownames(allZ), rownames(tumDE)), "padj"]))
                                           & tumDE[match(rownames(allZ), rownames(tumDE)), "log2FoldChange"] > 0, "#b2182b",
                                           ifelse((tumDE[match(rownames(allZ), rownames(tumDE)), "padj"] >= 0.05 | is.na(tumDE[match(rownames(allZ), rownames(tumDE)), "padj"]))
                                                  & tumDE[match(rownames(allZ), rownames(tumDE)), "log2FoldChange"] > 0, "#fddbc7",
                                                  ifelse(is.na(tumDE[match(rownames(allZ), rownames(tumDE)), "padj"]) & tumDE[match(rownames(allZ), rownames(tumDE)), "log2FoldChange"] > 0, "#fddbc7",
                                                         ifelse((tumDE[match(rownames(allZ), rownames(tumDE)), "padj"] < 0.05 & !is.na(tumDE[match(rownames(allZ), rownames(tumDE)), "padj"]))
                                                                & tumDE[match(rownames(allZ), rownames(tumDE)), "log2FoldChange"] < 0, "#2166ac", "#d1e5f0"))))),
    border = TRUE,
    axis_param = list(at = seq(-2, 2, by = 2)), ylim = c(-2.5, 2.5), width = unit(3.5, "cm")
  ),
  Concur = tumConcur,
  col = list(Concur = c("True" = "#2166ac", "False" = "#b2182b")),
  show_legend = FALSE,
  annotation_label = c("KOvsWt PSCs", "KOvsWt Tumour", "Concur"),
)



##Custom Legends ----
legends <- list(
  Legend(col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), title = "Per Gene\nz Scores",
         grid_height  = unit(0.7, "cm"), grid_width = unit(0.7, "cm"), labels_gp = gpar(fontsize = 11), title_gp = gpar(fontsize = 11)),
  Legend(labels  = c("p<0.05 & +l2fc", "p>0.05 & +l2fc", "p>0.05 & -l2fc", "p<0.05 & -l2fc"), title = "\nDifferential\nExpression", legend_gp = gpar(fill = c("#b2182b" ,"#fddbc7", "#d1e5f0" , "#2166ac")),
         grid_height  = unit(0.7, "cm"), grid_width = unit(0.7, "cm"), labels_gp = gpar(fontsize = 11), title_gp = gpar(fontsize = 11)),
  Legend(labels  = c("True", "False"), title = "\nConcur", legend_gp = gpar(fill = c("#2166ac" ,"#b2182b")),
         grid_height  = unit(0.7, "cm"), grid_width = unit(0.7, "cm"), labels_gp = gpar(fontsize = 11), title_gp = gpar(fontsize = 11))
)



##Final Heatmap Plot ----
baseHeatmap <- Heatmap(
  name = "per gene z scores",
  as.matrix(allZ),

  row_order = rowOrder,
  row_names_side = "left",
  cluster_rows = FALSE,
  row_labels = labels,
  row_names_gp = gpar(fontface = labelFace),

  column_split = heatmapSplit,
  column_title = c("KO PSCs", "KO Tumour", "Wt PSC", "Wt Tumour"),
  right_annotation = rowAnnot,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_column_names = FALSE,
  cluster_column_slices = FALSE,

  show_heatmap_legend = FALSE
)
heatmap <- draw(baseHeatmap, annotation_legend_list = legends)
heatmap
