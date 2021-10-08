# Aims:
# - Plot a heatmap of matrix index data obtained from Oliver's lab.


##Preparation ----
#Required packages are loaded
library("ComplexHeatmap")
library("circlize")

#Required data is read into the environment
miNorm <- read.csv("./Data/Oliver Matrix Index Normalised Counts.csv", stringsAsFactors = FALSE, row.names = 1)



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



##Scaled Heatmap Matrix ----
miZ <- calculateZScore(miNorm)



## Heatmap Metadata ----
miSplit <- data.frame(
  row.names = rownames(miZ),
  split = c("Malignant", "Malignant", "Malignant",
            "Malignant", "Malignant", "Malignant",
            "Protective", "Protective", "Protective",
            "Protective", "Protective", "Protective",
            "Protective", "Protective", "Protective",
            "Protective", "Protective", "Protective",
            "Protective", "Protective", "Protective",
            "Protective"
))

miLabs <- data.frame(
  row.names = colnames(miZ),
  names = c("Wt", "Wt", "Wt", "Wt", "Wt", "KO", "KO", "KO", "KO", "KO", "KO")
)



##Heatmap Custom Legends ----
miLegend <- Legend(
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  title = "Per Gene\nz Scores",
  grid_height  = unit(0.7, "cm"),
  grid_width = unit(0.7, "cm"),
  labels_gp = gpar(fontsize = 11),
  title_gp = gpar(fontsize = 11)
)



##Plotting the Heatmap
miBaseHeatmap <- Heatmap(
  miZ,
  name = "Per Gene Z-Scores",
  row_split = miSplit,
  cluster_row_slices = FALSE,
  clustering_method_columns = "complete",
  clustering_distance_columns = "euclidean",
  column_labels = miLabs$names,
  show_heatmap_legend = FALSE
)

miHeatmap <- draw(miBaseHeatmap, annotation_legend_list = miLegend)
