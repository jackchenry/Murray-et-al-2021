# Aims:
# - Calculate the PKN2KO Signature score for the mouse pancreatic tumours to confirm signature
# - Save data for later plotting in graphpad



## Preparation ----
#Required data are loaded
#This is the normalised count data from the pancreatic tumours (Figure S5)
tumNorm <- read.csv("./Data/InVivo Tumour Normalised Count Matrix.csv", stringsAsFactors = FALSE, row.names = 1)
#This is the PKN2-KO signature from the paper (Figure 6D)
PKN2KOSig <- read.csv("./Data/PKN2KO Signature.csv", stringsAsFactors = FALSE)$mouse_ensembl_gene_id



##Calculate Zscore function ----
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



##Calculating the score ----
#The important genes are extracted
tumPKN2Sig <- tumNorm[PKN2KOSig, ]
tumPKN2SigZ <- calculateZScore(tumPKN2Sig)
finalScore <- colSums(tumPKN2SigZ)



#Saving data ----
write.csv(finalScore, "InVivo Tumour PKN2KOSig Score.csv")