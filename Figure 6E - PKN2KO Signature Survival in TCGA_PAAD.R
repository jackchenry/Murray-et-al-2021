# Aims:
# - Plot a Kaplan Meier graph comparing low and high PKN2KO signature in TCGA-PAAD data.



## Preparation ----
#Required packages are loaded
library("TCGAbiolinks")
library("dplyr")
library("survival")
library("survminer")

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

#This is the PKN2-KO signature from the paper
PKN2KOSig <- read.csv("./Data/PKN2KO Signature.csv", stringsAsFactors = FALSE)

#This is normalised TCGA-PAAD RNAseq count data
PAAD <- read.csv(
  "./../common/NormalisedCountMatrices/TCGA-PAAD-Normalised-Primary_Tumor.csv",
  stringsAsFactors = FALSE, row.names = 1
)



##Scoring patients for the PKN2KO Signature
#PKN2KO Signature genes are selected from the TCGA-PAAD count data
PAADSig <- PAAD[PKN2KOSig$human_ensembl_gene_id, ]
#Gene-wise z-scores are calculated for each signature gene
PAADSigZ <- calculateZScore(PAADSig)
#Patients are scored by the sum of the signature genes
PKN2KOScore <- data.frame(score = colSums(PAADSigZ))
#Patient column is made that will match with the clinic data information
PKN2KOScore$patient <- gsub(".", "-", substr(rownames(PKN2KOScore), 1, 12), fixed = TRUE)
row.names(PKN2KOScore) <- NULL



##Downloading the clinic data from TCGA-biolinks
#Creating a clinic query for TCGA-PAAD
clinicQuery <- GDCquery_clinic(
  project = "TCGA-PAAD",
  type = "clinical"
)
#Creates a dataframe containing patient ID, survival status and time
clinicData <- data.frame(
  patient = clinicQuery$submitter_id,
  status = case_when(
    clinicQuery$vital_status == "Alive" ~ 0,
    clinicQuery$vital_status == "Dead" ~ 1),
  time = case_when(
    clinicQuery$vital_status == "Alive" ~ clinicQuery$days_to_last_follow_up,
    clinicQuery$vital_status == "Dead" ~ clinicQuery$days_to_death),
  stage = clinicQuery$ajcc_pathologic_stage
  grade = clinicQuery$ajcc_pathologic_t)



##Preparing the survival data
#Clinical data is merged with PKN2KO signature score data
survData <- merge(clinicData, PKN2KOScore, by = "patient")
#Optimal cuts are determined
cuts <- surv_cutpoint(data = survData, time = "time", event = "status", variables = "score")
fitData <- surv_categorize(cuts)



##Fitting the survival data
fit <- survfit(Surv(time, status) ~ score, data = fitData)



##Plotting the KM curve
survPlot <- ggsurvplot(
  fit,
  data = fitData,
  risk.table = T,
  pval = T,
  pval.method = T,
  conf.int = F,
  xlab = "Time (Days)",
  palette = c("red", "blue"),
  break.time.by = 500,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
survPlot



##Saving annotated data
survData$split <- ifelse(survData$score > cuts$cutpoint$cutpoint, "High", "Low")
write.csv(survData, "./Data/TCGA-PAAD PKN2KO Signature Data.csv", row.names = FALSE)
