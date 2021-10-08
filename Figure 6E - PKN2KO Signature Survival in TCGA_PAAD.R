# Aims:
# - Plot a Kaplan Meier graph comparing low and high PKN2KO signature in TCGA-PAAD data.


## Preparation ----
#Required packages are loaded
library("TCGAbiolinks")
library("biomaRt")
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
PKN2KOSig <- c("Col6a3", "Fmod", "Mmp28", "Prelp", "Serping1", "Col4a1", "Gpc1", "Megf10", "Itga7", "Serpinb8", "Pcdh7")



##Downloading the clinic data from TCGA-biolinks
#Creating a clinic query for TCGA-PAAD
clinicQuery <- GDCquery_clinic(
  project = "TCGA-PAAD",
  type = "clinical"
)
#Creates a dataframe containg patient ID, survival status and time
clinicData <- data.frame(
  row.names = clinicQuery$submitter_id,
  status = case_when(
    clinicQuery$vital_status == "Alive" ~ 0,
    clinicQuery$vital_status == "Dead" ~ 1),
  time = case_when(
    clinicQuery$vital_status == "Alive" ~ clinicQuery$days_to_last_follow_up,
    clinicQuery$vital_status == "Dead" ~ clinicQuery$days_to_death))













cuts <- surv_cutpoint(data = basalPAADScoresClinic, time = "time", event = "status", variables = "Score")
basalPAADScoresClinic <- surv_categorize(cuts)

fit <- survfit(Surv(time, status) ~ Score, data = basalPAADScoresClinic)

survPlot <- ggsurvplot(
  fit,
  data = basalPAADScoresClinic,
  risk.table = T,
  pval = T, #p,
  pval.method = T,
  conf.int = F,
  xlab = "Time (Days)",
  palette = c("red", "blue"),
  break.time.by = 500,
  ggtheme = theme_bw()
  #surv.median.line = "hv",
  #ncensor.plot = TRUE,
)
