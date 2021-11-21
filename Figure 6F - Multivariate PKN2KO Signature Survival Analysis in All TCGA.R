# Aims:
# - Perform multivariate survival analysis on all TCGA projects using the PKN2KO Signature split
# - Plot a forest plot of this survival analysis

#Multivariates (where available):
# - Age --> >60 vs <60
# - Stage --> low stage vs high stage
# - Grade --> low grade vs high grade
# - Pharmaceutical Treatment --> treatment vs no treatment
# - Radiation Treatment --> treatment vs no treatment



## Preparation ----
#Required packages are loaded
library("TCGAbiolinks") #Used to download treatment data
library("survival") #Used for the cox proportional hazards model
library("ggplot2") #Used to plot the final plot

#Clinical data read into the environment
rawClinicData <- read.csv("./../common/PatientData/clinicalDataLiuetal.csv", stringsAsFactors = FALSE)

#This is the PKN2-KO signature from the paper (figure 6D)
PKN2KOSig <- read.csv("./Data/PKN2KO Signature.csv", stringsAsFactors = FALSE)$human_ensembl_gene_id

#A list of TCGA projects
TCGAProjects <- levels(factor(rawClinicData$type))
#Remove
TCGAProjects <- TCGAProjects[TCGAProjects != c("LAML")]

#An empty table where the results will be stored
coxTable <- data.frame(row.names = "Initialisation", "HR" = 1 , "lowerCI" = 1, "upperCI" = 1, "p" =  1, "Wald" = 1)



##Readings ----
#Some readings from the data that were required to tidy the data
levels(factor(rawClinicData$ajcc_pathologic_tumor_stage))
levels(factor(rawClinicData$histological_grade))



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



#Calculate PKN2KO score function
#A function to read in count data and calculate the PKN2KO score for each patient
calculatePKN2KOScore <- function(signature, project){
  #Reads the normalised count matrix into the environment and extracts the genes of interest
  counts <- read.csv(paste0("./../common/NormalisedCountMatrices/TCGA-", project, "-Normalised-Primary_Tumor.csv"), stringsAsFactors = FALSE, row.names = 1)
  counts <- counts[PKN2KOSig, ]

  #Calculates the per-gene z-scores for the counts
  zScores <- calculateZScore(counts)

  #Calculates the PKN2KO signature scores
  PKN2KOScore <- data.frame("score" = colSums(zScores))
  #Renames patients to match clinical data
  PKN2KOScore$patient <- gsub(".", "-", substr(rownames(PKN2KOScore), 1, 12), fixed = TRUE)
  return(PKN2KOScore)
}



##Create clinic data function ----
#A function to create a table of patient data for the survival analysis
createClinic <- function(rawClinicData, project){
  #Extracts the clinic data for the project from the Liuetal resource
  projectClinic <- rawClinicalData[rawClinicalData$type == project, ]

  #Requires a TCGA-biolinks request to acquire treatment data
  clinicQuery <- GDCquery_clinic(
    project = paste0("TCGA-", project),
    type = "clinical"
  )
  row.names(clinicQuery) <- clinicQuery$submitter_id

  #Orders treatment information so that it is the same as the Liuetal resource
  treatment <- data.frame(
    "patient" = projectClinic$bcr_patient_barcode,
    "pharmaceutical" = clinicQuery[projectClinic$bcr_patient_barcode, "treatments_pharmaceutical_treatment_or_therapy"],
    "radiation" = clinicQuery[projectClinic$bcr_patient_barcode, "treatments_radiation_treatment_or_therapy"]
  )

  #Creates the dataframe with the final information to be returned. Tidies up all the information into categories to compare.
  #Uses a lot of case_when() to reduce the categories down to 2.
  clinicData <- data.frame(
    "patient" = projectClinic$bcr_patient_barcode,
    "status" = as.numeric(projectClinic$OS),
    "time" = as.numeric(projectClinic$OS.time),
    "age" = factor(ifelse(projectClinic$age_at_initial_pathologic_diagnosis < 60, "Low", "High")),
    "stage" = factor(case_when(
      #Low stage
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IS" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage 0" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage I" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IA" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IB" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage II" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IIA" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IIB" ~ "Stage I/II",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IIC" ~ "Stage I/II",
      #High stage
      projectClinic$ajcc_pathologic_tumor_stage == "Stage III" ~ "Stage III/IV",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IIIA" ~ "Stage III/IV",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IIIB" ~ "Stage III/IV",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IIIC" ~ "Stage III/IV",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IV" ~ "Stage III/IV",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IVA" ~ "Stage III/IV",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IVB" ~ "Stage III/IV",
      projectClinic$ajcc_pathologic_tumor_stage == "Stage IVC" ~ "Stage III/IV",
      #Unknown data
      projectClinic$ajcc_pathologic_tumor_stage == "[Not Available]" ~ NA_character_,
      projectClinic$ajcc_pathologic_tumor_stage == "[Not Applicable]" ~ NA_character_,
      projectClinic$ajcc_pathologic_tumor_stage == "[Discrepancy]" ~ NA_character_,
      projectClinic$ajcc_pathologic_tumor_stage == "[Unknown]" ~ NA_character_,
      projectClinic$ajcc_pathologic_tumor_stage == "Stage X" ~ NA_character_
    )),
    "grade" = factor(case_when(
      #Low/Intermediate grade
      projectClinic$histological_grade == "G1" ~ "G1/2",
      projectClinic$histological_grade == "G2" ~ "G1/2",
      projectClinic$histological_grade == "Low Grade" ~ "G1/2",
      #High grade
      projectClinic$histological_grade == "G3" ~ "G3/4",
      projectClinic$histological_grade == "G4" ~ "G3/4",
      projectClinic$histological_grade == "High Grade" ~ "G3/4",
      #Unknown data
      projectClinic$histological_grade == "GX" ~ NA_character_,
      projectClinic$histological_grade == "[Discrepancy]" ~ NA_character_,
      projectClinic$histological_grade == "[Not Available]" ~ NA_character_,
      projectClinic$histological_grade == "[Unknown]" ~ NA_character_,
      projectClinic$histological_grade == "GB" ~ NA_character_
    )),
    "pharmaceutical" = factor(case_when(
      treatment$pharmaceutical == "yes" ~ "yes",
      treatment$pharmaceutical == "no" ~ "no",
      treatment$pharmaceutical == "not reported" ~ NA_character_,
      TRUE ~ NA_character_
    )),
    "radiation" = factor(case_when(
      treatment$radiation == "yes" ~ "yes",
      treatment$radiation == "no" ~ "no",
      treatment$radiation == "not reported" ~ NA_character_,
      TRUE ~ NA_character_
    ))
  )

  return(clinicData)
}



##Main loop ----
#A loop to go through each of the TCGA projects and collect the results for each
for(project in TCGAProjects){
  #Gives an idea of where the loop is currently at
  print(project)

  #Calculates the PKN2KO score for the patients
  lScore <- calculatePKN2KOScore(PKN2KOSig, project)

  #Prepares the clinical data for the project
  lClinic <- createClinic(rawClinicalData, project)

  #Removes any score data that is not in the clinical data and vice versa
  lScore <- lScore[lScore$patient %in% lClinic$patient, ]
  lClinic <- lClinic[lClinic$patient %in% lScore$patient, ]

  #Merges the clinical and score data together
  lScoreClinic <- merge(lScore, lClinic, by = "patient")

  #Calculates and applies the optimum cutpoint
  lCuts <- surv_cutpoint(data = lScoreClinic, time = "time", event = "status", variables = "score")
  lScoreClinic$scoreGroup <- factor(ifelse(lScoreClinic$score < lCuts$cutpoint$cutpoint, "low", "high"), levels = c("low", "high"))

  #A dynamic formula is created depending on what clinical information is available
  if(!all(is.na(lScoreClinic$stage))) lStage <- " + stage" else lStage <- NULL
  if(!all(is.na(lScoreClinic$grade))) lGrade <- " + grade" else lGrade <- NULL
  if(!all(is.na(lScoreClinic$pharmaceutical)) & length(levels(factor(lScoreClinic$pharmaceutical))) > 1) lPharmaceutical <- " + pharmaceutical" else lPharmaceutical <- NULL
  if(!all(is.na(lScoreClinic$radiation)) & length(levels(factor(lScoreClinic$radiation))) > 1) lRadiation <- " + radiation" else lRadiation <- NULL
  lFormula <- paste0("Surv(time, status) ~ scoreGroup + age", lStage , lGrade, lPharmaceutical, lRadiation)
  print(paste0(lFormula))

  #The formula created is used in the coxph calculations on the score and clinical data
  lCoxph <- coxph(as.formula(lFormula), data = lScoreClinic)
  lSummary <- summary(lCoxph)
  print(lCoxph)

  #Adds the required information into a dataframe
  lCoxTable <- data.frame(
    row.names = project,
   "HR" = lSummary$coef[1, "exp(coef)"] ,
   "lowerCI" = lSummary$conf.int[1 ,"lower .95"],
   "upperCI" = lSummary$conf.int[1 ,"upper .95"],
   "p" =  lSummary$coef[1, "Pr(>|z|)"],
   "Wald" = lSummary$coef[1, "z"]
  )

  #Combines the overall dataframe with new row
  coxTable <- rbind(coxTable, lCoxTable)
}



##Plot data preparation ----
#Removes the initialisation row
coxTable <- coxTable[rownames(coxTable) != "Initialisation", ]
#Makes a copy of the table so the code can be run again without creating problems
coxTableCopy <- coxTable

#Significance p labels are calculated
p <- case_when(coxTable$p >= 0.05 ~ "      ",
               coxTable$p < 0.0001 ~ " ****",
               coxTable$p < 0.001 ~ " *** ",
               coxTable$p < 0.01 ~ " **  ",
               coxTable$p < 0.05 ~ " *   ")

#Pastes the project name an significance together
coxTable$project <- paste0(rownames(coxTableCopy), p)

#If the range values exceed the maximum then it creates an arrow
coxTable$upperArrows <- ifelse(coxTableCopy$upperCI >= 10, 10, NA)
coxTable$lowerArrows <- ifelse(coxTableCopy$lowerCI <= 0.1,  0.1, NA)

#Sets range values to the maximum
coxTable[coxTableCopy$upperCI >= 10, "upperCI"] <- 10
coxTable[coxTableCopy$lowerCI <= 0.1, "lowerCI"] <- 0.1



##Final plot ----
#Plotted using ggplot2
hrPlot <- ggplot(coxTable, aes(x = HR, y = reorder(project, Wald), xmin = lowerCI, xmax = upperCI)) +
  geom_point(size = 3) +
  geom_linerange(size = 0.8) +
  labs(x = "Hazard Ratio", y = "TCGA Project") +
  geom_vline(xintercept = 1, linetype = 2, size = 1) +
  theme_bw() +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  geom_segment(aes(x = lowerArrows+0.01, xend = lowerArrows, y = project, yend = project), linejoin = "mitre", lineend = "butt", size = 1, arrow = arrow(length = unit(0.032, "npc"))) +
  geom_segment(aes(x = upperArrows-0.01, xend = upperArrows, y = project, yend = project), linejoin = "mitre", lineend = "butt", size = 1, arrow = arrow(length = unit(0.032, "npc"))) +
  theme(axis.text.y = element_text(size = 11, face = "bold"))
hrPlot
