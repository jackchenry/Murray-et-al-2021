# Aims:
# - Perform individual univariate survival analysis on all significant genes from the Qiaseq panel in TCGA-PAAD
# - Plot a forest plot of this survival analysis



## Preparation ----
#Required packages are loaded
library("biomaRt") #Used to convert mouse to human ensebl ids
library("survival") #Used for the cox proportional hazards model
library("ggplot2") #Used to plot the final plot

#Clinical data read into the environment
rawClinicData <- read.csv("./../common/PatientData/clinicalDataLiuetal.csv", stringsAsFactors = FALSE, row.names = 1)

#The normalised TCGA-PAAD counts are read into the environment
paadNorm <- read.csv("./../common/NormalisedCountMatrices/TCGA-PAAD-Normalised-Primary_Tumor.csv", stringsAsFactors = FALSE, row.names = 1)

#PSC differential expression results from Figure 2. Required for the list of genes.
pscDE <- read.csv("./Data/PSCs Basal DE Results.csv", stringsAsFactors = FALSE, row.names = 1)
#Only the significant genes are kept
pscDE <- pscDE[!is.na(pscDE$padj) & pscDE$padj < 0.05, ]

#The row order of the
rowOrder <- read.csv("./Data/Tumour PSC Heatmap Row Order", stringsAsFactors = FALSE)$x

#An empty table where the results will be stored
coxTable <- data.frame(row.names = "Initialisation", "HR" = 1 , "lowerCI" = 1, "upperCI" = 1, "p" =  1, "Wald" = 1)



##Clinic Data ----
rawClinicData <- rawClinicData[rawClinicData$type == "PAAD", c("bcr_patient_barcode", "OS", "OS.time")]
clinicData <- data.frame(
  row.names = rawClinicData$bcr_patient_barcode,
  "status" = as.numeric(rawClinicData$OS),
  "time" = as.numeric(rawClinicData$OS.time)
)



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



##Converting mouse ensemble IDs into human
#Mouse and human ensemble marts
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Uses getLDS to convert from mouse to human
conversionTable <- getLDS(
  attributes = c("ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = row.names(pscDE),
  mart = mouse,
  attributesL = c("ensembl_gene_id", "hgnc_symbol"),
  martL = human,
  uniqueRows=T
)
#Rename columns to make sense
colnames(conversionTable) <- c("mouse_ensembl_gene_id", "human_ensembl_gene_id", "hgnc_symbol")
#Reorder table so that it matches PSC data
conversionTable <- conversionTable[match(row.names(pscDE), conversionTable$mouse_ensembl_gene_id), ]



##Create zScores for each gene ----
pscPAADNorm <- paadNorm[conversionTable$human_ensembl_gene_id, ]
pscPAADZ <- calculateZScore(pscPAADNorm)
#Removes extra sample information so just left with patient identifier
colnames(pscPAADZ) <- gsub(".", "-", substr(colnames(pscPAADZ), 1, 12), fixed = TRUE)
#Re-orders the data so that it matches the heatmap
pscPAADZ <- pscPAADZ[match(conversionTable[match(rowOrder, conversionTable$mouse_ensembl_gene_id), "human_ensembl_gene_id"], row.names(pscPAADZ)), ]



##Main Loop ----
#A loop to go through each of the genes in the panel and collect the survival results for each
#Using apply would be faster but loops are more understandable
for(gene in row.names(pscPAADZ)){
  #Extracts the HGCN name for the gene
  lName <- conversionTable[match(gene, conversionTable$human_ensembl_gene_id), "hgnc_symbol"]

  #Extracts the zScores for that gene
  lScore <- data.frame("score" = pscPAADZ[gene, ])

  #Merges the clinic data and score data together
  lClinicScore <- merge(clinicData, lScore, by = "row.names")
  rownames(lClinicScore) <- lClinicScore$Row.names
  lClinicScore$Row.names <- NULL

  #The best cutpoint for this gene is calculated and applied
  lCuts <- surv_cutpoint(data = lClinicScore, time = "time", event = "status", variables = "score")
  lClinicScore <- surv_categorize(lCuts)
  lClinicScore$score <- factor(lClinicScore$score, levels = c("low", "high"))

  #Univariate Coxph is calculated
  lCoxph <- summary(coxph(Surv(time, status) ~ score, data = lClinicScore))
  lCoxTable <- data.frame(row.names = lName, "HR" = lCoxph$coef[1, "exp(coef)"] ,
                               "lowerCI" = lCoxph$conf.int[1 ,"lower .95"],
                               "upperCI" = lCoxph$conf.int[1 ,"upper .95"],
                               "p" =  lCoxph$coef[1, "Pr(>|z|)"],
                               "Wald" = lCoxph$coef[1, "z"])
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
#To preserve the order
coxTable$project <- factor(coxTable$project, levels = rev(coxTable$project))

#If the range values exceed the maximum then it creates an arrow
coxTable$upperArrows <- ifelse(coxTableCopy$upperCI >= 10, 10, NA)
coxTable$lowerArrows <- ifelse(coxTableCopy$lowerCI <= 0.1,  0.1, NA)

#Sets range values to the maximum
coxTable[coxTableCopy$upperCI >= 10, "upperCI"] <- 10
coxTable[coxTableCopy$lowerCI <= 0.1, "lowerCI"] <- 0.1



##Final plot ----
#Plotted using ggplot2
hrPlot <- ggplot(coxTable, aes(x = HR, y = project, xmin = lowerCI, xmax = upperCI)) +
  geom_point(size = 3) +
  geom_linerange(size = 0.8) +
  labs(title = "Univarite Cox Proportional Hazards\nper Gene in TCGA-PAAD", subtitle = "Best KM Split", x = "Hazard Ratio", y = "Gene") +
  geom_vline(xintercept = 1, linetype = 2, size = 1) +
  theme_bw() +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  geom_segment(aes(x = lowerArrows+0.01, xend = lowerArrows, y = project, yend = project), linejoin = "mitre", lineend = "butt", size = 1, arrow = arrow(length = unit(0.032, "npc"))) +
  geom_segment(aes(x = upperArrows-0.01, xend = upperArrows, y = project, yend = project), linejoin = "mitre", lineend = "butt", size = 1, arrow = arrow(length = unit(0.032, "npc"))) +
  theme(axis.text.y = element_text(size = 11, face = "bold"))
hrPlot
