# Aims:
# - Identify whether the PKN2KO signature enriches for any of the pancreatic cancer subtypes
# - Save data for later plotting in GraphPad
# - Identify whether the PKN2KO signature predicts survival within any of the pancreatic cancer subtypes

#Pancreatic cancer subtypes:
# - Moffitt - Basal/Classical
# - Collisson - Classical/QM/Exocrine
# - Bailey - Progenitor/Squamous/ADEX/Immunogenic



##Preparation ----
#The required packages are loaded
library("survival") #Used to perform the Cox proportional analysis
library("ggplot2") #Used to plot the forest plot
#Required data loaded into the environment
#This is TCGA-PAAD patients labelled with their PKN2KO signature score created in Figure 6E
PAADSig <- read.csv("./Data/TCGA-PAAD PKN2KO Signature Data.csv", stringsAsFactors = FALSE)
#This is labelled TCGA-PAAD clinical data downloaded from the pancreatic cancer expression database (PED)
PEDData <- read.csv("./Data/PED Bioinformatics Portal.csv", stringsAsFactors = FALSE)



##Preparing the data ----
#Removing the extra information from the patient identifier
PEDData$patient <- substr(PEDData$patient, 1, 12)
#Only keeping the columns of interest
PEDData <- PEDData[, c("patient", "Moffit_subtype", "Collisson_subtype", "Bailey_subtype")]

#Merging the PAAD PKN2KO signature with PED labels
PEDPKN2 <- merge(PAADSig, PEDData, by = "patient")

#Factoring the PKN2KO signature splits
PEDPKN2$split <- factor(PEDPKN2$split, levels = c("Low", "High"))



##Calculating pancreatic cancer counts ----
#Moffitt
moffit <- data.frame(
  "Study" = "Moffit",
  "Condition" = factor(c("Basal", "Basal", "Classical", "Classical", "Undefined", "Undefined")),
  "Score" = factor(c("Low", "High", "Low", "High","Low", "High")),
  "Count" = c(sum(
    PEDPKN2$Moffit_subtype  == "Basal" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Moffit_subtype  == "Basal" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Moffit_subtype  == "Classical" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Moffit_subtype  == "Classical" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Moffit_subtype  == "" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Moffit_subtype  == "" & PEDPKN2$split == "High")
  )
)
#Collisson
collisson <- data.frame(
  "Study" = "Collisson",
  "Condition" = factor(c("Classical", "Classical", "QM", "QM", "Exocrine", "Exocrine", "Undefined", "Undefined")),
  "Score" = factor(c("Low", "High", "Low", "High","Low", "High", "Low", "High")),
  "Count" = c(
    sum(PEDPKN2$Collisson_subtype  == "Classical" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Collisson_subtype  == "Classical" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Collisson_subtype  == "QM" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Collisson_subtype  == "QM" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Collisson_subtype  == "Exocrine" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Collisson_subtype  == "Exocrine" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Collisson_subtype  == "" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Collisson_subtype  == "" & PEDPKN2$split == "High")
  )
)
#Bailey
bailey <- data.frame(
  "Study" = "Bailey",
  "Condition" = factor(c("Progenitor", "Progenitor", "Squamous", "Squamous", "ADEX", "ADEX", "Immunogenic", "Immunogenic", "Undefined", "Undefined")),
  "Score" = factor(c("Low", "High", "Low", "High","Low", "High", "Low", "High", "Low", "High")),
  "Count" = c(
    sum(PEDPKN2$Bailey_subtype == "Progenitor" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Bailey_subtype == "Progenitor" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Bailey_subtype == "Squamous" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Bailey_subtype == "Squamous" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Bailey_subtype == "ADEX" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Bailey_subtype == "ADEX" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Bailey_subtype == "Immunogenic" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Bailey_subtype == "Immunogenic" & PEDPKN2$split == "High"),
    sum(PEDPKN2$Bailey_subtype == "" & PEDPKN2$split == "Low"),
    sum(PEDPKN2$Bailey_subtype == "" & PEDPKN2$split == "High")
  )
)
#Combines all the data frames together
moffitCollissonBailey <- rbind(moffit, collisson, bailey)



##Saving the data ----
#This data will be later plotted in GraphPad Prism
write.csv(moffitCollissonBailey, "./Data/MoffitCollissonBaileyData.csv")



##Preparing data for CoxPH survival analysis ----
#Extracts each overlapping patient groups and puts them into a list
coxMCB <- list(
  "Moffit Basal" = PEDPKN2[PEDPKN2$Moffit_subtype == "Basal", ],
  "Moffit Classical" = PEDPKN2[PEDPKN2$Moffit_subtype == "Classical", ],
  "Collisson Classical" = PEDPKN2[PEDPKN2$Collisson_subtype  == "Classical", ],
  "Collisson Exocrine" = PEDPKN2[PEDPKN2$Collisson_subtype  == "Exocrine", ],
  "Collisson QM" = PEDPKN2[PEDPKN2$Collisson_subtype  == "QM", ],
  "Bailey ADEX" = PEDPKN2[PEDPKN2$Bailey_subtype == "ADEX", ],
  "Bailey Immunogenic" = PEDPKN2[PEDPKN2$Bailey_subtype == "Immunogenic", ],
  "Bailey Progenitor" = PEDPKN2[PEDPKN2$Bailey_subtype == "Progenitor", ],
  "Bailey Squamous" = PEDPKN2[PEDPKN2$Bailey_subtype == "Squamous", ]
)



##CoxPH survival analysis ----
#Cox table initialisation
mcbCoxTable <- data.frame(row.names = "Initialisation", "HR" = 1 , "lowerCI" = 1, "upperCI" = 1, "p" =  1, "Wald" = 1)
#This will be for easy indexing of the name
i <- 1
for(type in coxMCB){
  #Performs univariate survival analysis
  lCoxph <- summary(coxph(Surv(time, status) ~ split, data = type))
  print(lCoxph)

  #Builds a dataframe with the information of interest
  lCoxphDF <- data.frame(
    row.names = names(coxMCB[i]),
    "HR" = lCoxph$coef[1, "exp(coef)"],
    "lowerCI" = lCoxph$conf.int[1 ,"lower .95"],
    "upperCI" = lCoxph$conf.int[1 ,"upper .95"],
    "p" =  lCoxph$coef[1, "Pr(>|z|)"],
    "Wald" = lCoxph$coef[1, "z"]
  )
  #Combines the dataframe with the previous ones
  mcbCoxTable <- rbind(mcbCoxTable, lCoxphDF)
  i <- i + 1
}



##Preparing forest plot data ----
#Removes the initialisation row
mcbCoxTable <- mcbCoxTable[rownames(mcbCoxTable) != "Initialisation", ]
#Makes a copy of the table so the code can be run again without creating problems
mcbCoxTableCopy <- mcbCoxTable

#Significance p labels are calculated
p <- case_when(mcbCoxTable$p >= 0.05 ~ "      ",
               mcbCoxTable$p < 0.0001 ~ " ****",
               mcbCoxTable$p < 0.001 ~ " *** ",
               mcbCoxTable$p < 0.01 ~ " **  ",
               mcbCoxTable$p < 0.05 ~ " *   ")

#Pastes the project name an significance together
mcbCoxTable$project <- paste0(substr(rownames(mcbCoxTableCopy), 1, 100), p)

#If the range values exceed the maximum then it creates an arrow
mcbCoxTable$upperArrows <- ifelse(mcbCoxTableCopy$upperCI >= 10, 10, NA)
mcbCoxTable$lowerArrows <- ifelse(mcbCoxTableCopy$lowerCI <= 0.1,  0.1, NA)

#Sets range values to the maximum
mcbCoxTable[mcbCoxTable$upperCI >= 10, "upperCI"] <- 10
mcbCoxTable[mcbCoxTable$lowerCI <= 0.1, "lowerCI"] <- 0.1



##Final forest plot ----
#Plotted using ggplot2
mcbHRPlot <- ggplot(mcbCoxTable, aes(x = HR, y = reorder(project, Wald), xmin = lowerCI, xmax = upperCI)) +
  geom_point(size = 3) +
  geom_linerange(size = 0.8) +
  labs(title = "Univariate Cox Proportional Hazards", subtitle = "Moffit/Collison/Bailey Subtypes", x = "Hazard Ratio", y = "Subtypes") +
  geom_vline(xintercept = 1, linetype = 2, size = 1) +
  theme_bw() +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  geom_segment(aes(x = lowerArrows+0.01, xend = lowerArrows, y = project, yend = project), linejoin = "mitre", lineend = "butt", size = 1, arrow = arrow(length = unit(0.032, "npc"))) +
  geom_segment(aes(x = upperArrows-0.01, xend = upperArrows, y = project, yend = project), linejoin = "mitre", lineend = "butt", size = 1, arrow = arrow(length = unit(0.032, "npc"))) +
  theme(axis.text.y = element_text(size = 11, face = "bold"))
mcbHRPlot
