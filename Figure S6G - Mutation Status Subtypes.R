# Aims:
# - Identify whether the PKN2KO signature enriches for any of the most commonly mutated genes

#Mutated genes:
# - KRAS
# - CDKN2A
# - SMAD4
# - TP53



##Preparation ----
#Required data loaded into the environment
#This is TCGA-PAAD patients labelled with their PKN2KO signature score created in Figure 6E
PAADSig <- read.csv("./Data/TCGA-PAAD PKN2KO Signature Data.csv", stringsAsFactors = FALSE)
#This is TCGA-PAAD patient data with labelled mutation status
mutations <- read.csv("./Data/TCGA-PAAD_Mutation_Status.csv", stringsAsFactors = FALSE)



##Data preparation ----
#Labels the PKN2KO signature data with whether or not the patient is mutated for each gene
PAADSig$KRAS <- ifelse(PAADSig$patient %in% mutations$KRAS, "Mutated", "Normal")
PAADSig$TP53 <- ifelse(PAADSig$patient %in% mutations$TP53, "Mutated", "Normal")
PAADSig$SMAD4 <- ifelse(PAADSig$patient %in% mutations$SMAD4, "Mutated", "Normal")
PAADSig$CDKN2A <- ifelse(PAADSig$patient %in% mutations$CDKN2A, "Mutated", "Normal")

#Creates a dataframe with the sum of the number of patients within each group
mutationCount <- data.frame(
  "Gene" = c("KRAS", "KRAS", "KRAS", "KRAS", "TP53", "TP53", "TP53", "TP53", "SMAD4", "SMAD4", "SMAD4", "SMAD4", "CDKN2A", "CDKN2A", "CDKN2A", "CDKN2A"),
  "Condition" = factor(c("Normal", "Normal", "Mutated", "Mutated",
                        "Normal", "Normal", "Mutated", "Mutated",
                        "Normal", "Normal", "Mutated", "Mutated",
                        "Normal", "Normal", "Mutated", "Mutated")),
  "Score" = factor(c("Low", "High", "Low", "High","Low", "High", "Low", "High", "Low", "High", "Low", "High","Low", "High", "Low", "High")),
  "Count" = c(
    sum(PAADSig$KRAS == "Normal" & PAADSig$split == "Low"),
    sum(PAADSig$KRAS == "Normal" & PAADSig$split == "High"),
    sum(PAADSig$KRAS == "Mutated" & PAADSig$split == "Low"),
    sum(PAADSig$KRAS == "Mutated" & PAADSig$split == "High"),
    sum(PAADSig$TP53 == "Normal" & PAADSig$split == "Low"),
    sum(PAADSig$TP53 == "Normal" & PAADSig$split == "High"),
    sum(PAADSig$TP53 == "Mutated" & PAADSig$split == "Low"),
    sum(PAADSig$TP53 == "Mutated" & PAADSig$split == "High"),
    sum(PAADSig$SMAD4 == "Normal" & PAADSig$split == "Low"),
    sum(PAADSig$SMAD4 == "Normal" & PAADSig$split == "High"),
    sum(PAADSig$SMAD4 == "Mutated" & PAADSig$split == "Low"),
    sum(PAADSig$SMAD4 == "Mutated" & PAADSig$split == "High"),
    sum(PAADSig$CDKN2A == "Normal" & PAADSig$split == "Low"),
    sum(PAADSig$CDKN2A == "Normal" & PAADSig$split == "High"),
    sum(PAADSig$CDKN2A == "Mutated" & PAADSig$split == "Low"),
    sum(PAADSig$CDKN2A == "Mutated" & PAADSig$split == "High")
  )
)



##Saving the data ----
#This data will be later plotted in GraphPad Prism
write.csv(mutationCount, "./Data/MutationStatusData.csv")
