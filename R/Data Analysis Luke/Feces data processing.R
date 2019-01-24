# setwd("Repositories/Eimeria_Lab/R/Data Analysis Luke/")

#load csv files#
Exp007a_design <- read.csv("../../data/2_designTables/Inf2a_Exp005.DESIGN.csv")
Exp007b_design <- read.csv("../../data/2_designTables/Inf2b_Exp005.DESIGN.csv")
E7aF <- read.csv("../../data/3_recordingTables/Exp007/Exp007a/Exp_007a feces.csv")
E7bF <- read.csv("../../data/3_recordingTables/Exp007/Exp007b/Exp_007b feces.csv")

# the columns we want to keep
col2keep <- c("Strain", "HybridStatus", "InfectionStrain", "EH_id")

Exp007a_design <- Exp007a_design[col2keep]
Exp007b_design <- Exp007b_design[col2keep]

# rename EH_id to EH_ID#
names(Exp007a_design)[names(Exp007a_design) == "EH_id"] <- "EH_ID"
names(Exp007b_design)[names(Exp007b_design) == "EH_id"] <- "EH_ID"

# let's make one big fat Expe007 design table
Expe007_design <- rbind(Exp007a_design, Exp007b_design)

# remove shit columns
E7aF <- E7aF[-grep(pattern = "X", x = names(E7aF))]
E7bF <- E7bF[-grep(pattern = "X", x = names(E7bF))]

# keep the batch information
E7aF$batch <- "october2018"
E7bF$batch <- "december2018"

# Make one big fat table Expe 7
Expe007_record <- rbind(E7aF, E7bF)

# Merge all
Expe007 <- merge(Expe007_design, Expe007_record)

