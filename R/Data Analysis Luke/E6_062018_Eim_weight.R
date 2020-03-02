# process E6 weight and oocysts
library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)

E6_weight <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E6_062018_Eim_full_RECORDweight.csv"))
E6_oocyst <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E6_062018_Eim_full_RECORDoocysts.csv"))
# replace "na" and "" with NA
E6_weight[E6_weight=="na"]<- NA
E6_weight[E6_weight==""]<- NA
E6_weight$fecweight <- as.numeric(as.character(E6_weight$fecweight))
# replace "n.a." with NA and convert oocyst squares to numbers
E6_oocyst[E6_oocyst=="n.a."]<- NA
E6_oocyst$oocyst_sq1 <- as.numeric(as.character(E6_oocyst$oocyst_sq1))
E6_oocyst$oocyst_sq2 <- as.numeric(as.character(E6_oocyst$oocyst_sq2))
E6_oocyst$oocyst_sq3 <- as.numeric(as.character(E6_oocyst$oocyst_sq3))
E6_oocyst$oocyst_sq4 <- as.numeric(as.character(E6_oocyst$oocyst_sq4))
# calculate oocyst mean and get rid of unecessary columns
E6_oocyst$oocyst_mean <- ((E6_oocyst$oocyst_sq1 
                             + E6_oocyst$oocyst_sq2 
                             + E6_oocyst$oocyst_sq3 
                             + E6_oocyst$oocyst_sq4) / 4) * 
  10000 * # because volume chamber
  E6_oocyst$dilution
# clean up and merge
E6_oocyst <- select(E6_oocyst, labels, oocyst_mean, OPG, Expe)
E6_weight <- select(E6_weight, EH_ID, labels, dpi, weight_dpi0, weightloss, weight, Expe, Eimeria, fecweight)
E6_WandO <- merge(E6_oocyst, E6_weight, all = T)

# calculate OPG
E6_WandO$OPG <- E6_WandO$oocyst_mean / E6_WandO$fecweight

# add mouse info and clean up before merge
E6_1a <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E5_1a_062018_Eim_DESIGN.csv"))
E6_1b <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E5_1b_062018_Eim_DESIGN.csv"))
E6_2a <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E5_2a_062018_Eim_DESIGN.csv"))
E6_2b <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E5_2b_062018_Eim_DESIGN.csv"))

E6_design <- rbind(E6_1a, E6_1b)
E6_design <- rbind(E6_design, E6_2a)
E6_design <- rbind(E6_design, E6_2b)

E6_design <- select(E6_design, Strain, HybridStatus, Batch, InfectionStrain, EH_id)
names(E6_design)[names(E6_design) == "EH_id"] <- "EH_ID"

E6_design$EH_ID <- as.character(E6_design$EH_ID)
E6_WandO$EH_ID <- as.character(E6_WandO$EH_ID)
E6_design$EH_ID <- gsub(" ", "", E6_design$EH_ID)
E6_WandO$EH_ID <- gsub(" ", "", E6_WandO$EH_ID)

E6_complete <- merge(E6_design, E6_WandO, all = T, by = "EH_ID")
# quick quality check, remove empty rows, rename columns to match other combine tables, etc.
E6_complete <- distinct(E6_complete)
E6_complete <- E6_complete[-c(1361:1382),]
E6_complete$EH_ID <- gsub("LM", "LM_", E6_complete$EH_ID)
names(E6_complete)[names(E6_complete) == "Eimeria"] <- "primary"
names(E6_complete)[names(E6_complete) == "Expe"] <- "EXP"
E6_complete$EXP <- as.character(E6_complete$EXP)
E6_complete[E6_complete=="Expe005_1a"] <- "5.1a"
E6_complete[E6_complete=="Expe005_2a"] <- "5.2a"
E6_complete[E6_complete=="Expe005_1b"] <- "5.1b"
E6_complete[E6_complete=="Expe005_2b"] <- "5.2b"
E6_complete$labels <- paste(E6_complete$EXP, E6_complete$labels, sep =  "")
E6_complete$challenge <- NA
E6_complete$weight <- as.numeric(as.character(E6_complete$weight))
E6_complete$Wchange <- (E6_complete$weight/E6_complete$weight_dpi0)*100


# write out
write.csv(E6_complete, "./Eimeria_Lab/data/3_recordingTables/E6_062018_Eim_WandO_complete.csv")
