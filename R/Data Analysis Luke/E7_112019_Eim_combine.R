# E7 combining script for weight, oocysts, qPCR, RT-qPCR, ELISA and hopefully FACS
library(httr)
library(RCurl)
library(dplyr)
library(Rmisc)
library(ggplot2)
library(tidyverse)
library(modelr)
library(ggrepel)

# load in weight and oocysts
E7_weightANDoocysts <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_Weight%26Oocyst_complete.csv"
E7_weightANDoocysts <- read.csv(text = getURL(E7_weightANDoocysts))
E7_weightANDoocysts$X <- NULL
E7_weightANDoocysts <- distinct(E7_weightANDoocysts)

################## at the moment one and the same
# load in qPCRs 
E7_qPCR <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT_and_qPCR_complete.csv"
E7_qPCR <- read.csv(text = getURL(E7_qPCR))
E7_qPCR$X <- NULL
E7_qPCR <- distinct(E7_qPCR)
# creating a mess when merging so let"s reduce it to only the necessary
E7_qPCR <- dplyr::select(E7_qPCR, EH_ID, delta, CXCR3, IL.12, IRG6, labels, Caecum, Eim_MC, infHistory)
E7_qPCR <- distinct(E7_qPCR)
# # load in RT-qPCRs
# P3_RT <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RT-qPCR_complete.csv"
# P3_RT <- read.csv(text = getURL(P3_RT))
# P3_RT$X <- NULL
###################################################

# # load in CEWE ELISA
E7_CEWE_ELISA <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_CEWE_ELISAs/E7_112018_Eim_CEWE_ELISA_complete.csv"
E7_CEWE_ELISA <- read.csv(text = getURL(E7_CEWE_ELISA))
E7_CEWE_ELISA$X <- NULL
colnames(E7_CEWE_ELISA)[1] <- "labels"
colnames(E7_CEWE_ELISA)[2] <- "IFNy_CEWE"
E7_CEWE_ELISA <- distinct(E7_CEWE_ELISA)

# load in FEC ELISA
# missing labels fix in E7 FEC ELISA script
E7_FEC_ELISA <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISAs_complete.csv"
E7_FEC_ELISA <- read.csv(text = getURL(E7_FEC_ELISA))
E7_FEC_ELISA$X <- NULL
colnames(E7_FEC_ELISA)[2] <- "IFNy_FEC"
colnames(E7_FEC_ELISA)[1] <- "labels"
E7_FEC_ELISA <- distinct(E7_FEC_ELISA)

# merge (E7 vs E7a and E7b makes mess, go back to E7_FEC_ELISA and rewrite labels (consult collection info table and boxes))
E7 <- merge(E7_weightANDoocysts, E7_FEC_ELISA, all = T)
E7 <- distinct(E7)

E7 <- merge(E7, E7_CEWE_ELISA, all.x = T)
E7 <- distinct(E7)

E7 <- merge(E7, E7_qPCR)
E7 <- distinct(E7)

write.csv(E7, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_COMPLETE.csv")
write.csv(E7, "~/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_COMPLETE.csv")

