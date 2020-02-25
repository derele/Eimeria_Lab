# P3 combining script for weight, oocysts, qPCR, RT-qPCR, ELISA and hopefully FACS
library(httr)
library(RCurl)
library(dplyr)
library(Rmisc)

# load in weight and oocysts
P3_weightANDoocysts <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_Weight%26Oocyst_complete.csv"
P3_weightANDoocysts <- read.csv(text = getURL(P3_weightANDoocysts))
P3_weightANDoocysts$X <- NULL

# load in qPCRs
P3_qPCR <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR_complete.csv"
P3_qPCR <- read.csv(text = getURL(P3_qPCR))
P3_qPCR$X <- NULL

# load in RT-qPCRs
P3_RT <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RTqPCR_complete.csv"
P3_RT <- read.csv(text = getURL(P3_RT))
P3_RT$X <- NULL

# load in CEWE ELISA (important to merge CEWE ELISAs with qPCR and RTqPCR to give them labels)
P3_CEWE_ELISA <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISAs_complete.csv"
P3_CEWE_ELISA <- read.csv(text = getURL(P3_CEWE_ELISA))
P3_CEWE_ELISA$X <- NULL
colnames(P3_CEWE_ELISA)[2] <- "IFNy_CEWE"

# load in FEC ELISA
P3_FEC_ELISA <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_FEC_ELISAs/P3_112019_Eim_FEC_ELISA1_complete.csv"
P3_FEC_ELISA <- read.csv(text = getURL(P3_FEC_ELISA))
P3_FEC_ELISA$X <- NULL
colnames(P3_FEC_ELISA)[2] <- "IFNy_FEC"

# load in qPCR MCs for Eimeria
P3_MC <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7%26P3_Eim_MCs.csv"), sep = ";")
P3_MC$X <- NULL
P3_MC$X.1 <- NULL

# start merging (important to merge CEWE ELISAs with qPCR and RTqPCR to give them labels)
P3 <- merge(P3_CEWE_ELISA, P3_qPCR)
P3 <- merge(P3, P3_RT)
P3 <- merge(P3, P3_MC)
names(P3_weightANDoocysts)[names(P3_weightANDoocysts) == "labels"] <- "label"
# P3bXAK is an unmeasured early sacrifice
P3 <- merge(P3, P3_weightANDoocysts, all = T)
names(P3_FEC_ELISA)[names(P3_FEC_ELISA) == "labels"] <- "label"
P3 <- merge(P3, P3_FEC_ELISA, all = T)

# up to date most complete P3 dataset
write.csv(P3, "./Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_COMPLETE.csv")
