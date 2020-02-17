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
P3_qPCR <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR1_clean.csv"
P3_qPCR <- read.csv(text = getURL(P3_qPCR))
P3_qPCR$X <- NULL

# load in RT-qPCRs
P3_RT <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RT-qPCR_complete.csv"
P3_RT <- read.csv(text = getURL(P3_RT))

# load in CEWE ELISA
P3_CEWE_ELISA <- ""
P3_CEWE_ELISA <- read.csv(text = getURL(P3_CEWE_ELISA))

# load in FEC ELISA
P3_FEC_ELISA <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_FEC_ELISAs/P3_112019_Eim_FEC_ELISA1_complete.csv"
P3_FEC_ELISA <- read.csv(text = getURL(P3_FEC_ELISA))
