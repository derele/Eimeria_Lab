# RT-qPCR_clean data processing for E1_012017, additional data from CEC
library(RCurl)
library(httr)
library(Rmisc)

#load in data from GitHub
RT-qPCRurl<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/20_05_2019_Samples_36_37_38_39/Luke_2019_05_21_Enas2017samples_36_37_38_39%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_36_37_38_39 <- read.csv(text = getURL(Samples_36_37_38_39url), sep = ";")
