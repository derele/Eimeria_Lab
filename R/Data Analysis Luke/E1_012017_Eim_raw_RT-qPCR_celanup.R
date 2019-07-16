# raw RT-qPCR cleanup for E1_012017, additional data from CEC
library(RCurl)
library(httr)
library(Rmisc)

#load in data from GitHub
Samples36_37_38_39url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/20_05_2019_Samples_36_37_38_39/Luke_2019_05_21_Enas2017samples_36_37_38_39%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples36_37_38_39 <- read.csv(text = getURL(Samples36_37_38_39url), sep = ";")

Samples_C_27_74_95url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/20_05_2019_Samples_C_27_74_95/admin_2019-05-20%2011-56-24_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_C_27_74_95 <- read.csv(text = getURL(Samples_C_27_74_95url), sep = ";")

Samples_40_41_43_46url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/20_05_2019_Samples_C_27_74_95/admin_2019-05-20%2011-56-24_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_40_41_43_46 <- read.csv(text = getURL(Samples_40_41_43_46url), sep = ";")