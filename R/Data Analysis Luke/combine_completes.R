library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)

P3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_COMPLETE.csv"))
E7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_complete.csv"))

P3 <- P3[,order(colnames(P3))]
P3$AVG <- NULL
P3$batch <- NULL
P3$X <- NULL
P3$oocyst_1 <- NULL 
P3$oocyst_2 <- NULL 
P3$oocyst_3 <- NULL 
P3$oocyst_4 <- NULL 
P3$N.oocyst <- NULL

E7 <- E7[,order(colnames(E7))]
E7$average <- NULL
E7$comment <- NULL
E7$X <- NULL
E7$oocyst_1 <- NULL 
E7$oocyst_2 <- NULL 
E7$oocyst_3 <- NULL 
E7$oocyst_4 <- NULL 
E7$totalOocysts <- NULL
E7$volume_PBS_mL <- NULL

colnames(E7)[1] <- "Eim_MC"
colnames(E7)[7] <- "faeces_weight"
colnames(P3)[13] <- "labels"
P3$HybridStatus <- NA
P3$Strain <- "SWISS"
colnames(P3)[21] <- "Wchange"
P3$totalOocysts <- NULL


complete <- rbind(P3, E7)

write.csv(complete, "./Eimeria_Lab/data/3_recordingTables/E7andP3_complete.csv")
