library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)

P3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_COMPLETE.csv"))
E7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_COMPLETE.csv"))
E6 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E6_062018_Eim_WandO_complete.csv"))


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
E7$Caecum <- NULL

E6$X <- NULL
E6 <- E6[,order(colnames(E6))]
E6$Batch <- NULL
E6$Eim_MC <- NA
E6$CXCR3 <- NA
E6$IL.12 <- NA
E6$IRG6 <- NA
E6$delta <- NA
colnames(E6)[5] <- "faeces_weight"
E6$IFNy_FEC <- NA #temporary
E6$IFNy_CEWE <- NA
E6$InfectionStrain <- NULL
E6$oocyst_mean <- NULL
E6$infHistory  <- NA

colnames(E7)[7] <- "faeces_weight"
colnames(P3)[13] <- "labels"
P3$HybridStatus <- NA
P3$Strain <- "SWISS"
P3$totalOocysts <- NULL
colnames(P3)[18] <- "Wchange"
P3$EXP <- "P3"
E7$EXP <- "E7"
E6$weightloss <- NULL
E7$primary <- NA


P3 <- P3[,order(colnames(P3))]
E7 <- E7[,order(colnames(E7))]


complete <- rbind(P3, E7)
complete <- rbind(complete, E6)
complete <- distinct(complete)

# let's see if the NAs in primary and challenge worked
# complete %>% transform(currentInf=ifelse(grepl("(P3|E7)a", labels), 
#                                          as.character(primary), 
#                                          as.character(challenge))) -> 
#   complete
# 
# complete %>% transform(isChallenge=ifelse(grepl("(P3|E7)a", labels), 
#                                          "primary", 
#                                          as.character(challenge))) -> 
#   complete

write.csv(complete, "./Eimeria_Lab/data/3_recordingTables/E7_P3_E6_complete.csv")

write.csv(complete, "~/Documents/Eimeria_Lab/data/3_recordingTables/E7_P3_E6_complete.csv")
