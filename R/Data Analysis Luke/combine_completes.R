library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)

P3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_COMPLETE.csv"))
E7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_complete.csv"))
E5 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E6_062018_Eim_WandO_complete.csv"))


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

E5$X <- NULL
E5 <- E5[,order(colnames(E5))]
E5$Batch <- NULL
E5$Eim_MC <- NA
E5$CXCR3 <- NA
E5$IL.12 <- NA
E5$IRG6 <- NA
E5$delta <- NA
colnames(E5)[5] <- "faeces_weight"
E5$IFNy_FEC <- NA #temporary
E5$IFNy_CEWE <- NA
E5$InfectionStrain <- NULL
E5$oocyst_mean <- NULL
E5$infHistory  <- NA


colnames(E7)[1] <- "Eim_MC"
colnames(E7)[7] <- "faeces_weight"
colnames(P3)[13] <- "labels"
P3$HybridStatus <- NA
P3$Strain <- "SWISS"
P3$totalOocysts <- NULL
colnames(P3)[18] <- "Wchange"
P3$EXP <- "P3"
E7$EXP <- "E7"
E5$weightloss <- NULL
E7$primary <- NA

complete <- rbind(P3, E7)
complete <- rbind(complete, E5)


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


