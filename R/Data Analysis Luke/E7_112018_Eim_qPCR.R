# E7 qPCR for infection intensity

library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(naniar)

# load in raw data MCs
E7EimMC <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7%26P3_Eim_MCs.csv"
E7EimMC <- read.csv(text = getURL(E7EimMC), sep = ";")
E7EimMC$X <- NULL
E7EimMC$X.1 <- NULL

# add intensity 
E7_inf <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_Anna_qPCR_DNA_ct_Zusammenfassung.csv"
E7_inf <- read.csv(text = getURL(E7_inf))
E7_inf$Ct.SYBR <- NULL
E7_inf$Pos <- NULL
E7_inf$Amount.SYBR..Copies. <- NULL
E7_inf$Amount.Mean.SYBR <- NULL
E7_inf$Amount.Dev..SYBR <- NULL
E7_inf <- distinct(E7_inf)
E7_inf <- E7_inf %>% 
  dcast(Name ~ Target.SYBR, value.var = "Ct.Mean.SYBR", fill = 0) %>% 
  mutate(delta = mouse - eimeria) %>% 
  dplyr::select(Name,delta)
names(E7_inf)[names(E7_inf) == "Name"] <- "EH_ID"
E7_inf <- merge(E7_inf, E7EimMC)
# write
write.csv(E7_inf, "../Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_qPCR.csv")
