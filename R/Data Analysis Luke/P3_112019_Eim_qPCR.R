# qPCRs for infection intensity
library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(naniar)
library(data.table)

######### add infection intensity
qPCR1 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR4.CSV"
qPCR1 <- read.csv(text = getURL(qPCR1))
# remove LM_340 because, added in qPCR3
qPCR1 <-qPCR1[!(qPCR1$Name=="LM_0340"),]

qPCR1.long <- dplyr::select(qPCR1, Name, Ct.Mean.SYBR, Target.SYBR)
qPCR1.long <- distinct(qPCR1.long)

# add qPCR 3 and use LM_340 from there
qPCR2 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR3.CSV"
qPCR2 <- read.csv(text = getURL(qPCR2))
qPCR2.long <- dplyr::select(qPCR2, Name, Ct.Mean.SYBR, Target.SYBR)
qPCR2.long <- distinct(qPCR2.long)

# make into deltas (mouse-eimeria)
qPCR1.long <- qPCR1.long %>% 
  dcast(Name ~ Target.SYBR, value.var = "Ct.Mean.SYBR", fill = 0) %>% 
  mutate(delta = mouse - eimeria) %>% 
  dplyr::select(Name,delta)
names(qPCR1.long)[names(qPCR1.long) == "Name"] <- "EH_ID"

qPCR2.long <- qPCR2.long %>% 
  dcast(Name ~ Target.SYBR, value.var = "Ct.Mean.SYBR", fill = 0) %>% 
  mutate(delta = mouse - eimeria) %>% 
  dplyr::select(Name,delta)
names(qPCR2.long)[names(qPCR2.long) == "Name"] <- "EH_ID"
# check MCs one more time but pick qPCR1 predominantly and add 332, 333 and 340 from qPCR2 for now)
qPCR2.long <- subset(qPCR2.long, EH_ID == c("LM_0332", "LM_0333", "LM_0340"))
# merge and write out
P3_qPCR <- rbind(qPCR1.long, qPCR2.long)
write.csv(P3_qPCR, "./Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR_complete.csv")

