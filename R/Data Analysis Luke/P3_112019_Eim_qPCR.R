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
qPCR1 <- qPCR1[-c(85:105), ]
qPCR1.long <- dplyr::select(qPCR1, Name, Ct.Mean.SYBR, Target.SYBR)
qPCR1.long <- distinct(qPCR1.long)

qPCR1.long <- qPCR1.long %>% 
  dcast(Name ~ Target.SYBR, value.var = "Ct.Mean.SYBR", fill = 0) %>% 
  mutate(delta = mouse - eimeria) %>% 
  dplyr::select(Name,delta)
names(qPCR1.long)[names(qPCR1.long) == "Name"] <- "EH_ID"
qPCR1.long$dpi <- 8

#### write out (will become merge later)
write.csv(qPCR1.long, "./Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR1_clean.csv")
