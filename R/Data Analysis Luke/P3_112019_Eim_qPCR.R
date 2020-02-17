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
P3 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR1.CSV"
P3 <- read.csv(text = getURL(P3))
P3.long <- dplyr::select(P3, Name, Ct.Mean.SYBR, Target.SYBR)
P3.long <- distinct(P3.long)

P3.long <- P3.long %>% 
  dcast(Name ~ Target.SYBR, value.var = "Ct.Mean.SYBR", fill = 0) %>% 
  mutate(delta = eimeria - mouse) %>% 
  dplyr::select(Name,delta)
names(P3.long)[names(P3.long) == "Name"] <- "EH_ID"

ggplot(P3.long, aes(x = delta, y = EH_ID)) + 
  geom_point()
