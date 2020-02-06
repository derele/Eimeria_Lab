# P3 RT-qPCR and infection intensity qPCR analysis

library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(naniar)

# load in raw

RT1 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR1.CSV"
RT1 <- read.csv(text = getURL(RT1))

