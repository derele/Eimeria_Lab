#E7_112018_Eim_RT-qPCR
library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)

RT <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCR1.csv"
RT <- read.csv(text = getURL(RT))
