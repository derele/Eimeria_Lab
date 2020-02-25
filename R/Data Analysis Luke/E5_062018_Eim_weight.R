# process E6 weight and oocysts
library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)

E5_weight <- read.csv(text = getURL("https://github.com/derele/Eimeria_Lab/blob/master/data/3_recordingTables/E6_062018_Eim_full_RECORDweight.csv"))
E5_oocyst <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E6_062018_Eim_full_RECORDoocysts.csv"))
