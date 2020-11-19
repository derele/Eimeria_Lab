# example extract E139 and E64 (E. ferrisi) weightloss + oocyst shedding data from the Heitlinger group repositories
# looking for primary infections only, separate by experiment and date for sure

# load in necessary libraries
library(httr)
library(tidyverse)
library(Rmisc) #masks quite few dplyr objects use dplyr:: where necessary 
library(RCurl)

# load in data from github and extract E64 and E139 primary infections
# have to scavenge manually atm
# consult this table for overview of strains used in experiments and passaging:
# https://docs.google.com/spreadsheets/d/1ejzceB7z5nZh0V_kQUBDVbW4AQEWwRG2FO5MV-oOX1g/edit?usp=sharing
# use as many "coplete" (or likewise labeled) tables as possible
# Capital letter denotes Passaging or Experiment (P,E), number denotes order and small letter denotes primary
# or challenge infection. In this case we look at "a" (primary).
# If batch is present in columns, those serve the same purpose

# From the table we can see that E3, E4, E8, P1, P2 and P4 contain E139
# E64 is in E1, E2, E3, E4, E5, E7, E10, P1, P2, P3 and P4 

############ load in  and process (down to EH_ID, labels, OPG, dpi, weight loss and Eimeria strain) P3 and P4
CLS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_Weight%26Oocyst_complete.csv"))
# split take out only batch "a" for primary + E139 and E64
CLS <- subset(CLS, CLS$batch == "a")
CLSE64 <- subset(CLS, CLS$primary == "E64")
CLSE139 <- subset(CLS, CLS$primary == "E139")
CLS <- rbind(CLSE139, CLSE64)
# keep only necessary columns for this task
CLS <- select(CLS, EH_ID, labels, totalOocysts, AVG, dpi, OPG, weight, wloss, primary)

############ load in P2
P2 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P2_052019_Eim_record.csv"))

############ load in P1
P1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P1_112018_Eim_oocyst_only.csv"))

############ load in E2
read.csv(text = getURL(""))

############ load in E3

