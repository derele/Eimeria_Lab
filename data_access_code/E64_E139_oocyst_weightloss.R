# example extract E139 and E64 (E. ferrisi) weightloss + oocyst shedding data from the Heitlinger group repositories
# looking for primary infections only, separate by experiment and date for sure

# load in necessary libraries
library(httr)
library(tidyverse)
library(Rmisc) #masks quite few dplyr objects use dplyr:: where necessary 
library(RCurl)

# load in data from github and extract E64 and E139 primary infections
# have to scavenge manually atm
# use as many "coplete" (or likewise labeled) tables as possible
# Capital letter denotes Passaging or Experiment (P,E), number denotes order and small letter denotes primary
# or challenge infection. In this case we look at "a" (primary).

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

