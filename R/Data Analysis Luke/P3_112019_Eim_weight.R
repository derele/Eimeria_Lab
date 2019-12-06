# P3 (16 SWISS mice) infection experiment (Eflab, E64, E139, E88, UNI)
library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
library(ggplot2)
library(tidyr)


# load in raw data
P3a_record <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3a_112019_Eim_Record.csv"
P3a_record <- read.csv(text = getURL(P3a_record))

P3b_record <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3b_112019_Eim_Record.csv"
P3b_record <- read.csv(text = getURL(P3b_record))
P3b_record$X <- NULL

# load in design
P3_design <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/P3_112019_Eim_design.csv"
P3_design <- read.csv(text = getURL(P3_design))
P3_design <- select(P3_design, EH_ID, primary, challenge)
P3_design <- unique(P3_design)

# mak infection history
P3_design$infHistory <- paste(P3_design$primary, P3_design$challenge, sep =  ":") 
# merge to acquire infection history (also fix any num/factor class discrepancies and data frame columns)
P3a_record <- merge(P3a_record, P3_design, by = "EH_ID")
P3b_record <- merge(P3b_record, P3_design, by = "EH_ID")
P3b_record$wloss <- as.numeric(as.character(P3b_record$wloss))

# graph to see overall state
ggplot(P3a_record, aes(x = dpi, y = wloss)) +
  geom_point() +
  geom_smooth() + 
  facet_wrap("primary")

ggplot(P3b_record, aes(x = dpi, y = wloss, color = infHistory)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("infHistory")


