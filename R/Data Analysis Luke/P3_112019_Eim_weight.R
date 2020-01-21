# P3 (16 SWISS mice) infection experiment (Eflab, E64, E139, E88, UNI)
library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
library(ggplot2)
library(tidyr)


# load in raw data (big mistake, at least one sample has duplicate label across batches. 
# (Physical samples are separated howevers so only needs fixing here)
P3a_record <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3a_112019_Eim_Record.csv"
P3a_record <- read.csv(text = getURL(P3a_record))
P3a_record$batch <- "a"
P3a_record$day_change <- NULL
P3a_record$labels <- sub("^", "1", P3a_record$labels)

P3b_record <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3b_112019_Eim_Record.csv"
P3b_record <- read.csv(text = getURL(P3b_record))
P3b_record$X <- NULL
P3b_record$batch <- "b"
P3b_record$day_change <- NULL
P3b_record$labels <- sub("^", "2", P3b_record$labels)

# load in design
P3_design <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/P3_112019_Eim_design.csv"
P3_design <- read.csv(text = getURL(P3_design))
P3_design <- select(P3_design, EH_ID, primary, challenge)
P3_design <- unique(P3_design)

# mak infection history
P3_design$infHistory <- paste(P3_design$primary, P3_design$challenge, sep =  ":") 
# merge to acquire infection history (also fix any num/factor class discrepancies and data frame columns)
P3a_record <- merge(P3a_record, P3_design, by = "EH_ID", all = T)
P3b_record <- merge(P3b_record, P3_design, by = "EH_ID", all = T)

# graph to see overall state
ggplot(P3a_record, aes(x = dpi, y = wloss)) +
  geom_point() +
  geom_smooth() + 
  facet_wrap("primary")

ggplot(P3b_record, aes(x = dpi, y = wloss)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("challenge")

ggplot(P3b_record, aes(x = dpi, y = wloss)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("primary")

# add oocyst data
oocysts <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_oocysts.csv"
oocysts <- read.csv(text = getURL(oocysts))
oocysts$X <- NULL 

# make overall table with labels
P3_record <- rbind(P3a_record, P3b_record)
P3_record_full <- merge(P3_record, oocysts, by = "labels", all = T)

# add Eim species column
P3_record_full <- P3_record_full %>% mutate(Eim_sp = case_when(
  infHistory == "Eflab:E64" ~ "Efal:Efer",
  infHistory == "E64:E88" ~ "Efer:Efal",
  infHistory == "E139:E88" ~ "Efer:Efal",
  infHistory == "E88:E64" ~ "Efal:Efer",
  
  infHistory == "E88:UNI" ~ "Efal",
  infHistory == "Eflab:UNI" ~ "Efal",
  infHistory == "UNI:E88" ~ "Efal",
  infHistory == "E88:E88" ~ "Efal",
  infHistory == "Eflab:E88" ~ "Efal",
  
  infHistory == "UNI:E64" ~ "Efer",
  infHistory == "E64:E64" ~ "Efer",
  infHistory == "E139:UNI" ~ "Efer",
  infHistory == "E139:E64" ~ "Efer",
  infHistory == "E64:UNI" ~ "Efer",
  
  infHistory == "UNI:UNI" ~ "uni"))

# calculate OPG
P3_record_full$OPG <- P3_record_full$AVG / P3_record_full$faeces_weight

#let's see what we made
ggplot(P3_record_full, aes(x = dpi, y = AVG, color = challenge)) +
  geom_point() +
  facet_wrap("primary")

ggplot(P3_record_full, aes(x = dpi, y = OPG, color = challenge)) +
  geom_point() +
  facet_wrap("primary")

ggplot(P3_record_full, aes(x = dpi, y = OPG, color = infHistory)) +
  geom_point() +
  geom_line() +
  facet_wrap("batch")

ggplot(P3_record_full, aes(x = primary, y = OPG, color = Eim_sp))+ 
  geom_boxplot() +
  geom_point() + 
  facet_wrap("challenge")


# graph like Anna's for comparison
