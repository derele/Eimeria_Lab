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
  facet_wrap("primary") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 primary infection weighloss by parasite strain")

ggplot(P3b_record, aes(x = dpi, y = wloss)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("challenge") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 challenge infection weighloss by parasite strain")

ggplot(P3b_record, aes(x = dpi, y = wloss)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("infHistory") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 weightloss infhist")
# add oocyst data
oocysts <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_oocysts.csv"
oocysts <- read.csv(text = getURL(oocysts))
oocysts$X <- NULL 

# make overall table with labels
P3_record <- rbind(P3a_record, P3b_record)
P3_record_full <- merge(P3_record, oocysts, by = "labels")

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
P3_record_full$N.oocyst <- (P3_record_full$AVG * 10^4)/2

# create just primary and challenge oocyst and weightloss graphs
P3a <- merge(P3a_record, oocysts)

P3a$N.oocyst <- (P3a$AVG * 10^4)/2
P3a$OPG <- P3a$N.oocyst / P3a$faeces_weight

ggplot(P3a, aes(x = dpi, y = N.oocyst)) +
  geom_point() +
  facet_wrap("primary") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  ggtitle("P3 primary oocyst shedding by parasite strain")


P3b <- merge(P3b_record, oocysts)
P3b <- merge(P3b_record, oocysts)

P3b$N.oocyst <- (P3b$AVG * 10^4)/2
P3b$OPG <- P3b$N.oocyst / P3b$faeces_weight

ggplot(P3b, aes(x = dpi, y = N.oocyst, color = challenge)) +
  geom_point(size = 3) +
  facet_wrap("EH_ID") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  ggtitle("P3 challenge oocyst shedding by parasite strain")

ggplot(P3b, aes(x = dpi, y = N.oocyst)) +
  geom_point() +
  geom_line() +
  facet_wrap("infHistory") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 challenge oocyst shedding by parasite strain")


#let's see what we made


# makes no sense to graph this 
# ggplot(P3_record_full, aes(x = dpi, y = AVG, color = challenge, shape = batch)) +
#   geom_point(size = 3) +
#   facet_wrap("primary") +
#   theme(axis.text=element_text(size=12, face = "bold"), 
#         axis.title=element_text(size=14,face="bold"),
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold"))+
#   ggtitle("P3 oocyst shedding by parasite strain AVG")
P3_record_full$N.oocyst <- (P3_record_full$AVG * 10^4)/2

ggplot(P3_record_full, aes(x = dpi, y = OPG, color = challenge)) +
  geom_point() +
  facet_wrap("primary") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 oocyst shedding by parasite strain OPG")

ggplot(P3_record_full, aes(x = dpi, y = N.oocyst, color = batch)) +
  geom_point() +
  geom_line() +
  facet_wrap("infHistory") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 batch shedding")


ggplot(P3_record_full, aes(x = primary, y = N.oocyst, color = Eim_sp))+ 
  geom_boxplot() +
  geom_point() + 
  facet_wrap("challenge")


# add qPCRs
P3_qPCRs <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_qPCRs_complete.csv"
P3_qPCRs <- read.csv(text = getURL(P3_qPCRs))

P3_qPCRs <- merge(P3_qPCRs, P3_design, by = "EH_ID")

ggplot(P3_qPCRs, aes(x = delta, y = NE, color = infHistory)) +
  geom_point(size = 3) +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 gene expression vs delta")

ggplot(P3_qPCRs, aes(x = delta, y = NE, color = challenge)) +
  geom_point(size = 3) +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 gene expression vs delta")

P3b_q <- select(P3b, EH_ID, OPG, N.oocyst, wloss, dpi)
P3b_q <- merge(P3b_q, P3_qPCRs, by = "EH_ID")

ggplot(P3b_q, aes(x = N.oocyst, y = NE, color = challenge)) +
  geom_point(size = 3) +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("P3 gene expression vs delta")
