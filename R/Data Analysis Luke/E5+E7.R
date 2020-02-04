# E5 + E7 analysis script

library(Rmisc)
library(httr)
library(RCurl)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(data.table)

E7 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_complete.csv"
E7 <- read.csv(text = getURL(E7))
E7$EH_ID <- sub("0", "_0", E7$EH_ID)
E7$X <- NULL


E5 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E6_062018_Eim_complete(Wchange).csv"
E5 <- read.csv(text = getURL(E5))
E5 <- select(E5, EH_ID, labels, Expe, dilution_ml, OPG, dpi, weight, weight_dpi0, Wchange, fecweight, Eimeria, Strain, HybridStatus)
# rename columns to match (EH_ID, labels, Expe, dilution_ml, OPG, dpi, weight, weight_dpi0, Wchange, fecweight, Eimeria, 
# Strain, Hybridstatus) + keep only mice that went into reinfection
names(E5)[names(E5) == "Eimeria"] <- "primary"
names(E5)[names(E5) == "dilution_ml"] <- "volume_PBS_mL"
E5$labels <- sub("^", "0", E5$labels)
E5$OPG <- varhandle::unfactor(E5$OPG)
E5$OPG <- sub(",", ".", E5$OPG)
E5$OPG <- as.numeric(as.character(E5$OPG))

############ graph out OPG

ggplot(E5, aes(x = dpi, y = OPG)) +
  geom_point() +
  facet_wrap("primary", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Primary infection shedding")

ggplot(E7, aes(x = dpi, y = OPG)) +
  geom_point() +
  facet_wrap("challenge", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Challenge infection shedding")

ggplot(E7, aes(x = dpi, y = OPG)) +
  geom_point() +
  facet_wrap("infHistory", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Challenge infection shedding vs. infection history")

##### include gene expression and infection intensity
E7_exp <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT_and_qPCR_complete.csv"
E7_exp <- read.csv(text = getURL(E7_exp))
E7_exp$X <- NULL
E7 <- merge(E7, E7_exp, by = "EH_ID")

###### look at ferrisi - falciformis shedding
# on primary
ggplot(E5, aes(x = dpi, y = OPG)) +
  geom_point() +
  facet_wrap("primary", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Shedding of oocysts during primary infection (parasite)")
# on challenge
ggplot(E7, aes(x = dpi, y = OPG)) +
  geom_point() +
  facet_wrap("challenge", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Shedding of oocysts during challenge infection (parasite)")

######## look at shedding vs. intensity (only challenge)

ggplot(E7, aes(x = delta, y = OPG, color = challenge)) +
  geom_point() +
  facet_wrap("dpi", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Shedding of oocysts during challenge infection (intensity)")

