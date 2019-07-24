# RT-qPCR_clean data processing for E1_012017, additional data from CEC
library(RCurl)
library(httr)
library(Rmisc)
library(ggplot2)
library(dplyr)
library(compare)
library(viridis)
library(RColorBrewer)

#load in data from GitHub
RTqPCRurl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/b4efd8df335199ff9037634c5ba1d909a7d58baa/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv"
RTqPCR <- read.csv(text = getURL(RTqPCRurl))
# manual loading: 
RTqPCR <- read.csv(file = "./Eimeria_Lab/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv", sep = ";")
#change colnames to match standard
names(RTqPCR)[names(RTqPCR) == "Sample"] <- "EH_ID"
#average duplicates and add SD
AVG <- RTqPCR %>% group_by(Target, EH_ID) %>% 
  summarize(SD = sd(Cq.Mean),
            Cq.Mean = mean(Cq.Mean))

#add mouse data
InfectionURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E1_012017_Eim_Experiment_Table_raw_NMRI.csv"
Infection <- read.csv(text = getURL(InfectionURL))
#rename columns and merge
names(Infection)[names(Infection) == "mouseID"] <- "EH_ID"
names(Infection)[names(Infection) == "inf.strain"] <- "InfectionStrain"
AVG_INF <- merge(AVG, Infection, by = "EH_ID", all = TRUE)
#missing 80 values, compare
comparison <- compare(AVG,AVG_INF,allowAll=TRUE)
comparison$tM
#write out AVG_inf to correct based on lab book notes
write.csv(AVG_INF, file = "./AVG_INF.csv")

#subset for graphing 
Inf.groups <- subset(AVG_INF, 
                InfectionStrain%in%"Eflab"|
                InfectionStrain%in%"Efwild"|
                InfectionStrain%in%"EI70"|
                InfectionStrain%in%"EI64"|
                InfectionStrain%in%"Uninf")

ggplot(AVG_INF, aes(x = InfectionStrain, y = Cq.Mean, color = Target)) +
  geom_boxplot() +
  scale_color_viridis(option = "A")
