# RT-qPCR_clean data processing for E1_012017, additional data from CEC
library(RCurl)
library(httr)
library(Rmisc)
library(ggplot2)
library(dplyr)
library(compare)
library(viridis)
library(RColorBrewer)
library(ggsci)

#load in data from GitHub, doesn't work atm (more columns than column names)
RTqPCRurl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/b4efd8df335199ff9037634c5ba1d909a7d58baa/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv"
RTqPCR <- read.csv(text = getURL(RTqPCRurl))
# manual loading: Win
RTqPCR <- read.csv(file = "./Eimeria_Lab/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv", sep = ";")
# manual loading Deb laptop
RTqPCR <- read.csv(file = "./Documents/Eimeria_Lab/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv", sep = ";")
#change colnames to match standard
names(RTqPCR)[names(RTqPCR) == "Sample"] <- "EH_ID"
#average duplicates and add SD
AVG <- RTqPCR %>% group_by(Target, EH_ID) %>% 
  summarize(SD = sd(Cq.Mean),
            Cq.Mean = mean(Cq.Mean))

#add mouse data
InfectionURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E1_012017_Eim_Experiment_Table_raw_NMRI.csv"
Infection.design <- read.csv(text = getURL(InfectionURL))
#rename columns and merge
names(Infection)[names(Infection) == "mouseID"] <- "EH_ID"
names(Infection)[names(Infection) == "inf.strain"] <- "InfectionStrain"
AVG_INF <- merge(AVG, Infection, by = "EH_ID", all = TRUE)
#add Uninf to LM00C (otherwise gets removed by subset)
AVG_INF[715:724,"InfectionStrain"] <- "Uninf"
#remove rows for missing samples
AVG_INF <- AVG_INF[ !(AVG_INF$EH_ID %in% c("LM0021", "LM0033", "LM0035", "LM0052")), ]
#subset by infection 
Infection <- subset(AVG_INF, InfectionStrain%in%"Eflab"|InfectionStrain%in%"Efwild"|InfectionStrain%in%"EI70"|
                       InfectionStrain%in%"EI64"|InfectionStrain%in%"Uninf")
# subset by gene
Genes <- subset(AVG_INF, Target%in%"CDC42"|Target%in%"CXCL9"|Target%in%"IFN-g"|Target%in%"IL-10"|Target%in%"IL-12"|
                Target%in%"IL-6"|Target%in%"Ppia"|Target%in%"Ppib"|Target%in%"STAT6"|Target%in%"TGF-b")
# subset by dpi (dpi7 says "dip")
Dpi <- subset(AVG_INF,dpi.diss%in%"3dpi"|dpi.diss%in%"5dpi"|dpi.diss%in%"7dip"|dpi.diss%in%"9dpi"|dpi.diss%in%"11dpi")
#graph by subsets
#Infection
ggplot(Infection, aes(x = InfectionStrain, y = Cq.Mean, color =InfectionStrain)) +
  geom_boxplot()
#Gene
ggplot(Genes, aes(x = Target, y = Cq.Mean, color = Target)) +
  geom_boxplot()
#dpi
ggplot(Dpi, aes(x = dpi.diss, y = Cq.Mean, color = dpi.diss)) +
  geom_boxplot()
