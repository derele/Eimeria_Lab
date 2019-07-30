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
library(tidyverse)
library(reshape2)

#------------------------------ add and process RTqPCR data from GitHub-----------------------------------
RTqPCRurl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv"
RTqPCR <- read.csv(text = getURL(RTqPCRurl), sep = ",")
#change colnames and misnamed rows to match standard
names(RTqPCR)[names(RTqPCR) == "Sample"] <- "EH_ID"
RTqPCR[RTqPCR=="IFN-y"] <- "IFN-g"
# just averages and add SD
RTqPCR <- data.frame(RTqPCR %>% group_by(Target, EH_ID) %>% 
  summarize(SD = sd(Cq.Mean),
            Cq.Mean = mean(Cq.Mean)))
#convert columns to char + remove multiple classses
RTqPCRcharacters <- sapply(RTqPCR, is.factor)
RTqPCR[RTqPCRcharacters] <- lapply(RTqPCR[RTqPCRcharacters], as.character)
RTqPCR = as.data.frame(RTqPCR)


#------------------add and process design table---------------------------------------------------------
InfectionURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E1_012017_Eim_Experiment_Table_raw_NMRI.csv"
Infection.design <- read.csv(text = getURL(InfectionURL))
#rename columns and merge
names(Infection.design)[names(Infection.design) == "mouseID"] <- "EH_ID"
names(Infection.design)[names(Infection.design) == "inf.strain"] <- "InfectionStrain"

#--------------------------all.data--------------------------------------------
## wide dateset for merging in overall table (don't forget to subtract standards, see Emanuel's script)
all.data <- merge(RTqPCR, Infection.design, by = "EH_ID", all = TRUE)
#big NA introduction, replacing with 0s
all.data$SD[is.na(all.data$SD)] <- 0
#add Uninf to LM00C (otherwise gets removed)
all.data[715:724,"InfectionStrain"] <- "Uninf"
#add dpi0 to LM00C (won't work as factor)
all.data.characters <- sapply(all.data, is.factor)
all.data[all.data.characters] <- lapply(all.data[all.data.characters], as.character)
all.data = as.data.frame(all.data)
all.data[715:724,"dpi.diss"] <- "dpi0"
#remove NAs (missing samples from Infection design)
all.data <- na.omit(all.data)
#convert numeric to alow continuous scale
all.data$dpi.diss <- as.numeric(gsub("dpi|dip", "", all.data$dpi.diss))
# all.data[715:724, "dpi.diss"] <- as.factor(x = "dpi0")
#remove rows for missing samples
all.data <- all.data[ !(all.data$EH_ID %in% c("LM0021", "LM0033", "LM0035", "LM0052")), ]
#make dpi.diss numeric and a separate column
all.data$dpi <- as.numeric(gsub("dpi|dip", "", all.data$dpi.diss))

ggplot(subset(all.data, nchar(all.data$Target)>2), aes(dpi.diss, Cq.Mean, color=InfectionStrain)) +
  geom_jitter(width=0.2) +
  geom_smooth(se=FALSE) +
  scale_x_continuous(breaks=c(3, 5, 7, 9, 11),
                     labels=c("3dpi", "5dpi", "7dip", "9dpi", "11dpi")) +
  facet_wrap(~Target, scales="free_y", nrow=2)+
  scale_colour_brewer("infection\nisolate", palette = "Dark2") +
  scale_y_continuous("normalized mRNA expression")+
  theme_bw()
dev.off()

#----------------add and process infection intensity expression-----------------------------------#needs rewrite
RtissueURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Eimeria-NMRI_Relative%20quantification_clean.csv"
Rtissue <- read.csv(text = getURL(RtissueURL))
## only the means
RtMeans <- Rtissue[seq(1, nrow(Rtissue), by=2), c("Sample", "Cq.Mean", "Cq.Mean.1")]
#naming and upper case
names(RtMeans) <- c("EH_ID", "Mouse_gDNA", "Eimeria_mDNA")
RtMeans$EH_ID <- toupper(RtMeans$EH_ID)
## LM0065 was measured twice with the same outcome
RtMeans <- RtMeans[!duplicated(RtMeans$EH_ID),]
RtMeans <- merge(RtMeans, Infection.design, all=TRUE)
RtMeans$dpi.diss <- gsub("dip$", "dpi", RtMeans$dpi.diss)
RtMeans$dpi_count<- as.numeric(gsub("dpi", "", RtMeans$dpi.diss))

#--------------------------add and process oocyst data---------------------------------------------------
oocystsURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Clean_oocyst_data.csv"
oocysts <- read.csv(text = getURL(oocystsURL))

#--------------------------add and process mouse weight data---------------------------------------------
weightURL <- "https://raw.githubusercontent.com/derele/Jan2017Exp/master/Weight_Expe_Jan_2017.csv"
weight <- read.csv(text = getURL(weightURL))




##subset by infection 
Infection <- subset(all.data, InfectionStrain%in%"Eflab"|InfectionStrain%in%"Efwild"|InfectionStrain%in%"EI70"|
                      InfectionStrain%in%"EI64"|InfectionStrain%in%"Uninf")
# subset by gene
Genes <- subset(all.data, Target%in%"CDC42"|Target%in%"CXCL9"|Target%in%"IFN-g"|Target%in%"IL-10"|Target%in%"IL-12"|
                  Target%in%"IL-6"|Target%in%"Ppia"|Target%in%"Ppib"|Target%in%"STAT6"|Target%in%"TGF-b")
# subset by dpi (dpi7 says "dip")
Dpi <- subset(all.data,dpi.diss%in%"3dpi"|dpi.diss%in%"5dpi"|dpi.diss%in%"7dip"|dpi.diss%in%"9dpi"|dpi.diss%in%"11dpi")
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