## Experiment 007 Reinfection experiment Nov-Dec 2018
## Eimeria isolates : E64 (E. ferrisi), E88 (E. falciformis)

###########################################################################

## Infection graphs

###########################################################################

pathToData <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp007/"
class(pathToData)

# Import our data
Expe007a_oocysts <- read.csv(
  paste(pathToData, "Exp007a/Exp_007a_oocyst_counts.csv", sep = ""))
Expe007a_weight <- read.csv(
  paste(pathToData, "Exp007a/Exp_007a_feces.csv", sep = ""))
Expe007b_oocysts <- read.csv(
  paste(pathToData, "Exp007b/Exp_007b_oocyst_counts.csv", sep = ""))
Expe007b_weight <- read.csv(
  paste(pathToData, "Exp007b/Exp_007b_feces.csv", sep = ""))

pathToDesign <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/"
Expe007a_design <- read.csv(
  paste(pathToDesign, "Exp007a_design.csv", sep = ""))
Expe007b_design <- read.csv(
  paste(pathToDesign, "Exp007b_design.csv", sep = ""))

# First, merge the tables
Expe007a <- merge(Expe007a_oocysts, Expe007a_weight)
Expe007a <- merge(Expe007a, Expe007a_design)

Expe007b <- merge(Expe007b_oocysts, Expe007b_weight)
Expe007b <- merge(Expe007b, Expe007b_design)

Expe007 <- merge(Expe007a, Expe007b, all = TRUE)

InfectionHistory <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp007/Exp_007_infection_history_complete.csv")

Expe007 <- merge(Expe007, InfectionHistory, all=TRUE)

# To export
# write.csv(Expe007, 
#           file = "../../../../home/alice/Schreibtisch/GIT/Eimeria_Lab/data/3_recordingTables/Exp007/Expe007FullTable.csv", )

# NB. Let's not consider which parent is which, but make A_B mouse = B_A mouse
Expe007$Strain <- as.character(Expe007$Strain)
x <- strsplit(Expe007$Strain, "_")
y <- lapply(x, sort)
z <- unlist(lapply(y, FUN = function(x){paste(x, collapse="-")}))
Expe007$Strain <- z

############################################################################
# Then prepare the oocysts counts

Expe007$totalOocysts <- ((Expe007$oocyst_1 
                          + Expe007$oocyst_2 
                          + Expe007$oocyst_3 
                          + Expe007$oocyst_4) / 4) * 
  10000 * # because volume chamber
  Expe007$volume_PBS_mL

Expe007$OPG <- Expe007$totalOocysts / Expe007$fecweight 

###############################################################################

# Make plot cause we can
#install.packages("ggplot2")
library(ggplot2) # call the package
library(grid)
library(gridExtra)

############ Plot OPG/dpi ############

# explore
#table(Expe007$batch, Expe007$challenge)

#Plot OPG/day
ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = EH_ID)) +
  geom_point()+geom_line()+theme_bw() # pretty

#Plot OPG/day per strain
ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = Strain)) +
  geom_point()+geom_line()+theme_bw()

ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = Strain)) +
  geom_point()+geom_line()+theme_bw()+facet_grid(.~challenge)

#Plot OPG/dpi per hybrid status
ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = HybridStatus)) +
  geom_point()+geom_line()+theme_bw()

ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = HybridStatus)) +
  geom_point()+geom_line()+theme_bw()+facet_grid(.~HybridStatus)

ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = challenge)) +
  geom_point()+geom_line()+theme_bw()+facet_grid(.~HybridStatus)

ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = HybridStatus)) +
  geom_point()+geom_line()+theme_bw()+facet_grid(.~challenge)

#Plot OPG/day per Eimeria 
ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = challenge)) +
  geom_point()+geom_line()+theme_bw()

ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = challenge)) +
  geom_point()+geom_line()+theme_bw()+facet_grid(.~challenge)

ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = challenge)) +
  geom_point()+geom_line()+theme_bw()+facet_grid(.~Strain)

ggplot(Expe007, aes(x = dpi, y = OPG, group = EH_ID, col = HybridStatus)) +
  geom_point()+geom_line()+theme_bw()+facet_grid(.~Strain)


############ Plot weight/dpi ############

#Plot weight/dpi
ggplot(Expe007a, aes(x=dpi, y=weight, group=EH_ID, col= EH_ID)) +
  geom_point()+geom_line()+theme_bw()

#Plot weight/dpi per strain
ggplot(Expe007a, aes(x=dpi, y=weight, group=EH_ID, col= Strain)) +
  geom_point()+geom_line()+theme_bw()

#Plot weight/dpi per hybrid
ggplot(Expe007a, aes(x=dpi, y=weight, group=EH_ID, col= HybridStatus)) +
  geom_point()+geom_line()+theme_bw()

#Plot weight/dpi per Eimeria
ggplot(Expe007a, aes(x=dpi, y=weight, group=EH_ID, col= challenge)) +
  geom_point()+geom_line()+theme_bw()


############ Additional plots ############

#Plot weight change/dpi 
#this is relative weight change...
# ggplot(Expe007a, aes(x=dpi, y=Wchange, group=EH_ID, col= EH_ID)) + geom_point()+geom_line()+theme_bw()

# Make a plot with group = mean + 95%CI

#OPG mean per Strain
#StrainMean <- aggregate(Expe007[ , 10], list(Expe007$Strain), mean)


###########################################################################

## qPCR results

###########################################################################

#Say R where to search data
setwd("C:/Users/Anna/Documents/HU Berlin/AG HU/Bachelorarbeit/DNA qPCR")
#rm(list=ls())
getwd()

#Open and save data
qPCR <- read.table("qPCR_DNA_ct_Zusammenfassung.txt", sep="\t", header=T, dec=",")

###########################################################################

##### Clean table #####

library(dplyr)
library(tidyr)

#Delete columns that I don't need
#qPCR <- qPCR[-(1)]
#qPCR <- qPCR[-(4:7)] #delete everything between 4 and 7

qPCR <- qPCR[!names(qPCR) %in% c("Pos", "Ct.SYBR", "Amount.SYBR..Copies.", "Amount.Mean.SYBR", "Amount.Dev..SYBR")]

#Remove duplicate ct means per mouse (eimeria and mouse)
qPCR <- qPCR %>% distinct()


##### Calculate delta ct #####

#qPCR <- qPCR %>% spread(Target.SYBR)
#qPCR <- merge("Name","Target.SYBR", by = intersect("eimeria", "mouse"))

#Create new tables for Eimeria and mouse data
sumDataMouse <- qPCR[qPCR$Target.SYBR %in% "mouse",]
sumDataMouse$Ct.Mean.SYBR.mouse <- sumDataMouse$Ct.Mean.SYBR

sumDataEimeria <- qPCR[qPCR$Target.SYBR %in% "eimeria",]
sumDataEimeria$Ct.Mean.SYBR.eimeria <- sumDataEimeria$Ct.Mean.SYBR

#Merge tables and calculate delta ct
mergedData <- merge(sumDataEimeria, sumDataMouse, by = "Name")
mergedData$deltaCt_MminusE <- as.numeric(as.character(mergedData$Ct.Mean.SYBR.mouse)) -
  as.numeric(as.character(mergedData$Ct.Mean.SYBR.eimeria)) # DeltaCt MOUSE minus EIMERIA

#Remove unnecessary data from environment
rm(sumDataEimeria, sumDataMouse)

#Clean table
mergedData <- mergedData[!names(mergedData) %in% c("Ct.Mean.SYBR.x", "Ct.Dev..SYBR.x", "Target.SYBR.x", "Ct.Mean.SYBR.y", "Ct.Dev..SYBR.y", "Target.SYBR.y" )]

#Status pos/neg
mergedData$Status <- ifelse(mergedData$deltaCt_MminusE > -5 & mergedData$deltaCt_MminusE < 5, "positive", "negative")

# a little plot to be happy
library(ggplot2)
ggplot(mergedData, aes(x = Name, y = deltaCt_MminusE,
                       col = Status)) + geom_point() + theme_bw()


##### Merge Infection History with qPCR #####

library(dplyr)
library(tidyr)

#Rename Name to EH_ID
colnames(mergedData)[colnames(mergedData)=="Name"] <- "EH_ID"

#Clean table
InfectionHistory <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp007/Exp_007_infection_history_complete.csv")
InfectionHistory <- InfectionHistory %>% distinct()

#Merge + clean tables
qPCRInfectionHistory <- merge(mergedData, InfectionHistory, by = "EH_ID")
qPCRInfectionHistory <- qPCRInfectionHistory%>% select(labels,EH_ID,Ct.Mean.SYBR.eimeria, Ct.Mean.SYBR.mouse, deltaCt_MminusE, Status, primary, challenge, Ct.Dev..SYBR.x, Ct.Dev..SYBR.y)

#Plot because we can
library(ggplot2)
ggplot(qPCRInfectionHistory, aes(x = challenge, y = deltaCt_MminusE,
                       col = Status)) + geom_point() + theme_bw()

ggplot(qPCRInfectionHistory, aes(x = EH_ID, y = deltaCt_MminusE,
                       col = challenge)) + geom_point() + theme_bw()

#Plot challenge and primary infection comparison
ggplot(qPCRInfectionHistory, aes(x = EH_ID, y = deltaCt_MminusE,
                                 col = challenge)) + geom_point() + theme_bw()+facet_grid(.~primary)

ggplot(qPCRInfectionHistory, aes(x = challenge, y = deltaCt_MminusE,
                                 col = challenge)) + geom_point() + theme_bw()+facet_grid(.~primary)

###########################################################################

##### Contaminated mice (before Expe005) #####

Expe005 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp005_complete.csv", sep = ";")

#Open table with labels
setwd("C:/Users/Anna/Documents/HU Berlin/AG HU/Bachelorarbeit/PCR contaminated mice")
labelscontaminatedmice <- read.table("labels.txt", sep="\t", header=T, dec=",")

#Merge tables and just keep contaminated mice
contaminatedmice <- merge(labelscontaminatedmice, Expe005)

#Clean table, just keep relevant info
library(dplyr)
contaminatedmice <- contaminatedmice %>% select(labels, EH_ID, Expe, InfectionStrain, Strain, HybridStatus, contamination)

#Remove unnecessary data from environment
rm(Expe005, labelscontaminatedmice, InfectionHistory, qPCR, mergedData)

#Reorder columns
contaminatedmice <- contaminatedmice[c(1,2,3,4,7,5,6)]

#Change ' ' unknown contamination
levels(contaminatedmice$contamination)[levels(contaminatedmice$contamination)==''] <- 'u.c.'

#Merge + clean tables
InfectionTable <- merge(qPCRInfectionHistory, contaminatedmice, all=TRUE)
InfectionTable <- InfectionTable %>% select(EH_ID, deltaCt_MminusE, primary, challenge, contamination)
InfectionTable <- InfectionTable[-c(54:58), ]

#New column with Infection History (primary+challenge)
InfectionTable$InfectionHistory <- paste(InfectionTable$primary,InfectionTable$challenge)

#Plot
ggplot(InfectionTable, aes(x = EH_ID, y = deltaCt_MminusE,
                                col = contamination)) + geom_point() + theme_bw()+facet_grid(.~InfectionHistory)

ggplot(InfectionTable, aes(x = contamination, y = deltaCt_MminusE, col = contamination)) + geom_point() + theme_bw() + facet_grid(.~InfectionHistory)


##### Group by infection history #####

#New column infection history+contamination
InfectionTable$InfectionHistoryTotal <- paste(InfectionTable$contamination,InfectionTable$InfectionHistory)

#Plot
ggplot(InfectionTable, aes(x = EH_ID, y = deltaCt_MminusE,
                           col = contamination)) + geom_point() + theme_bw()+facet_grid(.~InfectionHistoryTotal)


###########################################################################

##### Combine infection history and strain info #####

#Merge and clean design table
Expe007_design <- merge(Expe007a_design, Expe007b_design, all=TRUE)
Expe007_design <- Expe007_design %>% select(EH_ID, HybridStatus, Strain, Genome)

#A_B mouse = B_A mouse
Expe007_design$Strain <- as.character(Expe007_design$Strain)
a <- strsplit(Expe007_design$Strain, "_")
b <- lapply(a, sort)
c <- unlist(lapply(b, FUN = function(a){paste(a, collapse="-")}))
Expe007_design$Strain <- c

#Merge with Infection History and clean
InfectionTable <- merge(InfectionTable, Expe007_design, all=TRUE)
InfectionTable <- InfectionTable[-c(54), ]
#write.csv(InfectionTable, file = "C:/Users/Anna/Documents/HU Berlin/AG HU/Bachelorarbeit/Infection_history.csv")

#Plot
ggplot(InfectionTable, aes(x = Strain, y = deltaCt_MminusE,
                           col = Strain)) + geom_point() + theme_bw()+facet_grid(.~InfectionHistoryTotal)

ggplot(InfectionTable, aes(x = InfectionHistoryTotal, y = deltaCt_MminusE,
                           col = InfectionHistoryTotal)) + geom_point() + theme_bw()+facet_grid(.~Strain)

#Just to check: Infection History without contamination
ggplot(InfectionTable, aes(x = Strain, y = deltaCt_MminusE,col = Strain)) + geom_point() + theme_bw()+facet_grid(.~InfectionHistory) + coord_flip()
#ggplot(InfectionTable, aes(x = InfectionHistory, y = deltaCt_MminusE,col = InfectionHistory)) + geom_point() + theme_bw()+facet_grid(.~Strain)

#Plot
ggplot(InfectionTable, aes(x = InfectionHistoryTotal, y = deltaCt_MminusE,
                           col = InfectionHistoryTotal)) + geom_point() + theme_bw()+facet_grid(InfectionHistoryTotal~Strain)

##### Graph with deltact per oocysts on dpi_x #####

alldata <- merge(Expe007, InfectionTable, all = TRUE)

##dpi7##
#Like this, its just picking the rows that are true for the statement -> dplyr has it 
alldata_dpi7 <- alldata[alldata$dpi==7,]


ggplot(alldata_dpi7, aes(x = OPG, y = deltaCt_MminusE,
                  , col = InfectionHistoryTotal)) + geom_point() + theme_bw()+facet_grid(.~InfectionHistoryTotal)+coord_flip()

ggplot(alldata_dpi7, aes(x = OPG, y = deltaCt_MminusE,
                         , col = InfectionHistoryTotal)) + geom_point() + theme_bw()+facet_grid(.~Strain)


##dpi6##
alldata_dpi6 <- alldata[alldata$dpi==6,]

ggplot(alldata_dpi6, aes(x = OPG, y = deltaCt_MminusE,
                         col = InfectionHistoryTotal)) + geom_point() + theme_bw()

##dpi8##

alldata_dpi8 <- alldata[alldata$dpi==8,]

ggplot(alldata_dpi8, aes(x = OPG, y = deltaCt_MminusE,
                         col = InfectionHistoryTotal)) + geom_point() + theme_bw()


##### Explore Infection History #####

#table1 <- data.frame(InfectionTable$InfectionHistoryTotal, check.names = TRUE)
#?data.frame
#str(table1)

table(InfectionTable$InfectionHistoryTotal)
table(InfectionTable$InfectionHistory)
table(InfectionTable$HybridStatus, InfectionTable$InfectionHistory)
table(InfectionTable$Strain, InfectionTable$InfectionHistory)
table(InfectionTable$Genome, InfectionTable$InfectionHistory)


##### delta ct and oocyst counts #####

#calculate oocyst mean + clean table
OPGMean <- aggregate(alldata[ , 15], list(alldata$EH_ID), mean)
colnames(OPGMean)[colnames(OPGMean)=="x"] <- "OPGMean"
colnames(OPGMean)[colnames(OPGMean)=="Group.1"] <- "EH_ID"
#Status
OPGMean$OPGStatus <- ifelse(OPGMean$OPGMean > 0 , "Oocysts shed", "No oocysts")

#merge this with alldata
alldata <- merge(alldata, OPGMean, all=TRUE)

#try out with qPCRInfectionHistory
qPCRInfectionHistory <- merge(qPCRInfectionHistory, OPGMean, by = "EH_ID")

ggplot(alldata, aes(x = EH_ID, y = deltaCt_MminusE,
                     col = OPGStatus)) + geom_point() + theme_bw()

ggplot(alldata, aes(x = OPGMean, y = deltaCt_MminusE,
                    col = InfectionHistory)) + geom_point() + theme_bw()

ggplot(alldata, aes(x = OPGMean, y = deltaCt_MminusE,
                    col = InfectionHistoryTotal)) + geom_point() + theme_bw()

