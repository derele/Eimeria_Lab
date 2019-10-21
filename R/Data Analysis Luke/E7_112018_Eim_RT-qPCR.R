#E7_112018_Eim_RT-qPCR
library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)
library(reshape2)
# load in raw tables
RT1 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_112018_Eim_RT-qPCR1/E7_112018_Eim_RT-qPCR1.csv"
RT1 <- read.csv(text = getURL(RT1))

RT2 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_112018_Eim_RT-qPCR2/E7_112018_Eim_RT-qPCR2.CSV"
RT2 <- read.csv(text = getURL(RT2))

RT3 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_112018_Eim_RT-qPCR3/E7_112018_Eim_RT-qPCR3.CSV"
RT3 <- read.csv(text = getURL(RT3))
# remove  extra columns
RT1$Ct.Mean.SYBR <- NULL
RT1$Ct.Dev..SYBR <- NULL
RT1$Amount.SYBR <- NULL
RT2$Amount.SYBR <- NULL
RT3$Amount.SYBR <- NULL
RT1$Pos <- NULL
RT2$Pos <- NULL
RT3$Pos <- NULL
# bind
RT <- bind_rows(RT1, RT2)
RT <- bind_rows(RT, RT3)
#check and remove negative controls
RT <- RT[!grepl("IRG6A", RT$Name),]
RT <- RT[!grepl("CXCR3", RT$Name),]
RT <- RT[!grepl("IL-12rb1", RT$Name),]
RT <- RT[!grepl("B-actin", RT$Name),]
RT <- RT[!grepl("GAPDH", RT$Name),]
# name columns to match other data sets (Mouse_ID, Target) + make RT.CT numeric
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
names(RT)[names(RT) == "Name"] <- "Mouse_ID"
RT$RT.Ct <- as.numeric(as.character(RT$RT.Ct))
# remove NAs, calculate averages + save long
RT <- na.omit(RT)
RT.wide <- reshape(RT[, c("Target", "Mouse_ID","RT.Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")
RT.long <- RT %>% dplyr::group_by(Mouse_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct))
RT.long <- data.frame(RT.long)
# set ref and target genes
refGenes <- c("RT.Ct.B-actin", "RT.Ct.GAPDH")
targetGenes <- c("RT.Ct.CXCR3", "RT.Ct.IL.12", "RT.Ct.IRG6")
# calculate ref genes in new column and subtract (sweep) from target genes
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE,select(., refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)
RT.wide[,2:4] <- sweep(RT.wide[2:4],1,refMean,'-')
RT.wide$refMean <- NULL
RT.wide$RT.Ct.B.actin <- NULL
RT.wide$RT.Ct.GAPDH <- NULL

# rename columns after melt + convert to long with NE values
names(RT.wide)[names(RT.wide) == "RT.Ct.CXCR3"] <- "CXCR3"
names(RT.wide)[names(RT.wide) == "RT.Ct.IL.12"] <- "IL-12"
names(RT.wide)[names(RT.wide) == "RT.Ct.IRG6"] <- "IRG6"

RT.long <- melt(RT.wide, id.vars = "Mouse_ID")
names(RT.long)[names(RT.long) == "variable"] <- "Target"
names(RT.long)[names(RT.long) == "value"] <- "NE"
# graph
ggplot(RT.long, aes(x = NE, y = Mouse_ID)) +
  geom_point() +
  facet_wrap("Target")
# load in complete mouse info 
complete <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_complete.csv"
complete <- read.csv(text = getURL(complete))
