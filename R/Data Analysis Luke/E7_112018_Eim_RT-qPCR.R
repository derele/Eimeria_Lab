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
#remove LM_0231 because we have a repeat in RT6
RT1 <- RT1[!grepl("LM_0231", RT1$Name),] 

RT2 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_112018_Eim_RT-qPCR2/E7_112018_Eim_RT-qPCR2.CSV"
RT2 <- read.csv(text = getURL(RT2))

RT3 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_112018_Eim_RT-qPCR3/E7_112018_Eim_RT-qPCR3.CSV"
RT3 <- read.csv(text = getURL(RT3))

RT4 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_112018_Eim_RT-qPCR4/E7_112018_Eim_RT-qPCR4.CSV"
RT4 <- read.csv(text = getURL(RT4))

RT5 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_112018_Eim_RT-qPCR5/E7_112018_Eim_RT-qPCR5.CSV"
RT5 <- read.csv(text = getURL(RT5))

RT6 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR6/E7_RT-qPCR6.CSV"
RT6 <- read.csv(text = getURL(RT6))

# remove  extra columns
RT1$Ct.Mean.SYBR <- NULL
RT1$Ct.Dev..SYBR <- NULL
RT1$Amount.SYBR <- NULL
RT2$Amount.SYBR <- NULL
RT3$Amount.SYBR <- NULL
RT4$Amount.SYBR <- NULL
RT5$Amount.SYBR <- NULL
RT6$Amount.SYBR <- NULL
RT1$Pos <- NULL
RT2$Pos <- NULL
RT3$Pos <- NULL
RT4$Pos <- NULL
RT5$Pos <- NULL
RT6$Pos <- NULL

# convert all Ct.SYBR to numbers
RT1$Ct.SYBR <- as.numeric(levels(RT1$Ct.SYBR))[RT1$Ct.SYBR]
RT2$Ct.SYBR <- as.numeric(levels(RT2$Ct.SYBR))[RT2$Ct.SYBR]
RT3$Ct.SYBR <- as.numeric(levels(RT3$Ct.SYBR))[RT3$Ct.SYBR]
RT4$Ct.SYBR <- as.numeric(levels(RT4$Ct.SYBR))[RT4$Ct.SYBR]
RT5$Ct.SYBR <- as.numeric(levels(RT5$Ct.SYBR))[RT5$Ct.SYBR]
RT6$Ct.SYBR <- as.numeric(levels(RT6$Ct.SYBR))[RT6$Ct.SYBR]

# bind
RT <- dplyr::bind_rows(RT1, RT2)
RT <- bind_rows(RT, RT3)
RT <- bind_rows(RT, RT4)
RT <- bind_rows(RT, RT5)
RT <- bind_rows(RT, RT6)
#remove negative controls
RT <- RT[!grepl("IRG6A", RT$Name),]
RT <- RT[!grepl("CXCR3", RT$Name),]
RT <- RT[!grepl("IL-12rb1", RT$Name),]
RT <- RT[!grepl("B-actin", RT$Name),]
RT <- RT[!grepl("GAPDH", RT$Name),]
# load in MC analysis of RT
RTMC <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR_MC.csv"
RTMC <- read.csv(sep = ";", text = getURL(RTMC))
names(RTMC)[names(RTMC) == "EH_ID"] <- "Name"
names(RTMC)[names(RTMC) == "Target"] <- "Target.SYBR"
RT <- merge(RT, RTMC)
#RT <- RT[!(RT$Caecum == "neg"),]
RT$Observer <- NULL
RT$Date <- NULL 
# name columns to match other data sets (Mouse_ID, Target) + make RT.CT numeric
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
names(RT)[names(RT) == "Name"] <- "Mouse_ID"
RT$RT.Ct <- as.numeric(as.character(RT$RT.Ct))
# remove NAs, calculate averages + save long
RT <- na.omit(RT)
RT.long <- RT %>% dplyr::group_by(Mouse_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct))
RT.long <- data.frame(RT.long)
RT.wide <- reshape(RT.long[, c("Target", "Mouse_ID","RT.Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")
# set ref and target genes
refGenes <- c("RT.Ct.B-actin", "RT.Ct.GAPDH")
targetGenes <- c("RT.Ct.CXCR3", "RT.Ct.IL-12", "RT.Ct.IRG6")
# calculate ref genes in new column and subtract targets from HKG average, create new columns
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE,select(., refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)
RT.wide$CXCR3 <- (RT.wide$refMean - RT.wide$RT.Ct.CXCR3)
RT.wide$IRG6 <- (RT.wide$refMean - RT.wide$RT.Ct.IRG6)
RT.wide$IL.12 <- (RT.wide$refMean - RT.wide$RT.Ct.IL.12)
# remove non normalized expressions
RT.wide$RT.Ct.CXCR3 <- NULL
RT.wide$RT.Ct.IRG6 <- NULL
RT.wide$RT.Ct.IL.12 <- NULL
RT.wide$RT.Ct.GAPDH <- NULL
RT.wide$RT.Ct.B.actin <- NULL
RT.wide$refMean <- NULL

# #not necessary
# RT.wide[,2:4] <- sweep(RT.wide[2:4],1,refMean,'-')
# RT.wide$refMean <- NULL
# RT.wide$RT.Ct.B.actin <- NULL
# RT.wide$RT.Ct.GAPDH <- NULL

# # rename columns after melt + convert to long with NE values
# names(RT.wide)[names(RT.wide) == "RT.Ct.CXCR3"] <- "CXCR3"
# names(RT.wide)[names(RT.wide) == "RT.Ct.IL.12"] <- "IL-12"
# names(RT.wide)[names(RT.wide) == "RT.Ct.IRG6"] <- "IRG6"

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
names(RT.wide)[names(RT.wide) == "Mouse_ID"] <- "EH_ID"
names(RT.long)[names(RT.long) == "Mouse_ID"] <- "EH_ID"
complete$X <- NULL
# split EH_ID name and make with "_"
complete$EH_ID <- gsub("LM", "LM_", complete$EH_ID)
complete <- merge(complete, RT.long, by = "EH_ID")
#quick basic subtract before merge (based on E7 melting curves)
E7EimMC <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_detection_MC.csv"
E7EimMC <- read.csv(text = getURL(E7EimMC))
complete <- merge(complete, E7EimMC, by = "EH_ID")
complete <- complete[!(complete$Caecum == "neg"),]
# add intensity 
E7_inf <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_Anna_qPCR_DNA_ct_Zusammenfassung.csv"
E7_inf <- read.csv(text = getURL(E7_inf))
E7_inf$Ct.SYBR <- NULL
E7_inf$Pos <- NULL
E7_inf$Amount.SYBR..Copies. <- NULL
E7_inf$Amount.Mean.SYBR <- NULL
E7_inf$Amount.Dev..SYBR <- NULL
E7_inf <- distinct(E7_inf)
E7_inf <- E7_inf %>% 
  dcast(Name ~ Target.SYBR, value.var = "Ct.Mean.SYBR", fill = 0) %>% 
  mutate(delta = mouse - eimeria) %>% 
  dplyr::select(Name,delta)
E7_inf <- E7_inf %>% tidyr::separate(Name, c("LM", "EH_ID"))
E7_inf$EH_ID <- sub("^", "LM", E7_inf$EH_ID )
E7_inf$LM <- NULL
E7_inf$EH_ID <- gsub("LM", "LM_", E7_inf$EH_ID)
complete <- merge(complete, E7_inf, by = "EH_ID")
# graph 
ggplot(complete, aes(x = delta, y = NE, color = Target)) +
  geom_point() +
  facet_wrap("infHistory") +
  geom_smooth(method = "lm")
ggplot(complete, aes(x = NE, y = delta, color = Target)) +
  geom_point() + 
  facet_wrap("HybridStatus") +
  geom_smooth(method = "lm")
# load in wild data for comparison
HZ18 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ18_RT-qPCR_RTlong.csv"
HZ18 <- read.csv(text = getURL(HZ18))
names(HZ18)[names(HZ18) == "inf"] <- "Caecum"
HZ18$Caecum <- as.character(HZ18$Caecum)
HZ18$Caecum[HZ18$Caecum == "TRUE"] <- "pos"
HZ18$Caecum[HZ18$Caecum == "FALSE"] <- "neg"
# graph to compare NE between infected and non-infected





ggplot(HZ18, aes(x = HI, y = NE, color = Caecum)) +
  geom_point() +
  facet_wrap("Target") +
  geom_smooth(method = "lm")


HZ18 <- HZ18[!(HZ18$Caecum == "neg"),]

#process to graph together
HZ18 <- HZ18[!(HZ18$Target=="GBP2"),]
HZ18 <- HZ18[!(HZ18$Target=="IL-6"),]
i <- sapply(HZ18, is.factor)
HZ18[i] <- lapply(HZ18[i], as.character)
HZ18$Target[HZ18$Target == "IL-12b"] <- "IL-12"
HZ18$inf <- NULL

E7 <- merge(RT.long, E7_inf)
E7 <- merge(E7, E7EimMC, by = "EH_ID")
E7 <- E7[!(E7$Caecum == "neg"),]
names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"
E7 %>% mutate_if(is.factor, as.character) -> E7
E7$Target[E7$Target == "IL.12"] <- "IL-12"

ggplot(HZ18, aes(x = NE, y = delta)) +
  geom_point() +
  facet_wrap("Target") +
  coord_flip()

ggplot(E7, aes(x = NE, y = delta)) +
  geom_point() + 
  facet_wrap("Target") + 
  coord_flip()
# combine in one table and distinguish as batches (rename Mouse_ID to EH_ID for sake of merging)

HZ18$batch <- "HZ18"
E7$batch <- "E7"
names(HZ18)[names(HZ18) == "Mouse_ID"] <- "EH_ID"
All <- bind_rows(HZ18, E7)
#remove dodgy outliers
All <- All[-c(127, 129), ]

ggplot(All, aes(x = NE, y = delta, color = batch)) +
  geom_point() +
  facet_wrap("Target") +
  coord_flip() +
  geom_smooth(method = "lm")

# subest for graphing
IL6
