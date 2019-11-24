#E7_112018_Eim_RT-qPCR
library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)
library(reshape2)
library(ggpubr)
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
require(dplyr)
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT.wide, refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)

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
# and subset to make smaller (reduce repeating points)
intensity <- select(complete, EH_ID, delta, infHistory, Caecum)
intensity <- distinct(intensity)
# substract ref genes individually to check whether they might influence the gene expression when averaged together
B_actin <- RT.wide
B_actin$CXCR3 <- (B_actin$RT.Ct.B.actin - B_actin$RT.Ct.CXCR3)
B_actin$IRG6 <- (B_actin$RT.Ct.B.actin - B_actin$RT.Ct.IRG6)
B_actin$IL12 <- (B_actin$RT.Ct.B.actin - B_actin$RT.Ct.IL.12)

GAPDH <- RT.wide
GAPDH$CXCR3 <- (GAPDH$RT.Ct.GAPDH - GAPDH$RT.Ct.CXCR3)
GAPDH$IRG6 <- (GAPDH$RT.Ct.GAPDH - GAPDH$RT.Ct.IRG6)
GAPDH$IL12 <- (GAPDH$RT.Ct.GAPDH - GAPDH$RT.Ct.IL.12)
#graph out to compare (B actin)
B_actin.long <- gather(B_actin, Target, NE, CXCR3:IL12, factor_key=TRUE)
B_actin.long$RT.Ct.B.actin <- NULL
B_actin.long$RT.Ct.CXCR3 <- NULL
B_actin.long$RT.Ct.GAPDH <- NULL
B_actin.long$RT.Ct.IL.12 <- NULL
B_actin.long$RT.Ct.IRG6 <- NULL
B_actin.long$refMean <- NULL
B_actin.long <- merge(B_actin.long, intensity)

ggplot(B_actin.long, aes(infHistory, NE)) +
  geom_jitter() +
  geom_boxplot() + 
  facet_wrap("Target",  scales = "free_y")
#graph out to compare (B actin)
GAPDH.long <- gather(GAPDH, Target, NE, CXCR3:IL12, factor_key=TRUE)
GAPDH.long$RT.Ct.B.actin <- NULL
GAPDH.long$RT.Ct.CXCR3 <- NULL
GAPDH.long$RT.Ct.GAPDH <- NULL
GAPDH.long$RT.Ct.IL.12 <- NULL
GAPDH.long$RT.Ct.IRG6 <- NULL
GAPDH.long$refMean <- NULL
GAPDH.long <- merge(GAPDH.long, intensity)

ggplot(GAPDH.long, aes(infHistory, NE)) +
  geom_jitter() +
  geom_boxplot() + 
  facet_wrap("Target",  scales = "free_y")

# continue with averaging refgenes and subtracting targets from them
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

# graph E7 infection intensity before deleting negs


ggplot(intensity, aes(x = delta, y = Caecum, color = infHistory)) +
  geom_jitter(size = 3) +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))
########################################################################
# delete negs
# complete <- complete[!(complete$Caecum == "neg"),]

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
HZ18 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ18_complete.csv"
HZ18 <- read.csv(text = getURL(HZ18))
names(HZ18)[names(HZ18) == "inf"] <- "Caecum"
HZ18$Caecum <- as.character(HZ18$Caecum)
HZ18$Caecum[HZ18$Caecum == "TRUE"] <- "pos"
HZ18$Caecum[HZ18$Caecum == "FALSE"] <- "neg"
names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"
# graph to compare NE between infected and non-infected
ggplot(HZ18, aes(x = HI, y = NE, color = Caecum)) +
  geom_point() +
  facet_wrap("Target") +
  geom_smooth(method = "lm")

ggplot(HZ18, aes(x = HI, y = NE, color = Eimeria.subspecies)) +
  geom_point() +
  facet_wrap("Target") +
  geom_smooth(method = "lm")

ggplot(HZ18, aes(x = delta, y = NE, color = Eimeria.subspecies)) +
  geom_point() +
  facet_wrap("Target") +
  geom_smooth(method = "lm")

###################################################################
# HZ18 <- HZ18[!(HZ18$Caecum == "neg"),]

#process to graph together
HZ18 <- HZ18[!(HZ18$Target=="GBP2"),]
HZ18 <- HZ18[!(HZ18$Target=="IL-6"),]
i <- sapply(HZ18, is.factor)
HZ18[i] <- lapply(HZ18[i], as.character)
HZ18$Target[HZ18$Target == "IL-12b"] <- "IL-12"
HZ18$inf <- NULL

E7 <- merge(RT.long, E7_inf)
E7 <- merge(E7, E7EimMC, by = "EH_ID")
############################################################
#remove negs
 E7 <- E7[!(E7$Caecum == "neg"),]
names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"
E7 %>% mutate_if(is.factor, as.character) -> E7
E7$Target[E7$Target == "IL.12"] <- "IL-12"


ggplot(HZ18, aes(x = HI, y = NE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap("Target")

ggplot(HZ18, aes(x = HI, y = delta)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap("Target")


# combine in one table and distinguish as batches (rename Mouse_ID to EH_ID for sake of merging)

HZ18$batch <- "HZ18"
E7$batch <- "E7"
complete$batch <- "E7"
names(HZ18)[names(HZ18) == "Mouse_ID"] <- "EH_ID"
All <- bind_rows(HZ18, E7)
#remove dodgy outliers
# All <- All[-c(127, 129), ]

ggplot(All, aes(x = NE, y = delta, color = batch)) +
  geom_point() +
  facet_wrap("Target") +
  coord_flip() +
  geom_smooth(method = "lm")

# just complete minus outliers
complete <-complete[!(complete$NE > 0),]
All <-All[!(All$NE > 0),]

# E7 graph only sample positive for Eimeria in the caecum
ggplot(data = subset(complete, !is.na(x = Target) & !(complete$Caecum == "neg")), aes(x = infHistory, y = NE, color = infHistory)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap("Target") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Gene expression of Eimeria positive mice")

# HZ18 graph only sample positive for Eimeria in the caecum
ggplot(data = subset(HZ18, !is.na(x = Target) & !(HZ18$Caecum == "neg")), aes(x = Eimeria.subspecies, y = NE, color = Eimeria.subspecies)) +
  geom_point() + 
  geom_boxplot() +
  facet_wrap("Target") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Gene expression of Eimeria positive mice")

# Laboratory gene expression
plot_complete <- select(complete, EH_ID, HybridStatus, batch, infHistory, Target, NE, Caecum, delta, primary, challenge)
plot_complete <- distinct(plot_complete)
plot_complete %>% mutate_if(is.factor, as.character) -> plot_complete
plot_complete$Target[plot_complete$Target == "IL.12"] <- "IL-12"
plot_complete <-plot_complete[!(plot_complete$NE > 0),]
plot_complete <-plot_complete[!(plot_complete$delta > 0),]
plot_complete[plot_complete == "inter subsp. hybrids"] <- "Outbreds"
plot_complete[plot_complete == "outbred hybrids"] <- "Hybrids"
plot_complete[plot_complete == "parental strains"] <- "Parentals"
plot_complete[plot_complete == "E64:E64"] <- "E. ferrisi : E. ferrisi"
plot_complete[plot_complete == "E64:E88"] <- "E. ferrisi : E. falciformis"
plot_complete[plot_complete == "E88:E64"] <- "E. falciformis : E. ferrisi"
plot_complete[plot_complete == "E88:E88"] <- "E. falciformis : E. falciformis"
names(plot_complete)[names(plot_complete) == "infHistory"] <- "Infection_History"
plot_complete <- na.omit(plot_complete)



ggplot(data = plot_complete, aes(x = Infection_History, y = NE)) +
  geom_boxplot() +
  facet_wrap("Target",  scales = "free_y") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.margin=unit(c(1,1,1,2),"cm"),
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("Laboratory mice gene expression")


# Laboratory infection intensity
ggplot(data = subset(plot_complete, !is.na(x = Target)), aes(x = Caecum, y = delta, color = Infection_History)) +
  geom_jitter(size = 2) +
  coord_flip() +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("Laboratory mice infection intensity")

# Hybrid Status effect (or lack thereof)
ggplot(data = subset(plot_complete, !is.na(x = Target)), aes(x = HybridStatus, y = NE, color = Target)) +
  geom_boxplot() +
  facet_wrap("Infection_History") +
  geom_jitter(position = position_jitterdodge()) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("Laboratory mice gene expression")

# Wild 2018 gene expression
plot_HZ18 <- select(HZ18, EH_ID, Target, HI, delta, Eimeria.subspecies, NE, Caecum, batch)
plot_HZ18 <- distinct(plot_HZ18)
plot_HZ18 <- plot_HZ18[!(plot_HZ18$NE > 0),]

ggplot(data = subset(plot_HZ18, !is.na(x = Target)), aes(x = Eimeria.subspecies, y = NE, color = Target)) +
  geom_boxplot() +
  geom_jitter(position = position_jitterdodge()) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("Wild mice gene expression")
# Wild 2018 infection intensity
ggplot(data = subset(plot_HZ18, !is.na(x = Target)), aes(x = NE, y = delta)) +
  geom_jitter() +
  geom_violin() +
  facet_wrap("Eimeria.subspecies") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("Wild mice infection intensity")

# Wild Hybrid effect
ggplot(data = subset(plot_HZ18, !is.na(x = Target)), aes(x = HI, y = NE, color = Target)) +
  geom_point() +
  facet_wrap("Eimeria.subspecies") +
  geom_smooth(method = "lm")
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("Wild mice infection intensity")


