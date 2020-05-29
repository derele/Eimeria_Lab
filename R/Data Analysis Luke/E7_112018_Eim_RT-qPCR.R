#E7_112018_Eim_RT-qPCR
library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(naniar)
library(tidyr)

# load in raw tables
RT1 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR1/E7_112018_Eim_RT-qPCR1.csv"
RT1 <- read.csv(text = getURL(RT1))
#remove LM_0231 because we have a repeat in RT6
RT1 <- RT1[!grepl("LM_0231", RT1$Name),] 

RT2 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR2/E7_RT-qPCR2.CSV"
RT2 <- read.csv(text = getURL(RT2))

RT3 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR3/E7_RT-qPCR3.CSV"
RT3 <- read.csv(text = getURL(RT3))

RT4 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR4/E7_112018_Eim_RT-qPCR4.CSV"
RT4 <- read.csv(text = getURL(RT4))

RT5 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR5/E7_112018_Eim_RT-qPCR5.CSV"
RT5 <- read.csv(text = getURL(RT5))

RT6 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR6/E7_RT-qPCR6.CSV"
RT6 <- read.csv(text = getURL(RT6))

RT7 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_RT-qPCRs/E7_RT-qPCR7/E7_RT-qPCR7.CSV"
RT7 <- read.csv(text = getURL(RT7))

# bind
RT1$Ct.Mean.SYBR <- NULL
RT1$Ct.Dev..SYBR <- NULL
RT7$Ct.Mean.SYBR <- NULL
RT7$Ct.Dev..SYBR <- NULL

RT <- rbind(RT1, RT2)
RT <- rbind(RT, RT3)
RT <- rbind(RT, RT4)
RT <- rbind(RT, RT5)
RT <- rbind(RT, RT6)
RT <- rbind(RT, RT7)

# remove  extra columns

RT$Amount.SYBR <- NULL
RT$Pos <- NULL

#remove negative controls
RT <- RT[!grepl("IRG6A", RT$Name),]
RT <- RT[!grepl("CXCR3", RT$Name),]
RT <- RT[!grepl("IL-12rb1", RT$Name),]
RT <- RT[!grepl("B-actin", RT$Name),]
RT <- RT[!grepl("GAPDH", RT$Name),]

# RT$Date <- NULL 
# name columns to match other data sets (Mouse_ID, Target) + make RT.CT numeric
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
names(RT)[names(RT) == "Name"] <- "Mouse_ID"
RT$RT.Ct <- as.numeric(as.character(RT$RT.Ct))
# remove NAs, calculate averages + save long
#RT <- na.omit(RT)
RT.long <- RT %>% dplyr::group_by(Mouse_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct, na.rm = T))
RT.long <- data.frame(RT.long)
RT.wide <- reshape(RT.long[, c("Target", "Mouse_ID","RT.Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")
# compare graphically becaus I'm just disabled like that
HKG1 <- dplyr::filter(RT.long, Target == "B-actin")
HKG2 <- dplyr::filter(RT.long, Target ==  "GAPDH")
HKG <- rbind(HKG1, HKG2)
ggplot(HKG, aes(x = Target, y = RT.Ct, color = Target)) +
  geom_point() +
  geom_boxplot() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank()) +
  ggtitle("HKG differences E7")
HKG$EXP <- "E7"
write.csv(HKG, "C:/Users/Luke Bednar/Eimeria_Lab/data/3_recordingTables/HKG_E7.csv")

# set ref and target genes
refGenes <- c("RT.Ct.B-actin", "RT.Ct.GAPDH")
targetGenes <- c("RT.Ct.CXCR3", "RT.Ct.IL-12", "RT.Ct.IRG6")
# calculate ref genes in new column and subtract targets from HKG average, create new columns
require(dplyr)
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT.wide, refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)

# continue with averaging refgenes and subtracting targets from them
RT.wide$CXCR3 <- (RT.wide$refMean - RT.wide$RT.Ct.CXCR3)
RT.wide$IRG6 <- (RT.wide$refMean - RT.wide$RT.Ct.IRG6)
RT.wide$IL.12 <- (RT.wide$refMean - RT.wide$RT.Ct.IL.12)
RT.wide[RT.wide=="NaN"]<-NA

RT.long <- reshape(RT.wide, 
            direction = "long",
            varying = list(names(RT.wide)[2:4]),
            v.names = "NE",
            times = c("CXCR3", "IRG6", "IL.12"),
            timevar = "Target",
            idvar = "Mouse_ID")
rownames(RT.long) <- NULL

# remove non normalized expressions
RT.wide$RT.Ct.CXCR3 <- NULL
RT.wide$RT.Ct.IRG6 <- NULL
RT.wide$RT.Ct.IL.12 <- NULL
RT.wide$RT.Ct.GAPDH <- NULL
RT.wide$RT.Ct.B.actin <- NULL
RT.wide$refMean <- NULL
RT.wide <- filter(RT.wide)
RT.wide <- distinct(RT.wide)
# write
write.csv(RT.wide, file = "../Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_RT-qPCR_wide.csv")
write.csv(RT.long, file = "../Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_RT-qPCR_long.csv")
