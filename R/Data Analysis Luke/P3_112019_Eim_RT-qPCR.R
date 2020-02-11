# P3 RT-qPCR and infection intensity qPCR analysis

library(httr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(RCurl)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(naniar)
library(data.table)

# load in raw

RT1 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR1.CSV"
RT1 <- read.csv(text = getURL(RT1))

RT2 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR3.CSV"
RT2 <- read.csv(text = getURL(RT2))

RT3 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR5.CSV"
RT3 <- read.csv(text = getURL(RT3))

RT4 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_RT-qPCRs/P3_112019_Eim_RT-qPCR6.CSV"
RT4 <- read.csv(text = getURL(RT4))


RT1$Ct.Mean.SYBR <- NULL
RT1$Ct.Dev..SYBR <- NULL
RT1$Amount.SYBR <- NULL
RT2$Amount.SYBR <- NULL
RT2$Ct.Mean.SYBR <- NULL
RT2$Ct.Dev..SYBR <- NULL
RT3$Amount.SYBR <- NULL
RT3$Ct.Mean.SYBR <- NULL
RT3$Ct.Dev..SYBR <- NULL
RT4$Amount.SYBR <- NULL

RT1$Pos <- NULL
RT2$Pos <- NULL
RT3$Pos <- NULL
RT4$Pos <- NULL

RT <- rbind(RT1, RT2)
RT <- rbind(RT, RT3)
RT <- rbind(RT, RT4)


names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
names(RT)[names(RT) == "Name"] <- "EH_ID"
RT <- RT %>% dplyr::group_by(EH_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct))
RT <- ungroup(RT)
RT.wide <- reshape2::dcast(RT, EH_ID ~ Target, value.var = "RT.Ct")

refGenes <- c("B-actin", "GAPDH")
targetGenes <- c("CXCR3", "IL-12", "IRG6")

RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT.wide, refGenes)))
RT.wide <- data.frame(RT.wide)
# Ref - Target = lower value (bigger minus means less expression)
RT.wide$CXCR3 <- (RT.wide$refMean - RT.wide$CXCR3)
RT.wide$IRG6 <- (RT.wide$refMean - RT.wide$IRG6)
RT.wide$IL.12 <- (RT.wide$refMean - RT.wide$IL.12)
RT.wide$refMean <- NULL

RT.long <- reshape(RT.wide, 
                   direction = "long",
                   varying = list(names(RT.wide)[2:4]),
                   v.names = "NE",
                   times = c("CXCR3", "IRG6", "IL.12"),
                   timevar = "Target",
                   idvar = "EH_ID")
rownames(RT.long) <- NULL
RT.wide$B.actin <- NULL
RT.wide$GAPDH <- NULL

ggplot(RT.long, aes(x = EH_ID, y = NE)) +
  geom_point(size = 3) +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  ggtitle("P3 gene expression")

######### add infection intensity
P3 <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_qPCRs/P3_112019_Eim_qPCR1.CSV"
P3 <- read.csv(text = getURL(P3))
P3.long <- select(P3, Name, Ct.Mean.SYBR, Target.SYBR)
P3.long <- distinct(P3.long)

P3.long <- P3.long %>% 
  dcast(Name ~ Target.SYBR, value.var = "Ct.Mean.SYBR", fill = 0) %>% 
  mutate(delta = mouse - eimeria) %>% 
  dplyr::select(Name,delta)
names(P3.long)[names(P3.long) == "Name"] <- "EH_ID"

P3_qPCRs <- merge(RT.long, P3.long, by = "EH_ID")

ggplot(P3_qPCRs, aes(x = delta, y = NE)) +
  geom_point(size = 3) +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  ggtitle("P3 gene expression vs delta")

####### write out
write.csv(P3_qPCRs, "./Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_qPCRs_complete.csv")

