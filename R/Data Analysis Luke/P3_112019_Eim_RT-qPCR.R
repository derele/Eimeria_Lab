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
RT1$Pos <- NULL
RT1$Ct.SYBR <- NULL
RT1$Amount.SYBR <- NULL
RT1$Ct.Dev..SYBR <- NULL
RT1 <- distinct(RT1)

names(RT1)[names(RT1) == "Target.SYBR"] <- "Target"
names(RT1)[names(RT1) == "Ct.Mean.SYBR"] <- "RT.Ct"
names(RT1)[names(RT1) == "Name"] <- "EH_ID"
RT1 <- RT1 %>% dplyr::group_by(EH_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct))
RT1 <- ungroup(RT1)
RT1.wide <- reshape2::dcast(RT1, EH_ID ~ Target, value.var = "RT.Ct")

refGenes <- c("B-actin", "GAPDH")
targetGenes <- c("CXCR3", "IL-12", "IRG6")

RT1.wide <- RT1.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT1.wide, refGenes)))
RT1.wide <- data.frame(RT1.wide)
RT1.wide$CXCR3 <- (RT1.wide$refMean - RT1.wide$CXCR3)
RT1.wide$IRG6 <- (RT1.wide$refMean - RT1.wide$IRG6)
RT1.wide$IL.12 <- (RT1.wide$refMean - RT1.wide$IL.12)
RT1.wide$refMean <- NULL

RT1.long <- reshape(RT1.wide, 
                   direction = "long",
                   varying = list(names(RT1.wide)[2:4]),
                   v.names = "NE",
                   times = c("CXCR3", "IRG6", "IL.12"),
                   timevar = "Target",
                   idvar = "EH_ID")
rownames(RT1.long) <- NULL
RT1.wide$B.actin <- NULL
RT1.wide$GAPDH <- NULL

ggplot(RT1.long, aes(x = EH_ID, y = NE)) +
  geom_point(size = 3) +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  ggtitle("P3 gene expression")

