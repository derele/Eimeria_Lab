library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(stats)
library(ggsignif)
library(ggpmisc)

# good graphing
# ggplot(HZgraph, 
#        aes(x = MC, y = NE, color = MC)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter() +
#   facet_wrap("Target", scales = "free") +
#   labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
#   theme(axis.text=element_text(size=12, face = "bold"),
#         title = element_text(size = 16, face = "bold"),
#         axis.title=element_text(size=14,face="bold"),
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold"))+
#   ggtitle("Gene expression in wild samples")

######
complete <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_P3_E6_complete.csv"))
# make negative MCs into NAs in a new column
complete$delta_clean <- complete$delta
complete <- mutate(complete, delta_clean = ifelse(Eim_MC == "neg", -30, delta_clean))
complete$dpi <- as.factor(complete$dpi)
complete$X <- NULL
########## make column with E. ferrisi and E. falciformis only
complete$Eimeria.p <- gsub("E64|E139", replacement = "E.ferrisi", complete$primary)
complete$Eimeria.p <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria.p))

complete$Eimeria.c <- gsub("E64|E139", replacement = "E.ferrisi", complete$challenge)
complete$Eimeria.c <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria.c))

########## gene expression
genes <- dplyr::select(complete, EH_ID, CXCR3, IL.12, IRG6)
Eim <- dplyr::select(complete, Eim_MC, EH_ID)
Eim <- dplyr::distinct(Eim)
Eim <- na.omit(Eim)
genes <- dplyr::distinct(genes)
genes <- merge(genes, Eim, by.y = "EH_ID", all = T)
genes <- dplyr::distinct(genes)
# tranform into long
# genes <- 
genes <- reshape2::melt(genes,
                        direction = "long",
                        varying = list(names(genes)[2:4]),
                        v.names = "NE",
                        na.rm = T, value.name = "NE", 
                        id.vars = c("EH_ID", "Eim_MC"))
names(genes)[names(genes) == "variable"] <- "Target"

############################### rope in FACS stuff for E7 too ###########################
FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_complete.csv"))
FACS$X <- NULL
FACS.long <- dplyr::select(FACS, EH_ID, Position, infHistory, CD4, Treg, Div_Treg, Treg17, Treg_prop, Th1, Div_Th1,
                           Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)
FACS.long <- dplyr::distinct(FACS.long)
# transform into long
FACS.long <- reshape2::melt(FACS.long,
             direction = "long",
             varying = list(names(FACS.long)[19:34]),
             v.names = "cell.pop",
             na.rm = T, value.name = "counts", 
             id.vars = c("EH_ID", "Position", "infHistory"))
FACS.long <- na.omit(FACS.long)
names(FACS.long)[names(FACS.long) == "variable"] <- "pop"

##### make a summary for comparing with wild (no Pos or Ant difference)
FACSt <- FACS.long %>% dplyr::group_by(EH_ID, pop, infHistory) %>% dplyr::summarise(counts = mean(counts, na.rm = T))

write.csv(FACSt, "~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACSt.csv")

############################### rope in FACS stuff for HZ19 too ###########################
FACSHZ <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_complete.csv"))
FACSHZ$X <- NULL
FACSHZ.long <- dplyr::select(FACSHZ, Mouse_ID, Position, CD4, Treg, Div_Treg, Treg17, Treg_prop, Th1, Div_Th1,
                           Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)
FACSHZ.long <- dplyr::distinct(FACSHZ.long)
# transform into long
FACSHZ.long <- reshape2::melt(FACSHZ.long,
                            direction = "long",
                            varying = list(names(FACSHZ.long)[3:17]),
                            v.names = "cell.pop",
                            na.rm = T, value.name = "counts", 
                            id.vars = c("Mouse_ID", "Position"))
FACSHZ.long <- na.omit(FACSHZ.long)
names(FACSHZ.long)[names(FACSHZ.long) == "variable"] <- "pop"

##### make a summary for comparing with wild (no Pos or Ant difference)
FACStHZ <- subset(FACSHZ.long, Position == "mLN")
FACStHZ <- FACSHZ.long %>% dplyr::group_by(Mouse_ID, pop) %>% dplyr::summarise(counts = mean(counts, na.rm = T))

write.csv(FACStHZ, "~/Documents/Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACSt.csv")
# add EXP columns to FACSt and FACStHZ, make dummy infHistory column and rename Mouse ID to EH ID for sake of merging
FACSt$EXP <- "lab"
FACStHZ$EXP <- "wild"
FACStHZ$infHistory <- NA
colnames(FACStHZ)[1] <- "EH_ID"
FACScombine <- rbind(FACSt, FACStHZ)

################################################################################################################################


#####################################################################################################################################

# merge with genes
FACScombine <- data.frame(FACScombine)
immuno <- merge(FACScombine, genes, all = T)

inf <- dplyr::select(complete, EH_ID, challenge)
inf <- distinct(inf)
inf <- na.omit(inf)
immuno <- full_join(immuno, inf, all.x = T)


# add IFN ELISA results``````````````````````
IFN <- dplyr::select(complete, EH_ID, IFNy_CEWE, dpi)
IFN <- subset(IFN, IFN$dpi == "8")
IFN <- dplyr::select(IFN, EH_ID, IFNy_CEWE)
IFN <- na.omit(IFN) 
IFN <- distinct(IFN)
# super desructive but will have to do for now
# IFN <- IFN[complete.cases(IFN[ , 2:3]),]
immuno <- merge(immuno, IFN, by = "EH_ID")
immuno <- distinct(immuno)


delta <- dplyr::select(complete, delta, delta_clean, EH_ID, Eim_MC)
delta <- na.omit(delta)
immuno <- merge(immuno, delta)
immuno <- distinct(immuno)

immuno$species <- gsub("E64|E139", replacement = "E.ferrisi", immuno$challenge)
immuno$species <- paste(gsub("E88|Eflab", replacement = "E.falciformis", immuno$species))


Wch.p <- subset(complete, Eimeria.p == "E.ferrisi")
Wch.p1 <- subset(complete, Eimeria.p == "E.falciformis")
Wch.p <- full_join(Wch.p, Wch.p1)

Wch.c <- subset(complete, Eimeria.c == "E.ferrisi")
Wch.c1 <- subset(complete, Eimeria.c == "E.falciformis")
Wch.c <- full_join(Wch.c, Wch.c1)

Wch <- dplyr::select(Wch.c, Eimeria.c, Eimeria.p, EH_ID)
genes <- merge(genes, Wch)
genes <- distinct(genes)
genes <- na.omit(genes)

#### add wild gene expression
HZgenes <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-18_gene_expression.csv"))
HZgenes$X <- NULL
colnames(HZgenes)[1] <- "EH_ID"
colnames(HZgenes)[6] <- "Eim_MC"
HZgenes$Eim_MC[HZgenes$Eim_MC == "TRUE"] <- "positive"
HZgenes$Eim_MC[HZgenes$Eim_MC == "FALSE"] <- "negative"
HZgenes$EXP <- "wild"
immuno$Eim_MC <- as.character(immuno$Eim_MC)
immuno$Eim_MC[immuno$Eim_MC == "pos"] <- "positive"
immuno$Eim_MC[immuno$Eim_MC == "neg"] <- "negative"
genes$Eim_MC <- as.character(genes$Eim_MC)
genes$Eim_MC[genes$Eim_MC == "pos"] <- "positive"
genes$Eim_MC[genes$Eim_MC == "neg"] <- "negative"

# god knows what happens here (investigate later)
immuno <- merge(HZgenes, immuno, all = T)
#### add IFN CEWE HZ19
IFN_HZ <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISAs_complete.csv"))
IFN_HZ$X <- NULL
colnames(IFN_HZ)[1] <- "EH_ID"
colnames(IFN_HZ)[2] <- "IFNy_CEWE"

IFNcomplete <- rbind(IFN, IFN_HZ)
WxL <- merge(IFNcomplete, FACScombine, all = T)
MC <- dplyr::select(delta, EH_ID, Eim_MC)
WxL <- merge(WxL, MC, all =T)

Wch.p <- subset(Wch.p, !Wch.p$dpi == 0)
Wch.p$OPG[is.na(Wch.p$OPG)] <- 0

lab <- subset(immuno, immuno$EXP == "lab")
wild <- subset(immuno, immuno$EXP == "wild")
lab$species[lab$species == "E.falciformis"] <- "uninfected"

FACS_IFN <- merge(FACScombine, IFNcomplete, all = T)
FACScombine <- merge(FACScombine, Eim, all =T)


# sig <- subset(FACScombine, subset = pop %in% c("CD4", "Div_Treg", "Treg17", "Th17", "Div_Th17", "CD8", "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8"))
# sig <- data.frame(sig)
# sig$pop <- as.character(sig$pop)
# # remove horrible outliers
# sig <- sig[-c(1174),]
# sig <- sig[-c(1173),]
# sig <- sig[-c(634),]
# 
# sig1 <- merge(sig, Wch.c, by = "EH_ID")
# Pos <- select(immuno, EH_ID, Position)
# sig1 <- distinct(sig1)
# sig1 <- merge(Pos, sig1, by = "EH_ID")
# sig1 <- distinct(sig1)
# 
# PosHZ <- select(FACSHZ, Mouse_ID, Position, CD4)
# colnames(PosHZ)[1] <- "EH_ID"
# FACStHZ <- full_join(PosHZ, FACStHZ)
# immuno <- distinct(immuno)

##################### pure graphing from here, any general code above ####################################################
# ggplot(immuno, aes(x = EXP, y = counts)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter() +
#   facet_wrap("pop", scales = "free") +
#   labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
#   theme(axis.text=element_text(size=12, face = "bold"),
#         title = element_text(size = 16, face = "bold"),
#         axis.title=element_text(size=14,face="bold"),
#         strip.text.x = element_text(size = 14, face = "bold"),
#         legend.text=element_text(size=12, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold"))+
#   ggtitle("Gene expression in wild samples")
# 

################## Wchange graphs ###########################################################################

# falciformis plot
ggplot(subset(Wch.p, Wch.p$Eimeria.p == "E.falciformis"), 
       aes(x = dpi, y = Wchange, group = Eimeria.p, color = Eimeria.p)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  labs(y="% Weight change", x = "days post infection", colour = "species") +
  #facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold.italic"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.1))

# complete plot
ggplot(Wch.p,
       aes(x = dpi, y = Wchange, color = Eimeria.p, group = Eimeria.p)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  labs(y="% Weight change", x = "days post infection", colour = "species") +
  #facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold.italic"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.1))
  
##########################################################################################################
##################### OPG graphs
# falciformis shedding

ggplot(subset(Wch.p, Wch.p$Eimeria.p == "E.falciformis"), 
       aes(x = dpi, y = OPG, color = Eimeria.p, group = Eimeria.p)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  labs(y="oocysts per gramm of feces", x = "days post infection", colour = "species") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold.italic"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.7))

# falciformis and ferrisi shedding
ggplot(Wch.p, 
       aes(x = dpi, y = OPG, color = Eimeria.p, group = Eimeria.p)) +
  geom_jitter(width = 0.2) +
  geom_smooth() +
  labs(y="oocysts per gramm of feces", x = "days post infection", colour = "species") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold.italic"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.7))

#####################################################################################################################
############ genes
# wild gene expression
ggplot(HZgenes,
       aes(x = Eim_MC, y = NE, colour = Eim_MC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =0.95) +
  facet_wrap("Target", scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2))

# lab gene expression
genes <- genes[!(genes$Target=="IRG6"& genes$NE > 19),]

ggplot(genes, aes(x = Eim_MC, y = NE, color = Eim_MC, group = Eim_MC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..),  size = 8, label.y.npc =0.95) +
  facet_wrap("Target", scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())


##################################################################################################
###### IFN CEWE vs cell populations

### don't have that yet
# # IFN effect on wild

ggplot(subset(FACS_IFN, (FACS_IFN$EXP == "wild")),
       aes(x = IFNy_CEWE , y = counts)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap("pop", scales = "free", drop = T) +
  labs(y="% of populations", x = "IFNy (ng/mL)", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))

# lab IFN CEWE vs populations
ggplot(lab,aes(x = IFNy_CEWE , y = counts, color = Eim_MC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap("pop", scales = "free", drop = T) +
  labs(y="% of populations", x = "IFNy (ng/mL)", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
               axis.text.x = element_text(angle = 45, hjust = 1),
               title = element_text(size = 16, face = "bold"),
               axis.title=element_text(size=14,face="bold"),
               strip.text.x = element_text(size = 14, face = "bold"),
               legend.text=element_text(size=12, face = "bold"),
               legend.title = element_text(size = 12, face = "bold")) 



#### genes vs IFN
ggplot(lab, aes(y = NE, x = IFNy_CEWE, color = Eim_MC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Target, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))


########## E7 vs HZ19 FACS
# overall
FACScombine <- FACScombine[!(FACScombine$pop=="Treg_prop"),]
FACScombine <- FACScombine[!(FACScombine$pop=="Treg"& FACScombine$counts > 40),]
FACScombine <- FACScombine[!(FACScombine$pop=="Th17"& FACScombine$counts > 20),]
FACScombine <- FACScombine[!(FACScombine$pop=="Treg17"& FACScombine$counts > 50),]


ggplot(subset(FACScombine, !is.na(FACScombine$EXP)),
       aes(x = EXP , y = counts, color = EXP)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "origin", colour = "origin") +
  scale_color_manual(values=c("red", "darkgreen"))+
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))

immunogenes <- dplyr::select(immuno, EXP, NE, Target, EH_ID)
immunogenes <- distinct(immunogenes)

ggplot(subset(immunogenes, !is.na(immunogenes$EXP)),
       aes(x = EXP , y = NE, color = EXP)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  facet_wrap("Target", scales = "free") +
  labs(y="delta", x = "origin", colour = "origin") +
  scale_color_manual(values=c("red", "darkgreen"))+
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))



############################################################ infected lab
ggplot(subset(FACScombine, FACScombine$EXP == "lab"),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infected", colour = "infected") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())
################## graph the flow chart style of populations
### lab CD4
ggplot(subset(lab, (lab$pop == "CD4")),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1)) +
  ggtitle("CD4+ in laboratory mice")
### lab Th1
ggplot(subset(lab, (lab$pop == "Th1")),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))+
  ggtitle("Th1 in laboratory mice")

### lab Dividing Th1
ggplot(subset(lab, (lab$pop == "Div_Th1")),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text.x = element_blank(),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1)) +
  ggtitle("Actively dividing Th1 in laboratory mice")

### CD4 IFN+
ggplot(subset(lab, (lab$pop == "IFNy_CD4")),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))

#### CD8
ggplot(subset(lab, (lab$pop == "CD8")),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1)) +
  ggtitle("CD8+ in laboratory mice")
#### activated CD8

ggplot(subset(lab, (lab$pop == "Act_CD8")),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1)) +
  ggtitle("T-bet+ (activated) CD8+ in laboratory mice")

ggplot(subset(lab, (lab$pop == "Div_Act_CD8")),
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))+
  ggtitle("Dividing T-bet+ (activated) CD8+ in laboratory mice")

# lab regulatory
reg <- filter(lab, lab$pop == "Treg")
reg1 <- filter(lab, lab$pop == "Div_Treg")
reg2 <- filter(lab, lab$pop == "Th17")
reg3 <- filter(lab, lab$pop == "Treg17")
reg4 <- filter(lab, lab$pop == "Div_Th17")
reg5 <- filter(lab, lab$pop == "IL17A_CD4")

reg <- rbind(reg, reg1)
reg <- rbind(reg,reg2)
reg <- rbind(reg,reg3)
reg <- rbind(reg,reg4)
reg <- rbind(reg, reg5)

reg <- reg[!(reg$pop=="Treg"& reg$counts > 20),]
reg <- reg[!(reg$pop=="Th17"& reg$counts > 3),]

ggplot(reg,
       aes(x = Eim_MC , y = counts, color = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =0.9) +
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "origin") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())

########################### lab
########## IFNy CEWE vs IFNy+ CD8 lab
Lab_IFNCD8 <- dplyr::select(FACScombine, EH_ID, pop, counts, EXP, Eim_MC)
Lab_IFNCD8 <- subset(Lab_IFNCD8, Lab_IFNCD8$EXP== "lab")
Lab_IFNCD8 <- subset(Lab_IFNCD8, Lab_IFNCD8$pop== "IFNy_CD8")
IFN_CD8 <- merge(Lab_IFNCD8, IFN)
IFN_CD8$Eim_MC <- as.character(IFN_CD8$Eim_MC)
IFN_CD8$Eim_MC[IFN_CD8$Eim_MC == "pos"] <- "positive"
IFN_CD8$Eim_MC[IFN_CD8$Eim_MC == "neg"] <- "negative"

ggplot(IFN_CD8,
       aes(x = IFNy_CEWE, y = counts, color = Eim_MC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "IFN-g (ng/mL)", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))
  # ggtitle("IFN-g+ abundance correlates to increased IFN-g+ CD8+ in laboratory mice")
########## IFNy CEWE vs IFNy+ CD4 lab
Lab_IFNCD4 <- dplyr::select(FACScombine, EH_ID, pop, counts, EXP, Eim_MC)
Lab_IFNCD4 <- subset(Lab_IFNCD4, Lab_IFNCD4$EXP== "lab")
Lab_IFNCD4 <- subset(Lab_IFNCD4, Lab_IFNCD4$pop== "IFNy_CD4")
IFN_CD4 <- merge(Lab_IFNCD4, IFN)
IFN_CD4$Eim_MC <- as.character(IFN_CD4$Eim_MC)
IFN_CD4$Eim_MC[IFN_CD4$Eim_MC == "pos"] <- "positive"
IFN_CD4$Eim_MC[IFN_CD4$Eim_MC == "neg"] <- "negative"

ggplot(IFN_CD4,
       aes(x = IFNy_CEWE, y = counts, color = Eim_MC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "IFN-g (ng/mL)", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))
##################################################################################################################################
######################################################################## IFNy CEWE vs IFNy+ CD4 wild
Wild_IFNCD4 <- dplyr::select(FACScombine, EH_ID, pop, counts, EXP)
Wild_IFNCD4 <- subset(Wild_IFNCD4, Wild_IFNCD4$EXP== "wild")
Wild_IFNCD4 <- subset(Wild_IFNCD4, Wild_IFNCD4$pop== "IFNy_CD4")
IFN_HZCD4 <- merge(Wild_IFNCD4, IFN_HZ)

ggplot(IFN_HZCD4,
       aes(x = IFNy_CEWE, y = counts)) +
  geom_point() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))+
  ggtitle("IFN-g+ CD4+ in wild mice")
########## IFNy CEWE vs IFNy+ CD8 wild
Wild_IFNCD8 <- dplyr::select(FACScombine, EH_ID, pop, counts, EXP)
Wild_IFNCD8 <- subset(Wild_IFNCD8, Wild_IFNCD8$EXP== "wild")
Wild_IFNCD8 <- subset(Wild_IFNCD8, Wild_IFNCD8$pop== "IFNy_CD8")
IFN_HZCD8 <- merge(Wild_IFNCD8, IFN_HZ)

ggplot(IFN_HZCD8,
       aes(x = IFNy_CEWE, y = counts)) +
  geom_point() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  #facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "infection", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.1))+
  ggtitle("IFN-g+ CD8+ in wild mice")

#################################################################################################################################
# let's look at genes
# remove outlier 0287
lab <- lab[!(lab$EH_ID=="LM_0287"),]
lab <- lab[!(lab$EH_ID=="LM_0294"),]
lab <- lab[!(lab$EH_ID=="LM_0275"),]

ggplot(lab, aes(x = delta, y = NE, color = Eim_MC)) +
  geom_point() +
  #stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_smooth(method = "lm") +
  facet_grid(Target~Eim_MC, scales = "free_y") +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.position = "none",
        legend.text=element_blank(),
        legend.title = element_blank())
        #legend.position = c(0.91, 0.1))+
  #ggtitle("delta NE lab")

ggplot(wild, aes(x = delta, y = NE, color = Eim_MC)) +
  geom_point() +
  #stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_smooth(method = "lm") +
  facet_grid(Target~Eim_MC, scales = "free_y") +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())
        #legend.position = c(0.91, 0.1))+
  #ggtitle("delta NE wild")

##################### check this on cell sorting
ggplot(lab, aes(x = delta, y = counts, color = Eim_MC)) +
  geom_point() +
  #stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_smooth(method = "lm") +
  facet_wrap(~pop, scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())
#legend.position = c(0.91, 0.1))+
#ggtitle("delta NE wild")

#### and genes
ggplot(subset(lab, lab$Target == "IL.12"), aes(x = NE, y = counts, color = Eim_MC)) +
  geom_point() +
  #stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_smooth(method = "lm") +
  facet_grid(Target~pop, scales = "free_x") +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())
#legend.position = c(0.91, 0.1))+
#ggtitle("delta NE wild")

################################## lab Th1 and Th17 IL-12 decline?
labTh1 <- filter(lab, pop == "Th1")
labTh1div <- filter(lab, pop == "Div_Th1")
labTh1 <- full_join(labTh1, labTh1div)
labTh17 <- filter(lab, pop == "Th17")
labTh17div <- filter(lab, pop == "Div_Th17")
labTh1 <- full_join(labTh1, labTh17div)
labTh1 <- full_join(labTh1, labTh17)
CXCR3 <- filter(labTh1, Target == "CXCR3")
IL.12 <- filter(labTh1, Target == "IL.12")
IRG6 <- filter(labTh1, Target == "IRG6")
CXCR3.1 <- filter(lab, Target == "CXCR3")
IL.12.1 <- filter(lab, Target =="IL.12")
wildIL12 <- filter(wild, Target == "IL.12")
IL12.1 <- filter(lab, Target == "IL.12")
IRG6.1 <- filter(lab, Target == "IRG6")
wildIRG6 <- filter(wild, Target == "IRG6")

ggplot(wildIRG6, aes(x = NE, y = counts)) +
  geom_point() +
  #stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_smooth(method = "lm") +
  facet_wrap(~pop, scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())

ggplot(IL.12, aes(x = NE, y = counts, color = Eim_MC)) +
  geom_point() +
  #stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_smooth(method = "lm") +
  facet_wrap(~pop, scales = "free_y") +
  labs(y="% of populations", x = "deltaCT = Target - HKG", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())

ggplot(IRG6.1, aes(x = NE, y = counts, color = Eim_MC)) +
  geom_point() +
  #stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_smooth(method = "lm") +
  facet_wrap(~pop, scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())
