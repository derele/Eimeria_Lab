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
complete$Eimeria_p <- gsub("E64|E139", replacement = "E.ferrisi", complete$primary)
complete$Eimeria_p <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria_p))

complete$Eimeria_c <- gsub("E64|E139", replacement = "E.ferrisi", complete$challenge)
complete$Eimeria_c <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria_c))
complete$EXP_type <- "lab"
########## start creating lab_immuno from complete to avoid NAs in wrong places
##### first delta then genes then FACS then para
### 1) delta
lab_delta <- dplyr::select(complete, EH_ID, delta, delta_clean, Eim_MC, EXP_type, challenge)
lab_delta <- dplyr::distinct(lab_delta)
### 2) genes
lab_genes <- dplyr::select(complete, EH_ID, CXCR3, IL.12, IRG6)
lab_genes <- dplyr::distinct(lab_genes)
#lab_genes <- subset(lab_genes, !is.na(lab_genes$CXCR3))
# tranform into long
lab_genes_long <- reshape2::melt(lab_genes,

                        direction = "long",
                        varying = list(names(lab_genes)[2:4]),
                        v.names = "NE",
                        na.rm = T, value.name = "NE", 
                        id.vars = c("EH_ID"))
names(lab_genes_long)[names(lab_genes_long) == "variable"] <- "Target"
### 3) FACS
lab_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_complete.csv"))
lab_FACS$X <- NULL
lab_FACS <- dplyr::select(lab_FACS, EH_ID, Position, infHistory, CD4, Treg, Div_Treg, Treg17, Treg_prop, Th1, Div_Th1,
                           Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)
lab_FACS <- dplyr::distinct(lab_FACS)
# transform into long
lab_FACS_long <- reshape2::melt(lab_FACS,
             direction = "long",
             varying = list(names(lab_FACS)[19:34]),
             v.names = "cell.pop",
             na.rm = T, value.name = "counts", 
             id.vars = c("EH_ID", "Position", "infHistory"))
names(lab_FACS_long)[names(lab_FACS_long) == "variable"] <- "pop"
# make a summary for comparing with wild (no Pos or Ant difference)
lab_FACS_long <- lab_FACS_long %>% dplyr::group_by(EH_ID, pop, infHistory) %>% dplyr::summarise(counts = mean(counts, na.rm = T))
# write.csv(lab_FACS_long, "~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_lab_FACS_long.csv")

### para data
#lab_para <- dplyr::select(complete, )

### finally join into lab_immuno

lab_immuno <- merge(lab_genes_long, lab_delta)
lab_immuno <- subset(lab_immuno, !is.na(lab_immuno$delta))
lab_immuno <- merge(lab_immuno, lab_FACS_long)


########## start creating wild_immuno 
##### first delta then genes then FACS then para




wild_FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_complete.csv"))
wild_FACS$X <- NULL
wild_FACS <- dplyr::select(wild_FACS, Mouse_ID, Position, CD4, Treg, Div_Treg, Treg17, Treg_prop, Th1, Div_Th1,
                           Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)
wild_FACS <- dplyr::distinct(wild_FACS)
# transform into long
wild_FACS_long <- reshape2::melt(wild_FACS,
                            direction = "long",
                            varying = list(names(wild_FACS)[3:17]),
                            v.names = "cell.pop",
                            na.rm = T, value.name = "counts", 
                            id.vars = c("Mouse_ID", "Position"))
wild_FACS_long <- na.omit(wild_FACS_long)
names(wild_FACS_long)[names(wild_FACS_long) == "variable"] <- "pop"

##### make a summary for comparing with wild (no Pos or Ant difference)
wild_FACS_long <- subset(wild_FACS_long, Position == "mLN")
wild_FACS_long <- wild_FACS_long %>% dplyr::group_by(Mouse_ID, pop) %>% dplyr::summarise(counts = mean(counts, na.rm = T))

# write.csv(wild_FACS_long, "~/Documents/Mouse_Eimeria_Databasing/data/Field_data/HZ19_FACSt.csv")
# add EXP columns, dummy infHistory column and rename Mouse ID to EH ID for sake of merging
lab_FACS_long$EXP_type <- "lab"
wild_FACS_long$EXP_type <- "wild"
wild_FACS_long$infHistory <- NA
colnames(wild_FACS_long)[1] <- "EH_ID"
wild_FACS_long <- data.frame(wild_FACS_long)
lab_FACS_long <- data.frame(lab_FACS_long)
FACS_long <- rbind(wild_FACS_long, lab_FACS_long)



# add IFN ELISA results``````````````````````
lab_IFN <- dplyr::select(complete, EH_ID, IFNy_CEWE, dpi)
lab_IFN <- subset(lab_IFN, !is.na(lab_IFN$IFNy_CEWE))
lab_IFN <- dplyr::select(lab_IFN, EH_ID, IFNy_CEWE)
#

lab_IFN$EXP_type <- "lab"
lab_immuno <- merge(lab_immuno, lab_IFN, all.x = T)
lab_immuno <- distinct(lab_immuno)
lab_immuno <- subset(lab_immuno, !is.na(lab_immuno$pop))

lab_delta <- dplyr::select(complete, delta, delta_clean, EH_ID, Eim_MC, Eimeria_c)
lab_delta <- na.omit(lab_delta)
lab_delta <- data.frame(lab_delta)
# 22.05.2020
immuno <- merge(immuno, lab_delta, all.x = T)
immuno <- distinct(immuno)

immuno$species <- gsub("E64|E139", replacement = "E.ferrisi", immuno$challenge)
immuno$species <- paste(gsub("E88|Eflab", replacement = "E.falciformis", immuno$species))

#### add wild gene expression
wild_genes_long <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-18_gene_expression.csv"))
wild_genes_long$X <- NULL
colnames(wild_genes_long)[1] <- "EH_ID"
colnames(wild_genes_long)[6] <- "Eim_MC"
wild_genes_long$Eim_MC[wild_genes_long$Eim_MC == "TRUE"] <- "positive"
wild_genes_long$Eim_MC[wild_genes_long$Eim_MC == "FALSE"] <- "negative"
wild_genes_long$EXP <- "wild"
immuno$Eim_MC <- as.character(immuno$Eim_MC)
immuno$Eim_MC[immuno$Eim_MC == "pos"] <- "positive"
immuno$Eim_MC[immuno$Eim_MC == "neg"] <- "negative"
lab_genes$Eim_MC <- as.character(lab_genes$Eim_MC)
lab_genes$Eim_MC[lab_genes$Eim_MC == "pos"] <- "positive"
lab_genes$Eim_MC[lab_genes$Eim_MC == "neg"] <- "negative"

# god knows what happens here (investigate later)
immuno <- merge(wild_genes_long, immuno, all.y = T)

#### add IFN CEWE HZ19
wild_IFN <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISAs_complete.csv"))
wild_IFN$X <- NULL
colnames(wild_IFN)[1] <- "EH_ID"
colnames(wild_IFN)[2] <- "IFNy_CEWE"
lab_IFN$EXP <- "E7"
lab_IFN$EXP_type <- "lab"
wild_IFN$EXP <- "HZ19"
wild_IFN$EXP_type <- "wild"

IFN <- rbind(lab_IFN, wild_IFN)
# why?
immuno <- merge(immuno, IFN, all = T)
immuno <- merge(immuno, FACS, all = T)

WxL <- merge(IFN, FACS, all = T)
WxL <- merge(WxL, Eim_MC, all =T)

Wch.p <- subset(Wch.p, !Wch.p$dpi == 0)
Wch.p$OPG[is.na(Wch.p$OPG)] <- 0

lab <- subset(immuno, immuno$EXP == "lab")
wild <- subset(immuno, immuno$EXP == "wild")
lab$species[lab$species == "E.falciformis"] <- "uninfected"

FACS_IFN <- merge(FACScombine, IFNcomplete, all = T)
FACScombine <- merge(FACScombine, Eim, all =T)

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


#####################################################################################################################
############ lab_genes
# wild gene expression
ggplot(HZlab_genes,
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
lab_genes <- lab_genes[!(lab_genes$Target=="IRG6"& lab_genes$NE > 19),]

ggplot(lab_genes, aes(x = Eim_MC, y = NE, color = Eim_MC, group = Eim_MC)) +
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



#### lab_genes vs IFN
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

immunolab_genes <- dplyr::select(immuno, EXP, NE, Target, EH_ID)
immunolab_genes <- distinct(immunolab_genes)

ggplot(subset(immunolab_genes, !is.na(immunolab_genes$EXP)),
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
# let's look at lab_genes
# remove outlier 0287
lab <- lab[!(lab$EH_ID=="LM_0287"),]
lab <- lab[!(lab$EH_ID=="LM_0294"),]
lab <- lab[!(lab$EH_ID=="LM_0275"),]

immuno_compare <- subset(immuno, !immuno$Target == "GBP2")
immuno_compare <- subset(immuno_compare, !immuno_compare$Target == "IL-6")
immuno_compare <- subset(immuno_compare, !is.na(immuno_compare$NE))
ggplot(immuno_compare, aes(x = delta, y = NE, color = EXP)) +
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

ggscatter(immuno, x = "delta", y = "NE", add = "reg.line", color = "EXP") +
  facet_grid(Eim_MC~Target)+
  stat_cor(label.x = -10, label.y = 2) +
  stat_regline_equation(label.x = -10, label.y = 4) + 
  ggtitle("")




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

#### and lab_genes
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

############ qPCR intensity
Efer <- filter(.data = delta, Eimeria.c == "E.ferrisi")

ggplot(Efer, aes(x = Eim_MC, y = delta, color = Eim_MC)) +
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
        legend.title = element_blank())

uni <- filter(.data = delta, Eimeria.c == "UNI")
uni1 <- filter(.data = delta, Eimeria.c == "E.falciformis")
uni <- rbind(uni, uni1)

ggplot(uni, aes(x = Eim_MC, y = delta, color = Eimeria.c)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())

delta2 <- rbind(uni, Efer)

ggplot(delta2, aes(x = Eim_MC, y = delta, color = Eimeria.c)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Eimeria.c)+
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  labs(y="deltaCT = Target - HKG", x = "deltaCT = Mouse - Eimeria", colour = "infection") +
  theme(title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_blank())

################### check IFN abundance vs cell counts
  
ggplot(IFNcomplete, aes(x = EXP, y = IFNy_CEWE)) +
  geom_boxplot() +
  geom_jitter(width = 0.2)

IFN_eim <- merge(IFN, Eim)

ggplot(IFN_eim, aes(x = Eim_MC, y = IFNy_CEWE, color = Eim_MC)) +
  geom_boxplot() +
  stat_compare_means(aes(label = ..p.signif..), size = 8, label.y.npc =1) +
  geom_jitter(width = 0.2)

# T17s vs IFN

IL17A_CD4 <- filter(lab, pop == "IL17A_CD4")
Th17 <- filter(lab, pop == "Th17")
Treg17 <- filter(lab, pop == "Treg17")
Div_Th17 <- filter(lab, pop == "Div_Th17")
T17s <- rbind(IL17A_CD4, Th17)
T17s <- rbind(T17s, Treg17)
T17s <- rbind(T17s, Div_Th17)

ggscatter(T17s, x = "IFNy_CEWE", y = "counts", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~pop)+
  stat_cor(label.x = 50, label.y = 60) +
  stat_regline_equation(label.x = 50, label.y = 63) + 
  ggtitle("IL-17 related cell populations")

# Cd4s vs IFN
CD4 <- filter(lab, pop == "CD4")
Treg <- filter(lab, pop == "Treg")
Div_Treg <- filter(lab, pop == "Div_Treg")
Th1 <- filter(lab, pop == "Th1")
Div_Th1 <- filter(lab, pop == "Div_Th1")
IFNy_CD4 <- filter(lab, pop == "IFNy_CD4")
CD4s <- rbind(CD4, Treg)
CD4s <- rbind(CD4s, Div_Treg)
CD4s <- rbind(CD4s, Th1)
CD4s <- rbind(CD4s, Div_Th1)
CD4s <- rbind(CD4s, IFNy_CD4)

ggscatter(CD4s, x = "IFNy_CEWE", y = "counts", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~pop)+
  stat_cor(label.x = 30, label.y = 80) +
  stat_regline_equation(label.x = 30, label.y = 83) + 
  ggtitle("CD4 related cell populations")

# CD8s vs IFN
CD8 <- filter(lab, pop == "CD8")
Act_CD8 <- filter(lab, pop == "Act_CD8")
Div_Act_CD8 <- filter(lab, pop == "Div_Act_CD8")
IFNy_CD8 <- filter(lab, pop == "IFNy_CD8")
CD8s <- rbind(CD8, Act_CD8)
CD8s <- rbind(CD8s, Div_Act_CD8)
CD8s <- rbind(CD8s, IFNy_CD8)

ggscatter(CD8s, x = "IFNy_CEWE", y = "counts", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~pop)+
  stat_cor(label.x = 30, label.y = 80) +
  stat_regline_equation(label.x = 30, label.y = 83) + 
  ggtitle("CD8+ related cell populations")

################### check IFN abundance vs lab_genes

lab_lab_genes <- filter(immuno, Target == c("CXCR3", "IL.12", "IRG6"))
lab_lab_genes <- filter(lab_lab_genes, EXP == "lab")

#CXCR3
ggscatter(subset(lab_lab_genes, lab_lab_genes$Target == "CXCR3"), x = "IFNy_CEWE", y = "NE", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~Target)+
  stat_cor(label.x = 50, label.y = 0) +
  stat_regline_equation(label.x = 50, label.y = 2) + 
  ggtitle("IFN-y and CXCR3")

#irg6
ggscatter(subset(lab_lab_genes, lab_lab_genes$Target == "IRG6"), x = "IFNy_CEWE", y = "NE", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~Target)+
  stat_cor(label.x = 50, label.y = 0) +
  stat_regline_equation(label.x = 50, label.y = 2) + 
  ggtitle("IFN-y and IRG6")

ggscatter(subset(lab_lab_genes, lab_lab_genes$Target == "IL.12"), x = "IFNy_CEWE", y = "NE", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~Target)+
  stat_cor(label.x = 50, label.y = 0) +
  stat_regline_equation(label.x = 50, label.y = 2) + 
  ggtitle("IFN-y and IL-12")

############# check gene expression vs cell counts

lab1 <- lab[!(lab$pop== "Treg_prop"),]
lab1 <- lab1[!(lab1$pop== "Div_Th1"),]
lab1 <- lab1[!(lab1$pop== "Div_Act_CD8"),]
lab1 <- lab1[!(lab1$pop== "Div_Treg"),]
lab1 <- lab1[!(lab1$pop== "Div_Th17"),]

ggscatter(lab1, x = "counts", y = "NE", add = "reg.line", color = "Eim_MC") +
  facet_grid(pop~Target)+
  stat_cor(label.x = 50, label.y = 30) +
  stat_regline_equation(label.x = 50, label.y = 34) + 
  ggtitle("lab_genes and populations")
# still too big, filter by lab_genes
IRG6_counts <- filter(lab1, Target == "IRG6")
CXCR3_counts <- filter(lab1, Target == "CXCR3")
IL12_counts <- filter(lab1, Target == "IL.12")

ggscatter(IRG6_counts, x = "counts", y = "NE", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~pop)+
  stat_cor(label.x = 0, label.y = 2) +
  stat_regline_equation(label.x = 0, label.y = 4) + 
  ggtitle("IRG6 and cell counts")

ggscatter(CXCR3_counts, x = "counts", y = "NE", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~pop)+
  stat_cor(label.x = 0, label.y = 0) +
  stat_regline_equation(label.x = 0, label.y = 2) + 
  ggtitle("CXCR3 and cell counts")

ggscatter(IL12_counts, x = "counts", y = "NE", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~pop)+
  stat_cor(label.x = 0, label.y = 0) +
  stat_regline_equation(label.x = 0, label.y = 2) + 
  ggtitle("IL-12 and cell counts")



immunopar <- filter(complete, dpi == 7:8)
immunopar <- select(immunopar, EH_ID, Wchange)
immunopar <- distinct(immunopar)
immunopar <- merge(immunopar, immuno, by.y = "EH_ID")
immunopar <- distinct(immunopar)

# check weightloss (weight change at day of dissection against cell counts)
ggscatter(subset(immunopar, immunopar$pop == "Treg17"), x = "Wchange", y = "counts", add = "reg.line", color = "Eim_MC") +
  facet_grid(Eim_MC~pop)+
  stat_cor(label.x = 80, label.y = 1.5) +
  stat_regline_equation(label.x = 80, label.y = 2) + 
  ggtitle("whatever")

# check weightloss vs IFN.
ggscatter(immunopar, y = "Wchange", x = "IFNy_CEWE", add = "reg.line", color = "Eim_MC") +
  facet_wrap(~Eim_MC) +
  stat_cor(label.x = 150, label.y = 110) +
  stat_regline_equation(label.x = 150, label.y = 120) + 
  ggtitle("IFN-y vs wchange")

IFNpoop <- select(complete, delta, dpi,IFNy_CEWE, IFNy_FEC, EH_ID, Wchange, weight, OPG, Eim_MC)
ggplot(IFNpoop, aes(x = delta, y = Wchange)) + 
  geom_point() +
  facet_wrap(~Eim_MC)

########## 
ggscatter(FACS, x = "Treg17", y = "IL17A_CD4", color = "Caecum") +
  facet_wrap(~Caecum)+
  stat_cor(label.x = 0, label.y = 3) +
  stat_regline_equation(label.x = 0, label.y = 3.1) + 
  ggtitle("whatever")
