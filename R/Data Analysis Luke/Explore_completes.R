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
complete <- mutate(complete, delta_clean = ifelse(Eim_MC == "neg", -20, delta_clean))
complete$dpi <- as.factor(complete$dpi)
complete$X <- NULL
########## make column with E. ferrisi and E. falciformis only
complete$Eimeria.p <- gsub("E64|E139", replacement = "E.ferrisi", complete$primary)
complete$Eimeria.p <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria.p))

complete$Eimeria.c <- gsub("E64|E139", replacement = "E.ferrisi", complete$challenge)
complete$Eimeria.c <- paste(gsub("E88|Eflab", replacement = "E.falciformis", complete$Eimeria.c))

########## gene expression
genes <- dplyr::select(complete, EH_ID, CXCR3, IL.12, IRG6, EXP)
Eim <- dplyr::select(complete, Eim_MC, EH_ID)
Eim <- dplyr::distinct(Eim)
genes <- dplyr::distinct(genes)
genes <- merge(genes, Eim, by.y = "EH_ID")
genes <- dplyr::distinct(genes)
# tranform into long
# genes <- 
genes <- reshape2::melt(genes,
                        direction = "long",
                        varying = list(names(genes)[2:4]),
                        v.names = "NE",
                        na.rm = T, value.name = "NE", 
                        id.vars = c("EH_ID", "EXP", "Eim_MC"))
genes <- na.omit(genes)
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
# merge with genes
immuno <- merge(FACS.long, genes, all = T)
inf <- dplyr::select(complete, EH_ID, challenge, primary, infHistory)
inf <- distinct(inf)
immuno <- merge(immuno, inf, all = T)
# add IFN ELISA results
IFN <- dplyr::select(complete, EH_ID, IFNy_FEC, IFNy_CEWE)
IFN <- distinct(IFN)
immuno <- merge(immuno, IFN, all = T)
immuno <- distinct(immuno)


delta <- select(complete, delta, delta_clean, EH_ID)
genes <- merge(delta, genes)
genes <- distinct(genes)


Wch.p <- subset(complete, Eimeria.p == "E.ferrisi")
Wch.p1 <- subset(complete, Eimeria.p == "E.falciformis")
Wch.p <- full_join(Wch.p, Wch.p1)

Wch.c <- subset(complete, Eimeria.c == "E.ferrisi")
Wch.c1 <- subset(complete, Eimeria.c == "E.falciformis")
Wch.c <- full_join(Wch.c, Wch.c1)

Wch.c <- select(Wch.c, Eimeria.c, Eimeria.p, EH_ID)
genes <- merge(genes, Wch.c)
genes <- distinct(genes)
genes <- genes[-c(122, 123, 3),]

immuno <- merge(immuno, Wch.c)
immuno <- distinct(immuno)
immuno1 <- immuno[-c(5908:5910),]


sig <- subset(FACScombine, subset = pop %in% c("CD4", "Div_Treg", "Treg17", "Th17", "Div_Th17", "CD8", "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8"))
sig <- data.frame(sig)
sig$pop <- as.character(sig$pop)
# remove horrible outliers
sig <- sig[-c(1174),]
sig <- sig[-c(1173),]
sig <- sig[-c(634),]

sig1 <- merge(sig, Wch.c, by = "EH_ID")
Pos <- select(immuno, EH_ID, Position)
sig1 <- distinct(sig1)
sig1 <- merge(Pos, sig1, by = "EH_ID")
sig1 <- distinct(sig1)

PosHZ <- select(FACSHZ, Mouse_ID, Position, CD4)
colnames(PosHZ)[1] <- "EH_ID"
FACStHZ <- full_join(PosHZ, FACStHZ)
immuno <- distinct(immuno)

##################### pure graphing from here, any general code above ####################################################
ggplot(subset(immuno, !is.na(immuno$Position)), aes(x = Position, y = counts, color = Position)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter() +
      facet_wrap("pop", scales = "free") +
      labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
      theme(axis.text=element_text(size=12, face = "bold"),
               title = element_text(size = 16, face = "bold"),
               axis.title=element_text(size=14,face="bold"),
               strip.text.x = element_text(size = 14, face = "bold"),
               legend.text=element_text(size=12, face = "bold"),
               legend.title = element_text(size = 12, face = "bold"))+
      ggtitle("Gene expression in wild samples")





ggplot(genes, 
       aes(y = NE, x = delta_clean, color = Eimeria.c)) +
  geom_point() +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("WchangeXdpi_by_EXP_primary")
################## Wchange graphs ###########################################################################
ggplot(subset(complete, !is.na(complete$Eimeria.p)), 
       aes(x = dpi, y = Wchange, color = Eimeria.p, group = Eimeria.p)) +
  geom_point() +
  geom_smooth() +
  #facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("WchangeXdpi_by_EXP_primary")



ggplot(Wch.p, 
       aes(x = dpi, y = Wchange, color = Eimeria.p, group = Eimeria.p)) +
  geom_point() +
  geom_smooth() +
  #facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Weight change per dpi by Eimeria strain")


ggplot(subset(complete, !is.na(complete$challenge)), 
       aes(x = dpi, y = Wchange, color = challenge, group = challenge)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("WchangeXdpi_by_EXP_challenge")



ggplot(Wch.c,
       aes(x = dpi, y = Wchange, color = Eimeria.c, group = Eimeria.c)) +
  geom_point() +
  geom_smooth() +
  #facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Weight change per dpi by Eimeria strain")
##########################################################################################################
##################### OPG graphs
ggplot(subset(complete, !is.na(complete$primary)), 
       aes(x = dpi, y = OPG, color = primary, group = primary)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("OPGXdpi_by_EXP_primary")

ggplot(Wch.p, 
       aes(x = dpi, y = OPG, color = Eimeria.p, group = Eimeria.p)) +
  geom_point() +
  geom_smooth() +
  #facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("OPG per dpi by Eimeria strain")



ggplot(subset(complete, !is.na(complete$challenge)), 
       aes(x = dpi, y = OPG, color = challenge, group = challenge)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("OPGXdpi_by_EXP_challenge")

#####################################################################################################################
############ genes


ggplot(subset(genes, Eim_MC == "pos"), 
       aes(x = Eimeria.c, y = NE, color = Eimeria.c)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("Target", scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "infection strain", colour = "infection strain") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Gene expression difference between Eimeria infections")


ggplot(subset(genes, EXP == "E7"| Eim_MC == "pos"), 
       aes(x = Eimeria.c, y = NE,  group = Eimeria.c)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("Target", scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "infection strain") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Eimeria positive mice gene expression")



ggplot(genes, 
       aes(x = Eim_MC, y = NE, color = Eim_MC, group = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(Target~EXP, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

ggplot(genes, 
       aes(x = Eim_MC, y = NE, color = Eim_MC, group = Eim_MC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("Target", scales = "free") +
  labs(y="deltaCT = Target - HKG", x = "infected", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Gene expression in laboratory samples")
###########################################################
ggplot(subset(complete, !is.na(complete$challenge)), 
       aes(x = challenge, y = CXCR3, color = EXP)) +
  geom_jitter() +
  geom_boxplot() +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("CXCR3_experiment_difference")

ggplot(subset(complete, !is.na(complete$challenge)), 
       aes(x = challenge, y = IRG6, color = EXP)) +
  geom_jitter() +
  geom_boxplot() +
  # facet_wrap("infHistory") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("IRG6_experiment_difference")





##################################################################################################
################### IFN_CEWE 
ggplot(subset(complete, !is.na(complete$IFNy_CEWE)), 
       aes(x = OPG, y = IFNy_CEWE, color = challenge)) +
  geom_jitter() +
  #geom_boxplot() +
  facet_wrap("EXP") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("IFNy_CEWE_X_OPG_in_challenge_EXP")

ggplot(subset(complete, !is.na(complete$IFNy_CEWE)), 
       aes(x = infHistory, y = IFNy_CEWE, color = challenge)) +
  geom_jitter() +
  geom_boxplot() +
  # facet_wrap("challenge", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold", angle = 45), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("IFNy_CEWE_X_OPG_in_challenge_EXP")


cor.test(~IFNy_FEC+IFNy_CEWE, subset(complete, dpi == 8 & EXP %in% "E7"))
### NS!

ggplot(subset(complete, EXP%in%"P3"), 
       aes(y = IFNy_CEWE, x = delta, color = challenge, shape=Eim_MC)) +
  geom_point() +
  facet_wrap(~primary)+
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

ggplot(complete, aes(y = OPG, x = dpi, color = EXP)) +
  geom_jitter(width=0.2) +
  facet_grid(primary~challenge, scales="free_y") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("OPG_x_dpi_in_EXP_by_primaryxchallenge")

ggplot(complete, aes(y = OPG, x = dpi, color = primary, group_by("EH_ID"))) +
  geom_jitter(width=0.2) +
  facet_wrap("EXP", scales="free_y") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("OPG_x_dpi_by_EXP")

ggplot(P3, aes(y = IFNy_CEWE, x = delta, color = OPG)) +
  geom_point() +
  facet_wrap("challenge") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

ggplot(E7, aes(y = IFNy_FEC, x = IFNy_CEWE , shape = Caecum)) +
  geom_point() +
  facet_wrap("infHistory") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

ggplot(P3, aes(y = IFNy_CEWE, x = OPG)) +
  geom_text_repel(aes(label=delta)) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")



ggplot(IRG6, aes(y = IRG6, x = delta)) +
  geom_point() +
  facet_wrap("challenge", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")
################################################### Enter FACS


ggplot(immuno1, aes(x = Eimeria.c, y = counts, color = Eimeria.c)) +
  geom_boxplot() +
  geom_jitter() +
  # ylim(2, -17) +
  facet_wrap("pop", scales = "free") +
  stat_compare_means(aes(label = ..p.signif..)) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

###### IFN CEWE vs cell populations

ggplot(immuno,
       aes(x = IFNy_CEWE , y = counts, color = Position)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Eimeria.c~pop, scales = "free") +
  labs(y="% of populations", x = "IFNy (ng/mL)", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Caecum IFNy effect on cell populations")

######## ferrisi only
ggplot(subset(immuno, challenge == "E64"), aes(y = counts, x = IFNy_CEWE, color = Eim_MC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  # ylim(2, -17) +
  facet_grid(Position~pop, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("IFNy Caecum vs cell populations")

######### falciformis only
ggplot(subset(immuno, challenge == "E88"), aes(y = counts, x = IFNy_CEWE, color = Eim_MC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  # ylim(2, -17) +
  facet_grid(Position~pop, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("IFNy Caecum vs cell populations")

########### genes vs cell populations
ggplot(immuno, aes(y = counts, x = NE, color = Eim_MC)) +
  geom_point() +
  xlim(2, -17) +
  facet_grid(Target~pop, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Gene expression vs cell populations")
############ ferrisi only

########## E7 vs HZ19
FACScombine <- merge(FACScombine, Wch.c, by = "EH_ID")
FACScombine <- distinct(FACScombine)
ggplot(FACScombine,
       aes(x = EXP , y = counts, color = EXP)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "population", colour = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("FACS comparison of wild and lab")

##### keep only significant ones





ggplot(sig,
       aes(x = EXP , y = counts, color = EXP)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "population", colour = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("FACS comparison of wild and lab")


ggplot(FACScombine,
       aes(x = Eimeria.c , y = counts)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "Eimeria species", colour = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Eimeria immune cell populations comparison")

ggplot(sig1,
       aes(x = Position , y = counts, color = EXP)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_grid("pop", scales = "free") +
  labs(y="% of populations", x = "Eimeria species", colour = "origin") +
  theme(axis.text=element_text(size=12, face = "bold"),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Eimeria immune cell populations comparison")





ggplot(subset(immuno, !is.na(immuno$counts)),
       aes(x = IFNy_CEWE , y = counts, color = Eimeria.c)) +
  geom_smooth(method = "lm") +
  geom_point() +
  
  #stat_compare_means(aes(label = ..p.signif..)) +
  
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "IFNy (ng/mL)", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Caecum IFNy effect on cell populations")

# wild

ggplot(FACStHZ, 
       aes(x = Position , y = counts)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(aes(label = ..p.signif..)) +
  facet_wrap("pop", scales = "free") +
  labs(y="% of populations", x = "IFNy (ng/mL)", colour = "infected") +
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Caecum IFNy effect on cell populations")

