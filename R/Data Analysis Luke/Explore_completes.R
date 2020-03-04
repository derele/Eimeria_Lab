library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(dplyr)
library(reshape2)

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

######## rope in FACS stuff for E7 too ###########################
FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_complete.csv"))
FACS$X <- NULL
FACS.long <- dplyr::select(FACS, EH_ID, Position, infHistory, CD4, Treg, Div_Treg, Treg17, Treg_prop, Th1, Div_Th1,
                           Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)
FACS.long <- dplyr::distinct(FACS)
# transform into long
FACS.long <- melt(FACS.long,
             direction = "long",
             varying = list(names(FACS.long)[4:18]),
             v.names = "cell.pop",
             na.rm = T, value.name = "counts", 
             id.vars = c("EH_ID", "Position", "infHistory"))
FACS.long <- na.omit(FACS.long)
names(FACS.long)[names(FACS.long) == "variable"] <- "pop"

##### make a summary for comparing with wild (no Pos or Ant difference)
FACSt <- FACS.long %>% dplyr::group_by(EH_ID, pop, infHistory) %>% dplyr::summarise(counts = mean(counts, na.rm = T))
write.csv(FACSt, "~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACSt.csv")
# merge with genes
immuno <- merge(FACS.long, genes)
inf <- dplyr::select(complete, EH_ID, challenge, primary, infHistory)
inf <- distinct(inf)
immuno <- merge(immuno, inf)
# add IFN ELISA results
IFN <- dplyr::select(complete, EH_ID, IFNy_FEC, IFNy_CEWE)
IFN <- distinct(IFN)
immuno <- merge(immuno, IFN)
immuno <- distinct(immuno)

##################### pure graphing from here, any general code above ####################################################
################## Wchange graphs ###########################################################################
ggplot(subset(complete, !is.na(complete$primary)), 
       aes(x = dpi, y = Wchange, color = primary, group = primary)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("EXP", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("WchangeXdpi_by_EXP_primary")

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

ggplot(subset(genes, EXP == "P3"), 
       aes(x = Eim_MC, y = NE, color = Eim_MC, group = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Experiments_gene_expression")

ggplot(subset(genes, EXP == "E7"), 
       aes(x = Eim_MC, y = NE, color = Eim_MC,  group = Eim_MC)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Experiments_gene_expression")

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

ggplot(P3, aes(wloss, IFNy_CEWE)) +
  geom_point()

P3_88 <- P3[P3$challenge == "E88",]
P3_64 <- P3[P3$challenge == "E64",]
P3_UNI <- P3[P3$challenge == "UNI",]

mod88 <- lm(wloss ~ IFNy_CEWE, data = P3_88)
mod64 <- lm(wloss ~ IFNy_CEWE, data = P3_64)
modUNI <- lm(wloss ~ IFNy_CEWE, data = P3_UNI)

ggplot(P3_88, aes(OPG, IFNy_CEWE)) +
  geom_abline(data = mod88) + 
  geom_point() 

ggplot(P3_64, aes(OPG, IFNy_CEWE)) +
  geom_abline(data = mod64) +
  geom_point()

ggplot(P3_UNI, aes(OPG, IFNy_CEWE)) +
  geom_abline(data = modUNI) +
  geom_point()

ggplot(P3, aes(y = IFNy_CEWE, x = OPG)) +
  geom_text_repel(aes(label=delta)) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

IRG6 <- E7 %>% drop_na(IRG6)

ggplot(IRG6, aes(y = IRG6, x = delta)) +
  geom_point() +
  facet_wrap("challenge", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")






ggplot(subset(immuno, challenge == "E64"), aes(y = IFNy_CEWE, x = counts, color = Position)) +
  geom_point() +
  # ylim(2, -17) +
  facet_wrap("pop", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")


###### IFN CEWE vs cell populations
ggplot(immuno, aes(y = counts, x = IFNy_CEWE, color = Eim_MC)) +
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
  geom_smooth(method = "lm") +
  xlim(2, -17) +
  facet_grid(Target~pop, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Gene expression vs cell populations")
############ ferrisi only


