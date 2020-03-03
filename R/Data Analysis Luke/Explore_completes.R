library(tidyverse)
library(ggplot2)
library(Rmisc)
library(httr)
library(RCurl)
library(dplyr)

complete <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_P3_E6_complete.csv"))
# make negative MCs into NAs in a new column
complete$delta_clean <- complete$delta
complete <- mutate(complete, delta_clean = ifelse(Eim_MC == "neg", -20, delta_clean))
complete$dpi <- as.factor(complete$dpi)
complete$X <- NULL

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
########## gene expression
genes <- dplyr::select(complete, EH_ID, CXCR3, IL.12, IRG6, EXP)
Eim <- dplyr::select(complete, Eim_MC, EH_ID)
Eim <- dplyr::distinct(Eim)
genes <- dplyr::distinct(genes)
genes <- merge(genes, Eim, by.y = "EH_ID")
genes <- dplyr::distinct(genes)
# tranform into long
# genes <- 
genes <- melt(genes,
     direction = "long",
     varying = list(names(genes)[2:4]),
     v.names = "NE",
     na.rm = T, value.name = "NE", 
     id.vars = c("EH_ID", "EXP", "Eim_MC"))
genes <- na.omit(genes)
names(genes)[names(genes) == "variable"] <- "Target"


ggplot(genes, 
       aes(x = Target, y = NE, color = Eim_MC)) +
  geom_jitter() +
  geom_boxplot() +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Experiments_gene_expression")

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

######## rope in FACS stuff for E7 too ###########################
FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_complete.csv"))
FACS <- dplyr::select(E7, EH_ID, ThCD4p, TcCD8p, Th1IFNgp_in_CD4p, Th17IL17Ap_in_CD4p, Tc1IFNgp_in_CD8p, Treg_Foxp3_in_CD4p,
                      Dividing_Ki67p_in_Foxp3p, RORgtp_in_Foxp3p, Th1Tbetp_in_CD4pFoxp3n, Dividing_Ki67p_in_Tbetp,
                      Th17RORgp_in_CD4pFoxp3n, Dividing_Ki67p_in_RORgtp, Position, infHistory)

FACS <- dplyr::distinct(FACS)




# tranform into long

FACS <- melt(FACS,
              direction = "long",
              varying = list(names(FACS)[2:13]),
              v.names = "cell.pop",
              na.rm = T, value.name = "counts", 
              id.vars = c("EH_ID", "Position", "infHistory"))
FACS <- na.omit(FACS)
names(FACS)[names(FACS) == "variable"] <- "pop"

##### make a summary for comparing with wild (no Pos or Ant difference)
FACSt <- FACS %>% dplyr::group_by(EH_ID, pop, infHistory) %>% dplyr::summarise(counts = mean(counts, na.rm = T))
write.csv(FACSt, "~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACSt.csv")
# merge with genes
immuno <- merge(FACS, genes)
inf <- dplyr::select(complete, EH_ID, challenge, primary, infHistory)
inf <- distinct(inf)
immuno <- merge(immuno, inf)
immuno1 <- merge(FACS, genes)
immuno1 <- merge(immuno1, inf)
############# CXCR3 graphs ###################################
# CXCR3 is a chemokine receptor that is highly expressed on effector T cells and plays an important role in 
# T cell trafficking and function. CXCR3 is rapidly induced on naïve cells following activation and preferentially 
# remains highly expressed on Th1-type CD4+ T cells and effector CD8+ T cells.

# let's name the populations better
levels(immuno$pop)
levels(immuno$pop) <- c("CD4+", "CTL", "Th1 IFNy+", "TH17 IL17+", "CTL IFNy+", "Treg", "X^ Treg", "Treg17", "Th1 T-bet", 
                        "X^ T-bet", "Th17", "X^ Th17/Treg17")
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

IFN <- dplyr::select(complete, EH_ID, IFNy_FEC, IFNy_CEWE)
immuno <- merge(immuno, IFN)
immuno <- distinct(immuno)

ggplot(immuno, aes(y = NE, x = IFNy_CEWE, color = Eim_MC)) +
  geom_point() +
  # ylim(2, -17) +
  facet_grid(Target~infHistory, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

####################################### try


