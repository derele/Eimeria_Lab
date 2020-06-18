# because of large number of conflicting packages, use "package::function" when necessary
library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(readxl)
library(dplyr)

# guess this one will have to be hard coded
FACSraw1 <- read_xlsx("~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_Lubomir_FACS_data.xlsx", sheet = 1)
FACSraw2 <- read_xlsx("~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_Lubomir_FACS_data.xlsx", sheet = 2)
FACSraw3 <- read_xlsx("~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_Lubomir_FACS_data.xlsx", sheet = 3)
FACSraw4 <- read_xlsx("~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_Lubomir_FACS_data.xlsx", sheet = 4)
#remove NA rows and rubbish
FACSraw1 <- FACSraw1[-c(52:53),]
FACSraw2 <- FACSraw2[-c(52:53),]
FACSraw3 <- FACSraw3[-c(28:35),]
FACSraw3 <- FACSraw3[-c(55:56),]
FACSraw4 <- FACSraw4[-c(55:60),]


# extract sample names and position 
FACSraw1$EH_ID <-gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM_02\\2", FACSraw1$Sample)
FACSraw1$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", FACSraw1$Sample)
FACSraw2$EH_ID <-gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM_02\\2", FACSraw2$Sample)
FACSraw2$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", FACSraw2$Sample)
FACSraw3$EH_ID <-gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM_02\\2", FACSraw3$Sample)
FACSraw3$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", FACSraw3$Sample)
FACSraw4$EH_ID <-gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM_02\\2", FACSraw4$Sample)
FACSraw4$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", FACSraw4$Sample)
# remove that strage Sample column
FACSraw1 <- FACSraw1[,-c(1)]
FACSraw2 <- FACSraw2[,-c(1)]
FACSraw3 <- FACSraw3[,-c(1)]
FACSraw4 <- FACSraw4[,-c(1)]
# combine into one and remove wrong sample (Hongwei said)
FACS1 <- full_join(FACSraw1, FACSraw2)
FACS2 <- full_join(FACSraw3, FACSraw4)
FACS <- full_join(FACS1,FACS2)
FACS <- FACS[!FACS$EH_ID%in%"LM_0293",] 

FACS <- FACS %>%
  group_by(EH_ID) %>%
  summarise_at(vars(-Position), funs(mean(., na.rm=TRUE)))

#####################################################################################################################################
#introduce parasitological data
E7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_COMPLETE.csv"))
E7 <- select(E7, EH_ID, labels, dpi, Strain, HybridStatus, primary, challenge, weight, Wchange, OPG, infHistory, IFNy_FEC,
             IFNy_CEWE, Caecum, delta, CXCR3, IRG6, IL.12)
#merge FACS with para data (maybe replace with less stringent join)
E7 <- merge(E7[E7$dpi%in%8,], FACS, by = "EH_ID")
#create R and .csv friendly column names
colnames(E7)[19]<- "CD4"
colnames(E7)[20]<- "Treg"
colnames(E7)[21]<- "Div_Treg"
colnames(E7)[22]<- "Treg17"
colnames(E7)[24]<- "Th1"
colnames(E7)[25]<- "Div_Th1"
colnames(E7)[26]<- "Th17"
colnames(E7)[27]<- "Div_Th17"
colnames(E7)[28]<- "CD8"
colnames(E7)[29]<- "Act_CD8"
colnames(E7)[30]<- "Div_Act_CD8"
colnames(E7)[31]<- "IFNy_CD4"
colnames(E7)[32]<- "IL17A_CD4"
colnames(E7)[33]<- "IFNy_CD8"
colnames(E7)[23]<- "Treg_prop"

write.csv(E7, "~/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACS_complete.csv")



####### graph out IFN CEWE and IFN pivotal cells #################################################################################

ggplot(subset(E7, !is.na(E7$IFNy_CEWE)), 
       aes(x = TcCD8p, y = IFNy_CEWE, color = infHistory)) +
  geom_jitter() +
  #geom_boxplot() +
  facet_wrap("Position") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("CTLs and IFNy CEWE")

ggplot(subset(E7, !is.na(E7$IFNy_CEWE)), 
       aes(x = Th1IFNgp_in_CD4p, y = IFNy_CEWE, color = infHistory)) +
  geom_jitter() +
  #geom_boxplot() +
  facet_wrap("Position") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Th1 IFN-y+ CD4+ and IFNy CEWE")

ggplot(subset(E7, !is.na(E7$IFNy_CEWE)), 
       aes(x = Tc1IFNgp_in_CD8p , y = IFNy_CEWE, color = infHistory)) +
  geom_jitter() +
  #geom_boxplot() +
  facet_wrap("Position") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("CTL IFN-y+ CD8+ and IFNy CEWE")


##select cell population names (now using .cells to calculate with actual cell populations)
facs.measure.cols <- c("ThCD4p", "TcCD8p", "Th1IFNgp_in_CD4p", "Th17IL17Ap_in_CD4p", 
                       "Tc1IFNgp_in_CD8p", "Treg_Foxp3_in_CD4p", "Dividing_Ki67p_in_Foxp3p", 
                       "RORgtp_in_Foxp3p", "ThCD4p_Foxp3n", "Th1Tbetp_in_CD4pFoxp3n", "Dividing_Ki67p_in_Tbetp", 
                       "Th17RORgp_in_CD4pFoxp3n", "Dividing_Ki67p_in_RORgtp")

#test for normality
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr") use for every cell line
ggdensity(E7$ThCD4p.cells, 
          main = "Density plot of ThCD4p cells",
          xlab = "population counts")

## tabulate  medians for different infection histories and anterior vs posterior
## create list of cell populations summaries infection strains
cell.medians <- lapply(facs.measure.cols, function (x){
  tapply(E7[, x], list(E7$infHistory, as.factor(E7$Position)), median)
})
names(cell.medians) <- facs.measure.cols
cell.medians

cell.means <- lapply(facs.measure.cols, function (x){
  tapply(E7[, x], list(E7$infHistory, as.factor(E7$Position)), mean)
})
names(cell.means) <- facs.measure.cols
cell.means

#set OPG NAs to 0, add  and clean intensity data
E7[["OPG"]][is.na(E7[["OPG"]])] <- 0
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
  mutate(delta = eimeria - mouse) %>% 
  dplyr::select(Name,delta)
colnames(E7_inf)[1]<- "EH_ID"
E7 <- merge(E7, E7_inf)

# graph with intensity
ggplot(E7, aes(delta, cell_counts, color = infHistory)) +
  geom_point() +
  facet_wrap("Strain")

# subtract POS from ANT (as ANT is the internal control)  (preserves direction)
E7clean <- dplyr::select(E7, EH_ID, dpi, Strain, HybridStatus, primary, challenge, weight, fecweight, OPG, cell_counts, 
                         ThCD4p, TcCD8p, Th1IFNgp_in_CD4p, Th17IL17Ap_in_CD4p, Tc1IFNgp_in_CD8p, Treg_Foxp3_in_CD4p,
                         Dividing_Ki67p_in_Foxp3p, RORgtp_in_Foxp3p, Th1Tbetp_in_CD4pFoxp3n, Dividing_Ki67p_in_Tbetp,
                         Th17RORgp_in_CD4pFoxp3n, Dividing_Ki67p_in_RORgtp, Position, infHistory, delta) 

###################################### BIG CHUNK TO BE AUTOMATED ################
E7Position <- E7clean %>% 
  dcast(EH_ID ~ Position, value.var = "cell_counts", fill = 0) %>% 
  mutate(POSminusANT = Posterior - Anterior) %>% 
  dplyr::select(EH_ID, POSminusANT)

E7ThCD4p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "ThCD4p", fill = 0) %>%
  mutate(ThCD4p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, ThCD4p)

E7TcCD8p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "TcCD8p", fill = 0) %>%
  mutate(TcCD8p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, TcCD8p)

E7Th1IFNgp_in_CD4p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Th1IFNgp_in_CD4p", fill = 0) %>%
  mutate(Th1IFNgp_in_CD4p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Th1IFNgp_in_CD4p)

E7Th17IL17Ap_in_CD4p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Th17IL17Ap_in_CD4p", fill = 0) %>%
  mutate(Th17IL17Ap_in_CD4p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Th17IL17Ap_in_CD4p)

E7Tc1IFNgp_in_CD8p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Tc1IFNgp_in_CD8p", fill = 0) %>%
  mutate(Tc1IFNgp_in_CD8p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Tc1IFNgp_in_CD8p)

E7Treg_Foxp3_in_CD4p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Treg_Foxp3_in_CD4p", fill = 0) %>%
  mutate(Treg_Foxp3_in_CD4p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Treg_Foxp3_in_CD4p)

E7Dividing_Ki67p_in_Foxp3p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Dividing_Ki67p_in_Foxp3p", fill = 0) %>%
  mutate(Dividing_Ki67p_in_Foxp3p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Dividing_Ki67p_in_Foxp3p)

E7RORgtp_in_Foxp3p <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "RORgtp_in_Foxp3p", fill = 0) %>%
  mutate(RORgtp_in_Foxp3p = Posterior - Anterior) %>%
  dplyr::select(EH_ID, RORgtp_in_Foxp3p)

# doesn't exist without conversion script
# E7ThCD4p_Foxp3n <- E7clean %>%
#   dcast(EH_ID ~ Position, value.var = "ThCD4p_Foxp3n", fill = 0) %>%
#   mutate(ThCD4p_Foxp3n = Posterior - Anterior) %>%
#   dplyr::select(EH_ID, ThCD4p_Foxp3n)

E7Th1Tbetp_in_CD4pFoxp3n <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Th1Tbetp_in_CD4pFoxp3n", fill = 0) %>%
  mutate(Th1Tbetp_in_CD4pFoxp3n = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Th1Tbetp_in_CD4pFoxp3n)

E7Dividing_Ki67p_in_Tbetp <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Dividing_Ki67p_in_Tbetp", fill = 0) %>%
  mutate(Dividing_Ki67p_in_Tbetp = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Dividing_Ki67p_in_Tbetp)

E7Th17RORgp_in_CD4pFoxp3n <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Th17RORgp_in_CD4pFoxp3n", fill = 0) %>%
  mutate(Th17RORgp_in_CD4pFoxp3n = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Th17RORgp_in_CD4pFoxp3n)

E7Dividing_Ki67p_in_RORgtp <- E7clean %>%
  dcast(EH_ID ~ Position, value.var = "Dividing_Ki67p_in_RORgtp", fill = 0) %>%
  mutate(Dividing_Ki67p_in_RORgtp = Posterior - Anterior) %>%
  dplyr::select(EH_ID, Dividing_Ki67p_in_RORgtp)

E7cells <- merge(E7Position, E7ThCD4p, by = "EH_ID")
E7cells <- merge(E7cells,E7TcCD8p, by = "EH_ID")
E7cells <- merge(E7cells,E7Th1IFNgp_in_CD4p, by = "EH_ID")
E7cells <- merge(E7cells,E7Th17IL17Ap_in_CD4p, by = "EH_ID")
E7cells <- merge(E7cells,E7Tc1IFNgp_in_CD8p, by = "EH_ID")
E7cells <- merge(E7cells,E7Treg_Foxp3_in_CD4p, by = "EH_ID")
E7cells <- merge(E7cells,E7Dividing_Ki67p_in_Foxp3p, by = "EH_ID")
E7cells <- merge(E7cells,E7RORgtp_in_Foxp3p, by = "EH_ID")
E7cells <- merge(E7cells,E7Th1Tbetp_in_CD4pFoxp3n, by = "EH_ID")
E7cells <- merge(E7cells,E7Dividing_Ki67p_in_Tbetp, by = "EH_ID")
E7cells <- merge(E7cells,E7Th17RORgp_in_CD4pFoxp3n, by = "EH_ID")
E7cells <- merge(E7cells,E7Dividing_Ki67p_in_RORgtp, by = "EH_ID")
###############################################################################################

# #weight to cell count ratio NO IDEA WHAT IM DOING
# E7$cell_counts <- as.numeric(E7$cell_counts)
# E7weight.lm <- lm(formula =  Wchange ~  weight * cell_counts, data = E7)
# lapply(E7weight.lm, summary)
# summary(E7weight.lm)
# interplot(m = E7weight.lm, var1 = "cell_counts", var2 = "weight")
# 
# E7Strain.lm <- lm(formula =  weight ~  Strain * cell_counts, data = E7)
# lapply(E7Strain.lm, summary)
# summary(E7Strain.lm)
# interplot(m = E7Strain.lm, var1 = "cell_counts", var2 = "Strain")

######################## gap for sake of graphing unrelated to previous sections
# try with Tc1IFNgp_in_CD8p, Th1Tbetp_in_CD4pFoxp3n, Th1IFNgp_in_CD4p, TcCD8p
ggplot(E7, 
       aes(x = Tc1IFNgp_in_CD8p, y = IRG6, color = infHistory)) +
  geom_jitter() +
  #geom_boxplot() +
  facet_wrap("challenge") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")


###########################

#check distribution infHistory
plotCells.inf <- function (col){
  ggplot(E7, aes(infHistory, get(col))) +
    geom_boxplot() +
    geom_jitter(width=0.2) +
    facet_wrap(~Position) +
    ggtitle(col)
}

facs_boxplots.inf <- lapply(facs.measure.cols, plotCells.inf)
names(facs_boxplots.inf) <-  facs.measure.cols

for(i in seq_along(facs_boxplots.inf)){
  pdf(paste0(names(facs_boxplots.inf)[[i]], ".inf.pdf"))
  plot(facs_boxplots.inf[[i]])
  dev.off()
}

## #check distribution Position
plotCells.position<- function (col){
  ggplot(E7, aes(Position, get(col))) +
    geom_boxplot() +
    geom_jitter(width=0.2) +
    facet_wrap(~infHistory) +
    ggtitle(col)
}

facs_boxplots.position <- lapply(facs.measure.cols, plotCells.position)
names(facs_boxplots.position) <-  facs.measure.cols

for(i in seq_along(facs_boxplots.position)){
  pdf(paste0(names(facs_boxplots.position)[[i]], ".position.pdf"))
  plot(facs_boxplots.position[[i]])
  dev.off()
}

### raw counts are modeled either as poisson or negative binomial in
### either case one could use the overall count (cell_counts) as
### "offset" to specify the "duration of observation" (normally
### offsets are used as a ratio, counto over time). I tried that, but
### then figured out that I know too little about how to interprete
### counts... expecially because the overall cell numbers are varying
### SO MUCH that this changes the results completely!!!

# model interaction of cell populations with primary and secondary infection + constant position direction (PRIMARY : SECONDARY + POSITION)
mods.l <- lapply(facs.measure.cols, function (x) {
    lm(get(x) ~ (primary * challenge) + Position,
        data=E7)
})
names(mods.l) <- facs.measure.cols
lapply(mods.l, summary)

for(i in seq_along(facs.measure.cols)){
    eff <- ggpredict(mods.l[[i]], terms=c("primary", "challenge", "Position"))
    plot <-  plot(eff, rawdata=TRUE) +
        scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
        ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
    pdf(paste0(facs.measure.cols[[i]], ".priXcha+pos.pdf"))
    print(plot)
    dev.off()
}

# model interaction of cell populations with primary, secondary infection and position (PRIMARY : SECONDARY : POSITION)
mods.i <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ primary * challenge * Position,
     data=E7)
})
names(mods.i) <- facs.measure.cols
lapply(mods.i, summary)

for(i in seq_along(facs.measure.cols)){
  eff <- ggpredict(mods.i[[i]], terms=c("primary", "challenge", "Position"))
  plot <-  plot(eff, rawdata=TRUE) +
    scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
    ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
  pdf(paste0(facs.measure.cols[[i]], ".priXchaXpos.pdf"))
  print(plot)
  dev.off()
}

# comparison of models
lapply(seq_along(mods.i), function(i) anova(mods.i[[i]], mods.l[[i]]))

#check the model when using HI categories as well

modsHY.l <- lapply(facs.measure.cols, function (x) {
    lm(get(x) ~ (primary * challenge) + Position + HybridStatus,
        data=E7)
})

names(modsHY.l) <- facs.measure.cols

lapply(modsHY.l, summary)

# comparison of all 3 models
lapply(seq_along(mods.i), function(i) anova(mods.i[[i]], mods.l[[i]], modsHY.l[[i]]))
## And WOW (I reall wrote the above A PRIORY, otherwise... mayor
## fishing excursion ;-)...), but Tc1IFNgp_in_CD8p are lower in
## HYBRIDS look at THIS!!
summary(modsHY.l[["Tc1IFNgp_in_CD8p"]])

WOW <- ggpredict(modsHY.l[["Tc1IFNgp_in_CD8p"]],
                 terms=c("primary", "challenge", "HybridStatus"))

summary(modsHY.l[["Tc1IFNgp_in_CD8p"]])

pdf("WINNER_Tc1IFNgp_in_CD8p.effects.pdf")
plot(WOW)
dev.off()

# graph hybrid status
for(i in seq_along(facs.measure.cols)){
  eff <- ggpredict(modsHY.l[[i]], terms=c("primary", "challenge", "HybridStatus"))
  plot <-  plot(eff, rawdata=TRUE) +
    scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
    ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
  pdf(paste0(facs.measure.cols[[i]], ".priXchaXposHYB.pdf"))
  print(plot)
  dev.off()
}
## Now... I fooled myself a bit to that enthusiasm, as I expected
## "HybridStatusoutbred hybrids" to be ... well ... hybrids. Turns out
## this are the within subspecies outbreds. Let's do some PostHoc
## comparison. 

summary(glht(modsHY.l[["Tc1IFNgp_in_CD8p"]], mcp(HybridStatus="Tukey")))

## nothing too shocking here, just that "outbred hybrids" have a trend
## towards lower cell proportions compared to "inter subsp. hybrids"

# ---------------------------------------------------------- Make connections between facets of models--------
E7$primary <- as.character(E7$primary)
E7$challenge <- as.character(E7$challenge)
E7$primary[E7$primary == "E64"] <- "E. ferrisi"
E7$primary[E7$primary == "E88"] <- "E. falciformis"
E7$challenge[E7$challenge == "E64"] <- "E. ferrisi"
E7$challenge[E7$challenge == "E88"] <- "E. falciformis"


mods.l <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ (primary * challenge) + Position,
     data=E7)
})
names(mods.l) <- facs.measure.cols
lapply(mods.l, summary)

for(i in seq_along(facs.measure.cols)){
  eff <- ggpredict(mods.l[[i]], terms=c("primary", "challenge", "Position"))
  plot <-  plot(eff, rawdata=TRUE) +
    theme(axis.text=element_text(size=12, face = "bold"), 
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          strip.text.x = element_text(size = 14, face = "bold"),
          legend.text=element_text(size=12, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 20, face = "bold")) +
    scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
    ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
  pdf(paste0(facs.measure.cols[[i]], ".priXcha+pos.pdf"))
  print(plot)
  dev.off()
}

