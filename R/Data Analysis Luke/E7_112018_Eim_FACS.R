# because of large number of conflicting packages, use "package::function" when necessary
library(httr)
library(RCurl)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(lattice)
library(data.table)
library(ggeffects)
library(multcomp)
library(fitdistrplus)
library(interplot)
library(reshape2)

#read in cell counts (FACS) data
cell.countsURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_cell_counts_processed.csv"
cell.counts <- read.csv(text = getURL(cell.countsURL))

#remove mouse 293 as it was mixed both posterior and anterior + remove everpresent X column
cell.counts$X = NULL
cell.counts <- cell.counts[!cell.counts$EH_ID%in%"LM0293",] 

#####################################################################################################################################
#introduce parasitological data
paraURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_complete.csv"
E7 <- read.csv(text=getURL(paraURL))
E7$X = NULL
#merge FACS with para data
E7 <- merge(E7[E7$dpi%in%8,], cell.counts, by = "EH_ID")
#include combined infection history
E7$infHistory <- E7$primary:E7$challenge

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
E7_inf <- E7_inf %>% tidyr::separate(Name, c("LM", "EH_ID"))
E7_inf$EH_ID <- sub("^", "LM", E7_inf$EH_ID )
E7_inf$LM <- NULL
E7 <- merge(E7, E7_inf, by = "EH_ID")

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


## Now... I fooled myself a bit to that enthusiasm, as I expected
## "HybridStatusoutbred hybrids" to be ... well ... hybrids. Turns out
## this are the within subspecies outbreds. Let's do some PostHoc
## comparison. 

summary(glht(modsHY.l[["Tc1IFNgp_in_CD8p"]], mcp(HybridStatus="Tukey")))

## nothing too shocking here, just that "outbred hybrids" have a trend
## towards lower cell proportions compared to "inter subsp. hybrids"

# ---------------------------------------------------------- Make connections between facets of models--------
