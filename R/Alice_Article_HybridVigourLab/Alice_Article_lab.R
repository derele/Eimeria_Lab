# code Vivian & Franci & Alice's experiments
# March 2019
## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85
source("myFunctions.R")
#### Load data ####
source("loadExpe001toExpe005.R")

############## IDEAS block
# align by day +1 -1 weight loss and shedding
# tolerance / resistance / host impact
# change def tol cf raseberg 2015
##############

################# Part 1. data preparation ################# 
# CHOICE MADE : remove Ploen mice & NMRI. Variations are too high.
ALL_Expe <- merge(ExpeDF_003_4, ExpeDF_005, all = T)
ALL_Expe <- makeMiceGenotypeAndIsolate(ALL_Expe)

# NB if NA for oocysts.per.tube, set to 0 if the following day as 0 oocysts
ALL_Expe[is.na(ALL_Expe$oocysts.per.tube), "dpi"] 
# mainly dpi 1,2 and 3 (not taken) Check if zero at dpi4
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 4]
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 3]
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 2]
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 1]
miceWithNoOoDPI4 <- ALL_Expe$EH_ID[ALL_Expe$dpi %in% 4 & ALL_Expe$oocysts.per.tube %in% 0]
ALL_Expe$oocysts.per.tube[ALL_Expe$EH_ID %in% miceWithNoOoDPI4 &
                            ALL_Expe$dpi %in% 3 & 
                            is.na(ALL_Expe$oocysts.per.tube)] <- 0
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 4]
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 3]
# So we are safe to do : 
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 2] <- 0
ALL_Expe$oocysts.per.tube[ALL_Expe$dpi %in% 1] <- 0

# what's left
miceDeadBeforeEnd <- c("LM0168", "LM0187", "LM0189", "LM0193")
ALL_Expe <- ALL_Expe[!ALL_Expe$EH_ID %in% miceDeadBeforeEnd,]

ALL_Expe[is.na(ALL_Expe$oocysts.per.tube), ] 
# mouse LM0202 has one peak day missing for oocyst collection cause diarrhea
summaryDF <- makeToleranceTable(ALL_Expe) #!!! when dpi several, first chosen

################# END data preparation ################# 

plotsALL <- makeIntermPlots(ALL_Expe)
plotsALL[[1]] + coord_cartesian(ylim = c(-10,20))
plotsALL[[2]]# + coord_cartesian(ylim = c(-10,20))

## Experimental design
expedesign <- summaryDF %>%
  group_by(infection_isolate, Mouse_genotype)%>%
  summarise(n = n()) %>% 
  data.frame()
expedesign

## mean of peak day by host and parasite
by_host_parasite_peak <- summaryDF %>% 
  group_by(infection_isolate) %>% 
  summarise(
    meanPeakParasite = mean(dpi_maxOocysts),
    medianPeakParasite = median(dpi_maxOocysts),
    sdPeakParasite = sd(dpi_maxOocysts),
    meanPeakHost = mean(dpi_minWeight),
    medianPeakHost = median(dpi_minWeight),
    sdPeakHost = sd(dpi_minWeight)
    ) %>% 
  data.frame()
by_host_parasite_peak
# infection_isolate meanPeakParasite medianPeakParasite sdPeakParasite
# 1    E.ferrisi (E139)         6.208333                  6      0.7210600
# 2     E.ferrisi (E64)         6.260000                  6      0.6327781
# 3 E.falciformis (E88)         8.227273                  8      0.4289320
# meanPeakHost medianPeakHost sdPeakHost
# 1     5.541667              5   3.106503
# 2     4.740000              5   2.593634
# 3     8.090909              9   3.544400

# Shift all according to median peak, and keep equal window
mydataShifted <- data.frame(Exp_ID = ALL_Expe$Exp_ID,
                            EH_ID = ALL_Expe$EH_ID,
                            infection_isolate = ALL_Expe$infection_isolate,
                            Eimeria_species = ALL_Expe$Eimeria_species,
                            Mouse_genotype = ALL_Expe$Mouse_genotype,
                            dpi = ALL_Expe$dpi,
                            startingWeight = ALL_Expe$startingWeight)

# shift all centered on dpi 8 (centered around E88 parasite peak)
## oocysts
O_E88 <- ALL_Expe[ALL_Expe$infection_isolate == "E.falciformis (E88)",
                  c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E64 <- ALL_Expe[ALL_Expe$infection_isolate == "E.ferrisi (E64)",
                  c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E64$dpi <- O_E64$dpi + 2 
O_E139 <- ALL_Expe[ALL_Expe$infection_isolate == "E.ferrisi (E139)", 
                   c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E139$dpi <- O_E139$dpi + 2 
oocystsShifted <- rbind(O_E88, O_E64, O_E139)

## weight
W_E88 <- ALL_Expe[ALL_Expe$infection_isolate == "E.falciformis (E88)",
                  c("dpi", "weight", "EH_ID")]
W_E88$dpi <- W_E88$dpi - 1 
W_E64 <- ALL_Expe[ALL_Expe$infection_isolate == "E.ferrisi (E64)",
                  c("dpi", "weight", "EH_ID")]
W_E64$dpi <- W_E64$dpi + 3 
W_E139 <- ALL_Expe[ALL_Expe$infection_isolate == "E.ferrisi (E139)", 
                   c("dpi", "weight", "EH_ID")]
W_E139$dpi <- W_E139$dpi + 3 

weightShifted <- rbind(W_E88, W_E64, W_E139)

mydataShifted <- merge(mydataShifted, weightShifted)
mydataShifted <- merge(mydataShifted, oocystsShifted)

## restrict window to get same info on all
table(mydataShifted$dpi, mydataShifted$fecweight == 0)

#### Keep WINDOWS OF INFECTION #### 
# peak = dpi 8. 3 days before + peak + 2 days after full mice
mydataShifted <- mydataShifted[mydataShifted$dpi %in% 5:10,]
table(mydataShifted$dpi)

# consider weightBegWind as weight beginning of window
A = mydataShifted[mydataShifted$dpi == 5, c("weight", "EH_ID")]
names(A)[1] = "startingWeight"
mydataShifted <- mydataShifted %>% 
  dplyr::rename(previousStartingweight = startingWeight)
mydataShifted <- merge(mydataShifted, A)

# Cumulative sum of our weight loss
mydataShifted = mydataShifted %>% 
  group_by(EH_ID) %>% 
  dplyr::arrange(dpi, .by_group = TRUE) %>%
  dplyr::mutate(relativeWeight = weight / startingWeight * 100)
mydataShifted = data.frame(mydataShifted)

mydataShifted$OPG <- mydataShifted$oocysts.per.tube / mydataShifted$fecweight

# Summary DF within infection window
summaryDF2 <- makeToleranceTable(mydataShifted) #!!! when dpi several, first chosen

# check if our peaks are set at 8
by_host_parasite_peak2 <- summaryDF2 %>% 
  group_by(infection_isolate) %>% 
  summarise(
    meanPeakParasite = mean(dpi_maxOocysts),
    medianPeakParasite = median(dpi_maxOocysts),
    sdPeakParasite = sd(dpi_maxOocysts),
    meanPeakHost = mean(dpi_minWeight),
    medianPeakHost = median(dpi_minWeight),
    sdPeakHost = sd(dpi_minWeight)
  ) %>% 
  data.frame()
by_host_parasite_peak2

## Statistical models along dpi

# OFFSET: 
#https://stats.stackexchange.com/questions/237963/how-to-formulate-the-offset-of-a-glm

### Which distribution?
hist(mydataShifted$weight, col = "red")
hist(mydataShifted$oocysts.per.tube, breaks = 100, col = "red") # clearly negbin here :) 

x <- as.numeric(na.omit(mydataShifted$oocysts.per.tube))
x <- x[x>0]
library(fitdistrplus)

plot(x, pch = 20)
plotdist(x, histo = TRUE, demp = TRUE)
descdist(x, discrete=FALSE, boot=500)

fit_no  <- fitdist(x, "norm")
fit_po  <- fitdist(x, "pois")
fit_ne <- fitdist(x, "nbinom")
summary(fit_no)
summary(fit_po)
summary(fit_ne)

par(mfrow=c(2,2))
plot.legend <- c("normal", "poisson", "negbin")
denscomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
cdfcomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
qqcomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
ppcomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
par(mfrow=c(1,1))

# The density plot and the CDF plot may be considered as the basic classical 
# goodness-of-fits plots. The Q-Q plot emphasizes the lack-of-fit at the 
# distribution tails while the P-P plot emphasizes the lack-of-fit at the 
# distribution center. The nbinom distribution could be prefered for its better
# description of the tail of the empirical distribution.
# The negative binomial distribution seems to describe parasite load well 

### Resistance
# Raberg 2007 results on plasmodium
# To test for variation in resistance among mouse strains, we performed
# an analysis of peak parasite density against mouse strain and parasite clone.
# Peak parasite density differed between mouse strains
# [F4,102 = 15.54, P < 0.0001] and parasite clones
# [F2,103 = 64.81, P < 0.0001], but there was no strainby-clone interaction 
# [F8,102 = 0.66, P= 0.73]. There was also a significant effect of experimental 
# block (c2 = 47.4, P< 0.0001), but no interactions between
# block and strain and/or clone (P > 0.25). Thus, as in previous studies 
# (20–22, 25), mouse strains differed in resistance, and parasite clones 
# differed in the infection intensity that they induced.

# negative binomial requires INTEGER 
range(mydataShifted$oocysts.per.tube[mydataShifted$oocysts.per.tube > 0], na.rm = T)
mydataShifted$oocysts.per.tube.tsd <- as.integer(mydataShifted$oocysts.per.tube / 1000)
summaryDF2$oocysts.per.tube.tsd <- as.integer(summaryDF2$oocysts.per.tube / 1000)

# use of random factors : we assume a different baseline response per day
dev.off()
ggplot(mydataShifted, 
       aes(x = dpi, y = oocysts.per.tube.tsd + 1, group = dpi)) +
  geom_boxplot(fill = "red") + scale_y_log10()

ggplot(mydataShifted, 
       aes(x = dpi, y = weight, group = dpi)) +
  geom_boxplot(fill = "red") 

nrow(mydataShifted[mydataShifted$fecweight != 0,])
nrow(mydataShifted)

## General definition of resistance = 1 per mouse

## FIRST, test Exp_ID effect
m1.P.full <- glm.nb(oocysts.per.tube ~ Eimeria_species * Mouse_genotype + Exp_ID + 
                      offset(log(fecweight)), data = summaryDF2)
m1.P <- glm.nb(oocysts.per.tube ~ Eimeria_species * Mouse_genotype +
                      offset(log(fecweight)), data = summaryDF2)
# Likelihood ratio test
anova(m1.P.full, m1.P, test="LRT")
# EXPERIMENTAL BLOCK NOT SIGNIFICANT

summary(m1.P)
# post-hoc tests
mytukey(m1.P)
# $infection_isolate
# combi        padj      diff
# 1    E.ferrisi (E64)-E.ferrisi (E139) 0.005325540  616681.1
# 2 E.falciformis (E88)-E.ferrisi (E64) 0.006762052 -619478.9

## More detailed = along infection
# test significance of interactions
# Original model
model1.P <- glmer.nb(oocysts.per.tube.tsd ~
                       Eimeria_species * Mouse_genotype +
                       (1|dpi) + offset(log(fecweight)),
                     data = mydataShifted)
# Model without interaction
model2.P <- glmer.nb(oocysts.per.tube.tsd ~ 
                       Eimeria_species + Mouse_genotype + 
                     (1|dpi) + offset(log(fecweight)), 
                   data = mydataShifted)

# Likelihood ratio test
anova(model1.P, model2.P, test="LRT")
## --> INTERACTIONS NOT SLIGHTLY
summary(model2.P)

# Post-hoc analysis
s1 <- summary(glht(model2.P, mcp(Eimeria_species="Tukey")))
s1
s2 <- summary(glht(model2.P, mcp(Mouse_genotype="Tukey")))
s2

### Pathogenicity 

## General definition of resistance = 1 per mouse
m1.H <- glm(weight ~ Eimeria_species * Mouse_genotype + 
                 offset(log(startingWeight)), 
               data = summaryDF2)
m2.H <- glm(weight ~ Eimeria_species + Mouse_genotype + 
                 offset(log(startingWeight)), 
               data = summaryDF2)
# Likelihood ratio test
anova(m1.H, m2.H, test="LRT")
## Interaction NOT significant
summary(m2.H)
# post-hoc tests
mytukey(m2.H)

## More detailed = along infection
# test significance of interactions
# Original model
model1.H <- lmer(weight ~
                   Eimeria_species * Mouse_genotype +
                   (1|dpi) + offset(log(startingWeight)),
                 data = mydataShifted)
# Model without interaction
model2.H <- lmer(weight ~ 
                   Eimeria_species + Mouse_genotype + 
                   (1|dpi) + offset(log(startingWeight)), 
                 data = mydataShifted)

# Likelihood ratio test
anova(model1.H, model2.H, test="LRT")
## --> INTERACTIONS SUPER SIGNIFICANT
summary(model1.H)

# Post-hoc analysis
s1.H <- summary(glht(model1.H, mcp(Eimeria_species="Tukey")))
s1.H
s2.H <- summary(glht(model1.H, mcp(Mouse_genotype="Tukey")))
s2.H

### Tolerance (slope) Raberg 2007 Raber 2008
# resistance = inverse of peak parasite density
summaryDF2$peakParasiteDensity <- summaryDF2$oocysts.per.tube / summaryDF2$fecweight
summaryDF2$resistance <- - summaryDF2$peakParasiteDensity
summaryDF2$minWeightRelative <- summaryDF2$weight / summaryDF2$startingWeight
summaryDF2$tolerance <- summaryDF2$minWeightRelative / summaryDF2$peakParasiteDensity

ggplot(summaryDF2, aes(x = resistance, y = tolerance)) +
  geom_point(aes(fill = Mouse_genotype), size=4, pch = 21, color = "black") +
  facet_grid(. ~ Eimeria_species) +
  xlab("Resistance = inverse of peak parasite density")+
  ylab("Tolerance = minimun relative weight / peak parasite density")+
  theme(axis.text.x = element_text(angle = 0))

summaryTolerance <- summarySE(summaryDF2, 
                              measurevar="tolerance",
                              groupvars=c("Mouse_genotype",
                                          "Eimeria_species"), 
                              na.rm = T)
summaryTolerance$ci[is.na(summaryTolerance$ci)] <- 0
summaryResistance <- summarySE(summaryDF2, 
                              measurevar="resistance",
                              groupvars=c("Mouse_genotype",
                                          "Eimeria_species"), 
                              na.rm = T)
summaryResistance$ci[is.na(summaryResistance$ci)] <- 0

library(psych)
names(summaryTolerance)[names(summaryTolerance) %in% "N"] <- "n"
names(summaryTolerance)[names(summaryTolerance) %in% "tolerance"] <- "mean"
names(summaryResistance)[names(summaryResistance) %in% "N"] <- "n"
names(summaryResistance)[names(summaryResistance) %in% "resistance"] <- "mean"
summaryTolerance$axis <- "tolerance"
summaryResistance$axis <- "resistance"
summaryTable <- rbind(summaryTolerance, summaryResistance)

error.crosses(summaryTable[summaryTable$axis %in% "tolerance",], 
              summaryTable[summaryTable$axis %in% "resistance",],
              labels= summaryTable$Mouse_genotype,
              xlab = "tolerance", ylab = "resistance",
              pch=16,cex=1)

# tolerance = slope of a regression of minimum weight against peak parasite density
# modelTol <- lm(weight ~ peakParasiteDensity * infection_isolate * Mouse_genotype +
#                  offset(log(startingWeight)), data = summaryDF2)
# # 
# t <- summary(modelTol)
# d <- data.frame(t$coefficients)
# 
# rowIndx <- grepl("peak", row.names(d)) & grepl("infection", row.names(d)) & grepl("Mouse", row.names(d))
# d[rowIndx, "Estimate"]


mydataShifted$tolFac <- mydataShifted$relativeWeight / mydataShifted$OPG
mydataShifted$tolFac[mydataShifted$tolFac %in% c(Inf, -Inf)] <- NA

summaryTolerance <- summarySE(mydataShifted, 
                              measurevar="tolFac",
                              groupvars=c("Mouse_genotype",
                                          "Eimeria_species", "dpi"), 
                              na.rm = T)
summaryTolerance$ci[is.na(summaryTolerance$ci)] <- 0

ggplot(summaryTolerance, aes(x = dpi, y = tolFac))+
  geom_errorbar(aes(ymin = tolFac - ci, ymax = tolFac + ci),
                col = "gray") +
  geom_line(aes(group = Mouse_genotype, col = Mouse_genotype), size = 2, alpha = 0.5) +
  geom_point(aes(fill = Mouse_genotype), size=4, pch = 21, color = "black") +
  facet_grid(. ~ Eimeria_species) +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_log10() +
  geom_vline(xintercept = 8) +
  theme(axis.text.x = element_text(angle = 0))








## WRONG if we use full model. Too many assumptions. weight before peak and after peak shall 
## vary, not link with oocysts. The points before and after the peak are NOT independant.
## we can't treat them as repeated measures. 
## Hence, ONE tolerance value per animal. Defined by Raberg 2007.

#### Test interactions
modT1 <- lm(weight ~ oocysts.per.tube + 
              Eimeria_species + Mouse_genotype +
              oocysts.per.tube : Eimeria_species +
              oocysts.per.tube : Mouse_genotype + 
              offset(log(startingWeight)), 
            data = summaryDF2)

modT2 <- lm(weight ~ Eimeria_species + 
              Eimeria_species + Mouse_genotype +
              oocysts.per.tube : Mouse_genotype + 
              offset(log(startingWeight)), 
            data = summaryDF2)

modT3 <- lm(weight ~ oocysts.per.tube + 
              Eimeria_species + Mouse_genotype +
              oocysts.per.tube : Eimeria_species + 
              offset(log(startingWeight)), 
            data = summaryDF2)

modT4 <- lm(weight ~ oocysts.per.tube + 
              Eimeria_species + Mouse_genotype +
              offset(log(startingWeight)), 
            data = summaryDF2)

anova(modT1, modT2, test="LRT") # NOT SIGNIF DIFF WHEN DROP INF ISOLATE
anova(modT1, modT3, test="LRT") # NOT SIGNIF DIFF WHEN DROP MOUSE STRAIN
anova(modT1, modT4, test="LRT") # NO SIGNIF DIFF WHEN DROP BOTH
# --> choose T4

summary(modT4)
sTH <- summary(glht(modT4, mcp(Mouse_genotype="Tukey")))
sTP <- summary(glht(modT2, mcp(Eimeria_species="Tukey")))
sTH
sTP

## Damn! How comes interactions are not significant??? That's what we want to compare!
