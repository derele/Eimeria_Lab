# code Vivian & Franci & Alice's experiments
# March 2019
## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85
source("myFunctions.R")
#### Load data ####
source("loadExpe001toExpe005.R")
theme_set(theme_bw() +  theme(text = element_text(size = 20)))

## Packages
list.of.packages <- c("parasiteLoad", "bbmle", "devtools", "optimx", # for bbmle it needs to be required(?)
                      "ggplot2", "VennDiagram","fitdistrplus", # evaluate distribution
                      "epiR", # Sterne's exact method
                      # "simpleboot", # BS
                      "ggmap", "gridExtra",# several plots in one panel
                      "wesanderson", # nice colors
                      # "cowplot",# several plots in one panel
                      "ggpubr")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(list.of.packages)

## Reinstall the package in case I updated it
devtools::install_github("alicebalard/parasiteLoad")
library(parasiteLoad)

## Structure article ##
# Verif: Infect ALL --> the. prevalence, link with Victor paper
# matandmet : yound mice, weigth stabilisation period before (test no more growing)
# 1. test ANTH/no ANTH (expe 3-4 vs 5 parents) -> hyp, diff if 1 or 2 parasites (interaction)
# 2. HV in expe 1 & 5. Older mice suffer more? Test age. Could be, closer to "real" conditions 
# 3. H or P drivers? Restistance/tolerance. Important, test effect of strating weight (if small mouse, more affected)

## Mat&met
# age, sex controlled by balanced design within expe
## Orga: a table with all mice and their measures by day +
## a summary table with for each mouse peak shedding, time to peak shedding, peak WL, tol as ratio, and infos

# Remove Expe_001 until further notice

## Keep ONLY first batch!!! Contamination in the second one...
ALL_Expe <- merge(ExpeDF_003_4, ExpeDF_005[ExpeDF_005$Batch %in% 1,], all = T)
ALL_Expe <- makeMiceGenotypeAndIsolate(ALL_Expe)

# NB some mice died before the end of expe, treat carefully
# miceDeadBeforeEnd <- c("LM0168", "LM0187", "LM0189", "LM0193")
# ExpeDF_003_4 <- ExpeDF_003_4[!ExpeDF_003_4$EH_ID %in% miceDeadBeforeEnd,]
# ExpeDF_005 <- ExpeDF_005[!ExpeDF_005$EH_ID %in% miceDeadBeforeEnd,]

################# Part 1. data preparation ################# 

# NB if NA for oocysts.per.tube, set to 0 if the following day as 0 oocysts
# if zero at dpi 4 but 1,2,3 not collected, fill the gaps (falciformis)
temp <- ALL_Expe$EH_ID[ALL_Expe$dpi %in% 4 & ALL_Expe$oocysts.per.tube %in% 0]
ALL_Expe$oocysts.per.tube[
  ALL_Expe$EH_ID %in% temp & ALL_Expe$dpi %in% 1:3 & is.na(ALL_Expe$oocysts.per.tube)] <- 0

# if zero at dpi 3 but 1,2 not collected, fill the gaps (ferrisi)
temp <- ALL_Expe$EH_ID[ALL_Expe$dpi %in% 3 & ALL_Expe$oocysts.per.tube %in% 0]
ALL_Expe$oocysts.per.tube[
  ALL_Expe$EH_ID %in% temp & ALL_Expe$dpi %in% 1:2 & is.na(ALL_Expe$oocysts.per.tube)] <- 0

explore <- ALL_Expe[is.na(ALL_Expe$oocysts.per.tube), ] 
# LM0168 died at dpi 10 (keep, peak present)
# LM0187 died at dpi8, no feces collected (careful)
# LM0189 died at dpi9 (careful
# LM0193 died at dpi8, no feces collected (careful)
# LM0202 died at dpi8, no feces collected (careful) + mouse LM0202 has one day missing for oocyst collection cause diarrhea
# In ExpeDF_001, a lot of mice (falci) died before the end of expe, careful
ALL_summary <- makeSummaryTable(ALL_Expe)

###### IMPACT OF AGE ON RESULTS ###### 
p1 <- ggplot(ALL_summary, aes(x=ageAtInfection, y =minWeightRelative,
                              col = Mouse_genotype)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_grid(Eimeria_species~., scales = "free") +
  scale_x_continuous(name="Age at infection (weeks)") +
  scale_y_continuous(name="Min weight retained (%)") +
  theme(legend.position = "none")
p2 <- ggplot(ALL_summary, aes(x=ageAtInfection, y =max.OPG,
                              col = Mouse_genotype)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_y_log10(name="Max oocysts per gram of feces") +
  scale_x_continuous(name="Age at infection (weeks)") +
  facet_grid(Eimeria_species~., scales = "free") +
  theme(legend.position = "none")
grid.arrange(p1, p2, ncol=2)

###### 

###### Part I. No complete resistance to both Eimeria falciformis and Eimeria ferrisi###### 
length(unique(ALL_Expe$EH_ID))
length(unique(ALL_summary$EH_ID))

table(ALL_summary$max.oocysts.per.tube != 0, ALL_summary$infection_isolate)

# we remove uninfected mouse LM0187
toRM <- ALL_summary$EH_ID[ALL_summary$max.oocysts.per.tube == 0]
ALL_Expe <- ALL_Expe[!ALL_Expe$EH_ID %in% toRM,]
ALL_summary <- ALL_summary[!ALL_summary$EH_ID %in% toRM,]
table(ALL_summary$max.oocysts.per.tube != 0, ALL_summary$infection_isolate)

## Idea: include time. See how to geek that.

## Window of infection:

## mean of peak day by host and parasite
by_host_parasite_peak <- ALL_summary %>% 
  group_by(infection_isolate) %>% 
  summarise(
    meanPeakParasite = mean(dpi_max.oocysts.per.tube),
    medianPeakParasite = median(dpi_max.oocysts.per.tube),
    sdPeakParasite = sd(dpi_max.oocysts.per.tube),
    meanPeakHost = mean(dpi_minWeight),
    medianPeakHost = median(dpi_minWeight),
    sdPeakHost = sd(dpi_minWeight)
  ) %>% 
  data.frame()
by_host_parasite_peak
# infection_isolate meanPeakParasite medianPeakParasite sdPeakParasite meanPeakHost medianPeakHost sdPeakHost
# 1    E.ferrisi (E139)         6.200000                  6      0.7071068        5.720              5   3.169122
# 2     E.ferrisi (E64)         6.260000                  6      0.6327781        4.740              5   2.593634
# 3 E.falciformis (E88)         8.208333                  8      0.5089774        8.125              9   3.391966

# Shift all according to median peak, and keep equal window
mydataShifted <- data.frame(Exp_ID = ALL_Expe$Exp_ID,
                            EH_ID = ALL_Expe$EH_ID,
                            HybridStatus = ALL_Expe$HybridStatus,
                            infection_isolate = ALL_Expe$infection_isolate,
                            Eimeria_species = ALL_Expe$Eimeria_species,
                            Mouse_genotype = ALL_Expe$Mouse_genotype,
                            dpi = ALL_Expe$dpi,
                            startingWeight = ALL_Expe$startingWeight,
                            ageAtInfection = ALL_Expe$ageAtInfection,
                            Sex = ALL_Expe$Sex)

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
summaryDF2 <- makeSummaryTable(mydataShifted) #!!! when dpi several, first chosen

# check if our peaks are set at 8
by_host_parasite_peak2 <- summaryDF2 %>% 
  group_by(infection_isolate) %>% 
  summarise(
    meanPeakParasite = mean(dpi_max.oocysts.per.tube),
    medianPeakParasite = median(dpi_max.oocysts.per.tube),
    sdPeakParasite = sd(dpi_max.oocysts.per.tube),
    meanPeakHost = mean(dpi_minWeight),
    medianPeakHost = median(dpi_minWeight),
    sdPeakHost = sd(dpi_minWeight)
  ) %>% 
  data.frame()
by_host_parasite_peak2

## Statistical models along dpi 
# !! compare NB with Poisson later
fitmodel1 <- glmer(OPG ~ Mouse_genotype * infection_isolate * poly(dpi) +
                        dpi|EH_ID, data = mydataShifted, family = "poisson")


glmer.nb(OPG ~ Mouse_genotype * infection_isolate * poly(dpi) +
           dpi|EH_ID, data = mydataShifted)




mydataShifted

fitmodel1 <- lmer(OPG ~ 1 + 1|Exp_ID, data = mydataShifted) # random effect 1 way anova
summary(fitmodel1)
# https://personality-project.org/r/tutorials/summerschool.14/rosseel_lme3.pdf
# the intra-expe coefficient rho is the proportion of variation explained by diff between expe.
# also correlation between the OPG of mice within same expe
d2 = as.numeric(VarCorr(fitmodel1)$Exp_ID)
s2 = attr(VarCorr(fitmodel1), "sc")^2
rho = d2/(d2+s2)
rho # cluster effect 0.5%




############
kruskal.test(max.OPG ~ infection_isolate, data = ALL_summary) 
# Kruskal-Wallis rank sum test
#
# data:  max.OPG by infection_isolate
# Kruskal-Wallis chi-squared = 3.4747, df = 2, p-value = 0.176
ggplot(ALL_summary, aes(x = infection_isolate, y = max.OPG)) +
  geom_boxplot() +
  geom_segment(aes(x = 1.05, y = 3e+06, xend = 1.95, yend = 3e+06)) +
  annotate(geom="text", x=1.5, y=3.1e+06, label="n.s.") +
  geom_segment(aes(x = 2.05, y = 3e+06, xend = 2.95, yend = 3e+06)) +
  annotate(geom="text", x=2.5, y=3.1e+06, label="n.s.")

kruskal.test(minWeightRelative ~ infection_isolate, data = ALL_summary) 
# Kruskal-Wallis rank sum test
# 
# data:  minWeightRelative by infection_isolate
# Kruskal-Wallis chi-squared = 4.8525, df = 2, p-value = 0.08837
kruskal.test(tolfacRelative ~ infection_isolate, data = ALL_summary) 
# Kruskal-Wallis rank sum test
# 
# data:  tolfacRelative by infection_isolate
# Kruskal-Wallis chi-squared = 3.3433, df = 2, p-value = 0.1879
ggplot(ALL_summary, aes(x = infection_isolate, y = tolfacRelative)) +
  geom_boxplot() +
  scale_y_log10() +
  geom_segment(aes(x = 1.05, y = 1e-03, xend = 1.95, yend = 1e-03)) +
  annotate(geom="text", x=1.5, y=1.1e-03, label="n.s.") +
  geom_segment(aes(x = 2.05, y = 1e-03, xend = 2.95, yend = 1e-03)) +
  annotate(geom="text", x=2.5, y=1.1e-03, label="n.s.")

###### Part II. Test HV in lab res/tol ###### 
F1DF_expe <- ALL_Expe[grep("F1", ALL_Expe$Mouse_genotype),]
F1DF_summary <- ALL_summary[grep("F1", ALL_summary$Mouse_genotype),]

table(F1DF_summary$group, F1DF_summary$infection_isolate)
F1DF_summary$group <- "parental"
F1DF_summary$group[grep("Mmm-Mmd", F1DF_summary$Mouse_genotype)] <- "hybrid"

# F1DF_summary$HI <- NA
# F1DF_summary$HI[grep("MMm", F1DF_summary$Mouse_genotype)] <- 1
# F1DF_summary$HI[grep("MMd", F1DF_summary$Mouse_genotype)] <- 0
# F1DF_summary$HI[grep("Mmm-Mmd", F1DF_summary$Mouse_genotype)] <- 0.5

ggplot(F1DF_summary, aes(x = group, y = max.OPG, fill = infection_isolate)) +
  geom_boxplot() + geom_point(pch = 21, size = 3)

wilcox.test(max.OPG ~ group, 
            data = F1DF_summary[F1DF_summary$infection_isolate == "E.ferrisi (E64)",]) 

wilcox.test(max.OPG ~ group, 
            data = F1DF_summary[F1DF_summary$infection_isolate == "E.falciformis (E88)",]) 


###### 


###### BANANA model ###### 
bananalabDF <- ALL_summary[grep("F1", ALL_summary$Mouse_genotype),]
bananalabDF$Eimeria_species <- as.factor(bananalabDF$Eimeria_species)

bananalabDF$max.OPG.million <- bananalabDF$max.OPG / 1e6


hist(bananalabDF$max.OPG.million) # close to normal-ish

range(bananalabDF$max.OPG.million)
median(bananalabDF$max.OPG.million)
sd(bananalabDF$max.OPG.million)

speparam <- c(L1start = 1,
              L1LB = 0,
              L1UB = 3,
              L2start = 1,
              L2LB = 0,
              L2UB = 3,
              alphaStart = 0, alphaLB = -10, alphaUB = 10,
              mysdStart = 1, mysdLB = 1.e-9, mysdUB = 10)

fitMaxOPG <- parasiteLoad::analyse(data = bananalabDF,
                                   response = "max.OPG",
                                   model = "normal",
                                   group = "Eimeria_species",
                                   myparamBounds = speparam)
fitMaxOPG$H1
# "Testing H1 no alpha vs alpha"
# dLL dDF     pvalue
# 1 5.4   1 0.00101085

bananaPlot(mod = fitMaxOPG$H1,
           data = bananalabDF,
           response = "max.OPG",
           group = "Eimeria_species") +
  scale_fill_manual(values = c("purple", "green")) +
  scale_color_manual(values = c("purple", "green")) +
  scale_y_log10()+
  theme_bw() + theme(text = element_text(size = 20))

## Min Weight percent

fitMinWL <- parasiteLoad::analyse(data = bananalabDF,
                                   response = "minWeight",
                                   model = "normal",
                                   group = "Eimeria_species")
fitMinWL$H1
#"Testing H1 no alpha vs alpha"
# dLL dDF pvalue
# 1 478834   1      0

bananaPlot(mod = fitMinWL$H1,
           data = bananalabDF,
           response = "minWeight",
           group = "Eimeria_species") +
  scale_fill_manual(values = c("purple", "green")) +
  scale_color_manual(values = c("purple", "green")) +
  theme_bw() + theme(text = element_text(size = 20))



###### 

###### resistance
ggplot(ALL_summary, aes(x = Mouse_genotype, y = max.OPG, fill = Mouse_genotype)) +
  geom_boxplot() + geom_jitter() +
  scale_y_log10() + facet_grid(ageCat~infection_isolate, scales = "free") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

mod1 <- lm(data = oldMiceSummary, 
          formula = max.oocysts.per.tube ~ infection_isolate * Mouse_genotype +
            offset(log(fecweight)))
mod2 <- lm(data = oldMiceSummary, 
          formula = max.oocysts.per.tube ~ infection_isolate + Mouse_genotype +
            offset(log(fecweight)))
mod3 <- lm(data = oldMiceSummary, 
           formula = max.oocysts.per.tube ~ infection_isolate +
             offset(log(fecweight)))
mod4 <- lm(data = oldMiceSummary, 
           formula = max.oocysts.per.tube ~ Mouse_genotype +
             offset(log(fecweight)))
anova(mod1, mod2)
anova(mod2, mod3)
anova(mod2, mod4)

## impact of host health
ggplot(ALL_summary, aes(x = Mouse_genotype, y = min.weight.retained.percent, fill = Mouse_genotype)) +
geom_boxplot() + geom_jitter() +
  facet_grid(ageCat~infection_isolate, scales = "free") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

mod <- lm(data = oldMiceSummary, 
          formula = minWeight ~ infection_isolate * Mouse_genotype +
            offset(log(startingWeight)))
summary(mod) # diff E.falci/E.ferrisi

## tolerance
ggplot(ALL_summary, aes(x = Mouse_genotype, y = tolfac, fill = Mouse_genotype)) +
  geom_boxplot() + geom_jitter() + scale_y_log10() +
  facet_grid(ageCat~infection_isolate, scales = "free") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

mod <- lm(data = oldMiceSummary, 
          formula = tolfac ~ infection_isolate * Mouse_genotype)
summary(mod)

## time to peak Par
ggplot(oldMiceSummary, aes(x = Mouse_genotype, y = dpi_max.oocysts.per.tube)) +
  geom_boxplot() + geom_jitter() +
  facet_grid(.~infection_isolate)

mod <- lm(data = oldMiceSummary, 
          formula = dpi_max.oocysts.per.tube ~ infection_isolate)
summary(mod) # diff E.falci/E.ferrisi

## time to peak Host
ggplot(oldMiceSummary, aes(x = Mouse_genotype, y = dpi_minWeight)) +
  geom_boxplot() + geom_jitter() +
  facet_grid(.~infection_isolate)

mod <- lm(data = oldMiceSummary, 
          formula = dpi_minWeight ~ infection_isolate)
summary(mod) # diff E.falci/E.ferrisi




p1 <- ggplot(ALL_summary, aes(x=ageAtInfection, y =min.weight.retained.percent)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Eimeria_species~., scales = "free") +
  scale_x_continuous(name="Age at infection (weeks)") +
  scale_y_continuous(name="Min weight retained (%)")
p2 <- ggplot(ALL_summary, aes(x=ageAtInfection, y =min.weight.retained.percent,
                              col = Mouse_genotype)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_grid(Eimeria_species~ageCat, scales = "free") +
  scale_x_continuous(name="Age at infection (weeks)") +
  scale_y_continuous(name="Min weight retained (%)") +
  theme(legend.position = "none")
grid.arrange(p1, p2, ncol=2)

p3 <- ggplot(ALL_summary, aes(x=ageAtInfection, y =max.OPG)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10(name="Max oocysts per gram of feces") +
  facet_grid(Eimeria_species~., scales = "free") +
  scale_x_continuous(name="Age at infection (weeks)")
p4 <- ggplot(ALL_summary, aes(x=ageAtInfection, y =max.OPG,
                              col = Mouse_genotype)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_y_log10(name="Max oocysts per gram of feces") +
  facet_grid(Eimeria_species~ageCat, scales = "free") +
  theme(legend.position = "none")
grid.arrange(p3, p4, ncol=2)

###### 

######## Previously, when tryed to decompose (too noisy in my opinion)
summaryDF <- makeToleranceTable(ALL_Expe) #!!! when dpi several, first chosen

############## IDEAS block
# align by day +1 -1 weight loss and shedding
# tolerance / resistance / host impact
# change def tol cf raseberg 2015
##############


## Check age/sex/batch possible confounding factors

# Infection date : 
ALL_Expe$Sex <- gsub(" ", "", ALL_Expe$Sex)

ggplot(ALL_Expe[!duplicated(ALL_Expe$EH_ID),],
       aes(x = Mouse_genotype, y = ageAtInfection, fill = Sex)) +
  geom_violin() + geom_point(pch = 21, size = 4, position = position_jitter(width = 0.01))

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

## To do : N; time to peak weight loss; time to peak oocyst; maxWL; peakOO
summaryDF$dpi_maxOocysts






# OFFSET: 
#https://stats.stackexchange.com/questions/237963/how-to-formulate-the-offset-of-a-glm

### Which distribution?
hist(mydataShifted$weight, col = "red")
hist(mydataShifted$oocysts.per.tube, breaks = 100, col = "red") # clearly negbin here :) 
# 
# x <- as.numeric(na.omit(mydataShifted$oocysts.per.tube))
# x <- x[x>0]
# library(fitdistrplus)
# 
# plot(x, pch = 20)
# plotdist(x, histo = TRUE, demp = TRUE)
# descdist(x, discrete=FALSE, boot=500)
# 
# fit_no  <- fitdist(x, "norm")
# fit_po  <- fitdist(x, "pois")
# fit_ne <- fitdist(x, "nbinom")
# summary(fit_no)
# summary(fit_po)
# summary(fit_ne)
# 
# par(mfrow=c(2,2))
# plot.legend <- c("normal", "poisson", "negbin")
# denscomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
# cdfcomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
# qqcomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
# ppcomp(list(fit_no, fit_po, fit_ne), legendtext = plot.legend)
# par(mfrow=c(1,1))

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
  geom_boxplot(fill = "red") + scale_y_log10() + facet_grid(.~Mouse_genotype)

ggplot(mydataShifted, 
       aes(x = dpi, y = weight, group = dpi)) +
  geom_boxplot(fill = "red") + facet_grid(.~Mouse_genotype)

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

mydataShifted$millionOPG <- mydataShifted$OPG/1000000

# relevel for Bu-Bu first (less tolerant)
mydataShifted$Mouse_genotype <- relevel(mydataShifted$Mouse_genotype, ref = "MMm_F0 (Pw-Pw)")

modT <- lmer(weight ~ Eimeria_species * Mouse_genotype * millionOPG +
                (1|EH_ID) + offset(log(startingWeight)),
              data = mydataShifted)


#https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001989

# NO TOLERANCE DIFFERENT SIGNIFICATIVELY
library(dplyr)
library(sjPlot)
library(sjmisc)
library(ggplot2)

theme_set(theme_sjplot())

plot_model(modT, type = "pred", terms = c( "millionOPG", "Mouse_genotype", "Eimeria_species")) 

tab_model(modT)


## Time serie

# Control phase dpi 5 to 8
controlDF <- mydataShifted[mydataShifted$dpi %in% 5:8,]

ggplot(controlDF, aes(x = dpi, y = weight, col = Mouse_genotype)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  # geom_smooth(aes(group = EH_ID)) +
  facet_grid(.~infection_isolate)

fmC1 <- lmer(weight ~ dpi * Mouse_genotype * Eimeria_species +
              (dpi | EH_ID) + offset(log(startingWeight)), controlDF)
fmC2 <- lmer(weight ~ dpi * Mouse_genotype + Eimeria_species +
              (dpi | EH_ID) + offset(log(startingWeight)), controlDF)

anova(fmC1, fmC2)
#The term (Days|Subject) generates a vector-valued random effect
# (intercept and slope) for each of the 18 levels of the Subject factor

summary(fmC1)








# Create temporary data frame:
snarc_coefs <- vector()
for (i in unique(controlDF$EH_ID)) {
  snarc_tmp <-
    controlDF[controlDF$EH_ID %in% i,]
  # Perform regression:
  reg_result <- lm(snarc_tmp$relativeWeight ~
                     snarc_tmp$dpi)
  # Get coefficient:
  tmp_coef <- coef(reg_result)
  # Store coefficient:
  snarc_coefs[i] <- tmp_coef[2]
}

controlDF <- merge(controlDF,
                   data_frame(EH_ID = names(snarc_coefs), slopeRelWeight = snarc_coefs))

unicontrolDF <- controlDF[!duplicated(controlDF$EH_ID),]

unicontrolDF$host <- NA
for(i in c("Mmm-Mmd_F1","MMd_F0", "MMm_F0", "MMm_F1",  "MMd_F1")){
  unicontrolDF$host[grep(i, unicontrolDF$Mouse_genotype)] <- i
  
}

unicontrolDF$host <- as.factor(unicontrolDF$host)
unicontrolDF$host <- factor(unicontrolDF$host, 
                            levels = c("MMd_F0", "MMd_F1", "Mmm-Mmd_F1", "MMm_F1", "MMm_F0"))

ggplot(unicontrolDF, aes(x= host, y = slopeRelWeight, 
                         fill = host )) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(Eimeria_species~.)

# Resilience phase dpi 8 to 10
resiDF <- mydataShifted[mydataShifted$dpi %in% 8:10,]

# Create temporary data frame:
snarc_coefs <- vector()
for (i in unique(resiDF$EH_ID)) {
  snarc_tmp <-
    resiDF[resiDF$EH_ID %in% i,]
  # Perform regression:
  reg_result <- lm(snarc_tmp$relativeWeight ~
                     snarc_tmp$dpi)
  # Get coefficient:
  tmp_coef <- coef(reg_result)
  # Store coefficient:
  snarc_coefs[i] <- tmp_coef[2]
}

resiDF <- merge(resiDF,
                   data_frame(EH_ID = names(snarc_coefs), slopeRelWeight = snarc_coefs))

uniresiDF <- resiDF[!duplicated(resiDF$EH_ID),]

uniresiDF$host <- NA
for(i in c("Mmm-Mmd_F1","MMd_F0", "MMm_F0", "MMm_F1",  "MMd_F1")){
  uniresiDF$host[grep(i, uniresiDF$Mouse_genotype)] <- i
  
}

uniresiDF$host <- as.factor(uniresiDF$host)
uniresiDF$host <- factor(uniresiDF$host, 
                            levels = c("MMd_F0", "MMd_F1", "Mmm-Mmd_F1", "MMm_F1", "MMm_F0"))

ggplot(uniresiDF, aes(x= host, y = slopeRelWeight, 
                         fill = host )) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(Eimeria_species~.)
