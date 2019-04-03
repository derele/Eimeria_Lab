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

# Part 1: general course of experiment

## Franci, preliminary experiment
plotsExpe1 <- makeIntermPlots(ExpeDF_001)
plotsExpe1[[1]]
plotsExpe1[[2]]
## Passaging in NMRI
plotsExpe2 <- makeIntermPlots(ExpeDF_002)
plotsExpe2[[1]]
plotsExpe2[[2]]
## Parental inbred
plotsExpe3_4 <- makeIntermPlots(ExpeDF_003_4)
plotsExpe3_4[[1]]
plotsExpe3_4[[2]]
## Full design (F0, F1 outbred, F1 hybrids)
plotsExpe5 <- makeIntermPlots(ExpeDF_005[ExpeDF_005$Batch == 1,])
plotsExpe5[[1]]+ coord_cartesian(ylim = c(-10,20))
plotsExpe5[[2]]

## ALL TOGETHER
# CHOICE MADE : remove Ploen mice & NMRI. Variations are too high.
ALL_Expe <- merge(ExpeDF_003_4, ExpeDF_005, all = T)

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

plotsALL <- makeIntermPlots(ALL_Expe)
plotsALL[[1]] + coord_cartesian(ylim = c(-10,20))
plotsALL[[2]]# + coord_cartesian(ylim = c(-10,20))

# Part 2. Which is, by experiment, dpimaxOPG / dpimaxweightloss, 
# by host strain, by parasite isolate
toleranceTable <- merge(tolerance_003_4, tolerance_005, all = T)
toleranceTable <- toleranceTable[!toleranceTable$EH_ID %in% miceDeadBeforeEnd,]

toleranceTable$diff_maxWL_maxOPG <- 
  toleranceTable$dpi_maxweightloss - toleranceTable$dpi_maxOPG

## Experimental design
expedesign <- toleranceTable %>% group_by(infection_isolate, Mouse_genotype)
expedesign <- expedesign %>% summarise(n = n())
expedesign <- expedesign %>% data.frame()
expedesign

## mean of peak day by host and parasite
by_host_parasite_peak <- toleranceTable %>% 
  group_by(infection_isolate) %>% 
  summarise(
    meanPeakParasite = mean(dpi_maxOPG),
    medianPeakParasite = median(dpi_maxOPG),
    sdPeakParasite = sd(dpi_maxOPG),
    meanPeakHost = mean(dpi_maxweightloss),
    medianPeakHost = median(dpi_maxweightloss),
    sdPeakHost = sd(dpi_maxweightloss)
    ) %>% 
  data.frame()
by_host_parasite_peak
# infection_isolate meanPeakParasite medianPeakParasite
# 1    E.ferrisi (E139)         6.125000                  6
# 2     E.ferrisi (E64)         6.240000                  6
# 3 E.falciformis (E88)         8.090909                  8
# sdPeakParasite meanPeakHost medianPeakHost sdPeakHost
# 1      0.7408867     5.250000              5   3.110291
# 2      0.5174506     4.660000              5   2.615573
# 3      0.2942449     8.090909              9   3.544400

# Shift all according to median peak, and keep equal window
mydataShifted <- data.frame(Exp_ID = ALL_Expe$Exp_ID,
                            EH_ID = ALL_Expe$EH_ID,
                            infection_isolate = ALL_Expe$infection_isolate,
                            Mouse_strain = ALL_Expe$Mouse_strain,
                            dpi = ALL_Expe$dpi,
                            startingWeight = ALL_Expe$startingWeight)

# shift all centered on dpi 8 (centered around E88 parasite peak)
## oocysts
O_E88 <- ALL_Expe[ALL_Expe$infection_isolate == "E88", c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E64 <- ALL_Expe[ALL_Expe$infection_isolate == "E64", c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E64$dpi <- O_E64$dpi + 2 
O_E139 <- ALL_Expe[ALL_Expe$infection_isolate == "E139", c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E139$dpi <- O_E139$dpi + 2 

oocystsShifted <- rbind(O_E88, O_E64, O_E139)

## weight
W_E88 <- ALL_Expe[ALL_Expe$infection_isolate == "E88", c("dpi", "weight", "EH_ID")]
W_E88$dpi <- W_E88$dpi - 1 
W_E64 <- ALL_Expe[ALL_Expe$infection_isolate == "E64", c("dpi", "weight", "EH_ID")]
W_E64$dpi <- W_E64$dpi + 3 
W_E139 <- ALL_Expe[ALL_Expe$infection_isolate == "E139", c("dpi", "weight", "EH_ID")]
W_E139$dpi <- W_E139$dpi + 3 

weightShifted <- rbind(W_E88, W_E64, W_E139)

mydataShifted <- merge(mydataShifted, weightShifted)
mydataShifted <- merge(mydataShifted, oocystsShifted)

## restrict window to get same info on all
table(oocystsShifted$dpi)
table(weightShifted$dpi)
table(ALL_Expe$dpi)
table(mydataShifted$dpi)
table(mydataShifted$dpi, mydataShifted$fecweight == 0)

#### Keep WINDOWS OF INFECTION #### 
# peak = dpi 8. 3 days before + peak + 2 days after full mice
mydataShifted <- mydataShifted[mydataShifted$dpi %in% 5:10,]
table(mydataShifted$dpi)

# consider weightBegWind as weight beginning of window
A = mydataShifted[mydataShifted$dpi == 5, c("weight", "EH_ID")]
names(A)[1] = "weightBegWind"
mydataShifted = merge(mydataShifted, A)
# Cumulative sum of our weight loss
mydataShifted = mydataShifted %>% 
  group_by(EH_ID) %>% 
  dplyr::arrange(dpi, .by_group = TRUE) %>%
  dplyr::mutate(relativeWeight = weight / weightBegWind * 100)
mydataShifted = data.frame(mydataShifted)

mydataShifted$OPG <- mydataShifted$oocysts.per.tube / mydataShifted$fecweight

# Summary DF within infection window

### FROM THERE ON
All_summary_windows <- mydataShifted %>% 
  group_by(EH_ID) %>% 
  summarise(minWeight = min(weight),
            maxOo = max(OPG))
All_summary_windows
mydataShifted$weightBegWind

X <- mydataShifted %>% 
  group_by(EH_ID) %>%
  slice(which.min(weight)) 

X
X[c("EH_ID", "weight", "weightBegWind")]

  filter(weight == min(weight))# %>% # filter the data.frame to keep row where x is maximum
  select(dpi) # select column y


# dplyr::arrange(dpi, .by_group = TRUE) %>%
  dplyr::mutate(m = weight / weightBegWind * 100)

mydataShifted$weight

max.loss <- do.call("rbind", by(mydataShifted, mydataShifted$EH_ID, function (x){
  m.loss <- which(mydataShifted$weight == min(mydataShifted$weight, na.rm=TRUE))
  x[m.loss,]
}))
max.loss <- max.loss[!duplicated(max.loss$EH_ID),]
names(max.loss)[names(max.loss) %in% "dpi"] <- "dpi_maxweightloss"
names(max.loss)[names(max.loss) %in% "weight"] <- "minRelativeWeight"
max.loss



# create a table with tolerance factor, max.loss, max.OPG and sum.oocysts concatenated
# makeToleranceTable <- function(ExpeDF,
#                                dayVar = c("dpi", "mean_Neubauer", "OPG", "labels", "Neubauer1",
#                                           "Neubauer2", "Neubauer3", "Neubauer4", "dilution_ml", "weight",                   
#                                           "weightloss", "weightRelativeToInfection", "fecweight", 
#                                           "comment.x", "comment.y", "X.x", "relativeWeight", 
#                                           "X", "X.1", "X.2", "X.3", "X.4", "MM",
#                                           "oocysts.per.tube", "label")){
  # maximum weightloss
  max.loss = getMaxLoss(ExpeDF)
  max.loss = max.loss[!names(max.loss) %in% dayVar]
  # max and sum of oocysts
  max.OPG = getMaxOPG(ExpeDF)
  max.OPG = max.OPG[!names(max.OPG) %in% dayVar]
  # sum
  sum.oocysts = getSumOPG(ExpeDF)
  sum.oocysts = sum.oocysts[!names(sum.oocysts) %in% dayVar]
  # tolerance
  tolerance = merge(max.loss, max.OPG)
  tolerance = merge(tolerance, sum.oocysts)
  # maxweightloss
  tolerance$maxweightloss = 100 - tolerance$minRelativeWeight
  # REMOVE UNINFECTED ANIMALS
  tolerance <- tolerance[tolerance$maxOPG_inmillion != 0,]
  # http://science.sciencemag.org/content/318/5851/812
  # TO CORRECT LATER
  # !!!!!
  # resistance (inverse of peak parasite density) 
  # tolerance (slope of a regression of minimum weight retained against peak parasite load)
  tol = -(tolerance$maxweightloss / tolerance$maxOPG_inmillion)
  tolerance$toleranceFactor = tol + abs(min(tol))
  # erase useless level (Eflab)
  tolerance$infection_isolate <- droplevels(tolerance$infection_isolate)
  # NB. Let's not consider which parent is which, but make A_B mouse = B_A mouse
  # we don't have enough individuals to test this effect, and we are not interested in it anyway!
  x <- strsplit(tolerance$Mouse_strain, "_")
  y <- lapply(x, sort)
  z <- unlist(lapply(y, FUN = function(x){paste(x, collapse="-")}))
  tolerance$Mouse_genotype <- z
  ## Order the levels to be more clear in later plots (parents will be low and down 
  ## on the legend, hybrids in between...)
  tolerance$Mouse_genotype <- factor(tolerance$Mouse_genotype,
                                     levels = c("NMRI", 
                                                "WSB", "WP", "PWD1",
                                                "BUSNA-BUSNA", "PWD-PWD",
                                                "BUSNA-PWD",
                                                "BUSNA-STRA","PWD-SCHUNT",
                                                "SCHUNT-STRA",
                                                "SCHUNT-SCHUNT","STRA-STRA"),
                                     labels = c("NMRI", 
                                                "MMd_F0 (Ws-Ws)", "Mmm-Mmd_Hybrid (WP)", "MMm_F0 (Pw1-Pw1)",
                                                "MMm_F0 (Bu-Bu)",
                                                "MMm_F0 (Pw-Pw)",
                                                "MMm_F1 (Bu-Pw)",
                                                "Mmm-Mmd_F1 (Bu-St)",
                                                "Mmm-Mmd_F1 (Pw-Sc)",
                                                "MMd_F1 (Sc-St)",
                                                "MMd_F0 (Sc-Sc)",
                                                "MMd_F0 (St-St)"))
  tolerance$infection_isolate <- factor(tolerance$infection_isolate,
                                        levels = c("E139", "E64", "E88", "EfLab"),
                                        labels = c("E.ferrisi (E139)",
                                                   "E.ferrisi (E64)",
                                                   "E.falciformis (E88)",
                                                   "E.falciformis (EfLab)"))
  return(tolerance)
# }





## Statistical models along dpi

# OFFSET: 
#https://stats.stackexchange.com/questions/237963/how-to-formulate-the-offset-of-a-glm

### Which distribution?
hist(mydataShifted$weight)
hist(mydataShifted$oocysts.per.tube, breaks = 100) # clearly negbin here :) 

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
library(lme4)

# negative binomial requires INTEGER 
range(mydataShifted$oocysts.per.tube[mydataShifted$oocysts.per.tube > 0], na.rm = T)
mydataShifted$oocysts.per.tube.tsd <- as.integer(mydataShifted$oocysts.per.tube / 1000)

# use of random factors : we assume a different baseline response per day
ggplot(mydataShifted, aes(x = dpi, y = oocysts.per.tube +1, group = dpi)) +
  geom_boxplot() + scale_y_log10()

ggplot(mydataShifted, aes(x = dpi, y = weight, group = dpi)) +
  geom_boxplot() 

# test significance of interactions
# Original model
model1 <- glmer.nb(oocysts.per.tube.tsd ~ infection_isolate * Mouse_strain + 
                     (1|dpi) + offset(log(fecweight)),
                   data = mydataShifted[mydataShifted$fecweight != 0,])

# Model without an interaction
model2 <- glmer.nb(oocysts.per.tube.tsd ~ infection_isolate + Mouse_strain + 
                     (1|dpi) + offset(log(fecweight)), 
                   data = mydataShifted[mydataShifted$fecweight != 0,])

# Likelihood ratio test
anova(model1, model2, test="LRT")
## --> INTERACTIONS SLIGHTLY SIGNIFICANT
summary(model1)

library(sjPlot)
# I don't get that plot
# sjPlot::plot_model(model1, sort.est = T, show.values = TRUE, value.offset = .3)

# Post-hoc analysis
if(!require(multcomp)){install.packages("multcomp")}
library(multcomp)
s1 <- summary(glht(model1, mcp(infection_isolate="Tukey")))
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# E139 - E64 == 0  -1.0111     0.2997  -3.373  0.00207 **
#   E88 - E64 == 0    0.4042     0.4589   0.881  0.64587   
# E88 - E139 == 0   1.4153     0.4801   2.948  0.00860 **
s2 <- summary(glht(model1, mcp(Mouse_strain="Tukey")))
# only significant : 
# SCHUNT_SCHUNT - PWD_PWD == 0     -1.05978    0.27068  -3.915    <0.01 **

# and summary (sum of oocysts during window)


### Host suffering

# test significance of interactions
# Original model
mydataShifted[is.na(mydataShifted$weight),]

model1.w <- lmer(weight ~ infection_isolate * Mouse_strain + 
                   (1|dpi) + offset(log(weightBegWind)), 
                 data = mydataShifted)
# Model without an interaction
model2.w <- lmer(weight ~ infection_isolate + Mouse_strain + 
                    (1|dpi) + offset(log(weightBegWind)), 
                  data = mydataShifted)
# Likelihood ratio test
anova(model1.w, model2.w, test="LRT")

## --> INTERACTIONS SIGNIFICANT +++ p-value 4.194e-10 ***
summary(model1.w)

# Post-hoc analysis
if(!require(multcomp)){install.packages("multcomp")}
library(multcomp)
s.P.w <- summary(glht(model1.w, mcp(infection_isolate="Tukey")))
s.P.w # NOPE
s.H.w <- summary(glht(model1.w, mcp(Mouse_strain="Tukey")))
s.H.w # LOADS OF DIFFERENCES!!

### SUMMARY OVER EXPE (to compare)


m1 = lm(response ~ group * infection_isolate,
        data=mydata)
mytukey = TukeyHSD.lm(m1)
mytukeyDF = lapply(mytukey, function(x) {
  df <- data.frame(x)
  data.frame(combi = rownames(df[df$p.adj < 0.05,]),
             padj = df$p.adj[df$p.adj < 0.05],
             diff = df$p.adj[df$p.adj < 0.05])})
return(list(model = m1,
            summary.model = summary(m1),
            mytukeyHSD = mytukeyDF))


### Tolerance (slope) Raberg 2007 Raber 2008

ggplot(mydataShifted, aes(x = oocysts.per.tube.tsd, y = weight, col = EH_ID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme(legend.position = "none")

## WRONG if we use full model. Too many assumptions. weight before peak and after peak shall 
## vary, not link with oocysts. The points before and after the peak are NOT independant.
## we can't treat them as repeated measures. 
## Hence, ONE tolerance value per animal. Defined by Raberg 2007.

#### Test interactions
modT1 <- lmer(weight ~ oocysts.per.tube.tsd + 
                infection_isolate + Mouse_strain +
                oocysts.per.tube.tsd : infection_isolate +
                oocysts.per.tube.tsd : Mouse_strain + 
                (1|oocysts.per.tube.tsd:EH_ID) + 
                offset(log(startingWeight)), 
              data = mydataShifted)

modT1

modT2 <- lmer(weight ~ oocysts.per.tube.tsd + 
                infection_isolate + Mouse_strain +
                oocysts.per.tube.tsd : Mouse_strain + 
                (1|EH_ID) + offset(log(startingWeight)), 
              data = mydataShifted)
modT3 <- lmer(weight ~ oocysts.per.tube.tsd + 
                infection_isolate + Mouse_strain +
                oocysts.per.tube.tsd : infection_isolate +
                (1|EH_ID) + offset(log(startingWeight)), 
              data = mydataShifted)
modT4 <- lmer(weight ~ oocysts.per.tube.tsd + 
                infection_isolate + Mouse_strain +
                (1|EH_ID) + offset(log(startingWeight)), 
              data = mydataShifted)
anova(modT1, modT2, test="LRT") # NOT SIGNIF DIFF WHEN DROP INF ISOLATE
# --> choose T2
anova(modT2, modT3, test="LRT") # SIGNIF --> MOUSE STRAIN INTERACTION
# --> choose T2
anova(modT2, modT4, test="LRT") # SIGNIF
# --> choose T2

summary(modT2)
summary(glht(modT2, mcp(Mouse_strain="Tukey")))
summary(glht(modT2, mcp(infection_isolate="Tukey")))

# How to plot data???
# devtools::install_github("sjPlot/devel")
# library(sjmisc)
# library(lme4)
# library(sjPlot)

library(lsmeans)
library(multcompView)

lsm <- lsmeans(modT2, ~ oocysts.per.tube.tsd + 
                 infection_isolate + Mouse_strain +
                 oocysts.per.tube.tsd : Mouse_strain, 
               at=list(startingWeight=1))

a <- cld(lsm, type="response", sort=TRUE)

plot(a)

plot(lsm, type = "response") 
median(mydataShifted$startingWeight)
by(data = mydataShifted, mydataShifted$Mouse_strain, 
   function(x) mean(x$weight, na.rm = T))

lsm

plot_model(modT2, type = "pred", 
           terms = c("oocysts.per.tube.tsd", "startingWeight","Mouse_strain"))

mydataShifted$startingWeight

weight ~ oocysts.per.tube.tsd + 
  infection_isolate + Mouse_strain +
  oocysts.per.tube.tsd : Mouse_strain + 
  (1|EH_ID) + offset(log(startingWeight)

# mydataShifted$percentWeight <- mydataShifted$weight / mydataShifted$startingWeight

# an example
ggplot(mydataShifted[mydataShifted$EH_ID == "mouse_7",], 
       aes(x = oocysts.per.tube.tsd, y = percentWeight,
           shape = Mouse_strain)) +
  geom_point() +
  geom_smooth(aes(col = Mouse_strain), method = "lm", se = F)

library(sjPlot)
sjPlot::plot_model(modT2, sort.est = T)


ggplot(mydataShifted, 
       aes(x = oocysts.per.tube.tsd, y = percentWeight,
           shape = Mouse_strain, col = Mouse_strain)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_shape_manual(values = 1:20) +
  facet_grid(.~infection_isolate) 


library(nlme)
fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1|Subject)

lmer(weight ~ oocysts.per.tube.tsd + 
       infection_isolate + Mouse_strain +
       oocysts.per.tube.tsd : Mouse_strain + 
       (1|EH_ID) + offset(log(startingWeight)), 
     data = mydataShifted)

newdat <- expand.grid(infection_isolate=unique(mydataShifted$infection_isolate),
                      Mouse_strain= unique(mydataShifted$Mouse_strain),
                      oocysts.per.tube.tsd=c(min(mydataShifted$oocysts.per.tube.tsd),
                            max(mydataShifted$oocysts.per.tube.tsd)))

library(ggplot2)
p <- ggplot(mydataShifted, aes(x=oocysts.per.tube.tsd, y=weight, colour=Mouse_strain)) +
  geom_point(size=3) +
  geom_line(aes(y=predict(modT2), group=EH_ID, size="EH_ID")) +
  geom_line(data=newdat, aes(y=predict(modT2, level=0, newdata=newdat), size="Population")) +
  scale_size_manual(name="Predictions", values=c("EH_ID"=0.5, "Population"=3)) +
  theme_bw(base_size=22) 
print(p)


## End statistical models along dpi


ggplot(toleranceTable, aes(y = dpi_maxOPG, x = infection_isolate))+
  geom_violin() +
  geom_jitter(aes(fill = Mouse_genotype),
              pch = 21, height = 0.1, width = 0.2, size = 4, alpha = 0.9) 

ggplot(toleranceTable, aes(y = dpi_maxweightloss, x = infection_isolate))+
  geom_violin() +
  geom_jitter(aes(fill = Mouse_genotype),
              pch = 21, height = 0.1, width = 0.2, size = 4, alpha = 0.9) 

## CCL --> here already, host genotype seems to not influence the parasite
## reproductive output; But the peak of weight loss seems highly influenced by host!

## 


diffDaysDf <- data.frame(table(toleranceTable$diff_maxWL_maxOPG, 
                               toleranceTable$infection_isolate))

names(diffDaysDf) <- c("diff_maxWL_maxOPG", 
                       "infection_isolate", "Freq")

pdf("figures/fig1.pdf", width = 10, height = 6)
ggplot(toleranceTable, aes(x = infection_isolate, 
                       y = diff_maxWL_maxOPG, fill = infection_isolate ))+
  geom_violin() +
  geom_dotplot(dotsize = 0.5, binaxis = "y", stackdir = "center", fill = "black") +
  scale_y_continuous(breaks = -9:6) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = -1, linetype = "dotted") +
  ylab("Days between shedding peak and maximum weight loss") +
  xlab("")
dev.off()
  
## Results --> Eferrisi max shedding peak can be shifted of -1 day,
## Efalci of +1 day
## Then we align and compare



## Separated ##

## Test of hybrid effect
tolerance_001$HybridStatus <- "inbred"
tolerance_001[tolerance_001$Mouse_strain == "WP", "HybridStatus"] <- "hybrids"

### 1.1 On weight loss
myboxPlot(tolerance = tolerance_001,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "HybridStatus", groupName = "Hybrid status") +
  scale_fill_manual(values = c("darkorchid", "yellow", "pink"))

myMod(mydata = tolerance_001,
      response = tolerance_001$maxweightloss, 
      group = tolerance_001$HybridStatus)

### 1.2 On max of oocysts
myboxPlot(tolerance = tolerance_001,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "HybridStatus", groupName = "Hybrid status") +
  scale_fill_manual(values = c("darkorchid", "yellow", "pink"))

myMod(mydata = tolerance_001,
      response = tolerance_001$maxOPG_inmillion, 
      group = tolerance_001$HybridStatus)

### 1.2.3 On tolerance
myboxPlot(tolerance = tolerance_001,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "HybridStatus", groupName = "Hybrid status") 

myMod(mydata = tolerance_001,
      response = tolerance_001$toleranceFactor, 
      group = tolerance_001$HybridStatus)

## Part 2. Strain effect

tolerance_001_E64 <- tolerance_001[tolerance_001$infection_isolate %in% "E.ferrisi (E64)",] 

### 1.1 On weight loss
pdf("fig1.pdf", width = 4, height = 5)
myboxPlot(tolerance = tolerance_001_E64,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "Mouse_genotype", groupName = "Mouse genotype")
dev.off()
# Mmm-Mmd_Hybrid (WP)-MMd_F0 (Ws-Ws) 3.753162e-05 3.753162e-05
# MMm_F0 (Pw-Pw)-MMd_F0 (Ws-Ws) 2.900146e-04 2.900146e-04

myMod(mydata = tolerance_001,
      response = tolerance_001$maxweightloss, 
      group = tolerance_001$Mouse_genotype)

### 1.2 On sum of oocysts shed along expe
pdf("fig2.pdf", width = 4, height = 5)
myboxPlot(tolerance = tolerance_001_E64,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "Mouse_genotype", groupName = "Mouse genotype")
dev.off()
# MMm_F0 (Pw-Pw)-Mmm-Mmd_Hybrid (WP) 0.04848581 0.04848581

myMod(mydata = tolerance_001,
      response = tolerance_001$maxOPG_inmillion, 
      group = tolerance_001$Mouse_genotype)

### 1.2.3 On tolerance
pdf("fig3.pdf", width = 4, height = 5)
myboxPlot(tolerance = tolerance_001_E64,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "Mouse_genotype", groupName = "Mouse genotype") 
dev.off()
#nothing significant
myMod(mydata = tolerance_001,
      response = tolerance_001$toleranceFactor, 
      group = tolerance_001$Mouse_genotype)

### Part 2: Vivian
### comparison tolerance for Vivian
tolerance_comparison <- merge(tolerance_003_4, tolerance_005, all = T)
# Remove individuals that DIDN'T shed oocysts at all! (no infection)
tolerance_comparison <- tolerance_comparison[tolerance_comparison$sum.opg != 0,]
# Let's remove double infections (Exp_005_2a and Exp_005_2b)
tolerance_comparison <- tolerance_comparison[
  !tolerance_comparison$Exp_ID %in% c("Exp_005_2a", "Exp_005_2b"),]

## Check first strain
ggplot(tolerance_comparison, aes(maxOPG_inmillion, maxweightloss)) +
  geom_point(size = 5, shape = 21, alpha = 0.6,
             aes(fill = Mouse_genotype)) +
  facet_grid(.~infection_isolate, scales = "free_x") +
  scale_x_log10() # x axis in log 10 to see the fold differences

#######################
# Statistical modelling
#######################

## Preliminary : can we merge E139 and E64 as E.ferrisi?? --> YES
# tolerance_comparisonFERRISI <- tolerance_comparison[tolerance_comparison$infection_isolate != "E.falciformis (E88)",]
# myMod(mydata = tolerance_comparisonFERRISI,
#       response = tolerance_comparisonFERRISI$maxweightloss, 
#       group = tolerance_comparisonFERRISI$Mouse_genotype)
# myMod(mydata = tolerance_comparisonFERRISI,
#       response = tolerance_comparisonFERRISI$maxOPG, 
#       group = tolerance_comparisonFERRISI$Mouse_genotype)
# myMod(mydata = tolerance_comparisonFERRISI,
#       response = tolerance_comparisonFERRISI$toleranceFactor, 
#       group = tolerance_comparisonFERRISI$Mouse_genotype)
## -> no difference between E139 and E64 at the genotype level for all 3 measurements

tolerance_comparison$infection_isolate <- 
  as.character(tolerance_comparison$infection_isolate)

tolerance_comparison[
  tolerance_comparison$infection_isolate %in% 
    c("E.ferrisi (E64)", "E.ferrisi (E139)"), "infection_isolate"] <- "E.ferrisi (E64|E139)"

table(tolerance_comparison$HybridStatus,
      tolerance_comparison$infection_isolate)

## Part 1. Hybrid effect

### 1.1 On weight loss
myboxPlot(tolerance = tolerance_comparison,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "HybridStatus", groupName = "Hybrid status") +
  scale_fill_manual(values = c("darkorchid", "yellow", "pink"))

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxweightloss, 
           group = tolerance_comparison$HybridStatus)

### 1.2 On max of oocysts
myboxPlot(tolerance = tolerance_comparison,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "HybridStatus", groupName = "Hybrid status")

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxOPG_inmillion, 
           group = tolerance_comparison$HybridStatus)

### 1.2.3 On tolerance
myboxPlot(tolerance = tolerance_comparison,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "HybridStatus", groupName = "Hybrid status") 

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$toleranceFactor, 
           group = tolerance_comparison$HybridStatus)

## Part 2. Strain effect

table(tolerance_comparison$Mouse_genotype,
      tolerance_comparison$infection_isolate)

### 1.1 On weight loss
myboxPlot(tolerance = tolerance_comparison,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "Mouse_genotype", groupName = "Mouse genotype")

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxweightloss, 
      group = tolerance_comparison$Mouse_genotype)

### 1.2 On sum of oocysts shed along expe
myboxPlot(tolerance = tolerance_comparison,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "Mouse_genotype", groupName = "Mouse genotype")

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxOPG_inmillion, 
      group = tolerance_comparison$Mouse_genotype)

### 1.2.3 On tolerance
myboxPlot(tolerance = tolerance_comparison,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "Mouse_genotype", groupName = "Mouse genotype") 
myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$toleranceFactor, 
           group = tolerance_comparison$Mouse_genotype)



