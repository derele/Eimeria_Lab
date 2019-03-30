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
# BEFORE cleaning Expe001, we remove it cause weird values, different from all others

# ALL_Expe <- merge(ExpeDF_001, ExpeDF_002, all = T)
ALL_Expe <- merge(ExpeDF_002, ExpeDF_003_4, all = T)
ALL_Expe <- merge(ALL_Expe, ExpeDF_005, all = T)

plotsALL <- makeIntermPlots(ALL_Expe)
plotsALL[[1]] + coord_cartesian(ylim = c(-10,20))
plotsALL[[2]]# + coord_cartesian(ylim = c(-10,20))

# Part 2. Which is, by experiment, dpimaxOPG / dpimaxweightloss, 
# by host strain, by parasite isolate

# BEFORE cleaning Expe001, we remove it cause weird values, different from all others
# toleranceTable <- merge(tolerance_001, tolerance_002, all = T)
toleranceTable <- merge(tolerance_002, tolerance_003_4, all = T)
toleranceTable <- merge(toleranceTable, tolerance_005, all = T)

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
# 1      E.ferrisi (E139)         6.066667                  6
# 2       E.ferrisi (E64)         6.142857                  6
# 3   E.falciformis (E88)         8.266667                  8
# 4 E.falciformis (EfLab)         8.705882                  9
# sdPeakParasite meanPeakHost medianPeakHost sdPeakHost
# 1      0.6914918     5.700000              5   3.302559
# 2      0.4900670     5.057143              5   2.564449
# 3      0.5208305     7.600000              9   3.529189
# 4      0.4696682     8.764706              9   2.537947

# Shift all according to median peak, and keep equal window
mydataShifted <- data.frame(Exp_ID = ALL_Expe$Exp_ID,
                            EH_ID = ALL_Expe$EH_ID,
                            infection_isolate = ALL_Expe$infection_isolate,
                            Mouse_strain = ALL_Expe$Mouse_strain,
                            dpi = ALL_Expe$dpi,
                            startingWeight = ALL_Expe$startingWeight)

# shift all centered on dpi 9 (centered around Eflab)

## weight
W_Eflab <- ALL_Expe[ALL_Expe$infection_isolate == "EfLab", c("dpi", "weight", "EH_ID")]
W_E88 <- ALL_Expe[ALL_Expe$infection_isolate == "E88", c("dpi", "weight", "EH_ID")]
W_E64 <- ALL_Expe[ALL_Expe$infection_isolate == "E64", c("dpi", "weight", "EH_ID")]
W_E64$dpi <- W_E64$dpi + 4 
W_E139 <- ALL_Expe[ALL_Expe$infection_isolate == "E139", c("dpi", "weight", "EH_ID")]
W_E139$dpi <- W_E139$dpi + 4 

weightShifted <- rbind(W_Eflab, W_E88, W_E64, W_E139)

mydataShifted <- merge(mydataShifted, weightShifted)

## oocysts (TO CORRECT FOR MISSING DATA FRANCI NB)
O_Eflab <- ALL_Expe[ALL_Expe$infection_isolate == "EfLab", c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E88 <- ALL_Expe[ALL_Expe$infection_isolate == "E88", c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E88$dpi <- O_E88$dpi + 1 
O_E64 <- ALL_Expe[ALL_Expe$infection_isolate == "E64", c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E64$dpi <- O_E64$dpi + 3 
O_E139 <- ALL_Expe[ALL_Expe$infection_isolate == "E139", c("dpi", "oocysts.per.tube", "fecweight", "EH_ID")]
O_E139$dpi <- O_E139$dpi + 3 

oocystsShifted <- rbind(O_Eflab, O_E88, O_E64, O_E139)

mydataShifted <- merge(mydataShifted, oocystsShifted)

## restrict window to get same info on all
table(oocystsShifted$dpi)
table(weightShifted$dpi)
table(ALL_Expe$dpi)
table(mydataShifted$dpi)

mydataShifted <- mydataShifted[mydataShifted$dpi >= 4,]

## Statistical models along dpi

# OFFSET: 
#https://stats.stackexchange.com/questions/237963/how-to-formulate-the-offset-of-a-glm

### Which distribution?
hist(mydataShifted$weight) # not normal if NMRI are included, bimodal. Can we approx Normal?
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

mydataShifted$oocysts.per.tube.million <- mydataShifted$oocysts.per.tube / 1e6

glmer.nb(oocysts.per.tube.million ~ infection_isolate * Mouse_strain +
           (1|dpi), data = mydataShifted[!is.na(mydataShifted$oocysts.per.tube.million),])

         # + offset(log(fecweight)), 



glmmer(data = mydataShifted, 
       oocysts ~ infection_isolate * Mouse_genotype + 1|dpi + offset(grams))

### Host suffering
glmmer(data = mydataShifted, 
       weight ~ infection_isolate * Mouse_genotype + 1|dpi + offset(originalWeight))

### Tolerance (slope) To think harder better faster stronger...
glmmer(data = mydataShifted, 
       weight ~ infection_isolate * Mouse_genotype + 
         (Mouse_genotype : oocysts) + (infection_isolate : oocysts) +
         1|dpi + offset(originalWeight))

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



