## Data preparation
source("functions.R")

library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

mytheme <- theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        strip.text = element_text(size=25))

# ExpeDF need at least the following: 
# "OPG", "EH_ID", "dpi", "weight", "EH_id", "infection_isolate", "Mouse_strain"  

########################### Exp001 : April 2017 PWD, WSD and hybrids infection 
# by E.falciformis (Eflab) and E.ferrisi (E64)
ExpeDF_001 <- read.csv("../data/3_recordingTables/Exp001_May2017_crossing_infection.csv")

# Remove Eflab, mice suffered too much and had to be sacrificed earlier
# ExpeDF_001 <- ExpeDF_001[ExpeDF_001$infection_isolate %in% "EI64",]
# ExpeDF_001$infection_isolate <- droplevels(ExpeDF_001$infection_isolate)

# Add transect info
ExpeDF_001$transect <- "Commercial mice strains"

# Add mice subspecies info
ExpeDF_001$Mouse_subspecies <- "F1 hybrids"
ExpeDF_001$Mouse_subspecies[ExpeDF_001$Mouse_strain == "PWD"] <- "M.m.musculus"
ExpeDF_001$Mouse_subspecies[ExpeDF_001$Mouse_strain == "WSB"] <- "M.m.domesticus"

ExpeDF_001$Mouse_subspecies <- factor(ExpeDF_001$Mouse_subspecies, 
                                      levels = c("M.m.domesticus", "F1 hybrids", "M.m.musculus"))

# Mouse_strain: West should always be left 
ExpeDF_001$Mouse_strain <- factor(ExpeDF_001$Mouse_strain, 
                                  levels = c("WSB", "WP", "PWD"),
                                  labels = c("M.m.domesticus \n(WSB)", 
                                             "F1 hybrids \n(WP)", 
                                             "M.m.musculus \n(PWD)"))

# Keep only mice that survived 11 days
# survivors <- names(table(ExpeDF_001$EH_ID))[table(ExpeDF_001$EH_ID) %in% 12]
# ExpeDF_001 <- ExpeDF_001[ExpeDF_001$EH_ID %in% survivors,]
ExpeDF_001 <- calculateWeightLoss(ExpeDF_001, infectionDay = 1)

ExpeDF_001$infection_isolate

########################### Pass001: Nov 2017, passaging 4 isolates (some missing data)
# (Eflab, E88, E139, E64) in NMRI. 2 mice per cage. Only OPG recorded
PassDF_001 <- read.csv("../data/3_recordingTables/passaging_extra/Pass001_oocystsonly_Nov2017_Passaging_4Eimeria.csv")
PassDF_001$weightRelativeToInfection <- NA

########################### Exp002: March 2018 NMRI infected with 4 strains
ExpeDF_002 <- read.csv("../data/3_recordingTables/Exp002_March2018_NMRI_4strains_RECORDweightAndOocysts.csv")
ExpeDF_002$Mouse_strain <- "NMRI"
ExpeDF_002$EH_ID <- paste0("mouse_", ExpeDF_002$EH_ID)

ExpeDF_002$Mouse_subspecies = ExpeDF_002$Mouse_strain

# Calculate weight loss
ExpeDF_002 <- calculateWeightLoss(ExpeDF_002)

# Calculate OPG
ExpeDF_002 <- calculateOPG(ExpeDF_002)

ExpeDF_002$infection_isolate

# WRONG OPG just for plotting, for dpi 11 (we forgot to weight the feces...)
ExpeDF_002$OPG_temp[
  !is.na(ExpeDF_002$mean_Neubauer) & is.na(ExpeDF_002$OPG)] = 0

ExpeDF_002$OPG_temp[
  !is.na(ExpeDF_002$mean_Neubauer) & is.na(ExpeDF_002$OPG) & ExpeDF_002$mean_Neubauer > 0] = 1

ExpeDF_002$OPG[is.na(ExpeDF_002$OPG)] = ExpeDF_002$OPG_temp[is.na(ExpeDF_002$OPG)]

########################### Exp003 : May 2018 batch 1
oo <- read.csv("../data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDoocysts.csv")
we <- read.csv("../data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDweight.csv")
design <- read.csv("../data/2_designTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv")

ExpeDF_003 <- merge(oo, we, all = T)
ExpeDF_003 <- merge(ExpeDF_003, design, by = "EH_ID", all = T)
rm(design, oo, we)

# Keep for dpi 0 to 11
ExpeDF_003 <- ExpeDF_003[ExpeDF_003$dpi %in% 0:11, ]# remove stabilisation period

# correct abherante value
ExpeDF_003[ExpeDF_003$EH_ID %in% "LM0137" & ExpeDF_003$weight %in% 17.6, "weight"] <- NA

# Add transect
ExpeDF_003$transect <- "Commercial mice strains"
ExpeDF_003$transect[ExpeDF_003$Mouse_strain %in% c("BUSNA", "STRA")] <- "HMHZ"

# Add mice subspecies
ExpeDF_003$Mouse_subspecies <- "M.m.domesticus"
ExpeDF_003$Mouse_subspecies[ExpeDF_003$Mouse_strain %in% c("BUSNA", "PWD")] <- "M.m.musculus"

# Mouse_strain: West should always be left 
ExpeDF_003$Mouse_strain <- factor(ExpeDF_003$Mouse_strain,
                                  levels = c("STRA", "BUSNA", "SCHUNT", "PWD"),
                                  labels = c("M.m.domesticus \n(STRA)", 
                                             "M.m.musculus \n(BUSNA)", 
                                             "M.m.domesticus \n(SCHUNT)",
                                             "M.m.musculus \n(PWD)"))

# Calculate weight loss
ExpeDF_003 <- calculateWeightLoss(ExpeDF_003)

#Calculate OPG NOT DONE YET ;)
names(ExpeDF_003)[names(ExpeDF_003) %in% paste0("oocyst_sq", 1:4)] <- paste0("Neubauer", 1:4)
names(ExpeDF_003)[names(ExpeDF_003) %in% "dilution"] <- "dilution_ml"

ExpeDF_003 <- calculateOPG(ExpeDF_003)

ExpeDF_003$infection_isolate

########################### Exp004 : May 2018 batch 2
oo <- read.csv("../data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDoocysts.csv")
we <- read.csv("../data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDweight.csv")
design <- read.csv("../data/2_designTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_DESIGN.csv")

ExpeDF_004 <- merge(oo, we, all = T)
ExpeDF_004 <- merge(ExpeDF_004, design, by = "EH_ID", all = T)
rm(design, oo, we)

# Correct name
names(ExpeDF_004)[names(ExpeDF_004) == "Strain"] <- "Mouse_strain"

# Add transect
ExpeDF_004$transect <- "Commercial mice strains"
ExpeDF_004$transect[ExpeDF_004$Mouse_strain %in% c("BUSNA", "STRA")] <- "HMHZ"

# Add mice subspecies
ExpeDF_004$Mouse_subspecies <- "M.m.domesticus"
ExpeDF_004$Mouse_subspecies[ExpeDF_004$Mouse_strain %in% c("BUSNA", "PWD")] <- "M.m.musculus"

# Mouse_strain: West should always be left 
ExpeDF_004$Mouse_strain <- factor(ExpeDF_004$Mouse_strain,
                                  levels = c("STRA", "BUSNA", "SCHUNT", "PWD"),
                                  labels = c("M.m.domesticus \n(STRA)", 
                                             "M.m.musculus \n(BUSNA)", 
                                             "M.m.domesticus \n(SCHUNT)",
                                             "M.m.musculus \n(PWD)"))

# Calculate weight loss
ExpeDF_004 <- calculateWeightLoss(ExpeDF_004) 

# Calculate OPG NOT DONE YET ;)
# ExpeDF_004 <- calculateOPG(ExpeDF_004)

# Remove animals that died before the end of the experiment
tab <- table(ExpeDF_004$EH_ID[!is.na(ExpeDF_004$weight)])

ExpeDF_004 <- ExpeDF_004[!ExpeDF_004$EH_ID %in% names(tab)[tab < 12],]


########################### Exp004 : May 2018 batch 2
oo <- read.csv("../data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDoocysts.csv")
we <- read.csv("../data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDweight.csv")
design <- read.csv("../data/2_designTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_DESIGN.csv")

ExpeDF_004 <- merge(oo, we, all = T)
ExpeDF_004 <- merge(ExpeDF_004, design, by = "EH_ID", all = T)
rm(design, oo, we)

# Correct name
names(ExpeDF_004)[names(ExpeDF_004) == "Strain"] <- "Mouse_strain"

# Add transect
ExpeDF_004$transect <- "Commercial mice strains"
ExpeDF_004$transect[ExpeDF_004$Mouse_strain %in% c("BUSNA", "STRA")] <- "HMHZ"

# Add mice subspecies
ExpeDF_004$Mouse_subspecies <- "M.m.domesticus"
ExpeDF_004$Mouse_subspecies[ExpeDF_004$Mouse_strain %in% c("BUSNA", "PWD")] <- "M.m.musculus"

# Mouse_strain: West should always be left 
ExpeDF_004$Mouse_strain <- factor(ExpeDF_004$Mouse_strain,
                                  levels = c("STRA", "BUSNA", "SCHUNT", "PWD"),
                                  labels = c("M.m.domesticus \n(STRA)", 
                                             "M.m.musculus \n(BUSNA)", 
                                             "M.m.domesticus \n(SCHUNT)",
                                             "M.m.musculus \n(PWD)"))

# Calculate weight loss
ExpeDF_004 <- calculateWeightLoss(ExpeDF_004) 

# Calculate OPG NOT DONE YET ;)
# ExpeDF_004 <- calculateOPG(ExpeDF_004)

# Remove animals that died before the end of the experiment
tab <- table(ExpeDF_004$EH_ID[!is.na(ExpeDF_004$weight)])

ExpeDF_004 <- ExpeDF_004[!ExpeDF_004$EH_ID %in% names(tab)[tab < 12],]

########################### Exp005 : Full expe including hybrids, July 2018 
oo <- read.csv("../data/3_recordingTables/Exp005_full_RECORDoocysts.csv", na.strings = c("NA", " "))
we <- read.csv("../data/3_recordingTables/Exp005_full_RECORDweight.csv", na.strings = c("NA", " "))
design <- rbind(read.csv("../data/2_designTables/Inf1a_Exp005.DESIGN.csv", na.strings = c("NA", " ")),
                read.csv("../data/2_designTables/Inf2a_Exp005.DESIGN.csv", na.strings = c("NA", " ")),
                read.csv("../data/2_designTables/Inf1b_Exp005.DESIGN.csv", na.strings = c("NA", " ")),
                read.csv("../data/2_designTables/Inf2b_Exp005.DESIGN.csv", na.strings = c("NA", " ")))

# Correct error
names(design)[names(design) == "EH_id"] <- "EH_ID"
# Correct space error
design$EH_ID <- gsub(" ", "", design$EH_ID)
design$HybridStatus <- gsub(" ", "", design$HybridStatus)

ExpeDF_005 <- merge(oo, we, by = c("labels", "Expe"))

ExpeDF_005 <- merge(ExpeDF_005, design, by = "EH_ID", all.x = T)

## Correct error space
ExpeDF_005$Strain <- gsub(" ", "", ExpeDF_005$Strain)

ExpeDF_005$infection_isolate <- ExpeDF_005$Eimeria

# Correct error non numeric
ExpeDF_005$weight <- as.numeric(as.character(ExpeDF_005$weight))
ExpeDF_005$weight_dpi0 <- as.numeric(as.character(ExpeDF_005$weight_dpi0))

# Correct weight loss
ExpeDF_005$weightloss <- ExpeDF_005$weight_dpi0 - ExpeDF_005$weight

# Calculate weight Normalized to dpi0
ExpeDF_005$weightNormalized <- ExpeDF_005$weight / ExpeDF_005$weight_dpi0 * 100

# Add mice subspecies
ExpeDF_005$Mouse_subspecies <- ExpeDF_005$HybridStatus

# tapply(as.numeric(ExpeDF_005$weightloss), list(ExpeDF_005$dpi, ExpeDF_005$HybridStatus), FUN = function(x)mean(x, na.rm = T))

# Mouse_strain: West should always be left 
ExpeDF_005$Mouse_strain <- factor(as.factor(ExpeDF_005$Strain),
                                  levels = c("BUSNA_BUSNA","BUSNA_PWD","BUSNA_STRA","PWD_BUSNA",
                                             "PWD_PWD","PWD_SCHUNT","SCHUNT_PWD","SCHUNT_SCHUNT",
                                             "SCHUNT_STRA","STRA_BUSNA","STRA_SCHUNT","STRA_STRA"),
                                  labels = c("M.m.musculus P", "M.m.musculus F1",
                                             "Hybrid", "M.m.musculus F1",
                                             "M.m.musculus P","Hybrid",
                                             "Hybrid", "M.m.domesticus P",
                                             "M.m.domesticus F1","Hybrid",
                                             "M.m.domesticus F1","M.m.domesticus P"))

# Age at infection
ExpeDF_005$ageAtInfection[ExpeDF_005$Batch == 1] <- round(ExpeDF_005$ageAtdpi0expe1a)
ExpeDF_005$ageAtInfection[ExpeDF_005$Batch == 2] <- round(ExpeDF_005$ageAtdpi0expe1a +2)

########## Exclude potential covariates ########
dfcov <- ExpeDF_005[ExpeDF_005$dpi == 0,]

# age
ggplot(dfcov, aes(y = ageAtInfection, x = HybridStatus)) +
  geom_boxplot(col = "grey") +
  geom_point(position=position_dodge(width=0.5), 
             aes(group=EH_ID, col = HybridStatus), 
             size = 3) + mytheme

# original weight
ggplot(dfcov, aes(y = weight_dpi0, x = HybridStatus)) +
  geom_boxplot(col = "grey") +
  geom_point(position=position_dodge(width=0.5), 
             aes(group=EH_ID, col = HybridStatus), 
             size = 3) + mytheme

########## Plot weight ########
ExpeDF_005$OPG_plot[is.na(ExpeDF_005$OPG)] = "na"
ExpeDF_005$OPG_plot[!is.na(ExpeDF_005$OPG) & ExpeDF_005$OPG > 0] = "positive"
ExpeDF_005$OPG_plot[!is.na(ExpeDF_005$OPG) & ExpeDF_005$OPG == 0] = "negative"
ExpeDF_005$OPG_plot = as.factor(ExpeDF_005$OPG_plot)
  
# Enter Eimeria species
ExpeDF_005$Eimeria_species[ExpeDF_005$infection_isolate %in% c("E88", "Eflab", "EfLab")] = "E.falciformis"
ExpeDF_005$Eimeria_species[ExpeDF_005$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
  
ExpeDF_005[is.na(ExpeDF_005$Eimeria_species),]

ggplot(ExpeDF_005, aes(x=dpi, y=weight, colour=HybridStatus)) + 
  geom_line(aes(group=EH_ID)) +
  geom_point(aes(pch = Expe), size = 5) +
  mytheme +
  facet_grid(. ~ Eimeria_species) +
  theme(strip.text.y = element_text(size = 15))

mysum005 <- summarySE(ExpeDF_005[!is.na(ExpeDF_005$weightloss),], measurevar="weightloss", 
                      groupvars=c("Eimeria_species", "HybridStatus","dpi"))

pd <- position_dodge(0.2) # move them .05 to the left and right

plot005we <- ggplot(mysum005, aes(x = dpi, y = -weightloss, colour=HybridStatus)) + 
  geom_errorbar(aes(ymin=-weightloss-ci, ymax=-weightloss+ci), 
                width=1, position=pd) +
  geom_line(position=pd, size = 2) +
  geom_point(position=pd) +
  mytheme +
  facet_grid(. ~ Eimeria_species) +
  theme(strip.text.y = element_text(size = 15))

plot005we

## And if we take the same starting weight?
mysum005

ExpeDF_005$HybridStatus <- as.factor(ExpeDF_005$HybridStatus)

########## Stats weight ########
library(lme4)
library(lmerTest)
# Add age and sex
myfitE88 <- lmer(weight ~ HybridStatus * dpi + 
                (1|EH_ID) + (1|ageAtInfection) + (1|Sex), ExpeDF_005[ExpeDF_005$infection_isolate =="E88",])
anova(myfitE88)

# Plot
library(effects)
plot(Effect(c("dpi", "HybridStatus"),myfitE88))

# 
# myfitE64 <- glmer(weight ~ HybridStatus * dpi + 
#                 (1|EH_ID) + (1|ageAtInfection) + (1|Sex) + (1|dpi), 
#                 family = poisson, 
#                 data = ExpeDF_005[ExpeDF_005$infection_isolate =="E64",])

myfitE64 <- glmer(weight ~ HybridStatus * dpi + 
                    (1|EH_ID), 
                  family = poisson, 
                  data = ExpeDF_005[ExpeDF_005$infection_isolate =="E64",])
anova(myfitE64)

# Plot
library(effects)
plot(Effect(c("dpi", "HybridStatus"),myfitE64))


#http://lme4.r-forge.r-project.org/slides/2011-01-11-Madison/6NLMMH.pdf

# https://uknowledge.uky.edu/cgi/viewcontent.cgi?article=1006&context=epb_etds
# Non linear mixed effect model!!!
library(nlme)

# Repeated measurements on a group of individuals: to account for within subject
# as well as between subject variation simultaneously hierarchical modeling approach is necessary 
# To estimate biologically meaningful parameters in a longitudinal design, mixed-effects model =
# within and between subjects variation

Ex005Grouped_ferrisi <- groupedData(weight ~ dpi|EH_ID, 
                                    ExpeDF_005[ExpeDF_005$Eimeria_species == "E.ferrisi",])

plot(Ex005Grouped_ferrisi)

# Model fitting for the normal-errors NLMM was accomplished with the ‘nlme’
# package in R (R core team 2012). Initial parameter estimates for ‘nlme’ were set to their
# true values from the data generating model.

# "All models were fit using R version 2.15 interfaced with JAGS through the
# R2JAGS package (R core team 2012, Plummer 2003)"

# https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.3.0.tar.gz/download

# sudo aptitude install r-cran-rjags

library(rjags)
# install.packages("R2jags") 
library(R2jags)
R2jags::
# Observed subject-specific core body temperature profiles for eight stallions
# experimentally challenged with the KY-84 strain of EAV 

# Estimated subject-specific febrile response functions

# Composite residual plot

# myModel <- function() {
#   
#   for (i in 1:n) {
#     
#     Response[i] ~ dnorm(f[i], pow(s.e,-2))
#     
#     f[i] <- step(
#       p[Subject[i]]-DPI[i])*
#       (B[Subject[i]]+I[Subject[i]]*exp(-pow(DPI[i]-p[Subject[i]],2)/(2*pow(l[Subject[i]],2)))) +
#       step(DPI[i]-p[Subject[i]])*(B[Subject[i]]+I[Subject[i]]*exp(-pow(DPI[i]-p[Subject[i]],2)/(2*pow(r[Subject[i]],2))))
#     res[i] <- Response[i] - f[i]
#   } 
#   
  
  #example
  
  N <- 1000
  x <- rnorm(N, 0, 5)
  
#   write.table(x,
#               file = 'example1.data',
#               row.names = FALSE,
#               col.names = FALSE)
# # In every model specification file, you have to start out by telling JAGS that you’re specifying a model.
 # Then we , which are meant to be constant across the loop. We tell JAGS that mu is distributed normally with mean 0 and standard deviation 100. This is meant to serve as a non-informative prior, since our data set was designed to have all measurements substantially below 100. Then we specify tau in a slightly round-about way. We say that tau is a deterministic function (hence the deterministic <- instead of the distributional ~) of sigma, after raising sigma to the -2 power. Then we say that sigma has a uniform prior over the interval [0,100].

  # http://www.johnmyleswhite.com/notebook/2010/08/20/using-jags-in-r-with-the-rjags-package/
  
  model {
    for (i in 1:N) { # set up the model for every single data point using a for loop
      x[i] ~ dnorm(mu, tau) # response x[i] is distributed normally with mean mu and precision (reciprocal of the variance) tau
    }
    # specify our priors for mu and tau, constant across the loop
    # mu is distributed normally with mean 0 and standard deviation 100. non-informative prior, since our data set was designed to have all measurements substantially below 100
    mu ~ dnorm(0, .0001) 
    # tau is a deterministic function (= non random response) of sigma, after raising sigma to the -2 power
    tau <- pow(sigma, -2) 
    #  Then we say that sigma has a uniform prior over the interval [0,100].
    sigma ~ dunif(0, 100)
  }
  
  
  
  myModel <- function() {
    
    for (i in 1:n) {
      Response[i] ~ dnorm(f[i], pow(s.e,-2))
      f[i] <- step(p[Subject[i]]-DPI[i])*(B[Subject[i]]+I[Subject[i]]*exp(-
                                                                            pow(DPI[i]-p[Subject[i]],2)/(2*pow(l[Subject[i]],2)))) +
        step(DPI[i]-p[Subject[i]])*(B[Subject[i]]+I[Subject[i]]*exp(-
                                                                      pow(DPI[i]-p[Subject[i]],2)/(2*pow(r[Subject[i]],2))))
      res[i] <- Response[i] - f[i]
    }
    
    # Hyperparameters for mean, standard deviation, baseline and intensity
    for (j in 1:8) {
      B[j] ~ dnorm(m.B,pow(s.B,-2))
      I[j] ~ dnorm(m.I,pow(s.I,-2))
      p[j] ~ dnorm(m.p,pow(s.p,-2))
      l[j] ~ dlnorm(m.log_l,pow(s.log_l,-2))
      r[j] ~ dlnorm(m.log_r,pow(s.log_r,-2))
    }
    
    # Hyperpriors
    m.B ~ dunif(0,5)
    m.I ~ dunif(-10,0)
    m.p ~ dunif(0,42)
    m.log_l ~ dunif(0,5)
    m.log_r ~ dunif(0,5)
    
    s.e ~ dunif(0,100)
    s.B ~ dunif(0,6)
    s.I ~ dunif(0,6)
    s.p ~ dunif(0,22)
    s.log_l ~ dunif(0,2)
    s.log_r ~ dunif(0,2)
    
    # FWHM
    HWHM.o <- 2.355*0.5*exp(m.log_l)
    HWHM.r <- 2.355*0.5*exp(m.log_r)
    FWHM <- HWHM.l + HWHM.r
  } 
  
  pow(3, -2)
  
  
  
  
  # modeling single-peaked, longitudinal EI data that incorporates recent developments in
  # nonlinear hierarchical models and Bayesian statistics.
  
  # Model a post-challenge incubation, then an onset phase, then recovery (return to basal state) = pattern of
  # single-peaked response
  
  # WEIGHT: nonlinear mixed model (NLMM) for a symmetric infection response variable. We
  # employ a standard NLMM assuming normally distributed errors and a Gaussian mean
  # response function. The parameters of the model correspond directly to biologically
  # meaningful properties of the infection response, including baseline, peak intensity, time
  # to peak and spread. 
   
  # OOCYST LOAD: For several reasons, a normal-errors model is not appropriate for viral load. We propose and
  # illustrate a Bayesian NLHM with the individual responses at each time point modeled as
  # a Poisson random variable with the means across time points related through a Tricube
  # mean response function. 
  
  # WEIGHT
  We modeled the weight responses using the
  NLHM in (2.1)-(2.7) with 𝑔 ∈ {1,2}
  
  
  # We model the individual response profiles with the modified Gaussian function (1.6).
  # Level 1 (model of intra-individual variability)
  # 𝑦𝑖𝑗~𝑁(𝜇 , 𝜎2) (2.1)
  # 𝜇𝑖𝑗 = 𝐵𝑖 + 𝐼𝑖𝑒𝑥𝑝   2𝑠𝑖
  # )
  # 2
  # ] (2.2)
  # Bi: baseline response level for individual 𝑖
  # Pi: time to peak response
  # Ii:peak response intensity𝑖
  # si: can be interpreted as a measure of
  # response duration.
  # 