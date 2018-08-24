
## Expe 001 on Rmd document already

## TO DO : put all on this document

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

# first, plot the mean of weight along time to remove the variations during prepatent and post-infection period

mysum003 <- summarySE(ExpeDF_003[!is.na(ExpeDF_003$weight),], measurevar="weight", 
                      groupvars=c("infection_isolate", "Mouse_strain","dpi"))

ggplot(mysum003, aes(dpi, weight)) +
  geom_point() +
  facet_grid(Mouse_strain~infection_isolate, scales = "free") +
  geom_line() +
  scale_x_continuous(breaks = 0:11)+
  geom_vline(xintercept = 3) +
  geom_vline(xintercept = 9) +
  theme_bw() +
  ggtitle("Expe003, weight along infection")

# Stats 

# Calculate weight loss
ExpeDF_003 <- calculateWeightLoss(ExpeDF_003, infectionDay = 3)

#Calculate OPG NOT DONE YET ;)
names(ExpeDF_003)[names(ExpeDF_003) %in% paste0("oocyst_sq", 1:4)] <- paste0("Neubauer", 1:4)
names(ExpeDF_003)[names(ExpeDF_003) %in% "dilution"] <- "dilution_ml"

ExpeDF_003 <- calculateOPG(ExpeDF_003)

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

  
  
  # ###############################
  # # Input data recent:
  # ExpeDF <- ExpeDF_005
  # ###############################
  # 
  # ##### my theme #####
  # theme_alice <- theme_bw() +
  #   theme(plot.title = element_text(size = 20, face = "bold"),
  #         plot.subtitle = element_text(size = 20),
  #         legend.text = element_text(size=20),
  #         legend.title = element_blank(),
  #         text = element_text(size = 20))
  # ##### end my theme #####
  # 
  # # Violin plots of the total sum of oocysts collected during 11 days: 
  # sum.oocysts <- do.call("rbind", by(ExpeDF, ExpeDF$EH_ID, function (x){
  #   x$sum.oo <- sum(x$OPG, na.rm=TRUE)
  #   x
  # }))
  # 
  # ggplot(unique(sum.oocysts[c("sum.oo", "Mouse_strain", "infection_isolate")]), 
  #        aes(Mouse_strain, sum.oo)) +
  #   ggtitle("Sum of oocysts shed during the experiment") + 
  #   geom_boxplot(color = "black")+
  #   facet_wrap(~infection_isolate) +
  #   geom_jitter(width=0.1, size=7, alpha = 0.5, pch = 21, aes(fill = Mouse_strain)) +
  #   # scale_fill_manual(values = c("blue", "purple", "red")) +
  #   labs(y = "Total number of oocyst shed", x = "Mouse strain") +
  #   scale_y_continuous(labels = scientific) +
  #   theme_alice +
  #   theme(legend.position="none")
  # 
  # summary(lm(sum.oo~Mouse_strain, 
  #            data = sum.oocysts[sum.oocysts$infection_isolate == "EI64",]))
  # 