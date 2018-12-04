# Expe_005
# July 2018
# FULL experiment (parents, intra specific and inter species hybrids)
# BUSNA, STRA, SCHUNT, PWD
# infection with Eferrisi and Efalciformis (E64 and E88)
source("functions4InfectionExperiments.R")

##############Cleaning##############
oo <- read.csv("../data/3_recordingTables/Exp005_full_RECORDoocysts.csv", na.strings = c("NA", " ", "n.a."))
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
ExpeDF_005$Parent1 <- gsub(" ", "", ExpeDF_005$Parent1)
ExpeDF_005$Parent2 <- gsub(" ", "", ExpeDF_005$Parent2)
# Correct error non numeric
ExpeDF_005$weight <- as.numeric(as.character(ExpeDF_005$weight))
ExpeDF_005$weight_dpi0 <- as.numeric(as.character(ExpeDF_005$weight_dpi0))
# Correct weight loss
ExpeDF_005$weightloss <- ExpeDF_005$weight_dpi0 - ExpeDF_005$weight

# Add mice subspecies info
ExpeDF_005$Mouse_subspecies <- factor(as.factor(ExpeDF_005$Strain),
                                      levels = c("BUSNA_BUSNA","BUSNA_PWD","BUSNA_STRA","PWD_BUSNA",
                                                 "PWD_PWD","PWD_SCHUNT","SCHUNT_PWD","SCHUNT_SCHUNT",
                                                 "SCHUNT_STRA","STRA_BUSNA","STRA_SCHUNT","STRA_STRA"),
                                      labels = c("M.m.musculus P", "M.m.musculus F1",
                                                 "Hybrid", "M.m.musculus F1",
                                                 "M.m.musculus P","Hybrid",
                                                 "Hybrid", "M.m.domesticus P",
                                                 "M.m.domesticus F1","Hybrid",
                                                 "M.m.domesticus F1","M.m.domesticus P"))
# # Mouse_strain: West should always be left
ExpeDF_005$Mouse_strain <- factor(as.factor(ExpeDF_005$Strain),
                                  levels = c("BUSNA_BUSNA","BUSNA_PWD",
                                             "BUSNA_STRA","PWD_BUSNA",
                                             "PWD_PWD","PWD_SCHUNT",
                                             "SCHUNT_PWD","SCHUNT_SCHUNT",
                                             "SCHUNT_STRA","STRA_BUSNA",
                                             "STRA_SCHUNT","STRA_STRA"),
                                  labels = c("M.m.musculus P \n(BUSNA_BUSNA)", "M.m.musculus F1 \n(BUSNA_PWD)",
                                             "Hybrid \n(BUSNA_STRA)", "M.m.musculus F1 \n(PWD_BUSNA)",
                                             "M.m.musculus P \n(PWD_PWD)","Hybrid \n(PWD_SCHUNT)",
                                             "Hybrid \n(SCHUNT_PWD)", "M.m.domesticus P \n(SCHUNT_SCHUNT)",
                                             "M.m.domesticus F1 \n(SCHUNT_STRA)","Hybrid \n(STRA_BUSNA)",
                                             "M.m.domesticus F1 \n(STRA_SCHUNT)","M.m.domesticus P \n(STRA_STRA)"))

# Correct hybrid status
ExpeDF_005$HybridStatus <- factor(ExpeDF_005$HybridStatus,
                                  levels = c("intersubsp.hybrids", "outbredhybrids", "parentalstrains"),
                                  labels = c("hybrids", "outbred", "parents"))

# Add infection isolate
ExpeDF_005$infection_isolate <- ExpeDF_005$Eimeria

# Enter Eimeria species
ExpeDF_005$Eimeria_species[ExpeDF_005$infection_isolate %in% c("E88", "Eflab", "EfLab")] <- "E.falciformis"
ExpeDF_005$Eimeria_species[ExpeDF_005$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"

# # Age at infection
# ExpeDF_005$ageAtInfection[ExpeDF_005$Batch == 1] <- round(ExpeDF_005$ageAtdpi0expe1a)
# ExpeDF_005$ageAtInfection[ExpeDF_005$Batch == 2] <- round(ExpeDF_005$ageAtdpi0expe1a +2)

# calculate OPG 
# rename 
names(ExpeDF_005) <- gsub("oocyst_sq", "Neubauer", names(ExpeDF_005))
names(ExpeDF_005) <- gsub("dilution", "dilution_ml", names(ExpeDF_005))
ExpeDF_005$dilution_ml <- as.numeric(as.character(ExpeDF_005$dilution_ml))
ExpeDF_005$fecweight <- as.numeric(as.character(ExpeDF_005$fecweight))
# Replace NA by 0
ExpeDF_005$fecweight[is.na(ExpeDF_005$fecweight)] <- 0

ExpeDF_005 <- calculateOPG(ExpeDF = ExpeDF_005)

# Add OPG status
ExpeDF_005$OPGstatus[is.na(ExpeDF_005$OPG)] = "na"
ExpeDF_005$OPGstatus[!is.na(ExpeDF_005$OPG) & ExpeDF_005$OPG > 0] = "positive"
ExpeDF_005$OPGstatus[!is.na(ExpeDF_005$OPG) & ExpeDF_005$OPG == 0] = "negative"
ExpeDF_005$OPGstatus = as.factor(ExpeDF_005$OPGstatus)

# ERROR one fecal sample has Inf as value for OPG (cause 0g as fecal weight)
ExpeDF_005[ExpeDF_005$OPG %in% Inf, "OPG"] <- NA
# calculate relative weight loss compared to dpi0
ExpeDF_005 <- calculateWeightLoss(ExpeDF_005, startingDay = 0)
ExpeDF_005$relativeWeight <- as.numeric(as.character(ExpeDF_005$relativeWeight))

# Export
write.csv2(file = "../data/3_recordingTables/Exp005_complete.csv",
           x = ExpeDF_005, row.names = F, fileEncoding = "UTF-8")

## REMOVE SECOND BATCH. CROSS CONTAMINATION!!!!
ExpeDF_005 <- ExpeDF_005[ExpeDF_005$Batch==1,]

## Add info for lm
for (strain in c("BUSNA", "STRA", "SCHUNT", "PWD")){
  ExpeDF_005[[strain]] <- lengths(
    regmatches(ExpeDF_005$Strain, gregexpr(strain, ExpeDF_005$Strain)))
}

# Split in 2 infection batches (Efalciformis and Eferrisi are studied separetely)
ExpeDF_005_E88 <- ExpeDF_005[ExpeDF_005$infection_isolate %in% "E88",]
ExpeDF_005_E64 <- ExpeDF_005[ExpeDF_005$infection_isolate %in% "E64",]

# calculate summary statistics on weight loss
summaryWeight_005_E64 <- summarySE(ExpeDF_005_E64, measurevar = "relativeWeight",
                                   groupvars=c("HybridStatus", "dpi", "Expe"), na.rm = T)
summaryWeight_005_E64$ci[is.na(summaryWeight_005_E64$ci)] <- 0

summaryWeight_005_E88 <- summarySE(ExpeDF_005_E88, measurevar = "relativeWeight",
                                   groupvars=c("HybridStatus", "dpi", "Expe"), na.rm = T)
summaryWeight_005_E88$ci[is.na(summaryWeight_005_E88$ci)] <- 0

# calculate summary statistics on oocysts shedding
summaryOocysts_005_E64 <-summarySE(ExpeDF_005_E64, measurevar="OPG",
                                   groupvars=c("HybridStatus", "dpi", "Expe"), na.rm = T)
summaryOocysts_005_E64$ci[is.na(summaryOocysts_005_E64$ci)] <- 0

summaryOocysts_005_E88 <-summarySE(ExpeDF_005_E88, measurevar="OPG",
                                   groupvars=c("HybridStatus", "dpi", "Expe"), na.rm = T)
summaryOocysts_005_E88$ci[is.na(summaryOocysts_005_E88$ci)] <- 0

############## End cleaning##############
df <- ExpeDF_005 %>%
  group_by(dpi, infection_isolate, HybridStatus) %>%
  dplyr::summarise(mean.oo = mean(Noocysts, na.rm = T),
                   sd.oo = sd(Noocysts, na.rm = TRUE),
                   n.oo = n()) %>%
  dplyr::mutate(se.oo = sd.oo / sqrt(n.oo),
                lower.ci.oo = mean.oo - qt(1 - (0.05 / 2), n.oo - 1) * se.oo,
                upper.ci.oo = mean.oo + qt(1 - (0.05 / 2), n.oo - 1) * se.oo)

ggplot(df, aes(x = dpi, y = mean.oo))+
  geom_errorbar(aes(ymin = lower.ci.oo, ymax = upper.ci.oo, 
                    col = HybridStatus),
                alpha = 0.5) +
  geom_line(aes(group = HybridStatus, col = HybridStatus), size = 2) +
  geom_point(aes(fill = HybridStatus), size=4, pch = 21, color = "black") +
  mytheme +
  facet_grid(. ~ infection_isolate) +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "Nbr oocysts shed")
  
df <- ExpeDF_005 %>%
  group_by(dpi, infection_isolate, HybridStatus) %>%
  dplyr::summarise(mean.weight = mean(relativeWeight, na.rm = T),
                   sd.weight = sd(relativeWeight, na.rm = TRUE),
                   n.weight = n()) %>%
  dplyr::mutate(se.weight = sd.weight / sqrt(n.weight),
                lower.ci.weight = mean.weight - qt(1 - (0.05 / 2), 
                                                   n.weight - 1) * se.weight,
                upper.ci.weight = mean.weight + qt(1 - (0.05 / 2), 
                                                   n.weight - 1) * se.weight)
ggplot(df, aes(x = dpi, y = mean.weight))+
  geom_errorbar(aes(ymin = lower.ci.weight, ymax = upper.ci.weight, 
                    col = HybridStatus),
                alpha = 0.5) +
  geom_line(aes(group = HybridStatus, col = HybridStatus), size = 2) +
  geom_point(aes(fill = HybridStatus), size=4, pch = 21, color = "black") +
  mytheme +
  facet_grid(. ~ infection_isolate) +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "relative weight")

# find maximum weight loss during symptomatic period (cf shedding plot)
summary_005 <- ExpeDF_005 %>% 
  group_by(EH_ID) %>%
  filter(relativeWeight == min(relativeWeight)) 
summary_005 <- summary_005[!duplicated(summary_005$EH_ID),]

# calculate sum oocysts shed along experiment per animal
summary_005 <- merge(summary_005, 
                     ExpeDF_005 %>% 
                         group_by(EH_ID) %>%
                         dplyr::summarize(sum.oo = sum(Noocysts, na.rm = T)))

############# End data preparation for analysis #############

ggplot(summary_005, aes(x=sum.oo, y = relativeWeight)) +
  geom_point(aes(col = HybridStatus)) +
  facet_grid(.~infection_isolate)+
  theme_bw()

# For statistical analysis, we compare the maximum relative weight loss between the different groups. We limit our analysis to the period : dpi3 to dpi8 (symptomatic period for E64 strain).
plot1 <- ggplot(summary_005, aes(HybridStatus,100-relativeWeight)) +
  geom_boxplot(color = "black")+
  ggtitle("Maximum weigth loss during the infection") +
  facet_wrap(~infection_isolate) +
  geom_jitter(width=0.1, size=7, pch = 21,
              color = "black", aes(fill = Mouse_strain), alpha = 0.8) +
  labs(y= "Maximum weigth loss relative to weight at infection", x= "Mouse strain")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot1

# And now, strain by strain
plot_list <- list()
strains <- c("STRA", "BUSNA", "SCHUNT", "PWD")
for (i in 1:4){
  df <- summary_005[grep(strains[i], summary_005$Strain),]
  g <- ggplot(df, aes(HybridStatus,100-relativeWeight)) +
    geom_boxplot(color = "black")+
    ggtitle("Maximum weigth loss during the infection") +
    facet_wrap(~infection_isolate) +
    geom_jitter(width=0.1, size=7, pch = 21,
                color = "black", aes(fill = Mouse_strain), alpha = 0.8) +
    labs(y= "Maximum weigth loss relative to weight at infection", x= "Mouse strain")+
    theme_bw() +
    ggtitle(strains[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") 
  plot_list[[i]] <- g
}
do.call(grid.arrange, plot_list)

# kruskal.test(relativeWeight ~ Mouse_strain, data = summary_003_4)
# pairwise.wilcox.test(summary_003_4$relativeWeight, summary_003_4$Mouse_strain,
#                      p.adjust.method = "BH")

plot2 <-ggplot(summary_005, aes(HybridStatus, sum.oo)) +
  geom_boxplot(color = "black")+
  ggtitle("Sum of oocysts shed along the experiment") +
  facet_wrap(~infection_isolate) +
  geom_jitter(width=0.1, size=7, pch = 21,
              color = "black", aes(fill = Mouse_strain), alpha = 0.8) +
  labs(y= "Sum of oocysts shed along experiment", x= "Mouse strain")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") 
plot2

# And now, strain by strain
plot_list <- list()
strains <- c("STRA", "BUSNA", "SCHUNT", "PWD")
for (i in 1:4){
  df <- summary_005[grep(strains[i], summary_005$Strain),]
  g <- ggplot(df, aes(HybridStatus, sum.oo, color=infection_isolate)) +
    geom_boxplot(color = "black")+
    ggtitle("Sum of oocysts shed along the experiment") +
    facet_wrap(~infection_isolate) +
    geom_jitter(width=0.1, size=7, pch = 21,
                color = "black", aes(fill = Mouse_strain), alpha = 0.8) +
    labs(y= "Sum of oocysts shed along experiment", x= "Mouse strain")+
    theme_bw() +
    ggtitle(strains[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") 
  plot_list[[i]] <- g
}
do.call(grid.arrange, plot_list)
# kruskal.test(`sum(Noocysts)` ~ Mouse_strain, data = summary_003_4)
# pairwise.wilcox.test(summary_003_4$`sum(Noocysts)`, summary_003_4$Mouse_strain,
#                      p.adjust.method = "BH")

# grid.arrange(plot1, plot2, ncol=2)

# * Add variable for each 4 parents and test the linear relationships for each of these variables set to 0 (copy of DNA), 1 (copy of DNA) (2 we can remove as we want outbred vs hybrids) + another variable HybridStatus : hybrid or outbred. + mixed effect (1|EH_ID, 1|Expe)
# 
# * Depend on the angle, but could be really interesting to quantify this for each mouse strain (outbreeding effet + hybrid effect) and show that it is highly strain specific. The focus on the article could be on that.

# Question 1. is there a outbred effect and a hybrid effect in general?
# summary_005$sum.oo
# summary_005$relativeWeight
boxplot(relativeWeight ~ HybridStatus*infection_isolate,
        col=c("white","lightgray"),summary_005)

mod.HS <- lmer(data = summary_005, 
            formula = relativeWeight ~ HybridStatus + 1|infection_isolate,REML=FALSE)
mod.HS.null <- lmer(data = summary_005, 
               formula = relativeWeight ~ 1|infection_isolate,REML=FALSE)

anova(mod.HS, mod.HS.null)

boxplot(sum.oo ~ HybridStatus*infection_isolate,
        col=c("white","lightgray"),summary_005)

mod.HS <- lmer(data = summary_005, 
               formula = sum.oo ~ HybridStatus + 1|infection_isolate,REML=FALSE)
mod.HS.null <- lmer(data = summary_005, 
                    formula = sum.oo ~ 1|infection_isolate,REML=FALSE)

anova(mod.HS, mod.HS.null)

## 4 different blocs, overlaps !!
# strain effect
plot_list <- list()
strains <- c("STRA", "BUSNA", "PWD", "SCHUNT")
for (i in 1:4){
  df <- summary_005[grep(strains[i], summary_005$Strain),]
  g <- ggplot(df, aes(x = HybridStatus, col = infection_isolate, y = relativeWeight)) +
    geom_boxplot() +
    geom_point(position=position_jitterdodge()) +
    ggtitle(strains[i]) +
    theme_bw() +
    theme(legend.position="top")
  plot_list[[i]] <- g
}
do.call(grid.arrange, plot_list)

# strain effect
plot_list <- list()
strains <- c("STRA", "BUSNA", "PWD", "SCHUNT")
for (i in 1:4){
  df <- summary_005[grep(strains[i], summary_005$Strain),]
  g <- ggplot(df, aes(x = HybridStatus, col = infection_isolate, y = sum.oo)) +
    geom_boxplot() +
    geom_point(position=position_jitterdodge()) +
    ggtitle(strains[i]) +
    theme_bw() +
    theme(legend.position="top")
  plot_list[[i]] <- g
}
do.call(grid.arrange, plot_list)

