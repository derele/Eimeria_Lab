# Expe_003 & Expe_004
# April-May 2018
# first batch Parental strains (F0) BUSNA, STRA, SCHUNT, PWD
# infection with Eferrisi (E64 and E139) [2 batches]
source("functions4InfectionExperiments.R")

##############Cleaning##############
oo <- read.csv("../data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDoocysts.csv")
we <- read.csv("../data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDweight.csv")
design <- read.csv("../data/2_designTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv")
ExpeDF_003 <- merge(oo, we, all = T)
ExpeDF_003 <- merge(ExpeDF_003, design, by = "EH_ID", all = T)

oo <- read.csv("../data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDoocysts.csv")
we <- read.csv("../data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDweight.csv")
design <- read.csv("../data/2_designTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_DESIGN.csv")
ExpeDF_004 <- merge(oo, we, all = T)
ExpeDF_004 <- merge(ExpeDF_004, design, by = "EH_ID", all = T)
rm(design, oo, we)
ExpeDF_004$Mouse_strain <- ExpeDF_004$Strain
ExpeDF_004$Exp_ID <- "ExpeDF_004"

ExpeDF_003_4 <- merge(ExpeDF_003, ExpeDF_004, all = TRUE)

# Keep for dpi 0 to 11
ExpeDF_003_4 <- ExpeDF_003_4[ExpeDF_003_4$dpi %in% 0:11, ]# remove stabilisation period

# correct abherante value
ExpeDF_003_4[ExpeDF_003_4$EH_ID %in% "LM0137" & ExpeDF_003_4$weight %in% 17.6, "weight"] <- NA

# Add mice subspecies
ExpeDF_003_4$Mouse_subspecies <- "M.m.domesticus"
ExpeDF_003_4$Mouse_subspecies[ExpeDF_003_4$Mouse_strain %in% c("BUSNA", "PWD")] <- "M.m.musculus"

# Mouse_strain: West should always be left 
ExpeDF_003_4$Mouse_strain <- factor(ExpeDF_003_4$Mouse_strain,
                                    levels = c("STRA", "BUSNA", "SCHUNT", "PWD"),
                                    labels = c("M.m.domesticus \n(STRA)", 
                                               "M.m.musculus \n(BUSNA)", 
                                               "M.m.domesticus \n(SCHUNT)",
                                               "M.m.musculus \n(PWD)"))
# Enter Eimeria species
ExpeDF_003_4$Eimeria_species[ExpeDF_003_4$infection_isolate %in% c("E88", "Eflab", "EfLab")] = "E.falciformis"
ExpeDF_003_4$Eimeria_species[ExpeDF_003_4$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"

# Eimeria 64 "west" should always be left
ExpeDF_003_4$infection_isolate <- factor(ExpeDF_003_4$infection_isolate,
                                         levels = c("E64", "E139"))

# calculate OPG 
names(ExpeDF_003_4) <- gsub("oocyst_sq", "Neubauer", names(ExpeDF_003_4))
names(ExpeDF_003_4) <- gsub("dilution", "dilution_ml", names(ExpeDF_003_4))
ExpeDF_003_4$fecweight <- as.numeric(as.character(ExpeDF_003_4$fecweight))
ExpeDF_003_4 <- calculateOPG(ExpeDF = ExpeDF_003_4)

### until not counted, remove Expe4
ExpeDF_003_4 <- ExpeDF_003_4[!is.na(ExpeDF_003_4$OPG),]

##############End cleaning##############

############# Data preparation for analysis #############

# calculate relative weight loss compared to dpi0 (previous:3)
ExpeDF_003_4 <- calculateWeightLoss(ExpeDF_003_4, startingDay = 0)

# Add OPG status
ExpeDF_003_4$OPGstatus[is.na(ExpeDF_003_4$OPG)] = "na"
ExpeDF_003_4$OPGstatus[!is.na(ExpeDF_003_4$OPG) & ExpeDF_003_4$OPG > 0] = "positive"
ExpeDF_003_4$OPGstatus[!is.na(ExpeDF_003_4$OPG) & ExpeDF_003_4$OPG == 0] = "negative"
ExpeDF_003_4$OPGstatus = as.factor(ExpeDF_003_4$OPGstatus)

# calculate summary statistics on weight loss
summaryWeight_003_4 <- summarySE(ExpeDF_003_4, measurevar = "relativeWeight",
                                 groupvars=c("Mouse_strain", "Mouse_subspecies", "Eimeria_species",
                                             "infection_isolate", "dpi"), na.rm = T)
summaryWeight_003_4$ci[is.na(summaryWeight_003_4$ci)] <- 0

# calculate summary statistics on oocysts shedding
summaryOocysts_003_4 <-summarySE(ExpeDF_003_4, measurevar="OPG", 
                                 groupvars=c("Mouse_strain", "Mouse_subspecies", "Eimeria_species",
                                             "infection_isolate", "dpi"), na.rm = T)
summaryOocysts_003_4$ci[is.na(summaryOocysts_003_4$ci)] <- 0

# find maximum weight loss during symptomatic period (cf shedding plot)
summary_003_4 <- ExpeDF_003_4[ExpeDF_003_4$dpi %in% 4:8,] %>% 
  group_by(EH_ID) %>%
  filter(relativeWeight == min(relativeWeight)) 
# keep day of biggest shedding if equal weight
summary_003_4 <- summary_003_4[order(-summary_003_4[,"OPG"]),] 
summary_003_4 <- summary_003_4[!duplicated(summary_003_4$EH_ID),]

# calculate sum oocysts shed along experiment per animal
summary_003_4 <- merge(summary_003_4, ExpeDF_003_4 %>% 
        group_by(EH_ID) %>%
        dplyr::summarize(sum(Noocysts)))
      
############# End data preparation for analysis #############

ggplot(summary_003_4, aes(x=`sum(Noocysts)`, y = relativeWeight)) +
  geom_point(aes(col = Mouse_strain)) +
  facet_grid(.~infection_isolate)+
  theme_bw()

## 1. Weight loss
ggplot(ExpeDF_003_4, aes(x = dpi, y = weight, fill = OPGstatus)) +
  geom_line(aes(group = EH_ID, col = Mouse_strain), size = 2, alpha = 0.5) +
  geom_point(size=4, pch = 21, color = "black")+
  scale_fill_manual(values = c("lightgrey", "black", "yellow")) +
  mytheme +
  facet_grid(. ~ infection_isolate, scales = "free") +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "Weight (g)") +
  theme(strip.text.y = element_text(size = 15)) +
  # theme(legend.position="none") +
  ggtitle("Weight along experiment per indiidual", subtitle = "black : oocysts detected in feces")

ggplot(ExpeDF_003_4, aes(x = dpi, y = relativeWeight, fill = OPGstatus)) +
  geom_line(aes(group = EH_ID, col = Mouse_strain), size = 2, alpha = 0.5) +
  geom_point(size=4, pch = 21, color = "black")+
  scale_fill_manual(values = c("lightgrey", "black", "yellow")) +
  mytheme +
  facet_grid(. ~ infection_isolate, scales = "free_y", space = "free") +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "Weight relative to infection (%)") +
  theme(strip.text.y = element_text(size = 15)) +
  theme(legend.position="none") +
  ggtitle("Relative weight along experiment compared to dpi0", 
          subtitle = "black : oocysts detected in feces")


ggplot(summaryWeight_003_4, aes(x = dpi, y = relativeWeight))+
  geom_errorbar(aes(ymin = relativeWeight - ci, ymax = relativeWeight + ci, col = Mouse_strain),
                alpha = 0.5) +
  geom_line(aes(group = Mouse_strain, col = Mouse_strain), size = 2) +
  geom_point(aes(fill = Mouse_strain), size=4, pch = 21, color = "black") +
  mytheme +
  facet_grid(. ~ infection_isolate) +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "Weight relative to dpi3 (%)") +
  coord_cartesian(ylim = c(85,110)) +
  # theme(legend.position="none") +
  ggtitle("Mean relative weight along experiment compared to infection")

# For statistical analysis, we compare the maximum relative weight loss between the different groups. We limit our analysis to the period : dpi3 to dpi8 (symptomatic period for E64 strain).
plot1 <- ggplot(summary_003_4, aes(Mouse_strain,100-relativeWeight, color=infection_isolate)) +
  geom_boxplot(color = "black")+
  ggtitle("Maximum weigth loss") +
  facet_wrap(~infection_isolate) +
  geom_jitter(width=0.1, size=7, pch = 21,
              color = "black", aes(fill = Mouse_strain), alpha = 0.8) +
  labs(y= "Maximum weigth loss relative to weight at infection", x= "Mouse strain")+
  theme_bw() +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") 
plot1

kruskal.test(relativeWeight ~ Mouse_strain, data = summary_003_4)
pairwise.wilcox.test(summary_003_4$relativeWeight, summary_003_4$Mouse_strain,
                     p.adjust.method = "BH")

## 2. Parasite shedding
ggplot(ExpeDF_003_4, aes(x = dpi, y = OPG, fill = infection_isolate))+
  # geom_rect(aes(xmin = 3, xmax = 8, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = .5) +
  geom_line(aes(group = EH_ID, col = Mouse_strain), alpha = 0.5, size = 2) +
  geom_point(size=3, aes(color = Mouse_strain))+
  # scale_color_manual(values = c("blue", "purple", "red")) +
  facet_grid(Mouse_strain ~ infection_isolate) +
  mytheme +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  ylab("Oocysts shedding along infection (OPG)") +
  ggtitle("Oocyst shedding along experiment") +
  
  theme(legend.position="none") 

ggplot(summaryOocysts_003_4, aes(x = dpi, y = OPG))+
  geom_errorbar(aes(ymin = OPG - ci, ymax = OPG + ci, col = Mouse_strain), alpha = 0.7) +
  geom_line(aes(group = Mouse_strain, col = Mouse_strain), size = 2) +
  geom_point(aes(fill = Mouse_strain), size=4, pch = 21, color = "black") +
  mytheme +
  facet_grid(Mouse_strain ~ infection_isolate) +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "Oocysts shedding along infection (OPG)") +
  scale_y_log10() +
  theme(legend.position="none") +
  ggtitle("Mean oocysts shedding along experiment")

plot2 <-ggplot(summary_003_4, aes(Mouse_strain, `sum(Noocysts)`, color=infection_isolate)) +
  geom_boxplot(color = "black")+
  ggtitle("Sum of oocysts shed along the experiment") +
  facet_wrap(~infection_isolate) +
  geom_jitter(width=0.1, size=7, pch = 21,
              color = "black", aes(fill = Mouse_strain), alpha = 0.8) +
  labs(y= "Sum of oocysts shed along experiment", x= "Mouse strain")+
  theme_bw() +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") 
plot2
kruskal.test(`sum(Noocysts)` ~ Mouse_strain, data = summary_003_4)
pairwise.wilcox.test(summary_003_4$`sum(Noocysts)`, summary_003_4$Mouse_strain,
                     p.adjust.method = "BH")

grid.arrange(plot1, plot2, ncol=2)

ggplot(summary_003_4, aes(x = 100-relativeWeight, y = `sum(Noocysts)`,
                          group = Mouse_strain, col = Mouse_strain))+
  geom_point(size = 5) +
  facet_grid(.~infection_isolate) +
  labs(x= "Maximum weigth loss relative to weight at infection", 
       y= "Sum of oocysts shed along experiment")+
  theme_bw()

       