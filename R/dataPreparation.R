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
ExpeDF_001 <- ExpeDF_001[ExpeDF_001$infection_isolate %in% "EI64",]
ExpeDF_001$infection_isolate <- droplevels(ExpeDF_001$infection_isolate)

# Mouse_strain: West should always be left 
ExpeDF_001$Mouse_strain <- factor(ExpeDF_001$Mouse_strain, 
                                  levels = c("WSB", "WP", "PWD"),
                                  labels = c("M.m.domesticus \n(WSB)", 
                                             "F1 hybrids \n(WP)", 
                                             "M.m.musculus \n(PWD)"))
# Keep only mice that survived 11 days
survivors <- names(table(ExpeDF_001$EH_ID))[table(ExpeDF_001$EH_ID) %in% 12]
ExpeDF_001 <- ExpeDF_001[ExpeDF_001$EH_ID %in% survivors,]

########################### Pass001: Nov 2017, passaging 4 isolates
# (Eflab, E88, E139, E64) in NMRI. 2 mice per cage. Only OPG recorded
PassDF_002 <- read.csv("../data/3_recordingTables/passaging_extra/Pass001_oocystsonly_Nov2017_Passaging_4Eimeria.csv")
PassDF_002$weightRelativeToInfection <- NA

########################### Exp002: March 2018 NMRI infected with 4 strains
ExpeDF_002 <- read.csv("../data/3_recordingTables/Exp002_March2018_NMRI_4strains_RECORDweightAndOocysts.csv")
ExpeDF_002$Mouse_strain <- "NMRI"
ExpeDF_002$EH_ID <- paste0("mouse_", ExpeDF_002$EH_ID)

# Calculate weight loss
ExpeDF_002 <- calculateWeightLoss(ExpeDF = ExpeDF_002)

# Calculate OPG
ExpeDF_002 <- calculateOPG(ExpeDF_002)

########################### Exp003 : May 2018 batch 1
oo <- read.csv("../data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDoocysts.csv")
we <- read.csv("../data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDweight.csv")
design <- read.csv("../data/2_designTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv")

ExpeDF_003 <- merge(oo, we, all = T)
ExpeDF_003 <- merge(ExpeDF_003, design, by = "EH_ID", all = T)
rm(design, oo, we)

# correct abherante value
ExpeDF_003[ExpeDF_003$EH_ID %in% "LM0137" & ExpeDF_003$weight %in% 17.6, "weight"] <- NA

# Mouse_strain: West should always be left 
ExpeDF_003$Mouse_strain <- factor(ExpeDF_003$Mouse_strain,
                                  levels = c("STRA", "BUSNA", "SCHUNT", "PWD"),
                                  labels = c("M.m.domesticus \n(STRA)", 
                                             "M.m.musculus \n(BUSNA)", 
                                             "M.m.domesticus \n(SCHUNT)",
                                             "M.m.musculus \n(PWD)"))

# Calculate weight loss
ExpeDF_003 <- calculateWeightLoss(ExpeDF = ExpeDF_003)

# Calculate OPG NOT DONE YET ;)
# ExpeDF_003 <- calculateOPG(ExpeDF_003)
