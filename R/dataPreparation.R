## Data preparation
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

# functions used
calculateWeightLoss <- function(ExpeDF){
  A = ExpeDF[ExpeDF$dpi == 0, c("weight", "EH_ID")]
  names(A)[1] = "weightAtInfection"
  ExpeDF <- merge(ExpeDF, A)
  rm(A)
  ExpeDF$weightloss = ExpeDF$weightAtInfection - ExpeDF$weight
  ExpeDF$weightRelativeToInfection <- ExpeDF$weight /
    ExpeDF$weightAtInfection * 100
  return(ExpeDF)
}


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

########################### Exp002: March 2018 NMRI infected with
# ExpeDF <- read.csv("../data/3_recordingTables/March2018_NMRI_4strains_RECORDweightAndOocysts.csv")
# # general
# ExpeDF$Mouse_strain <- "NMRI"

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

# Calculate OPG
# ExpeDFMay2018batch1$mean_Neubauer <- (ExpeDFMay2018batch1$Neubauer1 + ExpeDFMay2018batch1$Neubauer2 + ExpeDFMay2018batch1$Neubauer3 + ExpeDFMay2018batch1$Neubauer4) / 4
# ExpeDFMay2018batch1$OPG <- ExpeDFMay2018batch1$mean_Neubauer * 10000 / ExpeDFMay2018batch1$dilution_ml / ExpeDFMay2018batch1$fecweight
# ExpeDFMay2018batch1$oocystsTotal <- ExpeDFMay2018batch1$mean_Neubauer * 10000 / ExpeDFMay2018batch1$dilution_ml
