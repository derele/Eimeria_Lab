# load all datasets Alice

# Expe_001, March 2017, Francisca's experiment. infection with E64 and Eflab
# ExpeDF need at least the following:
# "OPG", "EH_ID", "dpi", "weight", "EH_ID", "infection_isolate", "Mouse_strain"
ExpeDF_001 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp001_May2017_crossing_infection.csv")

ExpeDF_001$Exp_ID <- "Exp_001" 

ExpeDF_001$Infection_date <- as.character(ExpeDF_001$Infection_date)
ExpeDF_001$Infection_date[ExpeDF_001$Infection_date == "11_May_17"] <- "2017-05-11 CET"
ExpeDF_001$Infection_date[ExpeDF_001$Infection_date == "8_May_17"] <- "2017-05-08 CET"

ExpeDF_001$Born <- as.POSIXct(ExpeDF_001$born, format = "%d.%m.%Y")

ExpeDF_001$ageAtInfection <- difftime(ExpeDF_001$Infection_date, ExpeDF_001$Born,
                                      units = "weeks")

# correct names
names(ExpeDF_001)[names(ExpeDF_001) %in% "strain"] <- "Mouse_strain"
names(ExpeDF_001)[names(ExpeDF_001) %in% "Inf_strain"]  <- "infection_isolate"
names(ExpeDF_001)[names(ExpeDF_001) %in% "EH_id"]  <- "EH_ID"
names(ExpeDF_001)[names(ExpeDF_001) %in% "oocysts.per.g"]  <- "OPG"
names(ExpeDF_001)[names(ExpeDF_001) %in% "rel.weight"]  <- "relativeWeight"
names(ExpeDF_001)[names(ExpeDF_001) %in% "fec.weight"]  <- "fecweight"

# Add mice subspecies info
ExpeDF_001$Mouse_subspecies[ExpeDF_001$Mouse_strain == "PWD"] <- "M.m.musculus"
ExpeDF_001$Mouse_subspecies[ExpeDF_001$Mouse_strain == "WSB"] <- "M.m.domesticus"
ExpeDF_001$Mouse_subspecies <- factor(ExpeDF_001$Mouse_subspecies,
                                      levels = c("M.m.domesticus", "F1 hybrids", "M.m.musculus"))
# Enter Eimeria species
ExpeDF_001$Eimeria_species[ExpeDF_001$infection_isolate %in% c("E88", "Eflab", "EfLab")] = "E.falciformis"
ExpeDF_001$Eimeria_species[ExpeDF_001$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
# homogenize infection isolate
ExpeDF_001$infection_isolate <- factor(ExpeDF_001$infection_isolate,
                                       levels = c("Eflab", "EI64"),
                                       labels = c("EfLab", "E64"))
ExpeDF_001$Mouse_strain <- as.character(ExpeDF_001$Mouse_strain)
# Confident that there was a typo for LM0106 dpi4 --> set to 0
ExpeDF_001$OPG[ExpeDF_001$EH_ID %in% "LM0106" & ExpeDF_001$dpi %in% 4] <- 0

# weight loss NB compared to dpi1
ExpeDF_001 <- calculateWeightLoss(ExpeDF_001, startingDay = 1)
# # Set mice that died before the end of experiment at 20% weight loss (maximum) on last day WRONG
# ExpeDF_001[is.na(ExpeDF_001$relativeWeight) & ExpeDF_001$dpi != 0, "relativeWeight"] <- 80
# Add hybrid status
ExpeDF_001$HybridStatus <- "inbred"
ExpeDF_001[ExpeDF_001$Mouse_strain == "WP", "HybridStatus"] <- "hybrids"
# Call PWD --> PWD1 (other origin)
ExpeDF_001$Mouse_strain[ExpeDF_001$Mouse_strain == "PWD"] <- "PWD1"

ExpeDF_001$Sex <- NA
ExpeDF_001$Sex[ExpeDF_001$sex == "f"] <- "F"
ExpeDF_001$Sex[ExpeDF_001$sex == "m"] <- "M"

## Load info from passaging : EXPE_002 
ExpeDF_002 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp002_March2018_NMRI_4strains_RECORDweightAndOocysts.csv")
ExpeDF_002$Mouse_strain <- "NMRI"
ExpeDF_002$EH_ID <- paste0("mouse_", ExpeDF_002$EH_ID)
ExpeDF_002$Mouse_subspecies = ExpeDF_002$Mouse_strain
#rename colum "sex" to column "Sex"
names(ExpeDF_002)[names(ExpeDF_002) == "sex"] <- "Sex"
ExpeDF_002$Sex <- as.character(ExpeDF_002$Sex)
#rename "female"/"male" to "F" or "M"
ExpeDF_002$Sex[ExpeDF_002$Sex == "female"] <- "F" 
ExpeDF_002$Sex[ExpeDF_002$Sex == "male"] <- "M"
# create column "Hybridstatus"
ExpeDF_002$HybridStatus <- "inbred"
# Enter Eimeria species
ExpeDF_002$Eimeria_species[ExpeDF_002$infection_isolate %in% c("E88", "Eflab", "EfLab")] <- "E.falciformis"
ExpeDF_002$Eimeria_species[ExpeDF_002$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
# weight loss
ExpeDF_002 <- calculateWeightLoss(ExpeDF_002, startingDay = 0)
# oocysts
ExpeDF_002 <- calculateOPG(ExpeDF_002)

ExpeDF_002$Infection_date <- "19-03-18"

### EXPE_003 & 004
# read all information for EXPE_003 and merge to one BIG file
oo <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDoocysts.csv")
we <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_RECORDweight.csv")

design <- read.csv("../../data/2_designTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv")

# design <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/Exp003_April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv")
ExpeDF_003 <- merge(oo, we, all = T)
ExpeDF_003 <- merge(ExpeDF_003, design, by = "EH_ID", all = T)
ExpeDF_003$Born <- as.POSIXct(ExpeDF_003$born, format = "%d/%m/%y")

ExpeDF_003$Infection_date <- "2018-05-03 CET"
ExpeDF_003$ageAtInfection <- difftime(ExpeDF_003$Infection_date, ExpeDF_003$Born,
                                        units = "weeks")

# read all information for EXPE_004 and merge to one BIG file
oo <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDoocysts.csv")
we <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_RECORDweight.csv")
design <- read.csv("../../data/2_designTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_DESIGN.csv")
ExpeDF_004 <- merge(oo, we, all = T)
ExpeDF_004 <- merge(ExpeDF_004, design, by = "EH_ID", all = T)
ExpeDF_004$Born <- as.POSIXct(ExpeDF_004$Born, format = "%y-%m-%d")

ExpeDF_004$Infection_date <- "2018-06-03 CET"
ExpeDF_004$ageAtInfection <- difftime(ExpeDF_004$Infection_date, ExpeDF_004$Born,
                                      units = "weeks")
# rm dead
ExpeDF_004 <- ExpeDF_004[!ExpeDF_004$EH_ID %in% c("LM0155", "LM0160"),]

# rename column with "strain" for Expe_004 and merge BIG files into one for both experiments: ExpeDF_003_4
rm(design, oo, we)
ExpeDF_004$Mouse_strain <- ExpeDF_004$Strain
ExpeDF_004$Exp_ID <- "Exp_004"
names(ExpeDF_003)[names(ExpeDF_003) == "sex"] <- "Sex"
ExpeDF_003_4 <- merge(ExpeDF_003, ExpeDF_004, all = TRUE)
# # Keep for dpi 0 to 11
# ExpeDF_003_4 <- ExpeDF_003_4[ExpeDF_003_4$dpi %in% 0:11, ]# remove stabilisation period
# correct abherante value
ExpeDF_003_4[ExpeDF_003_4$EH_ID %in% "LM0137" & ExpeDF_003_4$weight %in% 17.6, "weight"] <- NA
# Add mouse subspecies
ExpeDF_003_4$Mouse_subspecies <- "M.m.domesticus"
ExpeDF_003_4$Mouse_subspecies[ExpeDF_003_4$Mouse_strain %in% c("BUSNA", "PWD")] <- "M.m.musculus"
# create column "HybriStauts"
ExpeDF_003_4$HybridStatus <- "inbred"
#rename mouse strains
ExpeDF_003_4$Mouse_strain <- as.character(ExpeDF_003_4$Mouse_strain)
ExpeDF_003_4$Mouse_strain[ExpeDF_003_4$Mouse_strain == "STRA"] <- "STRA_STRA"
ExpeDF_003_4$Mouse_strain[ExpeDF_003_4$Mouse_strain == "SCHUNT"] <- "SCHUNT_SCHUNT"
ExpeDF_003_4$Mouse_strain[ExpeDF_003_4$Mouse_strain == "BUSNA"] <- "BUSNA_BUSNA"
ExpeDF_003_4$Mouse_strain[ExpeDF_003_4$Mouse_strain == "PWD"] <- "PWD_PWD"
# Enter Eimeria species
ExpeDF_003_4$Eimeria_species[ExpeDF_003_4$infection_isolate %in% c("E88", "Eflab", "EfLab")] <- "E.falciformis"
ExpeDF_003_4$Eimeria_species[ExpeDF_003_4$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
# Eimeria 64 "west" should always be left
ExpeDF_003_4$infection_isolate <- factor(ExpeDF_003_4$infection_isolate,
                                         levels = c("E64", "E139"))
#weightloss
ExpeDF_003_4 <- calculateWeightLoss(ExpeDF_003_4, startingDay = 0)
# calculate OPG 
# rename 
names(ExpeDF_003_4) <- gsub("oocyst_sq", "Neubauer", names(ExpeDF_003_4))
names(ExpeDF_003_4) <- gsub("dilution", "dilution_ml", names(ExpeDF_003_4))
ExpeDF_003_4$fecweight <- as.numeric(as.character(ExpeDF_003_4$fecweight))
ExpeDF_003_4 <- calculateOPG(ExpeDF_003_4)
## NB LM0193 died before the peak, remove
ExpeDF_003_4 <- 
  ExpeDF_003_4[!ExpeDF_003_4$EH_ID %in% "LM0193",]

# Keep for dpi 0 to 11
ExpeDF_003_4 <- ExpeDF_003_4[ExpeDF_003_4$dpi %in% 0:11, ]# remove stabilisation period

#### Expe_005
# load all tables of all experiments, 
oo <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp005_full_RECORDoocysts.csv", na.strings = c("NA", " ", "n.a."))
we <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp005_full_RECORDweight.csv", na.strings = c("NA", " "))
design <- rbind(read.csv("../../data/2_designTables/Inf1a_Exp005.DESIGN.csv", na.strings = c("NA", " ")),
                read.csv("../../data/2_designTables/Inf2a_Exp005.DESIGN.csv", na.strings = c("NA", " ")),
                read.csv("../../data/2_designTables/Inf1b_Exp005.DESIGN.csv", na.strings = c("NA", " ")),
                read.csv("../../data/2_designTables/Inf2b_Exp005.DESIGN.csv", na.strings = c("NA", " ")))
# Correct error
names(design)[names(design) == "EH_id"] <- "EH_ID"
# Correct space error
design$EH_ID <- gsub(" ", "", design$EH_ID)
design$HybridStatus <- gsub(" ", "", design$HybridStatus)
ExpeDF_005 <- merge(oo, we, by = c("labels", "Expe"))
ExpeDF_005 <- merge(ExpeDF_005, design, by = "EH_ID", all.x = T)

ExpeDF_005$Born <- as.POSIXct(ExpeDF_005$Born, format = "%Y-%m-%d")

## Correct error space
ExpeDF_005$Strain <- gsub(" ", "", ExpeDF_005$Strain)

ExpeDF_005$Infection_date <- "2018-07-09 CET"

ExpeDF_005$ageAtInfection <- difftime(ExpeDF_005$Infection_date, ExpeDF_005$Born,
                                      units = "weeks")
# rename "Expe" -> "Exp_ID"
names(ExpeDF_005)[names(ExpeDF_005) == "Expe"] <- "Exp_ID"
names(ExpeDF_005)[names(ExpeDF_005) == "Strain"] <- "Mouse_strain"
#rename ExpeID
ExpeDF_005$Exp_ID <- as.character(ExpeDF_005$Exp_ID)
ExpeDF_005$Exp_ID[ExpeDF_005$Exp_ID == "Expe005_1a"] <- "Exp_005_1a"
ExpeDF_005$Exp_ID[ExpeDF_005$Exp_ID == "Expe005_1b"] <- "Exp_005_1b"
ExpeDF_005$Exp_ID[ExpeDF_005$Exp_ID == "Expe005_2a"] <- "Exp_005_2a"
ExpeDF_005$Exp_ID[ExpeDF_005$Exp_ID == "Expe005_2b"] <- "Exp_005_2b"
# Correct error non numeric
ExpeDF_005$weight <- as.numeric(as.character(ExpeDF_005$weight))
ExpeDF_005$weight_dpi0 <- as.numeric(as.character(ExpeDF_005$weight_dpi0))
# Correct weight loss
ExpeDF_005$weightloss <- ExpeDF_005$weight_dpi0 - ExpeDF_005$weight
# Add mouse subspecies info
ExpeDF_005$Mouse_subspecies <- factor(as.factor(ExpeDF_005$Mouse_strain),
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
ExpeDF_005$Full_info_mouse <- factor(as.factor(ExpeDF_005$Mouse_strain),
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
                                  labels = c("hybrids", "outbred", "inbred"))
# Add infection isolate
ExpeDF_005$infection_isolate <- ExpeDF_005$Eimeria
# Enter Eimeria species
ExpeDF_005$Eimeria_species[ExpeDF_005$infection_isolate %in% c("E88", "Eflab", "EfLab")] <- "E.falciformis"
ExpeDF_005$Eimeria_species[ExpeDF_005$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
# calculate OPG 
# rename 
names(ExpeDF_005) <- gsub("oocyst_sq", "Neubauer", names(ExpeDF_005))
names(ExpeDF_005) <- gsub("dilution", "dilution_ml", names(ExpeDF_005))
ExpeDF_005$dilution_ml <- as.numeric(as.character(ExpeDF_005$dilution_ml))
ExpeDF_005$fecweight <- as.numeric(as.character(ExpeDF_005$fecweight))
# Replace NA by 0
ExpeDF_005$fecweight[is.na(ExpeDF_005$fecweight)] <- 0
ExpeDF_005 <- calculateOPG(ExpeDF = ExpeDF_005)
# ERROR one fecal sample has Inf as value for OPG (cause 0g as fecal weight)
ExpeDF_005[ExpeDF_005$OPG %in% Inf, "OPG"] <- NA
# calculate relative weight loss 
ExpeDF_005 <- calculateWeightLoss(ExpeDF_005, startingDay = 0)
ExpeDF_005$relativeWeight <- as.numeric(as.character(ExpeDF_005$relativeWeight))

ExpeDF_005$Sex <- gsub(" ", "", ExpeDF_005$Sex)

saveRDS(ExpeDF_001, file = "ExpeDF_001_Alice.rds")
saveRDS(ExpeDF_002, file = "ExpeDF_002_Alice.rds")
saveRDS(ExpeDF_003_4, file = "ExpeDF_003_4_Alice.rds")
saveRDS(ExpeDF_005, file = "ExpeDF_005_Alice.rds")
