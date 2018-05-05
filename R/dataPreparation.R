## Data preparation

# ExpeDF <- read.csv("../data/3_recordingTables/March2018_NMRI_4strains_RECORDweightAndOocysts.csv")
# # general
# ExpeDF$Mouse_strain <- "NMRI"

########################### May 2018 batch 1
oo <- read.csv("../data/3_recordingTables/April2018_wildmice_Eferrisi_Firstbatch_RECORDoocysts.csv")
we <- read.csv("../data/3_recordingTables/April2018_wildmice_Eferrisi_Firstbatch_RECORDweight.csv")
design <- read.csv("../data/2_designTables/April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv")

ExpeDFMay2018batch1 <- merge(oo, we, all = T)
ExpeDFMay2018batch1 <- merge(ExpeDFMay2018batch1, design, by = "EH_ID", all = T)

# Calculate weight loss
A <- ExpeDFMay2018batch1[ExpeDFMay2018batch1$dpi == 0, c("weight", "EH_ID")]
names(A)[1] <- "weightAtInfection"
ExpeDFMay2018batch1 <- merge(ExpeDFMay2018batch1, A)

A <- ExpeDFMay2018batch1[ExpeDFMay2018batch1$dpi == -7, c("weight", "EH_ID")]
names(A)[1] <- "weightAtAnthelminthicTrt"
ExpeDFMay2018batch1 <- merge(ExpeDFMay2018batch1, A)

ExpeDFMay2018batch1$weightloss <- ExpeDFMay2018batch1$weightAtInfection - ExpeDFMay2018batch1$weight

ExpeDFMay2018batch1$weightRelativeToInfection <- ExpeDFMay2018batch1$weight /
  ExpeDFMay2018batch1$weightAtInfection * 100

ExpeDFMay2018batch1$weightRelativeToAnthelmTrtDay <- ExpeDFMay2018batch1$weight /
  ExpeDFMay2018batch1$weightAtAnthelminthicTrt * 100

# Calculate OPG
# ExpeDFMay2018batch1$mean_Neubauer <- (ExpeDFMay2018batch1$Neubauer1 + ExpeDFMay2018batch1$Neubauer2 + ExpeDFMay2018batch1$Neubauer3 + ExpeDFMay2018batch1$Neubauer4) / 4
# ExpeDFMay2018batch1$OPG <- ExpeDFMay2018batch1$mean_Neubauer * 10000 / ExpeDFMay2018batch1$dilution_ml / ExpeDFMay2018batch1$fecweight
# ExpeDFMay2018batch1$oocystsTotal <- ExpeDFMay2018batch1$mean_Neubauer * 10000 / ExpeDFMay2018batch1$dilution_ml

# remove junk
rm(A, design, oo, we)
