makeRecordTable <- function(designTable, myseed, ndays = 12){
  set.seed(myseed)
  # recordTable = data.frame(EH_ID = rep(designTable$EH_ID, ndays),
  #                          dpi = rep(0:(ndays-1), each=nrow(designTable)),
  #                          weight = "",
  #                          weight_dpi0 = "",
  #                          weightloss = "",
  #                          weightRelativeToInfection = "",
  #                          fecweight = "",
  #                          oocyst_sq1 = "",
  #                          oocyst_sq2 = "",
  #                          oocyst_sq3 = "",
  #                          oocyst_sq4 = "",
  #                          oocyst_mean = "",
  #                          dilution = "",
  #                          OPG = "")
  recordTable = data.frame(EH_ID = rep(designTable$EH_ID, 12),
                           dpi = rep(c(0,2,4,7,9,11,14,16,18,21,23,25), each=nrow(designTable)),
                           weight = "",
                           weight_dpi0 = "",
                           weightloss = "",
                           weightRelativeToInfection = "",
                           fecweight = "",
                           oocyst_sq1 = "",
                           oocyst_sq2 = "",
                           oocyst_sq3 = "",
                           oocyst_sq4 = "",
                           oocyst_mean = "",
                           dilution = "",
                           OPG = "")
  labels = sample(combn(LETTERS, 2, paste, collapse = ""), nrow(recordTable))
  #recordTable$labels = labels
  recordTable = cbind(labels = labels, recordTable)
  return(recordTable)
}

designTableExpe009 <- read.csv("../data/2_designTables/Exp009_CRYPTO_design.csv")
recordTableExpe009 <- makeRecordTable(designTableExpe009, 1234, ndays = 26)
# write.csv(recordTableExpe009, "../data/3_recordingTables/Exp009_CRYPTO.csv",
#           row.names = F, quote = F)
# designTableExpe008 <- read.csv("../data/2_designTables/Exp008_NMRI_DESIGN_jan2019.csv")
# recordTableExpe008 <- makeRecordTable(designTableExpe008, 1234)
# write.csv(recordTableExpe008, file = "../data/3_recordingTables/Exp008_Jan2019_NMRI_4isolatespassaging.csv", row.names = F)
# 
# # Experiment May 2018
# recordTable <- makeRecordTable("../data/2_designTables/April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv", myseed = 1344)
# 
# # 1. weight table
# write.csv(recordTable[names(recordTable) %in% c("labels", "EH_ID", "dpi", "weight", "weightloss",
#                                                 "weightRelativeToInfection", "fecweight")], 
#           file = "../data/3_recordingTables/April2018_wildmice_Eferrisi_Firstbatch_RECORDweight.csv", row.names = F)
# # 2. oocysts table
# write.csv(recordTable[!names(recordTable) %in% c("EH_ID", "dpi", "weight", "weightloss",
#                                                 "weightRelativeToInfection", "fecweight")], 
#           file = "../data/3_recordingTables/April2018_wildmice_Eferrisi_Firstbatch_RECORDoocysts.csv", row.names = F)

# Expe_005
## TO RERUN when I know which ones are dead

# For E. falciformis, we collect feces at dpi 0, 4:11
# For E. ferrisi, we collect feces at dpi 0, 3:10
tab1a = read.csv("../data/2_designTables/Inf1a_Exp005.DESIGN.csv")
tab1b = read.csv("../data/2_designTables/Inf1b_Exp005.DESIGN.csv")
tab2a = read.csv("../data/2_designTables/Inf2a_Exp005.DESIGN.csv")
tab2b = read.csv("../data/2_designTables/Inf2b_Exp005.DESIGN.csv")

correctID = function(x){
  x$EH_id = as.factor(gsub(" ", "", as.character(x$EH_id)))
  return(x)
}

tab1a = correctID(tab1a)
tab1b = correctID(tab1b)
tab2a = correctID(tab2a)
tab2b = correctID(tab2b)

designTable = rbind(tab1a, tab1b,tab2a, tab2b)

## LM0205 and LM0290 were dead before the transport
designTable = designTable[!designTable$EH_id %in% c("LM0205", "LM0290"),]

recordTable = makeRecordTable(designTable, 1234)

# Add Expe id
recordTable$Expe[recordTable$EH_ID %in% tab1a$EH_id] = "Expe005_1a"
recordTable$Expe[recordTable$EH_ID %in% tab1b$EH_id] = "Expe005_1b"
recordTable$Expe[recordTable$EH_ID %in% tab2a$EH_id] = "Expe005_2a"
recordTable$Expe[recordTable$EH_ID %in% tab2b$EH_id] = "Expe005_2b"

# Order by EH_id and dpi
recordTable = recordTable[order(recordTable$Expe, recordTable$dpi, recordTable$EH_ID),]

# 1. weight table
write.csv(recordTable[names(recordTable) %in% c("Expe", "labels", "EH_ID", "dpi", "weight", "weight_dpi0", "weightloss",
                                                "weightRelativeToInfection", "fecweight")],
          file = "../data/3_recordingTables/Exp005_full_RECORDweight.csv", row.names = F)
# 2. oocysts table
write.csv(recordTable[names(recordTable) %in% c("Expe", "labels", "oocyst_sq1","oocyst_sq2","oocyst_sq3","oocyst_sq4",
                                                "oocyst_mean","dilution","OPG")],
          file = "../data/3_recordingTables/Exp005_full_RECORDoocysts.csv", row.names = F)