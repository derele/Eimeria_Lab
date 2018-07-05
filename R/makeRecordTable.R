makeRecordTable <- function(designTable, myseed){
  set.seed(myseed)
  recordTable = data.frame(EH_ID = rep(designTable$EH_id, 12),
                           dpi = rep(0:11, each=nrow(designTable)),
                           weight = NA,
                           weight_dpi0 = NA,
                           weightloss = NA,
                           weightRelativeToInfection = NA,
                           fecweight = NA,
                           oocyst_sq1 = NA,
                           oocyst_sq2 = NA,
                           oocyst_sq3 = NA,
                           oocyst_sq4 = NA,
                           oocyst_mean = NA,
                           dilution = NA,
                           OPG = NA)
  labels = sample(combn(LETTERS, 3, paste, collapse = ""), nrow(recordTable))
  #recordTable$labels = labels
  
  recordTable = cbind(labels = labels, recordTable)
  
  return(recordTable)
}
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

designTable = rbind(tab1a, tab1b,tab2a, tab2b)

## TO RERUN when I know which ones are dead

recordTable = makeRecordTable(designTable, 1234)

# Add Expe id
recordTable$Expe[recordTable$EH_ID %in% tab1a$EH_id] = "Expe005_1a"
recordTable$Expe[recordTable$EH_ID %in% tab1b$EH_id] = "Expe005_1b"
recordTable$Expe[recordTable$EH_ID %in% tab2a$EH_id] = "Expe005_2a"
recordTable$Expe[recordTable$EH_ID %in% tab2b$EH_id] = "Expe005_2b"

# Order by EH_id and dpi
recordTable = recordTable[order(recordTable$dpi, recordTable$EH_ID),]

# 1. weight table
write.csv(recordTable[names(recordTable) %in% c("Expe", "labels", "EH_ID", "dpi", "weight", "weight_dpi0", "weightloss",
                                                "weightRelativeToInfection", "fecweight")],
          file = "../data/3_recordingTables/Exp005_full_RECORDweight.csv", row.names = F)
# 2. oocysts table

write.csv(recordTable[names(recordTable) %in% c("Expe", "labels", "oocyst_sq1","oocyst_sq2","oocyst_sq3","oocyst_sq4",
                                                "oocyst_mean","dilution","OPG")],
          file = "../data/3_recordingTables/Exp005_full_RECORDoocysts.csv", row.names = F)