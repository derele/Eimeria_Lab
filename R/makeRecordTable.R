makeRecordTable <- function(pathToDesignTable, myseed){
  set.seed(myseed)
  # Load information table
  designTable = read.csv(pathToDesignTable)
  
  recordTable = data.frame(EH_ID = rep(designTable$EH_id, 12),
                           dpi = rep(0:11, each=nrow(designTable)),
                           weight = NA,
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
  labels = sample(combn(LETTERS, 2, paste, collapse = ""), nrow(recordTable))
  #recordTable$labels = labels
  
  recordTable = cbind(labels = labels, recordTable)
  
  return(recordTable)
}

# Experiment May 2018
recordTable <- makeRecordTable("../data/2_designTables/April2018_wildmice_Eferrisi_Firstbatch_DESIGN.csv", myseed = 1344)

write.csv(recordTable, file = "../data/3_recordingTables/April2018_wildmice_Eferrisi_Firstbatch_RECORD.csv", row.names = F)
