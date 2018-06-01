makeDesignTable <- function(myseed, pathToInfoTable, firstEH_Id, Inf_strains_list){
  # Load information table
  infoTable = read.csv(pathToInfoTable)
  
  Nmice = nrow(infoTable)

  # Give EH_ids
  num = as.numeric(sub("LM", "", firstEH_Id))
  num = num + (1:(Nmice))
  EH_id = paste0("LM", sprintf("%04d", num))

  # Spread names randomly among mice
  set.seed(myseed)
  infoTable$EH_id <- sample(EH_id)

  # Spread infection equally among sex and Mouse_strain
  table(infoTable$Sex, infoTable$Strain)
  return(infoTable)
}

myDesignTable2 <- makeDesignTable(myseed = 1234 ,
                                  pathToInfoTable = "../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv",
                                  firstEH_Id = "LM0145")

# Separate equally between Mouse_strains
library(experiment)

expe <- randomize(data = myDesignTable2, group = c("E64", "E139"),
                   indx = myDesignTable2$EH_id, block = myDesignTable2$Strain)

trt <- data.frame(infection_isolate = expe$treatment)
trt$EH_id <- rownames(trt)
rownames(trt) <- NULL

# Create design table
designTable <- merge(trt, myDesignTable2)
print(table(designTable$Sex, designTable$Strain, designTable$infection_isolate))

write.csv(designTable, "../data/2_designTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_DESIGN.csv", row.names = F)