makeDesignTable <- function(myseed, pathToInfoTable, firstEH_Id, Inf_strains_list){
  # Load information table
  infoTable = read.csv(pathToInfoTable)
  
  Nmice = length(infoTable$original_mouse_id)
  
  # Give EH_ids
  num = as.numeric(sub("LM", "", firstEH_Id))
  num = num + (1:(Nmice))
  EH_id = paste0("LM", sprintf("%04d", num))
  
  # Spread names randomly among mice
  set.seed(myseed)
  infoTable$EH_id <- sample(EH_id)
  
  # Spread infection equally among sex and Mouse_strain
  table(infoTable$sex, infoTable$Mouse_strain)
  
  # Separate equally between Mouse_strains
  library(experiment)
  
  expe <- randomize(data = infoTable, group = Inf_strains_list,
                    indx = infoTable$EH_id, block = infoTable$Mouse_strain)
  
  trt <- data.frame(Eimeria_isolate = expe$treatment)
  trt$EH_id <- rownames(trt)
  rownames(trt) <- NULL
  
  # Create design table
  designTable <- merge(trt, infoTable)
  print(table(designTable$sex, designTable$Mouse_strain, designTable$Eimeria_isolate))
  
  return(designTable)
} 

# Experiment May 2018
myDesignTable1 <- makeDesignTable(myseed = 1234 , pathToInfoTable = "../data/1_informationTables/April2018_wildmice_Eferrisi_Firstbatch_INFO.csv", 
                                  firstEH_Id = "LM0122", Inf_strains_list = c("E64", "E139"))

write.csv(myDesignTable1, file = "../data/2_designTables/April2018_wildmice_Eferrisi_Firstbatch_DESIGN", row.names = F)
