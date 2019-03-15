# Load information table

###### data for Expe009
pathToInfoTable <- "../data/1_informationTables/Exp009_CRYPTO_march2019.csv"
lastEH_ID <- "LM0303"
infection_isolate <- c(rep("Cry1", 4), rep("Cry2", 4), rep("Cry3", 4), rep("Cry4", 4))
######

# infoTable = read.csv(pathToInfoTable)
# Nmice = nrow(infoTable)
# 
# # Give EH_IDs
# num = as.numeric(sub("LM", "", lastEH_ID))
# num = num + (1:(Nmice))
# EH_ID = paste0("LM", sprintf("%04d", num))
# 
# # Assign infection isolate
# designTable <- data.frame(infection_isolate = infection_isolate,
#                           EH_ID= EH_ID)
# 
# # Spread names randomly among mice
# infoTable$EH_ID <- sample(EH_ID)
# 
# # merge full stuff
# finaldesignTable <- merge(infoTable, designTable)
# 
# write.csv(finaldesignTable,
#           "../data/2_designTables/Exp009_CRYPTO_design.csv", 
#           row.names = F, quote = F )

# library(experiment)
# 
# expe <- randomize(data = myDesignTable2, group = c("E64", "E139"),
#                   indx = myDesignTable2$EH_ID, block = myDesignTable2$Strain)

# trt <- data.frame(infection_isolate = expe$treatment)
# trt$EH_ID <- rownames(trt)
# rownames(trt) <- NULL





