# makeDesignTable <- function(myseed, pathToInfoTable, firstEH_Id, Inf_strains_list){
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
# }

# raw <- read.csv("../data/1_informationTables/rawData/Exp005_WDS_21.6.18.csv")
# 
# raw$Parent1 = raw$Strain
# raw$Parent2 = raw$Sir
# raw$Strain = paste0(raw$Parent1, "_", raw$Parent2)
# 
# raw$HybridStatus[
#   raw$Strain %in% c("BUSNA_BUSNA","SCHUNT_SCHUNT", "STRA_STRA", "PWD_PWD")] = "parental strains"
# 
# raw$HybridStatus[
#   raw$Strain %in% c("BUSNA_PWD", "PWD_BUSNA", "STRA_SCHUNT", "SCHUNT_STRA")] = "outbred hybrids"
# 
# raw$HybridStatus[
#   raw$Strain %in% c("BUSNA_STRA", "STRA_BUSNA", "PWD_SCHUNT", "SCHUNT_PWD")] = "inter subsp. hybrids"
# 
# write.csv(raw, "../data/1_informationTables/Exp005_WDS_21.6.18.csv", row.names = F)

infoTable = read.csv("../data/1_informationTables/Exp005_WDS_21.6.18.csv")

Nmice = nrow(infoTable)

# Split in 4 groups based on age



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



myDesignTable3 <- makeDesignTable(myseed = 1234 ,
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
