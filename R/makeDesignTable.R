pathToInfoTable = "../data/1_informationTables/April2018_wildmice_Eferrisi_INFO.csv"
# How many mice do we have?
Nmice = 23
# Which one is the first EH_id to be used?
firstEH_Id = "LM0122"
# Which Eimeria isolates do we want to use?
Inf_strains_list = c("E64", "EfLab", "E88", "E139")

# makeDesignTable <- function(pathToInfoTable, Nmice, firstEH_Id, Inf_strains_list){
# Load information table
  infoTable <- read.csv(pathToInfoTable)
  
  # Give EH_ids
  num = as.numeric(sub("LM", "", firstEH_Id))
  num = num + (1:(Nmice))
  EH_id = paste0("LM", sprintf("%04d", num))
  
  # Spread infection equally among sex and Mouse_strain
#  }

makeInfectionExperimentTable(myDFName)


########### ************+
# How many days is the infection? 
lengthExpe <- 11

## makerecordtable
#  EH_id = rep(EH_id, lengthExpe+1)

# Add DPI
dpi = rep(0:lengthExpe, Nmice)
dpi = dpi[order(dpi)]

# Bind all columns in a DF
myDF <- cbind(EH_id, dpi)








  
  