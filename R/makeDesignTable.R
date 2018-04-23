pathToInfoTable = "../data/1_informationTables/April2018_wildmice_Eferrisi_Firstbatch_INFO.csv"

# Which one is the first EH_id to be used?
firstEH_Id = "LM0122"

# Which Eimeria isolates do we want to use?
Inf_strains_list = c("E64", "E139")

# makeDesignTable <- function(pathToInfoTable, Nmice, firstEH_Id, Inf_strains_list){
# Load information table
  infoTable = read.csv(pathToInfoTable)
  
  Nmice = length(infoTable$original_mouse_id)
  
  # Give EH_ids
  num = as.numeric(sub("LM", "", firstEH_Id))
  num = num + (1:(Nmice))
  EH_id = paste0("LM", sprintf("%04d", num))
  
  # Spread names randomly among mice
  set.seed(1567)
  infoTable$EH_id <- sample(EH_id)
  
  # Spread infection equally among sex and Mouse_strain
  table(infoTable$sex, infoTable$Mouse_strain)
  
  infoTable$original_mouse_id[infoTable$sex >  9]) 


infoTable$mygroups <- paste0(infoTable$sex, infoTable$Mouse_strain)

sample(infoTable$EH_id[infoTable$mygroups %in% "FPWD"],
       1.5)
       
data.frame(table(infoTable$mygroups)/2)

#AAAAAAAAAAAAAAAAAAAAA


library(dplyr)

G1 <- infoTable %>% group_by(`mygroups`) %>% 
  filter(row_number() %in% sample(seq_len(n()), n()/2))

table(G1$mygroups)
table(G1$Mouse_strain)
table(infoTable$Mouse_strain, infoTable$sex)
table(G1$Mouse_strain, G1$sex)

3/2

infoTable$Eimeria_isolate <- "E139"
infoTable$Eimeria_isolate[infoTable$EH_id %in% G1$EH_id] <- "E64"

table(infoTable$Eimeria_isolate)
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








  
  