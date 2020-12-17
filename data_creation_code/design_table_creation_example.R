# create design table
# Load information table
library(RCurl)
# load in initial dataset from GitHub (must be raw.)
infoTable = read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_creation_code/E10_112020_Eim_INFO.csv"))
# check last experiment and get highest EH_ID
lastEH_ID <- "LM0399"
# divide dataset into groups as desired
## e.g.: 100 mice would divide into c(rep("isolate1", 25), rep("isolate2", 25), 
##                                  rep("isolate3", 25), rep("uninfected", 25))
# for challenge infection plans this has to be spread further
primary_infection <- c(rep("E64", 11), rep("E88", 11), rep("UNI", 10))
challenge_infection <- c(rep("E64", 11), rep("E88", 11), rep("UNI", 10))
# number of mice
Nmice = nrow(infoTable)
#Give EH_IDs
num = as.numeric(sub("LM", "", lastEH_ID))
num = num + (1:(Nmice))
EH_ID = paste0("LM", sprintf("%04d", num))
#Assign infection isolate
designTable <- data.frame(primary_infection = primary_infection,
                          challenge_infection = challenge_infection,
                          EH_ID= EH_ID)
# Spread names and infections randomly among mice (restrict by total amount of mice)
infoTable$EH_ID <- sample(EH_ID, size = 32)
infoTable$primary_infection <- sample(primary_infection, size = 32)
infoTable$challenge_infection <- sample(challenge_infection, size = 32)
# merge
finaldesignTable <- merge(infoTable, designTable, all.x = T)
####### check necessary columns at https://github.com/derele/Eimeria_Lab
# add experiment column
finaldesignTable$experiment <- "E10"
# rename columns to match other design tables as stated in the repo
names(finaldesignTable)[names(finaldesignTable) == "Strain"] <- "mouse_strain"
finaldesignTable$infection_history <- paste(finaldesignTable$primary_infection, finaldesignTable$challenge_infection, sep = ":")


# write out
write.csv(finaldesignTable,
          "~/GitHub/Eimeria_Lab/data_creation_code/design_table_creation_example.csv",
          row.names = F, quote = F )
