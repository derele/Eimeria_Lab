# create oocyst table
## load libraries
library(httr)
library(RCurl)
# load in design table from GitHub
designTable <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E11_012021_Eim_DESIGN.csv")
## ndays = number of days the experiment is running for (12 standard)
## dpi = rep(c( dpis to be recorded))
makeOocystTable <- function(designTable, myseed, ndays = 8){
  set.seed(myseed)
  OocystTable = data.frame(EH_ID = rep(designTable$EH_ID, 8),
                           dpi = rep(c(1,2,3,4,5,6,7,8), each=nrow(designTable)),
                           oocyst_sq1 = "",
                           oocyst_sq2 = "",
                           oocyst_sq3 = "",
                           oocyst_sq4 = "",
                           oocyst_mean = "",
                           dilution = "",
                           OPG = "")
  labels = sample(combn(LETTERS, 3, paste, collapse = ""), nrow(OocystTable))
  OocystTable = cbind(labels = labels, OocystTable)
  return(OocystTable)
}

# load in design table for desired experiment and run through function
## make sure to use same ndays as you need
OocystTableExp <- makeOocystTable(designTable, myseed = 5678, ndays = 8)
# remove EH_ID and dpi to remove bias
OocystTableExp$EH_ID <- NULL
OocystTableExp$dpi <- NULL
# for clarity and avoiding mixups, paste experiment no. and batch denominator (e.g. E10a, E11b, etc.)
OocystTableExp$labels <- sub("^", "E11b", OocystTableExp$labels)
# write out the oocyst table
write.csv(OocystTableExp,"../GitHub/Eimeria_Lab/data/Experiment_results/E11b_012021_Eim_oocyst.csv")
