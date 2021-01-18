# create record table
## load libraries
library(httr)
library(RCurl)
# load in design table from GitHub
designTable <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E11_012021_Eim_DESIGN.csv")
## ndays = number of days the experiment is running for (11 standard for primary, 8 for chalenge)
## dpi = rep(c( dpis to be recorded))

makeRecordTable <- function(designTable, myseed, ndays = 8){
  set.seed(myseed)
  recordTable = data.frame(EH_ID = rep(designTable$EH_ID, 8),
                           dpi = rep(c(1,2,3,4,5,6,7,8), each=nrow(designTable)),
                           weight = "",
                           weight_dpi0 = "",
                           relative_weight = "",
                           feces_weight = "",
                           dpi_dissect = "")
  labels = sample(combn(LETTERS, 3, paste, collapse = ""), nrow(recordTable))
  recordTable = cbind(labels = labels, recordTable)
  return(recordTable)
}
# load in design table for desired experiment and run through function
## make sure to use same ndays as you need
recordTableExpeE11 <- makeRecordTable(designTable, myseed = 5678, ndays = 8)
# for clarity and avoiding mixups, paste experiment no. and batch denominator (e.g. E10a, E11b, etc.)
recordTableExpeE11$labels <- sub("^", "E11b", recordTableExpeE11$labels)
# write out the record table
write.csv(recordTableExpeE11,"../GitHub/Eimeria_Lab/data/Experiment_results/E11b_012021_Eim_record.csv",
          row.names = F, quote = F)
