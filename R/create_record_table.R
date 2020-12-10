# create record table
## load libraries
library(httr)
library(RCurl)
# load in design table from GitHub
designTable <- read.csv(text = getURL("path to raw design table"))
## ndays = number of days the experiment is running for (12 standard)
## dpi = rep(c( dpis to be recorded))

makeRecordTable <- function(designTable, myseed, ndays = 12){
  set.seed(myseed)
  recordTable = data.frame(EH_ID = rep(designTable$EH_ID, 12),
                           dpi = rep(c(1,2,3,4,5,6,7,8,9,10,11,12), each=nrow(designTable)),
                           weight = "",
                           weight_dpi0 = "",
                           relative_weight = "",
                           feces_weight = "",
                           dpi_dissect = "")
  labels = sample(combn(LETTERS, 3, paste, collapse = ""), nrow(recordTable))
  recordTable$labels = labels
  recordTable = cbind(labels = labels, recordTable)
  return(recordTable)
}
# load in design table for desired experiment and run through function
## make sure to use same ndays as you need
designTable <- read.csv("design table for your experiment")
recordTableExpe??? <- makeRecordTable(designTable, myseed = 1234, ndays = 12)
# write out the record table
write.csv(,"")