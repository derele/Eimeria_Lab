library(dplyr)
library(tidyr)

# functions used for data preparation 
calculateWeightLoss <- function(x){
  # define weight at infection
  A = x[x$dpi == 0, c("weight", "EH_ID")]
  names(A)[1] = "weightAtInfection"
  x = merge(x, A)
  # Cumulative sum of our weight loss
  x = x %>% 
    group_by(EH_ID) %>% 
    dplyr::arrange(dpi, .by_group = TRUE) %>%
    dplyr::mutate(weightNormal = weight / weightAtInfection) %>%
    dplyr::mutate(weightGainNormal = weightNormal - c(1, weightNormal[-length(weightNormal)])) %>%
    dplyr::mutate(csWeightGainNormal = cumsum(weightGainNormal)) 
  x = data.frame(x)
  return(x)
}

calculateOPG <- function(ExpeDF){
  ExpeDF$mean_Neubauer <- 
    (ExpeDF$Neubauer1 + ExpeDF$Neubauer2 + ExpeDF$Neubauer3 + ExpeDF$Neubauer4) / 4
  ExpeDF$OPG <- ExpeDF$mean_Neubauer * 10000 / ExpeDF$dilution_ml / ExpeDF$fecweight
  return(ExpeDF)
}

# For following before the infection the weight of mice
calculateWeightLossBeforeInf <- function(ExpeDF){
  A = ExpeDF[ExpeDF$dayFollowWeight == 0, c("weight", "original.label")]
  names(A)[1] = "weightAtStart"
  ExpeDF <- merge(ExpeDF, A)
  rm(A)
  ExpeDF$weightlossBeforeInf = ExpeDF$weightAtStart - ExpeDF$weight
  ExpeDF$weightRelativeToStart <- ExpeDF$weight / ExpeDF$weightAtStart * 100
  return(ExpeDF)
}

## Source: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

getCorrectedDates <- function(mytab = read.csv("../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv")){
  goodDates <- as.Date(mytab$Born)
  # # goodDates[which(is.na(goodDates))] <- as.Date(mytab$Born[is.na(goodDates)], format = "%d/%m/%Y")
  # # goodDates <- sub("2018-", replacement = "18-", goodDates)
  # goodDates <- sub("2017-", replacement = "17-", goodDates)
  return(data.frame(original = mytab$Born, 
                    corrected = goodDates))
}

getCorrectedDates()

# write.csv(correctedTable) once checked

# calculate age at infection
getAgeAtInfection <- function(mytab = read.csv("../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv"),
                              infectionDate = "18-06-05"){
  age <- difftime(infectionDate, mytab$Born, units = "weeks")
  return(age)
}

# mytab$ageAtInfection <- getAgeAtInfection()
# write.csv(mytab, "../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv", row.names = F)

# Histogrammes to visualise
# hist(as.numeric(getAgeAtInfection()), breaks = 50)