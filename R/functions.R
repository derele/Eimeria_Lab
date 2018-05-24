# functions used for data preparation 
calculateWeightLoss <- function(ExpeDF){
  A = ExpeDF[ExpeDF$dpi == 0, c("weight", "EH_ID")]
  names(A)[1] = "weightAtInfection"
  ExpeDF <- merge(ExpeDF, A)
  rm(A)
  ExpeDF$weightloss = ExpeDF$weightAtInfection - ExpeDF$weight
  ExpeDF$weightRelativeToInfection <- ExpeDF$weight /
    ExpeDF$weightAtInfection * 100
  return(ExpeDF)
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
