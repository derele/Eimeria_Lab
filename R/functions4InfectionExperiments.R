library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(plyr)

mytheme <- theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        strip.text = element_text(size=25))

calculateWeightLoss <- function(x, startingDay = 0){
  # define weight at infection
  A = x[x$dpi == startingDay, c("weight", "EH_ID")]
  names(A)[1] = "startingWeight"
  x = merge(x, A)
  # Cumulative sum of our weight loss
  x = x %>% 
    group_by(EH_ID) %>% 
    dplyr::arrange(dpi, .by_group = TRUE) %>%
    dplyr::mutate(relativeWeight = weight / startingWeight * 100)
  x = data.frame(x)
  return(x)
}

calculateOPG <- function(ExpeDF){
  ExpeDF$mean_Neubauer <- 
    (ExpeDF$Neubauer1 + ExpeDF$Neubauer2 + ExpeDF$Neubauer3 + ExpeDF$Neubauer4) / 4
  ExpeDF$Noocysts <- ExpeDF$mean_Neubauer * 10000 / ExpeDF$dilution_ml
  ExpeDF$OPG <- ExpeDF$Noocysts / ExpeDF$fecweight
  return(ExpeDF)
}

## Source: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
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

# calculate age at infection
getAgeAtInfection <- function(mytab = read.csv("../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv"),
                              infectionDate = "18-06-05"){
  age <- difftime(infectionDate, mytab$Born, units = "weeks")
  return(age)
}
