## Part 1: functions used for data preparation of each experiment
## Part 2: functions used for data analysis, plot, stats, etc.

library(dplyr)
library(tidyr)
library(plyr)

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

getCorrectedDates()

# calculate age at infection
getAgeAtInfection <- function(mytab = read.csv("../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv"),
                              infectionDate = "18-06-05"){
  age <- difftime(infectionDate, mytab$Born, units = "weeks")
  return(age)
}

####### Part 2: functions used for data analysis, plot, stats, etc. ####### 

## Plot comparison weightloss vs OPG, all day counfounded
plotCompOPGWeight <- function(ExpeDF, ylim = c(80,110)){
  
  # Enter Eimeria species
  ExpeDF$Eimeria_species[ExpeDF$infection_isolate %in% c("E88", "Eflab", "EfLab")] = "E.falciformis"
  ExpeDF$Eimeria_species[ExpeDF$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
  
  ggplot(ExpeDF, 
         aes(x = OPG+1, y = weightNormalized, group = Mouse_strain, fill = Mouse_strain)) +
    geom_point(pch = 21, size=3)  +
    facet_grid(Eimeria_species ~ infection_isolate, scales = "free_y", space = "free") +
    geom_smooth(method = "lm", col = "grey") +
    scale_y_continuous(name = "Weight relative to infection (%)") +
    scale_x_log10() +
    coord_cartesian(ylim = ylim) +
    mytheme
}