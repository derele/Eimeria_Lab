# Alice Balard
# useful functions for data analysis infection experiments
# March 2019

# Load libraries
if(!require(multcomp)){install.packages("multcomp")}
listLib <- c("ggplot2", "gridExtra", "reshape2", "scales", "lme4", 
             "lmerTest", "plyr", "dplyr", "tidyr", "reshape", "multcomp")
lapply(listLib, require, character.only = TRUE)

# Define a theme
mytheme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
theme_set(mytheme)

# To calculate means & 95%CI
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
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# calculate weightloss
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

# calculate OPG
calculateOPG <- function(ExpeDF){
  ExpeDF$mean_Neubauer <- 
    (ExpeDF$Neubauer1 + ExpeDF$Neubauer2 + ExpeDF$Neubauer3 + ExpeDF$Neubauer4) / 4
  # NB! Limit of detection = 1 oocysts
  ExpeDF$mean_Neubauer[ExpeDF$Neubauer1 + ExpeDF$Neubauer2 + ExpeDF$Neubauer3 + ExpeDF$Neubauer4 == 1] <- 0
  ExpeDF$oocysts.per.tube <- ExpeDF$mean_Neubauer * 10000 * ExpeDF$dilution_ml
  ExpeDF$OPG <- ExpeDF$oocysts.per.tube / ExpeDF$fecweight
  ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
  ExpeDF$oocysts.per.tube[ExpeDF$fecweight == 0 & ExpeDF$mean_Neubauer == 0] <- 0
  ExpeDF$OPG[ExpeDF$fecweight == 0 & ExpeDF$mean_Neubauer == 0] <- 0
  return(ExpeDF)
}

# calculate maximum realtive weight loss
getMaxLoss <- function(df){
  max.loss <- do.call("rbind", by(df, df$EH_ID, function (x){
    m.loss <- which(x$relativeWeight == min(x$relativeWeight, na.rm=TRUE))
    x[m.loss,]
  }))
  max.loss <- max.loss[!duplicated(max.loss$EH_ID),]
  names(max.loss)[names(max.loss) %in% "dpi"] <- "dpi_maxweightloss"
  names(max.loss)[names(max.loss) %in% "relativeWeight"] <- "minRelativeWeight"
  return(max.loss)
}

# calculate day with higher shedding peak
getMaxOPG <- function(df){
  max.opg <- do.call("rbind", by(df, df$EH_ID, function (x){
    m.opg <- which(x$OPG == max(x$OPG, na.rm=TRUE))
    x[m.opg,]
  }))
  max.opg <- max.opg[!duplicated(max.opg$EH_ID),]
  names(max.opg)[names(max.opg) %in% "dpi"] <- "dpi_maxOPG"
  names(max.opg)[names(max.opg) %in% "OPG"] <- "maxOPG"
  max.opg$maxOPG_inmillion = max.opg$maxOPG/1e6
  return(max.opg)
}

# to take with care if the animal is still infected when sacrificed...
getSumOPG <- function(df){
  all.sum.opg <- do.call("rbind", by(df, df$EH_ID, function (x){
    x$sum.opg <- sum(x$OPG, na.rm = TRUE)
    x
  }))
  sumOpg <- all.sum.opg[!duplicated(all.sum.opg$EH_ID),]
  sumOpg$sum.oocysts_inmillion = sumOpg$sum.opg/1e6
  return(sumOpg)
}

makeMiceGenotypeAndIsolate <- function(df){
  # NB. Let's not consider which parent is which, but make A_B mouse = B_A mouse
  # we don't have enough individuals to test this effect, and we are not interested in it anyway!
  df$Mouse_strain <- as.character(df$Mouse_strain)
  x <- strsplit(df$Mouse_strain, "_")
  y <- lapply(x, sort)
  z <- unlist(lapply(y, FUN = function(x){paste(x, collapse="-")}))
  df$Mouse_genotype <- z
  ## Order the levels to be more clear in later plots (parents will be low and down 
  ## on the legend, hybrids in between...)
  df$Mouse_genotype <- factor(df$Mouse_genotype,
                              levels = c("NMRI", 
                                         "WSB", "WP", "PWD1",
                                         "BUSNA-BUSNA", "PWD-PWD",
                                         "BUSNA-PWD",
                                         "BUSNA-STRA","PWD-SCHUNT",
                                         "SCHUNT-STRA",
                                         "SCHUNT-SCHUNT","STRA-STRA"),
                              labels = c("NMRI", 
                                         "MMd_F0 (Ws-Ws)", "Mmm-Mmd_Hybrid (WP)", "MMm_F0 (Pw1-Pw1)",
                                         "MMm_F0 (Bu-Bu)",
                                         "MMm_F0 (Pw-Pw)",
                                         "MMm_F1 (Bu-Pw)",
                                         "Mmm-Mmd_F1 (Bu-St)",
                                         "Mmm-Mmd_F1 (Pw-Sc)",
                                         "MMd_F1 (Sc-St)",
                                         "MMd_F0 (Sc-Sc)",
                                         "MMd_F0 (St-St)"))
  df$infection_isolate <- factor(df$infection_isolate,
                                 levels = c("E139", "E64", "E88", "EfLab"),
                                 labels = c("E.ferrisi (E139)",
                                            "E.ferrisi (E64)",
                                            "E.falciformis (E88)",
                                            "E.falciformis (EfLab)"))
  # erase useless level
  df$infection_isolate <- droplevels(df$infection_isolate)
  df$Mouse_genotype <- droplevels(df$Mouse_genotype)
  return(df)
}

# create a table with tolerance factor, max.loss, max.OPG and sum.oocysts concatenated
makeSummaryTable <- function(df){
  # minimun weight and associated original weight
  X <- as.data.frame(
    df %>% 
      dplyr::group_by(EH_ID) %>%
      dplyr::slice(which.min(weight)) %>%
      dplyr::select(EH_ID, weight, startingWeight, ageAtInfection, Sex,
             Mouse_genotype, Eimeria_species, infection_isolate, Exp_ID, dpi))
  # time to min host weight loss peak
  names(X)[names(X) %in% "dpi"] = "dpi_minWeight"
  names(X)[names(X) %in% "weight"] = "minWeight"
  X$minWeightRelative <- X$minWeight / X$startingWeight * 100
  # maximum oocysts and associated fecweight
  Y <- as.data.frame(
    df %>% 
      dplyr::group_by(EH_ID) %>%
      dplyr::slice(which.max(oocysts.per.tube)) %>%
      dplyr::select(EH_ID, oocysts.per.tube, fecweight, dpi))
  # time to parasite shedding peak
  names(Y)[names(Y) %in% "dpi"] = "dpi_max.oocysts.per.tube"
  names(Y)[names(Y) %in% "oocysts.per.tube"] = "max.oocysts.per.tube"
  Y$max.OPG <- Y$max.oocysts.per.tube / Y$fecweight
  # merge
  fullDF <- merge(X, Y)
  # tolerance factor (HF/PF = host min weight/parasite maximum shedding)
  fullDF$tolfacAbsolute <- fullDF$minWeight / fullDF$max.oocysts.per.tube
  fullDF$tolfacRelative <- fullDF$minWeightRelative / fullDF$max.OPG
  fullDF$ageAtInfection <- as.numeric(fullDF$ageAtInfection)
  # add age categorie (<25wo = young, >50we = old)
  fullDF$ageCat <- NA
  fullDF$ageCat[fullDF$ageAtInfection < 25] <- "young"
  fullDF$ageCat[fullDF$ageAtInfection > 50] <- "old"
  fullDF$ageCat <- relevel(as.factor(fullDF$ageCat), ref = "young")
  return(fullDF)
}

# Make groups boxplots
myboxPlot <- 
  function(tolerance, response, respName, group, groupName){
    ggplot(tolerance, 
           aes_string(x = group, y = response,
                      fill = group)) +
      geom_boxplot() +
      geom_jitter(position = position_jitter(0.1), alpha =0.5)+
      scale_x_discrete(name = groupName)+
      theme(legend.position = "none") +
      facet_grid(.~infection_isolate)
  }

#https://rdrr.io/github/ProjectMOSAIC/mosaic/src/R/Tukey.R
# Compute Tukey Honest Significant Differences
# Create a set of confidence intervals on the differences between the means of the levels of a factor with the specified family-wise probability of coverage. The intervals are based on the Studentized range statistic, Tukey's ‘Honest Significant Difference’ method.
TukeyHSD.lm <- function(x, which, ordered = FALSE, conf.level=0.95, ...) {
  stats::TukeyHSD( aov(x), which = which, ordered = ordered, conf.level = conf.level, ...)
}

mytukey <- function(m1){
  mytukey = TukeyHSD.lm(m1)
  mytukeyDF = lapply(mytukey, function(x) {
    df <- data.frame(x)
    data.frame(combi = rownames(df[df$p.adj < 0.05 & !is.na(df$p.adj),]),
               padj = df$p.adj[df$p.adj < 0.05 & !is.na(df$p.adj)],
               diff = df$diff[df$p.adj < 0.05 & !is.na(df$p.adj)])})
  return(mytukeyDF)
}

# linear modelling
myMod <- function(mydata, response, group){
  m1 = lm(response ~ group * infection_isolate,
          data=mydata)
  mytukey = mytukey(m1)
  return(list(model = m1,
              summary.model = summary(m1),
              mytukeyHSD = mytukeyDF))
}

## Plots to follow experiment course
makeIntermPlots <- function(df){
  df$relativeWeightLoss <- 100 - df$relativeWeight
  # calculate summary statistics on weight loss
  summaryWeight <- summarySE(df, measurevar = "relativeWeightLoss",
                             groupvars=c("Mouse_genotype",
                                         "Eimeria_species", "dpi"), na.rm = T)
  summaryWeight$ci[is.na(summaryWeight$ci)] <- 0
  
  # plot mean relative weightloss
  P1 <- ggplot(summaryWeight, aes(x = dpi, y = relativeWeightLoss))+
    geom_errorbar(aes(ymin = relativeWeightLoss - ci,
                      ymax = relativeWeightLoss + ci),
                  col = "gray")+
    geom_line(aes(group = Mouse_genotype, col = Mouse_genotype), size = 2, alpha = 0.5) +
    geom_point(aes(fill = Mouse_genotype), size=4, pch = 21, color = "black") +
    facet_grid(. ~ Eimeria_species, scales = "free_y", space = "free") +
    scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
    scale_y_continuous(name = "Weight loss relative to dpi0 (%)") +
    theme(axis.text.x = element_text(angle = 0))
  
  # calculate summary statistics on oocysts shedding
  summaryOocysts <-summarySE(df, measurevar="OPG",
                             groupvars=c("Mouse_genotype",
                                         "Eimeria_species", "dpi"), na.rm = T)
  summaryOocysts$ci[is.na(summaryOocysts$ci)] <- 0
  
  # mean parasite shedding
  P2 <- ggplot(summaryOocysts, aes(x = dpi, y = OPG +1))+
    geom_errorbar(aes(ymin = OPG - ci, ymax = OPG + ci),
                  col = "gray") +
    geom_line(aes(group = Mouse_genotype, col = Mouse_genotype), size = 2, alpha = 0.5) +
    geom_point(aes(fill = Mouse_genotype), size=4, pch = 21, color = "black") +
    facet_grid(. ~ Eimeria_species) +
    scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
    scale_y_log10(name = "Oocysts shedding along infection (OPG)",
                  breaks = c(1, 11, 101, 1001, 10001, 100001, 1000001, 10000001), 
                  labels = format(c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
                                  scientific = TRUE)) +
    theme(axis.text.x = element_text(angle = 0))
  return(list(P1, P2))
}

# Extra functions
getCorrectedDates <- function(mytab = read.csv("../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv")){
  goodDates <- as.Date(mytab$Born)
  # # goodDates[which(is.na(goodDates))] <- as.Date(mytab$Born[is.na(goodDates)], format = "%d/%m/%Y")
  # # goodDates <- sub("2018-", replacement = "18-", goodDates)
  # goodDates <- sub("2017-", replacement = "17-", goodDates)
  return(data.frame(original = mytab$Born, 
                    corrected = goodDates))
}
