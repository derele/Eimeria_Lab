## Part 1: functions used for data preparation of each experiment
## Part 2: functions used for data analysis, plot, stats, etc.

library(dplyr)
library(tidyr)
library(plyr)

####### Part 1 : functions used for data preparation ####### 
calculateWeightLoss <- function(x, infectionDay = 0){
  # define weight at infection
  A = x[x$dpi == infectionDay, c("weight", "EH_ID")]
  names(A)[1] = "weightAtInfection"
  x = merge(x, A)
  # Cumulative sum of our weight loss
  x = x %>% 
    group_by(EH_ID) %>% 
    dplyr::arrange(dpi, .by_group = TRUE) %>%
    dplyr::mutate(weightNormalized = weight / weightAtInfection * 100)
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

## Weight evolution along the infection
plotWeightAlongInf <- function(ExpeDF, ylim = c(85, 115)){
  # Code by OPG
  ExpeDF$OPG_plot[is.na(ExpeDF$OPG)] = "na"
  ExpeDF$OPG_plot[!is.na(ExpeDF$OPG) & ExpeDF$OPG > 0] = "positive"
  ExpeDF$OPG_plot[!is.na(ExpeDF$OPG) & ExpeDF$OPG == 0] = "negative"
  ExpeDF$OPG_plot = as.factor(ExpeDF$OPG_plot)

  # Enter Eimeria species
  ExpeDF$Eimeria_species[ExpeDF$infection_isolate %in% c("E88", "Eflab", "EfLab")] = "E.falciformis"
  ExpeDF$Eimeria_species[ExpeDF$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
  
  ggplot(ExpeDF, aes(x = dpi, y = weightNormalized, fill = OPG_plot)) +
    geom_line(aes(group = EH_ID, col = infection_isolate), size = 2, alpha = 0.5) +
    geom_point(size=4, pch = 21, color = "black")+
    scale_fill_manual(values = c("lightgrey", "black", "red")) +
    mytheme +
    facet_grid(Eimeria_species ~ Mouse_subspecies, scales = "free_y", space = "free") +
    scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
    scale_y_continuous(name = "Weight relative to infection (%)") +
    coord_cartesian(ylim = ylim)
}

# Mean + 95%CI
plotWeightAlongInfSUM <- function(ExpeDF, ylim = c(85, 115)){
 
  # Enter Eimeria species
  ExpeDF$Eimeria_species[ExpeDF$infection_isolate %in% c("E88", "Eflab")] = "E.falciformis"
  ExpeDF$Eimeria_species[ExpeDF$infection_isolate %in% c("E64", "EI64", "E139")] = "E.ferrisi"
  
  summaryWeight <- summarySE(ExpeDF, measurevar = "weightNormalized",
                             groupvars=c("Mouse_strain", "Mouse_subspecies", "Eimeria_species",
                                         "infection_isolate", "dpi"))
  summaryWeight$ci[is.na(summaryWeight$ci)] <- 0
  
  ggplot(summaryWeight, aes(x = dpi, y = weightNormalized, fill = infection_isolate))+
    geom_errorbar(aes(ymin = weightNormalized - ci,
                      ymax = weightNormalized + ci,
                      col = infection_isolate)) +
    geom_line(aes(group = infection_isolate, col = infection_isolate), size = 2, alpha = 0.5) +
    geom_point(aes(fill = infection_isolate), size=4, pch = 21, color = "black") +
    # scale_fill_manual(values = c("lightgrey", "black", "red")) +
    mytheme +
    facet_grid(Eimeria_species ~ Mouse_subspecies, scales = "free_y", space = "free") +
    scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
    scale_y_continuous(name = "Weight relative to infection (%)") +
    coord_cartesian(ylim = ylim)
}

## OPG evolution along the infection
plotOPGAlongInf <- function(ExpeDF){
  ggplot(ExpeDF, aes(x = dpi, y = OPG, fill = infection_isolate))+
    geom_line(aes(group = EH_ID, col = infection_isolate), alpha = 0.5) +
    geom_point(size=3, pch = 21, color = "black")+
    mytheme +
    facet_grid(Eimeria_species ~ Mouse_subspecies, scales = "free_y", space = "free") +
    scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
    ylab("Oocysts shedding along infection (OPG)")
}

# Mean + 95%CI
plotOPGAlongInfSUM <- function(ExpeDF){
  summaryOPG <- summarySE(ExpeDF, measurevar = "OPG",
                             groupvars=c("Mouse_strain", "Mouse_subspecies",
                                         "infection_isolate", "Eimeria_species", "dpi"))
  summaryOPG$ci[is.na(summaryOPG$ci)] <- 0
  
  ggplot(summaryOPG, aes(x = dpi, y = OPG +1, fill = infection_isolate))+
    geom_errorbar(aes(ymin = OPG +1 - ci,
                      ymax = OPG +1 + ci,
                      col = infection_isolate)) +
    geom_line(aes(group = infection_isolate), alpha = 0.5) +
    geom_point(aes(fill = infection_isolate), size=3, pch = 21, color = "black") +
    mytheme +
    ylab("Oocysts shedding along infection (OPG)") +
    facet_grid(Eimeria_species ~ Mouse_subspecies, scales = "free_y", space = "free") +
    scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") 
}

## Plot comparison weightloss vs OPG, all day counfounded
plotCompOPGWeight <- function(ExpeDF, ylim = c(80,110)){
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

# # maximum weight lost along experiment
# all.max.loss <- do.call("rbind", by(ExpeDF,
#                                 ExpeDF$EH_ID, function (x){
#   m.loss <- which(x$weightRelativeToInfection == min(x$weightRelativeToInfection,
#                                                      na.rm=TRUE))
#   x[m.loss,]
# }))
# 
# all.max.loss$weightLossRelativeToInfection <- 100 - all.max.loss$weightRelativeToInfection
# 
# # What is the mean of weight loss?
# max.loss <- all.max.loss[!duplicated(all.max.loss$EH_ID),]
# 
# tapply(max.loss$weightRelativeToInfection,
#        max.loss$infection_isolate:max.loss$Mouse_strain, 
#        mean)
# 
# ggplot(max.loss,
#        aes(Mouse_strain, weightLossRelativeToInfection)) +
#   geom_violin(color = "black")+
#   geom_jitter(width=0.1, size=7, alpha = 0.8,
#               pch = 21, aes(fill = Mouse_strain)) +
#   labs(y = "Maximum weight loss along infection", x = "Mouse strain") +
#   facet_grid(.~infection_isolate) +
#   mytheme +
#   theme(legend.position = "none")
# 
# if (length(levels(ExpeDF$infection_isolate)) < 2){
#   summary(glm(data = max.loss, weightLossRelativeToInfection ~ Mouse_strain))
# } else if (length(levels(ExpeDF$Mouse_strain)) < 2){
#   summary(glm(data = max.loss, weightLossRelativeToInfection ~ infection_isolate))
# } else {
#   summary(glm(data = max.loss, weightLossRelativeToInfection ~ Mouse_strain*infection_isolate))
# }
# ## using offsets and real weight instead of relative weight for modeling
# 
# # ?? Try Kolmogorv-smirnov approach to compare curves
# # https://stats.stackexchange.com/questions/129449/package-of-r-for-comparing-graphs-of-daily-activity-of-birds
# 
# # Calculate cumulative weight loss
# ggplot(ExpeDF, aes(dpi, ExpeDF$csWeightGainNormal, color = Mouse_strain)) +
#   stat_summary(fun.y = mean, geom = "line") +
#   stat_summary(fun.data = mean_se, geom = "pointrange") +
#   facet_grid(.~ infection_isolate) +
#   mytheme
# 
# # install.packages("lme4", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
# library(lme4)
# 
# # we study the two strains separetely, as we can't compare within (possible different inoculate)
# ExpeDF_E64 <- ExpeDF[ExpeDF$infection_isolate == "E64",]
# 
# ExpeDF_E139 <- ExpeDF[ExpeDF$infection_isolate == "E139",]
# 
# m <- lmer(csWeightGainNormal ~ dpi * Mouse_strain + 
#             (1 | EH_ID), data = ExpeDF_E64, REML = F)
# summary(m, corr = F)
# 
# ExpeDF_E64$m <- fitted(m)
# 
# ggplot(ExpeDF_E64, aes(dpi, csWeightGainNormal)) + 
#   facet_wrap(~ Mouse_strain) + theme_bw() + 
#   stat_summary(fun.y=mean, geom="point") +
#   stat_summary(aes(y=m), fun.y=mean, geom="line") +
#   labs(y="Fixation Proportion", x="Time since word onset (ms)") + 
#   scale_color_manual(values=c("red", "blue"))
# 
# # This fit is SHEIT
# 
# # https://www.sciencedirect.com/science/article/pii/S2468042717300234#bib14
# 
# # To which order should it fit?
# summary(ExpeDF)
# 
# # The first step is to create a third-order polynomial in the range of dpi
# ExpeDF$timebin <- ExpeDF$dpi + 1
# 
# t <- poly(unique(ExpeDF$timebin),3)
# 
# # The next step is to create variables ot1, ot2, ot3 corresponding to the orthogonal polynomial time terms and populate their values with the timeBin-appropriate orthogonal polynomial values:
# ExpeDF[ ,paste0("ot", 1:3)] <- t[ExpeDF$timebin, 1:3]
# 
# # Since this is a simple case with just one within-subjects fixed effect that has only two levels, we can skip to the full model and examine its parameter estimates:
# m.full <- lmer(weightNormal ~ (ot1+ot2+ot3) * Mouse_strain + 
#                      (ot1+ot2+ot3 | EH_ID) + 
#                      (ot1+ot2 | EH_ID:Mouse_strain), 
#                    control = lmerControl(optimizer="bobyqa"),
#                    data=ExpeDF, REML=F)
# 
# coef(summary(m.full))
# 
# # Notice that the parameter estimates do not have p-values. 
# # There are good reasons for that (see this FAQ for more information),
# # but this is cold comfort to most experimental psychologists, who need to report p-values. 
# # The quick and easy solution is to assume that, because we have relatively many
# # observations, the t-distribution converges to the z-distribution, so we can use a
# # normal approximation:
# coefs <- data.frame(coef(summary(m.full))) 
# coefs$p <- format.pval(2*(1-pnorm(abs(coefs$t.value))), digits=2, eps=0.0001) #also make the p-values a bit more readable
# coefs
# 
# ExpeDF$mfit <- fitted(m.full)
# 
# ggplot(ExpeDF, aes(dpi, weightNormal, color = Mouse_strain)) + 
#   theme_bw() + facet_wrap(~ Mouse_strain) +
#   stat_summary(fun.y=mean, geom="point") +
#   stat_summary(aes(y=mfit), fun.y=mean, geom="line") +
#   labs(y="Fixation Proportion", x="Time since word onset (ms)")
# 
# # fourth order?
# # The first step is to create a third-order polynomial in the range of dpi
# ExpeDF$timebin <- ExpeDF$dpi + 1
# 
# t <- poly(unique(ExpeDF$timebin),4)
# 
# # The next step is to create variables ot1, ot2, ot3 corresponding to the orthogonal polynomial time terms and populate their values with the timeBin-appropriate orthogonal polynomial values:
# ExpeDF[ ,paste0("ot", 1:4)] <- t[ExpeDF$timebin, 1:4]
# 
# # Since this is a simple case with just one within-subjects fixed effect that has only two levels, we can skip to the full model and examine its parameter estimates:
# m4 <- lmer(weightNormal ~ (ot1+ot2+ot3) * Mouse_strain + 
#                  (1+ot1+ot2+ot3+ot4 | EH_ID) + 
#                  (1+ot1+ot2 | EH_ID:Mouse_strain), 
#                control = lmerControl(optimizer="bobyqa"),
#                data=ExpeDF, REML=F)
# 
# coefs <- data.frame(coef(summary(m4))) 
# coefs$p <- format.pval(2*(1-pnorm(abs(coefs$t.value))), digits=2, eps=0.0001) #also make the p-values a bit more readable
# coefs
# 
# ExpeDF$mfit <- fitted(m4)
# 
# ggplot(ExpeDF, aes(dpi, weightNormal, color = Mouse_strain)) + 
#   theme_bw() + facet_wrap(~ Mouse_strain) +
#   stat_summary(fun.y=mean, geom="point") +
#   stat_summary(aes(y=mfit), fun.y=mean, geom="line") +
#   labs(y="Fixation Proportion", x="Time since word onset (ms)")
# 
# ###########################################
# # oocyst shedding evolution
# ggplot(ExpeDF, aes(x = dpi, y = OPG))+
#   geom_line(col = "black", aes(group = EH_ID), alpha = 0.5) +
#   geom_point(aes(fill = infection_isolate), 
#              size=3, pch = 21, color = "black", alpha = 0.78) +
#   mytheme +
#   facet_wrap(~Mouse_strain, scales = "free")+
#   scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)" )+
#   geom_smooth(aes(col = infection_isolate)) +
#   scale_y_continuous(labels = scientific)
#   
# # Mean + 95%CI
# summaryOocysts <-summarySE(ExpeDF, measurevar="OPG", 
#                            groupvars=c("Mouse_strain", "infection_isolate", "dpi"))
# 
# ggplot(summaryOocysts, aes(x = dpi, y = OPG))+
#   geom_errorbar(aes(ymin = OPG - ci,
#                     ymax = OPG + ci,
#                     col = infection_isolate)) +
#   geom_line(aes(group = infection_isolate)) +
#   geom_point(aes(fill = infection_isolate), 
#              size=3, pch = 21, color = "black") +
#   mytheme +
#   facet_wrap(~Mouse_strain)+
#   scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)")
# 
# # maximum OPG within expe
# all.max.shed <- do.call("rbind", 
#                         by(ExpeDF,
#                            ExpeDF$EH_ID, function (x){
#                              m.loss <- which(x$OPG == max(x$OPG, na.rm=TRUE))
#                              x[m.loss,]}))
# 
# max.shed <- all.max.shed[!duplicated(all.max.shed$EH_ID),]
# 
# ggplot(max.shed,
#        aes(Mouse_strain, OPG/1000000)) +
#   geom_violin(color = "black")+
#   geom_jitter(width=0.1, size=7, alpha = 0.8,
#               pch = 21, aes(fill = Mouse_strain)) +
#   labs(y = "Maximum OPG at shedding peak (millions)", x = "Mouse strain") +
#   facet_grid(.~infection_isolate) +
#   mytheme +
#   theme(legend.position = "none")
# 
# if (length(levels(ExpeDF$infection_isolate)) < 2){
#   summary(glm(data = max.shed, OPG ~ Mouse_strain))
# } else if (length(levels(ExpeDF$Mouse_strain)) < 2){
#   summary(glm(data = max.shed, OPG ~ infection_isolate))
# } else {
#   summary(glm(data = max.shed, OPG ~ Mouse_strain*infection_isolate))
# }
# 
# # Cumulative OPG along 11 days: 
# all.sum.oocysts <- do.call("rbind", by(ExpeDF, ExpeDF$EH_ID, function (x){
#   x$sum.oo <- sum(x$OPG, na.rm = TRUE)
#   x
# }))
# sum.oocysts <- all.sum.oocysts[!duplicated(all.sum.oocysts$EH_ID),]
# 
# ggplot(sum.oocysts,
#        aes(Mouse_strain, sum.oo/1000000)) +
#   geom_violin(color = "black")+
#   geom_jitter(width=0.1, size=7, alpha = 0.8,
#               pch = 21, aes(fill = Mouse_strain)) +
#   labs(y = "Cumulative OPG along infection (millions)", x = "Mouse strain") +
#   facet_grid(.~infection_isolate) +
#   mytheme +
#   theme(legend.position = "none")
# 
# if (length(levels(ExpeDF$infection_isolate)) < 2){
#   summary(glm(data = sum.oocysts, sum.oo ~ Mouse_strain))
# } else if (length(levels(ExpeDF$Mouse_strain)) < 2){
#   summary(glm(data = sum.oocysts, sum.oo ~ infection_isolate))
# } else {
#   summary(glm(data = sum.oocysts, sum.oo ~ Mouse_strain*infection_isolate))
# }
# 
# ############################
# # highest day of shedding vs highest weight loss?
# shedVsLossdpi <- rbind(data.frame(EH_ID = all.max.loss$EH_ID, 
#                                   what = "max weight loss",
#                                   dpi = all.max.loss$dpi),
#                        data.frame(EH_ID = all.max.shed$EH_ID, 
#                                   what = "max oocyst shedding",
#                                   dpi = all.max.shed$dpi))
# 
# shedVsLossdpi <- merge(shedVsLossdpi, 
#                        unique(ExpeDF[c("EH_ID", "infection_isolate", "Mouse_strain")]), 
#                        all = T)
# 
# 
# dirDF <- merge(data.frame(EH_ID = all.max.loss$EH_ID, 
#                           maxLoss = all.max.loss$dpi),
#                data.frame(EH_ID = all.max.shed$EH_ID, 
#                           maxShed = all.max.shed$dpi))
# 
# if (dirDF$maxShed - dirDF$maxLoss > 0){
#   dirDF$dir <- "darkred"
# } else if {dirDF$maxShed - dirDF$maxLoss < 0){
#   dirDF$dir <- "darkgreen"
# } else {
#   dirDF$dir <- "black"
# }
# 
# ## TBC...
# shedVsLossdpi
# 
# ggplot(shedVsLossdpi, aes(x = what,
#                           y = dpi)) +
#   geom_jitter(aes(fill = Mouse_strain), 
#               size=5, pch = 21, color = "black", alpha = .7, 
#               position = position_jitter(.1, .1)) +
#   geom_line(aes(group = EH_ID)) +
#   facet_grid(.~infection_isolate) +
#   mytheme +
#   scale_y_continuous(breaks = 1:11) +
#   theme(axis.text.x = element_text(angle = 90, size = 10),
#         axis.title.x = element_blank())

