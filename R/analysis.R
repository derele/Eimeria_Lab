source("dataPreparation.R")

###########################################
# Choose the experiment:

##### Expe_001
# March 2017, Francisca's experiment. infection with E64 and Eflab
ExpeDF <- ExpeDF_001

##### Expe_002
# March 2018 NMRI infected with 4 strains
ExpeDF <- ExpeDF_002

##### Expe_003
# April-May 2018, first batch Parental strains (F0) BUSNA, STRA, SCHUNT, PWD
# Infection with Eferrisi (E64 and E139)
ExpeDF <- ExpeDF_003

##### Expe_004
# June 2018, second batch Parental strains (F0) BUSNA, STRA, SCHUNT, PWD
# Infection with Eferrisi (E64 and E139)
ExpeDF <- ExpeDF_004

## Full F0s?
# ExpeDF <- rbind(ExpeDF_003[c("EH_ID", "dpi", "weight", "weightRelativeToInfection", "cumsumWeightLosRelToInfPercent", "infection_isolate", "Mouse_strain")],
#                 ExpeDF_004[c("EH_ID", "dpi", "weight", "weightRelativeToInfection", "cumsumWeightLosRelToInfPercent", "infection_isolate", "Mouse_strain")])
# ExpeDF <- ExpeDF[ExpeDF$dpi %in% 0:11, ]

###########################################
## Weight evolution along the infection

ggplot(ExpeDF, aes(x = dpi, y = weightNormal, fill = infection_isolate))+
  geom_line(aes(group = EH_ID), alpha = 0.3) +
  geom_point(size=3, pch = 21, color = "black")+
  geom_smooth(aes(col = infection_isolate), alpha = 0.3) +
  mytheme +
  facet_grid(~Mouse_strain, scales = "free_y", space = "free") +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)" ) 

# Mean + 95%CI
summaryWeight <-summarySE(ExpeDF, measurevar = "weightRelativeToInfection", 
                           groupvars=c("Mouse_strain", "infection_isolate", "dpi"))

# ggplot(summaryWeight, aes(x = dpi, y = weightRelativeToInfection))+
#   geom_errorbar(aes(ymin = weightRelativeToInfection - ci,
#                     ymax = weightRelativeToInfection + ci,
#                     col = infection_isolate)) +
#   geom_line(aes(group = infection_isolate)) +
#   geom_point(aes(fill = infection_isolate), 
#              size=3, pch = 21, color = "black") +
#   mytheme +
#   facet_wrap(~Mouse_strain)+
#   scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
#   scale_y_continuous(name = "Weight relative to infection (%)" )

# Now, using the CUMULATIVE weight loss
ggplot(ExpeDF, aes(x = dpi, y = csWeightGainNormal, fill = infection_isolate))+
  geom_smooth(aes(col = infection_isolate), alpha = 0.3) +
  geom_line(aes(group = EH_ID), col = "black", alpha = 0.5) +
  geom_point(size=3, pch = 21, color = "black")+
  mytheme +
  facet_grid(~Mouse_strain, scales = "free_y", space = "free") +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)" ) 

# maximum weight lost along experiment
all.max.loss <- do.call("rbind", by(ExpeDF,
                                ExpeDF$EH_ID, function (x){
  m.loss <- which(x$weightRelativeToInfection == min(x$weightRelativeToInfection,
                                                     na.rm=TRUE))
  x[m.loss,]
}))

all.max.loss$weightLossRelativeToInfection <- 100 - all.max.loss$weightRelativeToInfection

# What is the mean of weight loss?
max.loss <- all.max.loss[!duplicated(all.max.loss$EH_ID),]

tapply(max.loss$weightRelativeToInfection,
       max.loss$infection_isolate:max.loss$Mouse_strain, 
       mean)

ggplot(max.loss,
       aes(Mouse_strain, weightLossRelativeToInfection)) +
  geom_violin(color = "black")+
  geom_jitter(width=0.1, size=7, alpha = 0.8,
              pch = 21, aes(fill = Mouse_strain)) +
  labs(y = "Maximum weight loss along infection", x = "Mouse strain") +
  facet_grid(.~infection_isolate) +
  mytheme +
  theme(legend.position = "none")

if (length(levels(ExpeDF$infection_isolate)) < 2){
  summary(glm(data = max.loss, weightLossRelativeToInfection ~ Mouse_strain))
} else if (length(levels(ExpeDF$Mouse_strain)) < 2){
  summary(glm(data = max.loss, weightLossRelativeToInfection ~ infection_isolate))
} else {
  summary(glm(data = max.loss, weightLossRelativeToInfection ~ Mouse_strain*infection_isolate))
}
## using offsets and real weight instead of relative weight for modeling

# ?? Try Kolmogorv-smirnov approach to compare curves
# https://stats.stackexchange.com/questions/129449/package-of-r-for-comparing-graphs-of-daily-activity-of-birds

# Calculate cumulative weight loss
ggplot(ExpeDF, aes(dpi, ExpeDF$csWeightGainNormal, color = Mouse_strain)) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  facet_grid(.~ infection_isolate) +
  mytheme

# install.packages("lme4", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
library(lme4)

# we study the two strains separetely, as we can't compare within (possible different inoculate)
ExpeDF_E64 <- ExpeDF[ExpeDF$infection_isolate == "E64",]
ExpeDF_E139 <- ExpeDF[ExpeDF$infection_isolate == "E139",]

m <- lmer(csWeightGainNormal ~ dpi * Mouse_strain + 
            (1 | EH_ID), data = ExpeDF_E64, REML = F)
summary(m, corr = F)

ExpeDF_E64$m <- fitted(m)

ggplot(ExpeDF_E64, aes(dpi, csWeightGainNormal)) + 
  facet_wrap(~ Mouse_strain) + theme_bw() + 
  stat_summary(fun.y=mean, geom="point") +
  stat_summary(aes(y=m), fun.y=mean, geom="line") +
  labs(y="Fixation Proportion", x="Time since word onset (ms)") + 
  scale_color_manual(values=c("red", "blue"))

# This fit is SHEIT

# https://www.sciencedirect.com/science/article/pii/S2468042717300234#bib14

# To which order should it fit?
summary(ExpeDF)

# The first step is to create a third-order polynomial in the range of dpi
ExpeDF$timebin <- ExpeDF$dpi + 1

t <- poly(unique(ExpeDF$timebin),3)

# The next step is to create variables ot1, ot2, ot3 corresponding to the orthogonal polynomial time terms and populate their values with the timeBin-appropriate orthogonal polynomial values:
ExpeDF[ ,paste0("ot", 1:3)] <- t[ExpeDF$timebin, 1:3]

# Since this is a simple case with just one within-subjects fixed effect that has only two levels, we can skip to the full model and examine its parameter estimates:
m.full <- lmer(weightNormal ~ (ot1+ot2+ot3) * Mouse_strain + 
                     (ot1+ot2+ot3 | EH_ID) + 
                     (ot1+ot2 | EH_ID:Mouse_strain), 
                   control = lmerControl(optimizer="bobyqa"),
                   data=ExpeDF, REML=F)

coef(summary(m.full))

# Notice that the parameter estimates do not have p-values. 
# There are good reasons for that (see this FAQ for more information),
# but this is cold comfort to most experimental psychologists, who need to report p-values. 
# The quick and easy solution is to assume that, because we have relatively many
# observations, the t-distribution converges to the z-distribution, so we can use a
# normal approximation:
coefs <- data.frame(coef(summary(m.full))) 
coefs$p <- format.pval(2*(1-pnorm(abs(coefs$t.value))), digits=2, eps=0.0001) #also make the p-values a bit more readable
coefs

ExpeDF$mfit <- fitted(m.full)

ggplot(ExpeDF, aes(dpi, weightNormal, color = Mouse_strain)) + 
  theme_bw() + facet_wrap(~ Mouse_strain) +
  stat_summary(fun.y=mean, geom="point") +
  stat_summary(aes(y=mfit), fun.y=mean, geom="line") +
  labs(y="Fixation Proportion", x="Time since word onset (ms)")

# fourth order?
# The first step is to create a third-order polynomial in the range of dpi
ExpeDF$timebin <- ExpeDF$dpi + 1

t <- poly(unique(ExpeDF$timebin),4)

# The next step is to create variables ot1, ot2, ot3 corresponding to the orthogonal polynomial time terms and populate their values with the timeBin-appropriate orthogonal polynomial values:
ExpeDF[ ,paste0("ot", 1:4)] <- t[ExpeDF$timebin, 1:4]

# Since this is a simple case with just one within-subjects fixed effect that has only two levels, we can skip to the full model and examine its parameter estimates:
m4 <- lmer(weightNormal ~ (ot1+ot2+ot3) * Mouse_strain + 
                 (1+ot1+ot2+ot3+ot4 | EH_ID) + 
                 (1+ot1+ot2 | EH_ID:Mouse_strain), 
               control = lmerControl(optimizer="bobyqa"),
               data=ExpeDF, REML=F)

coefs <- data.frame(coef(summary(m4))) 
coefs$p <- format.pval(2*(1-pnorm(abs(coefs$t.value))), digits=2, eps=0.0001) #also make the p-values a bit more readable
coefs

ExpeDF$mfit <- fitted(m4)

ggplot(ExpeDF, aes(dpi, weightNormal, color = Mouse_strain)) + 
  theme_bw() + facet_wrap(~ Mouse_strain) +
  stat_summary(fun.y=mean, geom="point") +
  stat_summary(aes(y=mfit), fun.y=mean, geom="line") +
  labs(y="Fixation Proportion", x="Time since word onset (ms)")

###########################################
# oocyst shedding evolution
ggplot(ExpeDF, aes(x = dpi, y = OPG))+
  geom_line(col = "black", aes(group = EH_ID), alpha = 0.5) +
  geom_point(aes(fill = infection_isolate), 
             size=3, pch = 21, color = "black", alpha = 0.78) +
  mytheme +
  facet_wrap(~Mouse_strain*infection_isolate, scales = "free")+
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)" )+
  geom_smooth(aes(col = infection_isolate)) +
  scale_y_continuous(labels = scientific)
  
# Mean + 95%CI
summaryOocysts <-summarySE(ExpeDF, measurevar="OPG", 
                           groupvars=c("Mouse_strain", "infection_isolate", "dpi"))

ggplot(summaryOocysts, aes(x = dpi, y = OPG))+
  geom_errorbar(aes(ymin = OPG - ci,
                    ymax = OPG + ci,
                    col = infection_isolate)) +
  geom_line(aes(group = infection_isolate)) +
  geom_point(aes(fill = infection_isolate), 
             size=3, pch = 21, color = "black") +
  mytheme +
  facet_wrap(~Mouse_strain)+
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)")

# maximum OPG within expe
all.max.shed <- do.call("rbind", 
                        by(ExpeDF,
                           ExpeDF$EH_ID, function (x){
                             m.loss <- which(x$OPG == max(x$OPG, na.rm=TRUE))
                             x[m.loss,]}))

max.shed <- all.max.shed[!duplicated(all.max.shed$EH_ID),]

ggplot(max.shed,
       aes(Mouse_strain, OPG/1000000)) +
  geom_violin(color = "black")+
  geom_jitter(width=0.1, size=7, alpha = 0.8,
              pch = 21, aes(fill = Mouse_strain)) +
  labs(y = "Maximum OPG at shedding peak (millions)", x = "Mouse strain") +
  facet_grid(.~infection_isolate) +
  mytheme +
  theme(legend.position = "none")

if (length(levels(ExpeDF$infection_isolate)) < 2){
  summary(glm(data = max.shed, OPG ~ Mouse_strain))
} else if (length(levels(ExpeDF$Mouse_strain)) < 2){
  summary(glm(data = max.shed, OPG ~ infection_isolate))
} else {
  summary(glm(data = max.shed, OPG ~ Mouse_strain*infection_isolate))
}

# Cumulative OPG along 11 days: 
all.sum.oocysts <- do.call("rbind", by(ExpeDF, ExpeDF$EH_ID, function (x){
  x$sum.oo <- sum(x$OPG, na.rm = TRUE)
  x
}))
sum.oocysts <- all.sum.oocysts[!duplicated(all.sum.oocysts$EH_ID),]

ggplot(sum.oocysts,
       aes(Mouse_strain, sum.oo/1000000)) +
  geom_violin(color = "black")+
  geom_jitter(width=0.1, size=7, alpha = 0.8,
              pch = 21, aes(fill = Mouse_strain)) +
  labs(y = "Cumulative OPG along infection (millions)", x = "Mouse strain") +
  facet_grid(.~infection_isolate) +
  mytheme +
  theme(legend.position = "none")

if (length(levels(ExpeDF$infection_isolate)) < 2){
  summary(glm(data = sum.oocysts, sum.oo ~ Mouse_strain))
} else if (length(levels(ExpeDF$Mouse_strain)) < 2){
  summary(glm(data = sum.oocysts, sum.oo ~ infection_isolate))
} else {
  summary(glm(data = sum.oocysts, sum.oo ~ Mouse_strain*infection_isolate))
}

############################
# highest day of shedding vs highest weight loss?
shedVsLossdpi <- rbind(data.frame(EH_ID = all.max.loss$EH_ID, 
                                  what = "max weight loss",
                                  dpi = all.max.loss$dpi),
                       data.frame(EH_ID = all.max.shed$EH_ID, 
                                  what = "max oocyst shedding",
                                  dpi = all.max.shed$dpi))

shedVsLossdpi <- merge(shedVsLossdpi, 
                       unique(ExpeDF[c("EH_ID", "infection_isolate", "Mouse_strain")]), 
                       all = T)


dirDF <- merge(data.frame(EH_ID = all.max.loss$EH_ID, 
                          maxLoss = all.max.loss$dpi),
               data.frame(EH_ID = all.max.shed$EH_ID, 
                          maxShed = all.max.shed$dpi))

if (dirDF$maxShed - dirDF$maxLoss > 0){
  dirDF$dir <- "darkred"
} else if {dirDF$maxShed - dirDF$maxLoss < 0){
  dirDF$dir <- "darkgreen"
} else {
  dirDF$dir <- "black"
}

## TBC...
shedVsLossdpi

ggplot(shedVsLossdpi, aes(x = what,
                          y = dpi)) +
  geom_jitter(aes(fill = Mouse_strain), 
              size=5, pch = 21, color = "black", alpha = .7, 
              position = position_jitter(.1, .1)) +
  geom_line(aes(group = EH_ID)) +
  facet_grid(.~infection_isolate) +
  mytheme +
  scale_y_continuous(breaks = 1:11) +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.title.x = element_blank())
