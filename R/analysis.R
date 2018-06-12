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
## Weight evolution compared to dpi 0
ggplot(ExpeDF, aes(x = dpi, y = weightRelativeToInfection, fill = infection_isolate))+
  geom_smooth(aes(col = infection_isolate), alpha = 0.3) +
  geom_line(aes(group = EH_ID), col = "black", alpha = 0.5) +
  geom_point(size=3, pch = 21, color = "black")+
  mytheme +
  facet_grid(~Mouse_strain, scales = "free_y", space = "free") +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)" ) #+
  # geom_hline(yintercept = 80, col = "red")

# Mean + 95%CI
summaryWeight <-summarySE(ExpeDF, measurevar = "weightRelativeToInfection", 
                           groupvars=c("Mouse_strain", "infection_isolate", "dpi"))

ggplot(summaryWeight, aes(x = dpi, y = weightRelativeToInfection))+
  geom_errorbar(aes(ymin = weightRelativeToInfection - ci,
                    ymax = weightRelativeToInfection + ci,
                    col = infection_isolate)) +
  geom_line(aes(group = infection_isolate)) +
  geom_point(aes(fill = infection_isolate), 
             size=3, pch = 21, color = "black") +
  mytheme +
  facet_wrap(~Mouse_strain)+
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "Weight relative to infection (%)" )

# Now, using the CUMULATIVE weight loss
ggplot(ExpeDF, aes(x = dpi, y = cumsumWeightLosRelToInfPercent))+
#  geom_smooth(aes(col = infection_isolate), alpha = 0.3) +
  geom_line(aes(group = EH_ID), alpha = 0.5) +
 # scale_color_gradient(low = "pink", high = "black") +
  geom_point(size=3, pch = 21, color = "black")+
  mytheme +
#  facet_grid(~Mouse_strain, scales = "free_y", space = "free") +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)" ) #+
# geom_hline(yintercept = 80, col = "red")

hist(ExpeDF$ageAtInfection)


ggplot(ExpeDF, aes(x=dpi, y = weighlosstRelativeToInfection)) +
  geom_line(aes(group = EH_ID)) +
  geom_point(aes(col = infection_isolate)) +

ggplot(ExpeDF, aes(dpi, -cumsumWeightLosRelToInfPercent, color = Mouse_strain)) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "pointrange")






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

# Try Kolmogorv-smirnov approach to compare curves
# https://stats.stackexchange.com/questions/129449/package-of-r-for-comparing-graphs-of-daily-activity-of-birds

## 1. get cumulative weight of the mice
library(dplyr)

df <- ExpeDF_004 %>% 
  group_by(dpi, Mouse_strain, infection_isolate) %>%
  summarise(weightLossRelativeToInfection = sum(weightLossRelativeToInfection)) %>%
  mutate(csum = cumsum(weightLossRelativeToInfection))
df <- data.frame(df)
  
df$group <- paste0(df$Mouse_strain, df$infection_isolate)

ggplot(df, aes(x = dpi, y = csum, 
               group = group)) +
  geom_point(aes(fill = infection_isolate), pch =21, size = 5) + 
  geom_line(aes(color = Mouse_strain)) +
  mytheme

# not working!!

# Calculate cumulative weight loss
ExpeDF$

ExpeDF$group <- paste0(ExpeDF$Mouse_strain, ExpeDF$infection_isolate)
library(lme4)
m.0 <- lmer(ExpeDF$weightRelativeToInfection ~ dpi + group + (1 | EH_ID), data=ExpeDF, REML=F)


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
