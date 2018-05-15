source("dataPreparation.R")
source("summarySE.R")

###########################################
# Choose the experiment:

##### Expe_001
# March 2017, Francisca's experiment. infection with E64 and Eflab
ExpeDF <- ExpeDF_001

##### Expe_002

##### Expe_003
# April-May 2018, first batch. Parental strains (F0) BUSNA, STRA, SCHUNT, PWD
# Infection with Eferrisi (E64 and E139)
ExpeDF <- ExpeDF_003[ExpeDF_003$dpi %in% 0:11, ]# remove stabilisation period

###########################################
## Weight evolution compared to dpi 0
ggplot(ExpeDF, aes(x = dpi, y = weightRelativeToInfection))+
  geom_line(col = "black", aes(group = EH_ID), alpha = 0.5) +
  geom_point(aes(fill = infection_isolate), 
              size=3, pch = 21, color = "black") +
  mytheme +
  facet_wrap(~Mouse_strain)+
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)" )+
  geom_smooth(aes(col = infection_isolate)) 

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

if (length(levels(ExpeDF$infection_isolate)) <2){
  summary(glm(data = max.loss, weightLossRelativeToInfection ~ Mouse_strain))
} else {
  summary(glm(data = max.loss, weightLossRelativeToInfection ~ Mouse_strain*infection_isolate))
}
## using offsets and real weight instead of relative weight for modeling

###########################################
# oocyst shedding evolution
ggplot(ExpeDF, aes(x = dpi, y = OPG))+
  geom_line(col = "black", aes(group = EH_ID), alpha = 0.5) +
  geom_point(aes(fill = infection_isolate), 
             size=3, pch = 21, color = "black", alpha = 0.78) +
  mytheme +
  facet_wrap(~Mouse_strain)+
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
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.6, fill="grey")

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

if (length(levels(ExpeDF$infection_isolate)) <2){
  summary(glm(data = max.shed, OPG ~ Mouse_strain))
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

if (length(levels(ExpeDF$infection_isolate)) <2){
  summary(glm(data = sum.oocysts, sum.oo ~ Mouse_strain))
} else {
  summary(glm(data = sum.oocysts, sum.oo ~ Mouse_strain*infection_isolate))
}

############################
# highest day of shedding vs highest weight loss?
shedVsLossdpi <- merge(data.frame(EH_ID = all.max.loss$EH_ID, 
                                  dpi_maxLoss = all.max.loss$dpi),
                       data.frame(EH_ID = all.max.shed$EH_ID, 
                                  dpi_maxShed = all.max.shed$dpi), all = T)

shedVsLossdpi <- merge(shedVsLossdpi, 
                       unique(ExpeDF[c("EH_ID", "infection_isolate", "Mouse_strain")]), 
                       all = T)
shedVsLossdpi$diffMaxLossMaxShed <- shedVsLossdpi$dpi_maxLoss - shedVsLossdpi$dpi_maxShed

table(shedVsLossdpi$diffMaxLossMaxShed, shedVsLossdpi$infection_isolate, shedVsLossdpi$Mouse_strain)

ggplot(shedVsLossdpi, aes(x = diffMaxLossMaxShed,
                          y = EH_ID)) +
  geom_jitter(aes(fill = Mouse_strain), 
              size=5, pch = 21, color = "black",
              position = position_jitter(height = .1, width = 0)) +
  facet_grid(.~infection_isolate) +
  mytheme +
  annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           col = "red", alpha = .2) +
  annotate(geom = "rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "purple", alpha = .2) +
  annotate(geom = "text", x = -1.5, y = 1, label = "max weight loss \n before shedding peak")+
  annotate(geom = "text", x = 1.5, y = 1, label = "shedding peak \n before max weight loss")
  