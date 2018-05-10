library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

mytheme <- theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        strip.text = element_text(size=25))

########################### Data loading #######################################
# ExpeDF need at least the following: 
# "OPG", "EH_ID", "dpi", "weight", "EH_id", "infection_isolate"                   

source("dataPreparation.R")

ExpeDF <- ExpeDFMay2018batch1

## Weight evolution
plainWeight <- ggplot(ExpeDF, 
                      aes(x = dpi, y = weight))+
  geom_line(aes(col = EH_ID, group = EH_ID)) +
  geom_point(size=3, pch = 21, color = "black", aes(fill = EH_ID), alpha = 0.78) +
  mytheme +
  facet_wrap(~Mouse_strain, scales = "free_y")+
  scale_x_continuous(breaks = -7:11, name = "Day post infection (dpi)" )
plot(plainWeight)

# smooth by Mouse_strain
smoothplainWeight <- ggplot(ExpeDF, 
                      aes(x = dpi, y = weight))+
  geom_point(size=3, pch = 21, color = "black", aes(fill = Mouse_strain), alpha = 0.78) +
  geom_smooth(aes(col = Mouse_strain)) +
  mytheme +
  scale_x_continuous(breaks = -7:11, name = "Day post infection (dpi)" )
plot(smoothplainWeight)

# compared to dpi 0
plotWeight <- ggplot(ExpeDF, 
                     aes(x = dpi, y = weightRelativeToInfection))+
  geom_line(aes(col = EH_ID, group = EH_ID)) +
  geom_point(size=3, pch = 21, color = "black", aes(fill = EH_ID), alpha = 0.78) +
  mytheme +
  scale_x_continuous(breaks = -7:11, name = "Day post infection (dpi)" )
plot(plotWeight)

# If we have the weight at arrival
if (!is.na(ExpeDF$weightAtAnthelminthicTrt)){
  plotWeight2 <- ggplot(ExpeDFMay2018batch1, 
                        aes(x = dpi, y = weightRelativeToAnthelmTrtDay))+
#    geom_line(aes(col = EH_ID, group = EH_ID)) +
    geom_point(size=3, aes(col = Mouse_strain, pch = infection_isolate), alpha = 0.78) +
    mytheme +
    geom_smooth(aes(col = Mouse_strain, group = Mouse_strain, fill = Mouse_strain)) 
    plot(plotWeight2)
}


## PLOT mice strains:
# PlotOoFollow <- ggplot(ExpeDF, aes(x=dpi, y=OPG, group = Mouse_strain, col = Mouse_strain))+
#   geom_smooth(aes(fill = Mouse_strain), alpha = 0.2)+
#   ggtitle("Oocyst count at different days post infection (dpi)", 
#           subtitle = "Loess smoothing + 95% CI")+
#   scale_x_continuous(breaks = 0:11) +
#   facet_wrap(~infection_isolate,  scales="free_y") +
#   geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = Mouse_strain), alpha = 0.78) +
#   scale_y_continuous(labels = scientific) +
#   mytheme 
# plot(PlotOoFollow)
# 
# PlotOoTotFollow <- ggplot(ExpeDF, aes(x=as.factor(dpi), y=oocystsTotal, group = Mouse_strain, col = Mouse_strain))+
#   geom_smooth(aes(fill = Mouse_strain), alpha = 0.2)+
#   ggtitle("Oocyst count at different days post infection (dpi)", 
#           subtitle = "Loess smoothing + 95% CI")+
#   # scale_x_continuous(breaks = 0:11) +
#   facet_wrap(~infection_isolate,  scales="free_y") +
#   geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = Mouse_strain), alpha = 0.78) +
#   scale_y_continuous(labels = scientific) +
#   mytheme 
# plot(PlotOoTotFollow)



library(reshape2)

df <- na.omit(melt(ExpeDF, id.vars = c("EH_ID", "dpi"), measure.vars = c("weight")))

df <- df[df$variable == "weight",c("EH_ID", "dpi", "value")]

df1 <- df[df$dpi =="-7", ]
df2 <- df[df$dpi =="-1", ]
df <- merge(df1, df2, by = "EH_ID", all = T)

df$diff <- (df$value.x - df$value.y) / df$value.x * 100
mean(df$diff)


# Mean + 95%CI
source("summarySE.R")
summaryOocysts <-summarySE(ExpeDF, measurevar="oocystsTotal", 
                           groupvars=c("Mouse_strain","infection_isolate", "dpi"))

ggplot(summaryOocysts, aes(x=as.factor(dpi), y=oocystsTotal, group = Mouse_strain, col = Mouse_strain))+
  geom_errorbar(aes(ymin=oocystsTotal-ci, ymax=oocystsTotal+ci), colour="black", width=.1, 
                position=position_dodge(0.1)) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.1), size=3) +
  facet_wrap(~infection_isolate, scales="free_y") +
  scale_x_continuous(breaks = 0:11) +
  mytheme 

# Weight

PlotWeight <- ggplot(ExpeDF, aes(x=dpi, y=weightRelativeToInfection))+
  geom_smooth(aes(fill = Mouse_strain, col = infection_isolate), alpha = 0.2) +
  #scale_x_continuous(breaks = 0:11) +
  facet_wrap(~infection_isolate) +
  geom_point(pch = 21, color = "black", size = 3) +
  mytheme 
plot(PlotWeight)
# 
# ###############################
# # Input data recent:
# #ExpeDF <- read.csv("../data_clean/May2017_crossing_infection.csv")
# #ExpeDF <- read.csv(file = "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_clean/May2017_crossing_infection.csv")
# ###############################
# ## Part 1: West should always be left 
# # ExpeDF$strain <- factor(ExpeDF$strain, levels = c("WSB", "WP", "PWD"))
# 
# # Plot to follow the weight : if < 80% weight, mouse has to be sacrificed
# PlotWeightFollow <- ggplot(data=ExpeDF[ExpeDF$dpi > 0,], # from day 1
#                            aes(x=dpi, y=weightRelativeToInfection)) +
#   geom_line(aes(group = EH_ID, color = infection_isolate))+
#   geom_point(aes(color = infection_isolate), size = 3, alpha = 0.5)+
#   theme_bw()+
#   theme(legend.position=c(.15, .1), legend.title = element_text(size=20, face="bold"),
#         legend.text = element_text(size = 20),
#         axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold")) +
#   scale_x_continuous(breaks = 1:11, name = "Day post infection (dpi)" )
# 
# plot(PlotWeightFollow)
# 
# ## PLOT mice strains:
# PlotOoFollow <- ggplot(ExpeDF, aes(x=dpi, y=log10(OPG +1), group = Mouse_strain, col = Mouse_strain))+
#   geom_smooth(aes(fill = Mouse_strain), alpha = 0.2)+
#   ggtitle("Oocyst count at different days post infection (dpi)", 
#           subtitle = "Loess smoothing + 95% CI")+
#   scale_x_continuous(breaks = 0:11) +
#   facet_wrap(~infection_isolate) +
#   geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = Mouse_strain), alpha = 0.78) +
#   # scale_color_manual(values=c("blue", "purple", "red"),
#   #                    name="Mice genotypes",
#   #                    breaks=c("WSB", "WP", "PWD"),
#   #                    labels=c("M.m.domesticus (West)", "Hybrid", "M.m.musculus (East)"))+
#   # scale_fill_manual(values=c("blue", "purple", "red"),
#   #                   name="Mice genotypes",
#   #                   breaks=c("WSB", "WP", "PWD"),
#   #                   labels=c("M.m.domesticus (West)", "Hybrid", "M.m.musculus (East)"))+
#   # labs(y = "Oocyst per gram", x = "Day post infection (dpi)") +
#   scale_y_continuous(labels = scientific) +
#   mytheme 
# plot(PlotOoFollow)
# 
# PlotOoTotFollow <- ggplot(ExpeDF, aes(x=dpi, y=log10(oocystsTotal +1), group = Mouse_strain, col = Mouse_strain))+
#   geom_smooth(aes(fill = Mouse_strain), alpha = 0.2)+
#   ggtitle("Oocyst count at different days post infection (dpi)", 
#           subtitle = "Loess smoothing + 95% CI")+
#   scale_x_continuous(breaks = 0:11) +
#   facet_wrap(~infection_isolate) +
#   geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = Mouse_strain), alpha = 0.78) +
#   # scale_color_manual(values=c("blue", "purple", "red"),
#   #                    name="Mice genotypes",
#   #                    breaks=c("WSB", "WP", "PWD"),
#   #                    labels=c("M.m.domesticus (West)", "Hybrid", "M.m.musculus (East)"))+
#   # scale_fill_manual(values=c("blue", "purple", "red"),
#   #                   name="Mice genotypes",
#   #                   breaks=c("WSB", "WP", "PWD"),
#   #                   labels=c("M.m.domesticus (West)", "Hybrid", "M.m.musculus (East)"))+
#   # labs(y = "Oocyst per gram", x = "Day post infection (dpi)") +
#   scale_y_continuous(labels = scientific) +
#   mytheme 
# plot(PlotOoTotFollow)
# 
# 
# 
# ### To clean
# # Violin plots of the total sum of oocysts collected during 11 days: 
# sum.oocysts <- do.call("rbind", by(ExpeDF, ExpeDF$EH_ID, function (x){
#   x$sum.oo <- sum(x$OPG, na.rm=TRUE)
#   x
# }))
# 
# sum.oocysts <- sum.oocysts[!duplicated(sum.oocysts$EH_ID),]
# # NB: some mice died before!!
# 
# # sum.oocysts$strain <- factor(sum.oocysts$strain,
# #                              levels = c("WSB", "WP", "PWD"))
# # 
# # my_title <- expression(paste("Sum of oocysts shed during experimental ", 
# #                              italic("E. ferrisi"), " infection"))
# 
# PlotOoSum <- ggplot(sum.oocysts, #
#                     aes(Mouse_strain, sum.oo)) +
#   ggtitle(my_title) + 
#   geom_violin(color = "black")+
#   geom_jitter(width=0.1, size=7, alpha = 0.8,
#               pch = 21, aes(fill = Mouse_strain)) +
#   labs(y = "Oocyst count", x = "Mouse strain") +
#   facet_grid(.~infection_isolate) +
#   # scale_color_manual(values=c("blue", "purple", "red"))+
#   # scale_fill_manual(values=c("blue", "purple", "red"))+
#   theme_bw()
#   # theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
#   #       axis.text=element_text(size=20),
#   #       axis.title=element_text(size=20,face="bold"),
#   #       strip.text = element_text(size=25),
#   #       legend.position="none") +
#   # scale_x_discrete(labels=c("WSB" = expression(italic("M.m.domesticus")),
#   #                           "WP" = "hybrids",
#   #                           "PWD" = expression(italic("M.m.musculus"))))
# 
# #pdf(file="../figures/May2017_oocyst_sum.pdf", width=12, height=8)
# plot(PlotOoSum)
# #dev.off()
# 
# # mean and 95%CI
# library(Rmisc)
# aggregate(sum.oocysts$sum.oo, 
#           list(sum.oocysts$infection_isolate, sum.oocysts$Mouse_strain), 
#           CI, ci=0.95)
# 
# ## Some stats on the oocysts :
# sum.oocysts$sum.oo <- round(sum.oocysts$sum.oo, 0)
# levels(sum.oocysts$strain) <- c(0, 0.5, 1)
# sum.oocysts$strain
# sum.oocysts$strain <- as.numeric(as.character(sum.oocysts$strain))
# 
# kruskal.test(sum.oo ~ strain, data = sum.oocysts[sum.oocysts$Inf_strain == "EI64",])
# 
# # glm.hybrid::glm.hybrid(formula = sum.oo ~ strain, data = sum.oocysts, alpha.along = "strain")
# 
# # maximum weight lost before death
# max.loss <- do.call("rbind", by(ExpeDF, ExpeDF$EH_id, function (x){
#   m.loss <- which(x$rel.weight==min(x$rel.weight, na.rm=TRUE))
#   x[m.loss,]
# }))
# 
# max.loss$strain <- factor(max.loss$strain,
#                           levels = c("WSB", "WP", "PWD"))
# 
# max.loss$max_loss <- 100 - max.loss$rel.weight
# 
# table(max.loss$dpi, max.loss$strain)
# 
# table(max.loss$dpi, max.loss$Inf_strain)
# 
# table(max.loss$dpi, max.loss$Inf_strain, max.loss$strain)
# 
# tapply(max.loss$rel.weight, max.loss$Inf_strain:max.loss$strain, mean)
# 
# my_title <- expression(paste("Maximum weight loss during experimental ", 
#                              italic("E. ferrisi"), " infection"))
# 
# PlotWeightMax <- ggplot(max.loss[max.loss$Inf_strain == "EI64",],
#                         aes(strain, max_loss, color=Inf_strain)) +
#   ggtitle(my_title) +
#   geom_violin(color = "black")+
#   geom_jitter(width=0.1, size=7, pch = 21,
#               color = "black", aes(fill = strain), alpha = 0.8) +
#   labs(y= "Weight loss (%)", x= "Mouse strain")+
#   scale_color_manual(values=c("blue", "purple", "red"))+
#   scale_fill_manual(values=c("blue", "purple", "red"))+
#   theme_bw()+
#   theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
#         axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold"),
#         strip.text = element_text(size=25),
#         legend.position="none") +
#   scale_x_discrete(labels=c("WSB" = expression(italic("M.m.domesticus")),
#                             "WP" = "hybrids",
#                             "PWD" = expression(italic("M.m.musculus"))))
# 
# #pdf(file="../figures/May2017_weight_max.pdf", width=12, height=8)
# plot(PlotWeightMax)
# #dev.off()
# 
# summary(glm(rel.weight~strain + Inf_strain, data=max.loss))
# 
# summary(glm(rel.weight~strain, data = max.loss[max.loss$Inf_strain == "EI64",]))
# 
# kruskal.test(rel.weight ~ strain, data = max.loss[max.loss$Inf_strain == "EI64",])
# summary(aov(rel.weight ~ strain, data = max.loss[max.loss$Inf_strain == "EI64",]))
# 
# summary(max.loss$rel.weight[max.loss$Inf_strain == "EI64" & 
#                               max.loss$strain == "WSB"])
# 
# summary(max.loss$rel.weight[max.loss$Inf_strain == "EI64" & 
#                               max.loss$strain == "PWD"])
# 
# ## using offsets and real weight instead of relative weight for modeling