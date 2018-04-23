library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

###############################
# Input data recent:
ExpeDF <- read.csv("../data_clean/May2017_crossing_infection.csv")
#ExpeDF <- read.csv(file = "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_clean/May2017_crossing_infection.csv")
###############################

## Part 1: West should always be left 
ExpeDF$strain <- factor(ExpeDF$strain, levels = c("WSB", "WP", "PWD"))

###############################
## Part 2: how to follow up on the experiment:
# Plot to follow the weight : if < 80% weight, mouse has to be sacrificed
PlotWeightFollow <- ggplot(data=ExpeDF[ExpeDF$dpi > 0,], # from day 1
                           aes(x=dpi, y=rel.weight)) +
  geom_line(aes(group = EH_id, color = Inf_strain))+
  geom_point(aes(color=Inf_strain), size = 3, alpha = 0.5)+
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme_bw()+
  theme(legend.position=c(.15, .2), legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +
  scale_color_manual(values=c("#999999", "#E69F00"), 
                     name="Infection\nstrains")+
  scale_x_continuous(breaks = 1:11, name = "Day post infection (dpi)" )+ 
  scale_y_continuous(name = "Relative weight compared to dpi 1")

#pdf(file="../figures/May2017_weight_along.pdf", width=12, height=8)
plot(PlotWeightFollow)
#dev.off()

###############################
## Part 3: Plots:
## PLOT mice strains:
PlotOoFollow <- ggplot(ExpeDF, aes(x=dpi, y=oocysts.per.g, group = strain, col = strain))+
  geom_smooth(aes(fill = strain), alpha = 0.2)+
  ggtitle("Oocyst count at different days post infection (dpi)", 
          subtitle = "Loess smoothing + 95% CI")+
  scale_x_continuous(breaks = 0:11) +
  facet_wrap(~Inf_strain)+
  geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = strain), alpha = 0.78) +
  scale_color_manual(values=c("blue", "purple", "red"),  
                     name="Mice genotypes",
                     breaks=c("WSB", "WP", "PWD"),
                     labels=c("M.m.domesticus (West)", "Hybrid", "M.m.musculus (East)"))+
  scale_fill_manual(values=c("blue", "purple", "red"),  
                    name="Mice genotypes",
                    breaks=c("WSB", "WP", "PWD"),
                    labels=c("M.m.domesticus (West)", "Hybrid", "M.m.musculus (East)"))+
  labs(y = "Oocyst per gram", x = "Day post infection (dpi)") +
  scale_y_continuous(labels = scientific) +
  theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        legend.position=c(.2, .6), legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        strip.text = element_text(size=25))

#pdf(file="../figures/May2017_oocyst_along.pdf", width=12, height=8)
plot(PlotOoFollow)
#dev.off()

# Violin plots of the total sum of oocysts collected during 11 days: 
sum.oocysts <- do.call("rbind", by(ExpeDF, ExpeDF$EH_id, function (x){
  x$sum.oo <- sum(x$oocysts.per.g, na.rm=TRUE)
  x
}))

sum.oocysts <- sum.oocysts[!duplicated(sum.oocysts$EH_id),]
# NB: some mice died before!!

sum.oocysts$strain <- factor(sum.oocysts$strain,
                             levels = c("WSB", "WP", "PWD"))

my_title <- expression(paste("Sum of oocysts shed during experimental ", 
                             italic("E. ferrisi"), " infection"))

PlotOoSum <- ggplot(sum.oocysts[sum.oocysts$Inf_strain == "EI64",], #
                    aes(strain, sum.oo/1000000)) +
  ggtitle(my_title) + 
  geom_violin(color = "black")+
  geom_jitter(width=0.1, size=7, alpha = 0.8,
              pch = 21, aes(fill = strain)) +
  labs(y = "Oocyst count (millions)", x = "Mouse strain") +
  scale_color_manual(values=c("blue", "purple", "red"))+
  scale_fill_manual(values=c("blue", "purple", "red"))+
  theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        strip.text = element_text(size=25),
        legend.position="none") +
  scale_x_discrete(labels=c("WSB" = expression(italic("M.m.domesticus")),
                            "WP" = "hybrids",
                            "PWD" = expression(italic("M.m.musculus"))))

#pdf(file="../figures/May2017_oocyst_sum.pdf", width=12, height=8)
plot(PlotOoSum)
#dev.off()

# mean and 95%CI
data = sum.oocysts[sum.oocysts$Inf_strain == "EI64",]
aggregate(data$sum.oo, list(data$strain), mean)
library(Rmisc)
aggregate(data$sum.oo, list(data$strain), CI, ci=0.95)

## Some stats on the oocysts :
sum.oocysts$sum.oo <- round(sum.oocysts$sum.oo, 0)
levels(sum.oocysts$strain) <- c(0, 0.5, 1)
sum.oocysts$strain
sum.oocysts$strain <- as.numeric(as.character(sum.oocysts$strain))

kruskal.test(sum.oo ~ strain, data = sum.oocysts[sum.oocysts$Inf_strain == "EI64",])

# glm.hybrid::glm.hybrid(formula = sum.oo ~ strain, data = sum.oocysts, alpha.along = "strain")

# maximum weight lost before death
max.loss <- do.call("rbind", by(ExpeDF, ExpeDF$EH_id, function (x){
  m.loss <- which(x$rel.weight==min(x$rel.weight, na.rm=TRUE))
  x[m.loss,]
}))

max.loss$strain <- factor(max.loss$strain,
                          levels = c("WSB", "WP", "PWD"))

max.loss$max_loss <- 100 - max.loss$rel.weight

table(max.loss$dpi, max.loss$strain)

table(max.loss$dpi, max.loss$Inf_strain)

table(max.loss$dpi, max.loss$Inf_strain, max.loss$strain)

tapply(max.loss$rel.weight, max.loss$Inf_strain:max.loss$strain, mean)

my_title <- expression(paste("Maximum weight loss during experimental ", 
                             italic("E. ferrisi"), " infection"))

PlotWeightMax <- ggplot(max.loss[max.loss$Inf_strain == "EI64",],
                        aes(strain, max_loss, color=Inf_strain)) +
  ggtitle(my_title) +
  geom_violin(color = "black")+
  geom_jitter(width=0.1, size=7, pch = 21,
              color = "black", aes(fill = strain), alpha = 0.8) +
  labs(y= "Weight loss (%)", x= "Mouse strain")+
  scale_color_manual(values=c("blue", "purple", "red"))+
  scale_fill_manual(values=c("blue", "purple", "red"))+
  theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        strip.text = element_text(size=25),
        legend.position="none") +
  scale_x_discrete(labels=c("WSB" = expression(italic("M.m.domesticus")),
                            "WP" = "hybrids",
                            "PWD" = expression(italic("M.m.musculus"))))

#pdf(file="../figures/May2017_weight_max.pdf", width=12, height=8)
plot(PlotWeightMax)
#dev.off()

summary(glm(rel.weight~strain + Inf_strain, data=max.loss))

summary(glm(rel.weight~strain, data = max.loss[max.loss$Inf_strain == "EI64",]))

kruskal.test(rel.weight ~ strain, data = max.loss[max.loss$Inf_strain == "EI64",])
summary(aov(rel.weight ~ strain, data = max.loss[max.loss$Inf_strain == "EI64",]))

summary(max.loss$rel.weight[max.loss$Inf_strain == "EI64" & 
                              max.loss$strain == "WSB"])

summary(max.loss$rel.weight[max.loss$Inf_strain == "EI64" & 
                              max.loss$strain == "PWD"])

## using offsets and real weight instead of relative weight for modeling