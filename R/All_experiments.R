library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

# To keep all plot with the same theme:
theme_alice <- theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        text = element_text(size = 20))

###############################
# Input data recent:
ExpeDF <- read.csv(file = "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_clean/May2017_crossing_infection.csv")
###############################


## Part 1: West should always be left 
ExpeDF$strain <- factor(ExpeDF$strain, levels = c("WSB", "WP", "PWD"))


###############################
## Part 2: how to follow up on the experiment:
# Plot to follow the weight : if < 80% weight, mouse has to be sacrificed
PlotWeightFollow <- ggplot(data=ExpeDF,
                           aes(x=dpi, y=rel.weight,
                               group=EH_id, color=EH_id)) +
    geom_line()+
    geom_point()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    geom_hline(yintercept=80, linetype="dashed", color = "red", size=2)+
    geom_text(data=subset(ExpeDF, dpi == "9" | dpi == "11" ),
              aes(label=EH_id))+
    theme(legend.position="none")

pdf(file="figures/May2017_weight_along.pdf", width=12, height=8)
plot(PlotWeightFollow)
dev.off()


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
    scale_color_manual(values=c("blue", "purple", "red"))+
    scale_fill_manual(values=c("blue", "purple", "red"))+
  labs(y = "Oocyst per gram", x = "Day post infection (dpi)") +
  scale_y_continuous(labels = scientific) +
  theme_alice

pdf(file="./figures/May2017_oocyst_along.pdf", width=12, height=8)
plot(PlotOoFollow)
dev.off()

# Violin plots of the total sum of oocysts collected during 11 days: 
sum.oocysts <- do.call("rbind", by(ExpeDF, ExpeDF$EH_id, function (x){
  x$sum.oo <- sum(x$oocysts.per.g, na.rm=TRUE)
  x
}))

sum.oocysts <- sum.oocysts[!duplicated(sum.oocysts$EH_id),]
# NB: some mice died before!!

sum.oocysts$strain <- factor(sum.oocysts$strain,
                             levels = c("WSB", "WP", "PWD"))

PlotOoSum <- ggplot(sum.oocysts, aes(strain, sum.oo)) +
    ggtitle("Sum of oocysts shed during the experiment") + 
    geom_violin(color = "black")+
    facet_wrap(~Inf_strain) +
    geom_jitter(width=0.1, size=7, alpha = 0.8,
                pch = 21, aes(fill = strain)) +
    scale_color_manual(values=c("blue", "purple", "red"))+
    scale_fill_manual(values=c("blue", "purple", "red"))+
    labs(y = "Total number of oocyst shed", x = "Mouse strain")
    scale_y_continuous(labels = scientific) +
    theme_alice

pdf(file="./figures/May2017_oocyst_sum.pdf", width=12, height=8)
plot(PlotOoSum)
dev.off()

# maximum weight lost before death
max.loss <- do.call("rbind", by(ExpeDF, ExpeDF$EH_id, function (x){
  m.loss <- which(x$rel.weight==min(x$rel.weight, na.rm=TRUE))
  x[m.loss,]
}))

max.loss$strain <- factor(max.loss$strain,
                          levels = c("WSB", "WP", "PWD"))


table(max.loss$dpi, max.loss$strain)

table(max.loss$dpi, max.loss$Inf_strain)

table(max.loss$dpi, max.loss$Inf_strain, max.loss$strain)

tapply(max.loss$rel.weight, max.loss$Inf_strain:max.loss$strain, mean)

PlotWeightMax <- ggplot(max.loss, aes(strain, rel.weight, color=Inf_strain)) +
  ggtitle("Relative weight retained") + 
  geom_violin(color = "black")+
  facet_wrap(~Inf_strain) +
    geom_jitter(width=0.1, size=7, pch = 21,
                color = "black", aes(fill = strain), alpha = 0.8) +
    labs(y= "Minimum weigth retained relative to weight at infection",
         x= "Mouse strain")+
    scale_color_manual(values=c("blue", "purple", "red"))+
    scale_fill_manual(values=c("blue", "purple", "red"))+
  theme_alice

pdf(file="./figures/May2017_weight_max.pdf", width=12, height=8)
plot(PlotWeightMax)
dev.off()

summary(glm(rel.weight~strain + Inf_strain, data=max.loss))

summary(glm(rel.weight~strain, data = max.loss[max.loss$Inf_strain == "EI64",]))

## using offsets and real weight instead of relative weight for modeling
