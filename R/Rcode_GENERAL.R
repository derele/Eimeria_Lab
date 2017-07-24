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
ExpeDF <- read.csv(file = "../data_clean/First_crossing_infection.csv")
###############################
## Part 1: how to create a table (cf Emanuel)


###############################
## Part 2: how to follow up on the experiment:
# Plot to follow the weight : if < 80% weight, mouse has to be sacrificed
ggplot(data=ExpeDF,
       aes(x=dpi, y=rel.weight, group=EH_id, color=EH_id)) +
  geom_line()+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=2)+
  geom_text(data=subset(ExpeDF, dpi == "9" | dpi == "11" ),
            aes(label=EH_id))+
  theme(legend.position="none")

###############################
## Part 3: Plots:
## PLOT mice strains:
PlotOo <- ggplot(ExpeDF, aes(x=dpi, y=oocysts.per.g, group = strain, col = strain))+
  geom_smooth(aes(fill = strain), alpha = 0.2)+
  ggtitle("Oocyst count along the experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "Loess smoothing + 95% CI")+
  scale_x_continuous(breaks = 0:11) +
  facet_wrap(~Inf_strain)+
  geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = strain), alpha = 0.78) +
  labs(x = "Oocyst per gram", y = "Day post infection") +
  scale_y_continuous(labels = scientific) +
  theme_alice

#pdf(filename="../figures/FirstExpe_oocist_along.pdf")
#PlotOo
#dev.off()

# Violin plots of the total sum of oocysts collected during 11 days: 
sum.oocysts <- do.call("rbind", by(ExpeDF, ExpeDF$EH_id, function (x){
  x$sum.oo <- sum(x$oocysts.per.g, na.rm=TRUE)
  x
}))

sum.oocysts <- sum.oocysts[!duplicated(sum.oocysts$EH_id),]
# NB: some mice died before!!

plot_oo <- ggplot(sum.oocysts, aes(strain, sum.oo)) +
  ggtitle("Sum of oocysts shed during the experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "PWD : HMHZ-eastern-like mice; WSB : HMHZ-west-like mice; WP : HMHZ-hybrids-like mice \n EI64 : wild 'eastern' Eimeria, Eflab : lab 'western-like' Eimeria")+
  geom_violin(color = "black")+
  facet_wrap(~Inf_strain) +
  geom_jitter(width=0.1, size=7, alpha = 0.8, pch = 21, aes(fill = strain)) +
  scale_y_continuous(labels = scientific) +
  theme_alice

#pdf(filename="../figures/FirstExpe_violin_oocist.pdf")
#plot_oo
#dev.off()

# maximum weight lost before death
max.loss <- do.call("rbind", by(ExpeDF, ExpeDF$EH_id, function (x){
  m.loss <- which(x$rel.weight==min(x$rel.weight, na.rm=TRUE))
  x[m.loss,]
}))

table(max.loss$dpi, max.loss$strain)

table(max.loss$dpi, max.loss$Inf_strain)

table(max.loss$dpi, max.loss$Inf_strain, max.loss$strain)

tapply(max.loss$rel.weight, max.loss$Inf_strain:max.loss$strain, mean)

plot_weight <- ggplot(max.loss, aes(strain, rel.weight, color=Inf_strain)) +
  ggtitle("Relative weight lost during 11 days of experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "PWD : HMHZ-eastern-like mice; WSB : HMHZ-west-like mice; WP : HMHZ-hybrids-like mice \n EI64 : wild 'eastern' Eimeria, Eflab : lab 'western-like' Eimeria")+
  geom_violin(color = "black")+
  facet_wrap(~Inf_strain) +
  geom_jitter(width=0.1, size=7, pch = 21, color = "black", aes(fill = strain), alpha = 0.8) +
  theme_alice

#pdf(filename="../figures/FirstExpe_violin_weight.pdf")
#plot_weight
#dev.off()

summary(glm(rel.weight~strain + Inf_strain, data=max.loss))

summary(glm(rel.weight~strain, data = max.loss[max.loss$Inf_strain == "EI64",]))

## using offsets and real weight instead of relative weight for modeling