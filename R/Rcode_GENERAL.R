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

#### Input needed : 
## Oo_Df : a data frame of this format:
# EH_id dpi oocysts.per.g Inf_strain strain
# LM0096   6          0.00      Eflab    PWD
# LM0096   5          0.00      Eflab    PWD

## Weight_Df: a data frame of this format:
# LM_id   dpi Infection_date weight rel.weight mouse_id          name ear sex       born age age_weeks strain line  room.rack cage comment Inf_strain num.inf
# LM0096  dpi0      11_May_17  17.45   98.36528 50035522 PWD158.1.D.2F   0   f 16.04.2016 353        50    PWD  PWD WM1.12/5.1 2164              Eflab       9
# LM0096  dpi1      11_May_17  17.74  100.00000 50035522 PWD158.1.D.2F   0   f 16.04.2016 353        50    PWD  PWD WM1.12/5.1 2164              Eflab       9
# LM0096 dpi10      11_May_17     NA         NA 50035522 PWD158.1.D.2F   0   f 16.04.2016 353        50    PWD  PWD WM1.12/5.1 2164              Eflab       9

###############################
# Input data recent:
# Oo_Df <- read.csv("")
# Weight_Df <- read.csv("")
###############################

## Part 1: how to create a table (cf Emanuel)

## Part 2: how to follow up on the experiment:
# Plot to follow the weight : if < 80% weight, mouse has to be sacrificed
ggplot(data=Alldata,
       aes(x=dpi, y=rel.weight, group=LM_id, color=LM_id)) +
  geom_line()+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=2)+
  geom_text(data=subset(Alldata, dpi == "dpi9" | dpi == "dpi11" ),
            aes(label=LM_id))+
  theme(legend.position="none")

## Part 3: various plots:
## PLOT mice strains:
ggplot(Oo_Df, aes(x=dpi, y=oocysts.per.g, group = strain, col = strain))+
  geom_smooth(aes(fill = strain), alpha = 0.2)+
  ggtitle("Oocyst count along the experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "Loess smoothing + 95% CI")+
  scale_x_continuous(breaks = 0:11) +
  facet_wrap(~Inf_strain)+
  geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = strain), alpha = 0.78) +
  labs(x = "Oocyst per gram", y = "Day post infection") +
  scale_y_continuous(labels = scientific) +
  theme_alice

# Violin plots of the total sum of oocysts collected during 11 days: 
sum.oocysts <- do.call("rbind", by(Oo_Df, Oo_Df$EH_id, function (x){
  x$sum.oo <- sum(x$oocysts.per.g, na.rm=TRUE)
  x
}))

sum.oocysts <- sum.oocysts[!duplicated(sum.oocysts$EH_id),]
# NB: some mice died before!!

ggplot(sum.oocysts, aes(strain, sum.oo)) +
  ggtitle("Sum of oocysts shed during the experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "PWD : HMHZ-eastern-like mice; WSB : HMHZ-west-like mice; WP : HMHZ-hybrids-like mice \n EI64 : wild 'eastern' Eimeria, Eflab : lab 'western-like' Eimeria")+
  geom_violin(color = "black")+
  facet_wrap(~Inf_strain) +
  geom_jitter(width=0.1, size=7, alpha = 0.8, pch = 21, aes(fill = strain)) +
  scale_y_continuous(labels = scientific) +
  theme_alice

# maximum weight lost before death
max.loss <- do.call("rbind", by(Alldata, Alldata$LM_id, function (x){
  m.loss <- which(x$rel.weight==min(x$rel.weight, na.rm=TRUE))
  x[m.loss,]
}))

table(max.loss$dpi, max.loss$strain)

table(max.loss$dpi, max.loss$Inf_strain)

table(max.loss$dpi, max.loss$Inf_strain, max.loss$strain)

tapply(max.loss$rel.weight, max.loss$Inf_strain:max.loss$strain, mean)

ggplot(max.loss, aes(strain, rel.weight, color=Inf_strain)) +
  ggtitle("Relative weight lost during 11 days of experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "PWD : HMHZ-eastern-like mice; WSB : HMHZ-west-like mice; WP : HMHZ-hybrids-like mice \n EI64 : wild 'eastern' Eimeria, Eflab : lab 'western-like' Eimeria")+
  geom_violin(color = "black")+
  facet_wrap(~Inf_strain) +
  geom_jitter(width=0.1, size=7, pch = 21, color = "black", aes(fill = strain), alpha = 0.8) +
  theme_alice

summary(glm(rel.weight~strain + Inf_strain, data=max.loss))

summary(glm(rel.weight~strain, data = max.loss[max.loss$Inf_strain == "EI64",]))

## using offsets and real weight instead of relative weight for modeling

###################
## Plots weight and oocyst counts for a given strain (e.g. E64):
names(sum.oocysts)[1] <- "LM_id"
E64DF <- merge(max.loss, sum.oocysts, by = "LM_id")
E64DF <- E64DF[E64DF$Inf_strain.y == "EI64", ]
E64DF$weight.loss.max <- 100 - E64DF$rel.weight

p1 <- ggplot(E64DF, aes(strain.y, sum.oo)) +
  ggtitle("Sum of oocysts shed during the experiment") +
  geom_violin(color = "black") + 
  geom_jitter(width=0.1, size=7, alpha = 0.8, pch = 21, aes(fill = strain.y)) +
  scale_y_continuous(labels = scientific) +
  labs(x = "", y = "") + 
  theme_alice +
  theme(legend.position="none")

p2 <- ggplot(E64DF, aes(strain.y, weight.loss.max)) +
  ggtitle("Max weight loss reached during the experiment (%)") +
  geom_violin(color = "black") + 
  geom_jitter(width=0.1, size=7, alpha = 0.8, pch = 21, aes(fill = strain.y)) +
  labs(x = "", y = "") + 
  theme_alice

grid.arrange(p1, p2, ncol = 2)
