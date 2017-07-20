library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

#### Input needed : 
## Oo_Df
## Weight_Df

# To keep all plot with the same theme:
theme_alice <- theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        text = element_text(size = 20))

## PLOT mice strains:
ggplot(Oo_Df, aes(x=dpi, y=oocysts.per.g, group = strain, col = strain))+
  geom_smooth(aes(fill = strain), alpha = 0.2)+
  ggtitle("Oocyst count along the experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "Loess smoothing + 95% CI")+
  scale_x_continuous(breaks = 0:11) +
  facet_wrap(~Inf_strain)+
  geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = strain), alpha = 0.78) +
  labs(x = "Oocyst per gram", y = "Day post infection") +
  scale_fill_discrete(labels = c("Eastern mice", "Hybrids", "Western mice")) +
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