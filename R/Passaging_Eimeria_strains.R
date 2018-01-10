A <- read.csv("/home/alice/Schreibtisch/git_projects/Eimeria_Lab/data_clean/Passaging_Eimeria_strains.csv")

A$DPI <- factor(A$dpi, levels = levels(A$dpi)[c(4,5,6,7,8,9,1,2,3)])

library(ggplot2)
library(scales)


myplot <- ggplot(A, aes(x=DPI, y=OPG, group = Eimeria_strain, col = Eimeria_strain))+
  geom_smooth(aes(fill = Eimeria_strain), alpha = 0.2)+
  ggtitle("Oocyst count at different days post infection (dpi)", 
          subtitle = "Loess smoothing + 95% CI") +
  geom_jitter(width=0.1, size=5, 
              aes(fill = Eimeria_strain, shape = factor(cage_number), color = Eimeria_strain)) +
  labs(y = "Oocyst per gram", x = "Day post infection (dpi)") +
  scale_y_continuous(labels = scientific) +
  theme_classic()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        legend.position=c(.2, .6), legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.key = element_rect(size = 1),
        legend.key.size = unit(1.5, 'lines'),
        strip.text = element_text(size=25))

pdf(file="../figures/Dec2017_passaging_oocyst_along.pdf", width=12, height=8)
plot(myplot)
dev.off()
