source("dataPreparation.R")
batch1 <- ExpeDFMay2018batch1
batch2 <- read.csv("../data/3_recordingTables/Preliminary_April2018_wildmice_Eferrisi_second_RECORDweight .csv")

#############
library(ggplot2)
mytheme <- theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        strip.text = element_text(size=25))
#############


## Earliest weight we had
originalWeight <- batch1[batch1$dpi == -7, ]

table(originalWeight$weight, originalWeight$Mouse_strain)

by(data = originalWeight, 
   INDICES = list(originalWeight$sex, originalWeight$Mouse_strain), 
   FUN = function(x) mean(x$weight))
