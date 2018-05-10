sunFlo <- read.csv("../data/3_recordingTables/Preliminary_April2018_wildmice_Eferrisi_second_RECORDweight .csv")

library(ggplot2)
mytheme <- theme_bw()+
  theme(plot.title = element_text(size=22), plot.subtitle = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        strip.text = element_text(size=25))

sunFlo$prelim_labels <- as.character(sunFlo$prelim_labels)

ggplot(sunFlo) +
  geom_point(aes(x = date, y = weight)) +
  geom_line(aes(x = date, y = weight, group = prelim_labels)) +
  mytheme

begin <- sunFlo[sunFlo$date == "07/05/18", c("prelim_labels", "weight")]
plus3Days <- sunFlo[sunFlo$date == "10/05/18", c("prelim_labels", "weight")]

names(begin)[2] <- "weight_begin"      
names(plus3Days)[2] <- "weight_plus3Days"      

tot <- merge(begin, plus3Days, all = T)

tot$weightDiff <- tot$weight_plus3Days - tot$weight_begin

ggplot(tot) +
  geom_violin(aes(x = "Difference of weight", y = weightDiff))+
  mytheme
