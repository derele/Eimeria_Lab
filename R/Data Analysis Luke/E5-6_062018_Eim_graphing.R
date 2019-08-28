library(RCurl)
library(httr)
library(Rmisc)
library(magrittr)
library(ggplot2)

# load in data
E6_completeURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E6_062018_Eim_complete.csv"
E6_complete <- read.csv(text = getURL(E6_completeURL), sep = ";", dec = ",")

#convert number columns
cols = c(4, 5, 6, 7, 8, 10, 12, 13, 14, 15)
E6_complete[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

#graph weight
ggplot(E6_complete, aes(dpi, weight, color=Eimeria_species)) +
  # geom_jitter(width=0.2) +
  geom_smooth(se=FALSE) +
  # scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) +
  #                    labels=c("0" ,"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")) +
  facet_wrap(~HybridStatus, scales="free_y", nrow=2) +
  scale_colour_brewer("infection\nspecies", palette = "Dark2") +
  scale_y_continuous("weight (g)") +
  theme_bw()

# graph weight loss
ggplot(E6_complete, aes(dpi, weightloss, color=Eimeria_species)) +
  # geom_jitter(width=0.2) +
  geom_smooth(se=FALSE) +
  # scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) +
  #                    labels=c("0" ,"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")) +
  facet_wrap(~HybridStatus, nrow=2) +
  scale_colour_brewer("infection\nspecies", palette = "Dark2") +
  scale_y_continuous("weight (g)") +
  theme_bw()

# graph oocysts
ggplot(E6_complete, aes(dpi, OPG, color=Eimeria_species)) +
  # geom_jitter(width=0.2) +
  geom_smooth(se=FALSE) +
  # scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) +
  #                    labels=c("0" ,"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")) +
  facet_wrap(~HybridStatus, nrow=2) +
  scale_colour_brewer("infection\nspecies", palette = "Dark2") +
  scale_y_continuous("weight (g)") +
  theme_bw()


Cytokines <- ggarrange(CytokinesSP, CytokinesCE, nrow=2, common.legend = TRUE,
                       legend="right", labels = "AUTO")
dev.off()