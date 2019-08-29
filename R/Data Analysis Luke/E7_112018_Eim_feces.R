# setwd("Repositories/Eimeria_Lab/R/Data Analysis Luke/")
library(Rmisc)
library(httr)
library(RCurl)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(data.table)

E7a_design <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E7a_112018_Eim_design.csv"
E7b_design <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E7b_112018_Eim_design.csv"
E7aF <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7a_112018_Eim_feces.csv"
E7bF <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7b_112018_Eim_feces.csv"

E7a_design <- read.csv(text=getURL(E7a_design))
E7b_design <- read.csv(text=getURL(E7b_design))
E7aF <- read.csv(text=getURL(E7aF))
E7bF <- read.csv(text=getURL(E7bF))

#the columns we want to keep
col2keep <- c("Strain", "HybridStatus", "EH_ID")

E7a_design <- E7a_design[col2keep]
E7b_design <- E7b_design[col2keep]

# rename EH_id to EH_ID#
names(E7a_design)[names(E7a_design) == "EH_id"] <- "EH_ID"
names(E7b_design)[names(E7b_design) == "EH_id"] <- "EH_ID"

#add infection history
history_a <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E7a_112018_Eim_infection.history.csv"
history_a <- read.csv(text=getURL(history_a))
history_b <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E7b_112018_Eim_infection.history.csv"
history_b <- read.csv(text=getURL(history_b))
inf.history <- rbind(history_a, history_b) 

# let's make one big fat E7 design table (E88 = 31 entries, E64 = 38 entries)
E7_design <- rbind(E7a_design, E7b_design)
E7_design <- merge(E7_design, inf.history, by = "EH_ID")

# keep the batch information
E7aF$batch <- "october2018"
E7bF$batch <- "december2018"

# Make one big fat table Expe 7
E7_record <- rbind(E7aF, E7bF)

# Merge all, #
E7 <- merge(E7_design, E7_record)

#export HU
#write.csv(E7, "../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv", quote = FALSE)

#export IZW
#write.csv(Exp007, "../luke/Documents//Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv", quote = FALSE", quote = FALSE)

#export home (Win)
write.csv(E7, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv", quote = FALSE)

#include combined infection history
E7$infHistory <- E7$primary:E7$challenge

#distribution check
E7 <- E7[-c(447, 479),] #remove NAs
hist(E7$Wchange, col = grey)
summary(E7)
boxplot(Wchange ~ infHistory, data = E7)
boxplot(Wchange ~ HybridStatus, data = E7)
boxplot(Wchange ~ dpi, data = E7)
histogram(~Wchange | factor(infHistory), data = E7)

# #add SDs
# E7 <- data.frame(E7 %>% group_by(Wchange, EH_ID, primary, challenge, infHistory) %>% 
#                        summarize(SD = sd(Wchange),
#                                  Wchange.mean = mean(Wchange)))
# 
# #attach means of HybridStatus, infHistory, etc with data.table
# setDT(E7)[, WmeanHY := mean(Wchange), by = HybridStatus]
# setDT(E7)[, WmeanIH := mean(Wchange), by = infHistory]
# setDT(E7)[, WmeanD := mean(Wchange), by = dpi]


#plot weight and HS by infection history
ggplot(E7, aes(dpi, weight, color=HybridStatus)) +
  # geom_jitter(width=0.2) +
  geom_smooth(se=F) +
  # scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) +
  #                    labels=c("0" ,"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")) +
  facet_wrap(~infHistory, scales="free_y", nrow=2) +
  scale_colour_brewer("Hybrid\nStatus", palette = "Dark2") +
  scale_y_continuous("weight (g)") +
  theme_bw()

# plot weight and  by infection history
ggplot(E7, aes(dpi, weight, color=infHistory)) +
  # geom_jitter(width=0.2) +
  geom_smooth(se=FALSE) +
  # scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) +
  #                    labels=c("0" ,"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")) +
  facet_wrap(~HybridStatus, scales="free_y", nrow=2) +
  scale_colour_brewer("Infection\nHistory", palette = "Dark2") +
  scale_y_continuous("weight (g)") +
  theme_bw()
