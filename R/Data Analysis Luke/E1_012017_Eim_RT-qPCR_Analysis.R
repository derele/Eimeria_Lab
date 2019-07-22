# RT-qPCR_clean data processing for E1_012017, additional data from CEC
library(RCurl)
library(httr)
library(Rmisc)
library(ggplot2)

#load in data from GitHub
RTqPCRurl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/b4efd8df335199ff9037634c5ba1d909a7d58baa/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv"
RTqPCR <- read.csv(text = getURL(RTqPCRurl))
#RTqPCR <- read.csv(file = "./Documents/Eimeria_Lab/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv")

#very minor cleanup
RTqPCR$X <- NULL
RTqPCR$Cq <- NULL

#graph it out to know what I'm dealing with 
ggplot(data = RTqPCR, aes(x = Sample, y = Cq.Mean, color = Target)) +
  geom_boxplot()

