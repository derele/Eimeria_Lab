# Part of Eimeria_Lab cleanup
library(Rmisc)
library(httr)
library(RCurl)
library(tidyverse)

# load in all E5 design tables
E51a <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E5_1a_062018_Eim_DESIGN.csv"))
E51b <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E5_1b_062018_Eim_DESIGN.csv"))
E52a <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E5_2a_062018_Eim_DESIGN.csv"))
E52b <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E5_2b_062018_Eim_DESIGN.csv"))

E5 <- rbind(E51a, E51b)
E5 <- rbind(E5, E52a)
E5 <- rbind(E5, E52b)
colnames(E5)[32] <- "EH_ID"

write.csv(E5, "/Users/Luke Bednar/Eimeria_Lab/data/Experimental_design/E_062018_Eim_DESIGN.csv")
