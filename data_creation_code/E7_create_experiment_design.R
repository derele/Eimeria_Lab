# Create proper E7 design table

library(tidyverse)
library(Rmisc)
library(httr)
library(RCurl)

# load in infection history and design tables

E7inf <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E7_112018_Eim_infection.history.csv"))
E7a <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E7a_112018_Eim_design.csv"))
E7b <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E7b_112018_Eim_design.csv"))

E7 <- rbind(E7a, E7b) 
#strain is creating trouble as well as some other columns
E7inf$X <- NULL
E7inf$Strain <- NULL


E7 <- join(E7, E7inf)  

write.csv(E7, "/Users/Luke Bednar/Eimeria_Lab/data/Experimental_design/E7_112018_Eim_DESIGN.csv")
