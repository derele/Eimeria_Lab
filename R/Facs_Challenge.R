#load libraries
library(tidyr)
library(dplyr)
library(tidyverse)
library(purrr)

#load the challenge infection table 
Challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_p4_p3.csv")

#load the infection intensity data
P4_infection <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_CEWE_qPCR.csv") 
E57_infection <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_CEWE_qPCR.csv") %>%
  select(!X)
P3_infection <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_CEWE_qPCR.csv") %>%
  select(!X)

#load the facs data for P4, E57, E11
P4_FACS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_FACS.csv") %>%
  select(!X)

E57_FACS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E7_112018_Eim_FACS.csv") %>%
  select(!X)

E11_FACS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/E11_012021_Eim_FACS.csv") %>%
  select(!X) %>%
  select(!X.1)

#write a function to change LM_ to LM
#write a function to replace LM_ with LM
replace_LM <- function (x) {
  str_replace_all(x$EH_ID, "LM_", "LM") 
}

#replace LMs
P4_FACS$EH_ID <- replace_LM(P4_FACS)
E57_FACS$EH_ID <- replace_LM(E57_FACS)
E11_FACS$EH_ID <- replace_LM(E11_FACS)
P4_infection$EH_ID <- replace_LM(P4_infection)
E57_infection$EH_ID <- replace_LM(E57_infection)
P3_infection$EH_ID <- replace_LM(P3_infection)

#Eim_MC in some infection intensity column is a character and in others is logical
#change it so it is the same
P4_infection$Eim_MC <- as.character(P4_infection$Eim_MC)

#write a function to create a list of common column names
names_common <- function(x, y) {
  intersect(colnames(x), colnames(y))
}

#write a function to left join 
join_my_tables <- function (x, y) {
  x %>%
    left_join(y, by = c(names_common(x, y)))
}

#try to use the reduce to bind the tables together
#reduce works like apply, but you need to have a list for it to work on
Challenge_intensity <- list(Challenge, P4_infection, E57_infection, P3_infection) %>%
  reduce(join_my_tables)

Challenge_FACS_intensity <- list(Challenge_intensity, P4_FACS, E57_FACS, E11_FACS) %>%
  reduce(join_my_tables)
colnames(Challenge_FACS_intensity)

write.csv(Challenge_FACS_intensity, "~/Documents/GitHub/Eimeria_Lab/data_products/FACS_intensity.csv", row.names = FALSE)
