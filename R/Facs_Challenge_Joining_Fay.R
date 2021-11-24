#load libraries
library(tidyr)
library(dplyr)
library(tidyverse)
library(purrr)

#read the challnge infections
Challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infection_intensity.csv")

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

Challenge_FACS_intensity <- list(Challenge, P4_FACS, E57_FACS, E11_FACS) %>%
  reduce(join_my_tables)


write.csv(Challenge_FACS_intensity, "~/Documents/GitHub/Eimeria_Lab/data_products/Challenge_FACS.csv", row.names = FALSE)
