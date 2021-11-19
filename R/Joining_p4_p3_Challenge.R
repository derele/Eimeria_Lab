#load libraries
library(tidyr)
library(ggplot2)
library(dplyr)
library(stringr) #to manipulate strings
library(magrittr)
library(janitor)

#read challenge infections
Challenge_Infections <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")

#read experimental design P4
P4_Design <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P4_082020_Eim_DESIGN.csv")
P4_Design$EH_ID <- str_replace_all(P4_Design$EH_ID, "LM_", "LM") #replace string values in EH_ID Column so that table matches the challenge infections table

#read Experimental record P4
P4_record <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_record.csv") %>%
  select(!X)

P4_record$EH_ID <- str_replace_all(P4_record$EH_ID, "LM_", "LM") #replace string values in EH_ID Column so that table matches the challenge infections table

#write a function to create a list of common column names
names_common <- function(x, y) {
  intersect(colnames(x), colnames(y))
}

#write a function to left join 
join_my_tables <- function (x, y) {
  x %>%
    left_join(y, by = c(names_common(x, y)))
}

#write a function to replace LM_ with LM
replace_LM <- function (x) {
 str_replace_all(x$EH_ID, "LM_", "LM") 
}

#join P4_Design and P4_record
P4_Des_Rec <- join_my_tables(P4_record, P4_Design)

#add the oocyst counts
P4_oocysts <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P4_082020_Eim_oocyst.csv") %>%
  select(!X)

#join the oocyst counts to the P4_Des_Rec
P4_Experiment <- join_my_tables(P4_Des_Rec, P4_oocysts)

#which column names are common
names_common(P4_Experiment, Challenge_Infections)

#compare the columns of the two dataframes
compare_df_cols(P4_Experiment, Challenge_Infections)

#create column infection_history
P4_Experiment %>%
  mutate(infection_history = paste0(primary_infection, "_",
         challenge_infection)) -> P4_Experiment

#produce the column infection type
P4_Experiment %>%
  mutate(infection_type = case_when(
    P4_Experiment$infection == "primary" & primary_infection == "UNI" ~ paste0("UNI"),
    P4_Experiment$infection =="challenge" & challenge_infection == "UNI" ~ paste0("UNI"),
    P4_Experiment$infection == "primary" ~ paste0("primary_", primary_infection),
    P4_Experiment$infection == "challenge" ~ paste0("heterologous_", challenge_infection),
    TRUE ~ "other"
  )) -> P4_Experiment

#write the table 
write.csv(P4_Experiment, "/home/fay/Documents/GitHub/Eimeria-PhD-Project/Lab_mouse_eimeria/Products/P4_Experiment", row.names = FALSE)

#bind the P4 experiments to the challenge infections
Challenge_infections_with_p4 <- bind_rows(Challenge_Infections, P4_Experiment)

#write the combination table
write.csv(Challenge_infections_with_p4, "/home/fay/Documents/GitHub/Eimeria-PhD-Project/Lab_mouse_eimeria/Products/Challenge_p4", row.names = FALSE)

#continue now with the p3 data

#read experimental design P3
P3_Design <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/P3_112019_Eim_design.csv")

#replace the LM_ with LM
P3_Design$EH_ID <- replace_LM(P3_Design)

#read Experimental record P4
P3_record <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_record.csv") 
P3_record$EH_ID <- replace_LM(P3_record)

#join P3_Design and P3_record
P3_Des_Rec <- join_my_tables(P3_record, P3_Design)

#add the oocyst counts
P3_oocysts <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/P3_112019_Eim_oocyst.csv")

#join the oocyst counts to the P4_Des_Rec
P3_Experiment <- join_my_tables(P3_Des_Rec, P3_oocysts)

#which column names are common
names_common(P3_Experiment, Challenge_Infections)

#compare the columns of the two dataframes
compare_df_cols(P3_Experiment, Challenge_Infections)

#create column infection_history
P3_Experiment %>%
  mutate(infection_history = paste0(primary_infection, "_",
                                    challenge_infection)) -> P3_Experiment

#produce the column infection type
P3_Experiment %>%
  mutate(infection_type = case_when(
    P3_Experiment$infection == "primary" & primary_infection == "UNI" ~ paste0("UNI"),
    P3_Experiment$infection =="challenge" & challenge_infection == "UNI" ~ paste0("UNI"),
    P3_Experiment$infection == "primary" ~ paste0("primary_", primary_infection),
    P3_Experiment$infection == "challenge" ~ paste0("heterologous_", challenge_infection),
    TRUE ~ "other"
  )) -> P3_Experiment

#bind the P4 experiments to the challenge infections
#error oocyst_mean in one table as character and in the other as factor
P3_Experiment$oocyst_mean <- as.character(P3_Experiment$oocyst_mean)
Challenge_infections_with_p3 <- bind_rows(Challenge_infections_with_p4, P3_Experiment)

#write the combination table
write.csv(Challenge_infections_with_p3, "/home/fay/Documents/GitHub/Eimeria-PhD-Project/Lab_mouse_eimeria/Products/Challenge_p4_p3.csv", row.names = FALSE)
