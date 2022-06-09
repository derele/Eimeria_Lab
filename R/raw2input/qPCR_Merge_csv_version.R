library(tidyr)
library(dplyr)
library(janitor)
library(readr)

#Ms. Fay Webster's function was adpated to merge all qPCR files into one

setwd("/data/Experiment_results/Quant_Eimeria/qPCR_faeces/Results/Results_files/qPCR_csv")

#all individual qPCR.csv files form a list
list_faeces <- as.list(list.files())

#to all a function to perform itself on the list, it is transformed into a vector
list_names <- as.vector(unlist(list_faeces))

#function set by Fay Webster to deleted rows 1 to 24 and set a new column that 
#includes the file name i.e. date of qPCR experiment
read_qPCR_file <- function(x) {
  df1 = read_csv(x)
  filename <- colnames(df1[2])
  df1 <- df1 %>%
    filter(!row_number() %in% c(1:24))
  colnames(df1) <- df1[1, ]
  df1 <- df1 %>% filter(!row_number() %in% 1)
  df1 <- df1 %>% dplyr::mutate(plate = filename)
}

list_results <- lapply(list_names, read_qPCR_file)

#bind/merge all indivual csv files into one
df_results <- Reduce(rbind, list_results)

#removes duplicates
df_results <- unique(df_results)  

setwd("/data/Experiment_results/Quant_Eimeria/qPCR_faeces/Results/Results_files/")
write.csv(df_results, "qPCR_fecal_lab_merged.csv", row.names=FALSE)





