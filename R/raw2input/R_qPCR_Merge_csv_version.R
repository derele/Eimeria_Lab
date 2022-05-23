library(tidyr)
library(dplyr)
library(janitor)
library(readr)

setwd("/Users/vinuri/Documents/GitHub/Eimeria_Lab/data/Experiment_results/Quant_Eimeria/qPCR_faeces/Results/Results_files/qPCR_csv")
list_faeces <- as.list(list.files())
list_names <- as.vector(unlist(list_faeces))

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
df_results <- Reduce(rbind, list_results)
df_results <- unique(df_results)  

setwd("/Users/vinuri/Documents/GitHub/Eimeria_Lab/data/Experiment_results/Quant_Eimeria/qPCR_faeces/Results/Results_files/")
write.csv(df_results, "qPCR_faeces_lab_merged.csv", row.names=FALSE)





