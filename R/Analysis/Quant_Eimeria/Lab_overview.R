library(tidyr)
library(dplyr)

challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")


#create a new column in the challenge infections where we see the actual status of infection
#during the measurements
challenge <- challenge %>%
  mutate(status = case_when(
    infection == "primary" ~ primary_infection,
    infection == "challenge" ~ challenge_infection,
    TRUE ~ "NA"
  ))

#How many mice do we have in the challenge infections
length(unique(challenge$EH_ID)) #153 mice


#Summarizing the data
Summary_challenge <- challenge %>% 
  dplyr::group_by(EH_ID, status) %>%
  dplyr::summarize(n()) %>%
  dplyr::group_by()


#just discovered that the lab data challenge infections have misnamed values
#go back to the challenge and correct



#%>%
  count(status, infection)



