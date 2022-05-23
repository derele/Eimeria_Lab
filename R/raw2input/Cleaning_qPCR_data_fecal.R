library(tidyr)
library(dplyr)
library(tidyverse)


#### Reading the qPCR data that have already been merged 
## read the experimental design data 
ed <- read.csv("data_products/Quant2_E57.csv")

## let's filter out the mice infected with E. ferrisi (E64)
ed <- ed %>% 
  dplyr::filter(!primary_infection == "E64") %>%
  dplyr::select(-X)

## These qpcr reads have emerged from laboratory mice, from the faeces 
rq <- read.csv("data/Experiment_results/Quant_Eimeria/qPCR_faeces/Results/Results_files/qPCR_faeces_lab_merged.csv")

### Look at the data!
str(rq)
glimpse(rq)

# how many nas do we have?
sapply(rq, function(x) sum(is.na(x)))

### Cleaning:
# change the columns to match the Challenge infections df /quant eimeria df 
# preference would be to match the challenge infections 
rq <- rq %>% 
  rename(labels = Sample)

# compare the columns between the df
colnames(ed)
colnames(rq)
intersect(colnames(ed), colnames(rq))

rq <- rq %>%
  dplyr::mutate(labels = case_when(
    labels == "E57INR" ~ "E57aINR",
    TRUE ~ labels
   ))

# we have some labels that are incorrect, missing prefix E57a 
# we select the specific rows and add the prefix to the labels 
rq$labels[239:415] <- paste0('E57a', rq$labels[239:415])

# what happens with the nas in the labels? 


#https://www.google.com/search?q=data+cleaning+steps&oq=data+cleaning+steps&aqs=chrome..69i57j0i22i30l9.3172j0j7&sourceid=chrome&ie=UTF-8
# https://www.tableau.com/learn/articles/what-is-data-cleaning

#Step 1: Remove duplicate or irrelevant observations


# look at VIctor's script and see how he dealt with the triplicates 
# otherwise it will cause "issues" later in the merging, for example you will 
# get many more observations


#Step 2: Fix structural errors
# what about the columns sq in Victor's df
# look at what type is every column (use glimpse or typeof(), and then make
#changes if necessary) (numeric? , character?)

#Step 4: Handle missing data # if you wont handle then document 

#Step 5: Validate and QA
# join the two data frames- 
# example functions left_join /full_join / right_join
# https://www.youtube.com/watch?v=Yg-pNqzDuN4

#validate everything
# see if you created duplicated
# count 
