## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Lab

library(tidyr)
library(dplyr)
library(tidyverse)

#### Load the E88 data (include mice & flotation data)
E88 <- read.csv("data_products/Quant2_E57.csv")

### Load qPCR data (lab_fecal)
qPCR <- read.csv("data/Experiment_results/Quant_Eimeria/qPCR_faeces/Results/Results_files/qPCR_fecal_lab_merged.csv")

## Filter out the E64 (E. ferrisi) mice
E88 <- E88 %>% 
  dplyr::filter(!primary_infection == "E64") %>%
  dplyr::select(-X)

### adapt all data in files to allow merge into challenge_infection file:
# this includes changing column names of files to match that of challenge_infection file
# change column name 'Sample' to 'labels' of qPCR file
qPCR <- qPCR %>% 
  rename(labels = Sample)

#entries under column require the characters 'E57a' to allow cohesive merge
qPCR <- qPCR %>%
  dplyr::mutate(labels = case_when(
    labels == "E57INR" ~ "E57aINR",
    TRUE ~ labels
   ))
qPCR$labels[239:415] <- paste0('E57a', qPCR$labels[239:415])
qPCR$labels[196:198] <- paste0('E57a', qPCR$labels[196:198])

### adapt all data in files to allow consistency with Victors R script:
#change column name of qPCR from 'Ct' to 'Cq'.
#according to Bustin et al., 2009 (https://doi.org/10.1373/clinchem.2008.112797) 
#both terms refer to the same value
#line 95 of Victors data_prepartion R script
qPCR = rename(qPCR, c(Ct = Cq, Ct_mean = Cq.Mean, Sd_Ct = Cq.SD, Tm = Tm1))

##### this point onwards: logically replicate Victor's script 'data_preparation' to produce Figure 2 of his paper 

##Define numeric and factor variables (line 98)  
num.vars1 <- c("Ct", "Tm")
fac.vars1 <- c("labels", "plate")  

### Estimate mean Eimeria Tm for positive controls (lines 525 - 556) #........tbc

##### ask Professor the following : 
##found code to ' Estimate mean Eimeria Tm' in Victors script line 331
# What do we do with Tm2, Tm3,Tm4 - nothing found on Victor to clarify 
# Think of a way to organise the repeated qPCR

