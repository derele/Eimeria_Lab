library(tidyverse)
library(ggplot2)
library(dplyr)

#read challenge infections
Challenge_Infections_FACS <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/FACS_intensity.csv")



#facs
facs <- c("Th1", "Div_Act_CD8", "Div_Th1",  "IFNy_CD4", "CXCR3", "CD4", 
          "Th17", "IFNy_CD8", "IRG6", "Treg", "Div_Th17", "IL.12", "Div_Treg", 
          "CD8", "IFNy_FEC", "Treg_prop", "Treg17", "Act_CD8", "IFNy_CEWE", "IL17A_CD4")


#start plotting 

#first plot on infection intensity vs cells 
#let's see the differences between primary and challenge infections
#clean the data
Ch_th1 <- Challenge_Infections_FACS %>% 
  select(EH_ID, infection, delta, Th1, primary, mouse_strain, infection_type, primary_infection, 
         challenge_infection, challenge, infection_history, dpi, Eim_MC) %>%
  filter(!is.na(delta)) %>% #filter out the na's
  filter(!is.na(Th1)) %>% #filter out the na's
  filter(Eim_MC == TRUE) #filter for positive melting curve



ggplot(Ch_th1, aes(infection, Th1)) +
  geom_boxplot() +
  facet_wrap(~ infection_history)
        
 
  
