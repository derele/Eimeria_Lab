## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Lab
## Script adapted from Victor's 2_qPCR_data_preparation

library(ggpubr)
library(rcompanion)
library(dplyr)
library(gridExtra)
library(rstatix)
library(lmtest)
library(ggtext) 

##Load data
if(!exists("sample.data")){
  source("R/raw2input/Data_Preparation.R")
}
##Standard curves
data.std<- read.csv("data/Experiment_results/Quant_Eimeria/Eimeria_quantification_Std_Curve_data.csv")

##exclude E64 (E. ferrisi) data from standard curve data
data.std <- data.std %>% 
  dplyr::filter(!Parasite == "E_ferrisi") 

data.std%>%
  dplyr::mutate(Genome_copies= Oocyst_count*8)-> data.std

##Define numeric and factor variables 
num.vars <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", "Genome_copies")
fac.vars <- c("Well", "Sample.Name", "Detector", "Task",  "Std_series","Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction")  

## as.numeric alone will likely fail if stringsAsfactors is TRUE! 
data.std[, num.vars] <- apply(data.std[, num.vars], 2,
                              function (x) as.numeric(as.character(x)))
data.std[, fac.vars] <- apply(data.std[, fac.vars], 2, as.factor)

##Correct zero in NTC with not-detected results 
data.std$Ct[data.std$Ct == 0] <- NA
data.std$Ct_mean[data.std$Ct_mean == 0 & data.std$Task == "NTC"] <- NA
data.std$Sd_Ct[data.std$Sd_Ct == 0] <- NA

##Correct labels to have them homogeneous 
data.std$Sample.Name<- gsub(pattern = " ", replacement = "_", x = data.std$Sample.Name)

## Select just standards data
## Estimate the number of genome copies per ng of gDNA
data.std.lm<- subset(data.std, Task== "Standard") ## Select just data from standards 
data.std.lm %>% 
  dplyr::select(Sample.Name, Task, Ct, Cycler, Oocyst_count, Parasite, Genome_copies)%>%
  dplyr::mutate(Oocyst_DNA= Oocyst_count*(3.8E-4))%>% ##Estimation of DNA (ng) derived from Oocyst
  dplyr::mutate(DNA_PCR= Oocyst_DNA/30)%>% ##DNA (ng) in PCR considering 1uL from a stock of 30uL
  ##Considering that 1 ng of Eimeria gDNA is equivalent to 2.11E4 genome copies
  dplyr::mutate(Genome_copies_ngDNA= (2.11E4)*DNA_PCR)-> data.std.lm  

##Infection experiment: load qPCR data (lab_fecal) and clean 
data.inf.exp<-read.csv("data/Experiment_results/Quant_Eimeria/qPCR_faeces/Results/Results_files/qPCR_fecal_lab_merged.csv")

### adapt all data in files to allow merge into challenge_infection file:
# this includes changing column names of files to match that of challenge_infection file
# change column name 'Sample' to 'labels' of qPCR file
data.inf.exp <- data.inf.exp %>% 
  rename(labels = Sample)

#entries under column require the characters 'E57a' to allow cohesive merge
data.inf.exp$labels[239:415] <- paste0('E57a', data.inf.exp$labels[239:415])
data.inf.exp <- data.inf.exp %>%
  dplyr::mutate(labels = case_when(
    labels == "E57INR" ~ "E57aINR",
    labels == "E57CDE" ~ "E57aCDE", 
    labels == "E57CEW" ~ "E57aCEW", 
    labels == "E57EFU" ~ "E57aEFU",
    labels == "E57NTC" ~ "E57aNTC",
    labels == "E57aNTC1" ~ "E57aNTC",
    labels == "E57aNTC2" ~ "E57aNTC",
    labels == "E57aBHN2" ~ "E57aBHN",
    labels == "E57aFLV2" ~ "E57aFLV",
    labels == "E57aIJQ2" ~ "E57aIJQ",
    labels == "E57aRWY2" ~ "E57aRWY",
    labels == "NTC" ~ "E57aNTC",
    TRUE ~ labels
  ))

#unsuccessful qPCR replicates are deselected from data frame
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aABD" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aABL" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aAFR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aAJT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aAJT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aALU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aATW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCE" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCE" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBDY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBHN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBHN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBLN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBLX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBLX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBOR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBOR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBTZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCEG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCEG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_11042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCEG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCFG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCIW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCIW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCLT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCLT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCNS" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCQT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCQT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCUX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCUX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDMV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDTY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDTY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aAKO" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aBGU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aBGU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCFW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCJQ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCSX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDEK" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEIP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEIP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEIY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEKO" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEKO" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEQU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEQU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aESV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aESV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aETU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aETU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEKR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEMX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEUW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEMX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFJQ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFKM" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFKW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFRW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFRW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGHL" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGHU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGHU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGNZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGRV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHIS" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHJT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHQZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIJQ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIKZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIOY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIPT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_11042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIPY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aJOS" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKMX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKUV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKUV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLMW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLNU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLNU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLPU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLSZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLXZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aMQT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aMRT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_11042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aMRT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- filter(data.inf.exp, labels != "E57aNTC")
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aOPY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aOPY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aQVY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aQVY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRTV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRTV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aSUW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aTUY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aTUY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aUWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aUWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aUWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aVWX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aVWX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]

##E57aRSW is a peculiar case, check T1,T2,T3 values of all replicates to recall 
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aRSW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]

#to unsure the absence of qPCR replicates, the following code should produce no results
#the number of qPCR values for each label (E57aXXX) should be at a frequency of 3
#code displays labels containing more than 3 replicates 
Freq_list <- data.frame(table(data.inf.exp$labels))
Freq_list = Freq_list[Freq_list$Freq > 3,]
rm(Freq_list)

##apply Victors script logically
data.inf.exp <- rename (data.inf.exp, c(Ct = Cq, Tm = Tm1))

data.inf.exp%>%
  dplyr::mutate(Cycler= "BioRad")-> data.inf.exp #########CHECK IF CYCLER NEEDED----------------------<<<<<<<<<<

##Define numeric and factor variables 
num.vars3 <- c("Ct", "Tm")
fac.vars3 <- c("labels", "Task", "plate", "Cycler")  #########CHECK IF CYCLER NEEDED------------------<<<<<<<<<<
data.inf.exp[, num.vars3] <- apply(data.inf.exp[, num.vars3], 2,
                                   function (x) as.numeric(as.character(x)))
data.inf.exp[, fac.vars3] <- apply(data.inf.exp[, fac.vars3], 2, as.factor)

rm(fac.vars, num.vars, fac.vars3, num.vars3)

####### Standard curves 
set.seed(2020)
## Compute simple linear models from standards
## "Genome copies modeled by Ct"

##Ct modeled by Oocyst counts; data from different Cyclers   !!!!!!------------------> is this needed?
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Ct, color= Cycler)) +
  geom_smooth(method = "lm", se = T, aes(Oocyst_count, Ct, color=Cycler)) +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), 
              aes(size= 20,fill= Cycler, shape= Parasite), color= "black", alpha= 0.5)+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "b") -> tmp.fig.1

##Ct modeled by Gene counts; data from different Cyclers
data.std.lm%>%
  ggplot(aes(x = Ct, y = Genome_copies, color= Cycler)) +
  geom_smooth(method = "lm", se = F, size= 0.5) +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 15, fill= Cycler), color= "black", alpha= 0.25)+
  stat_cor(label.x = 25, label.y = c(8,7,6), 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+ # Add correlation coefficient
  stat_regline_equation(label.x = 25, label.y = c(8.5,7.5,6.5))+ # Add Regression equation lm log10(Genome_copies)~Ct+Cycler
  labs(tag = "a", y= "log 10 *Eimeria* genome copies")+ ##Make Eimeria in italics ;)
  theme_bw() +
  theme(text = element_text(size=20), legend.position= "none", axis.title.y = ggtext::element_markdown())+
  annotation_logticks(sides = "l")-> A

##Ct modeled by Oocyst_counts and extra predictors to be considered 
##Model 3: Ct modeled by oocyst count and cycler as predictors
lm.CtCyc<- lm(Ct~log10(Oocyst_count)+Cycler, data.std.lm)
rm(lm.CtCyc)
##Model 3 fit better the data... Cycler has major impact (confirm somehow our expectations)!

##Real standard curve##
##Genome copies modeled by Ct and extra predictors to be considered 
##Model 8: Genome copies modeled by Ct and cycler as predictors
lm.SCCyc<- lm(log10(Genome_copies)~Ct+Cycler, data.std.lm, na.action = na.exclude)
##Model 8 fit better the data... Cycler has major impact (again confirm expectations)!
##Linear model (Standard curve for the rest of experiments)
data.std.lm%>%
  ggplot(aes(x = Ct, y = Genome_copies)) +
  geom_smooth(method = "lm", se = T, color="black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "B)", y= "log 10 *Eimeria* genome copies")+
  theme_bw() +
  theme(text = element_text(size=20), legend.position= "none", axis.title.y = element_markdown())+
  annotation_logticks(sides = "l") -> tmp.fig.2

A +
  geom_smooth(method = "lm", se = T, color="black", size= 1.5)+
  geom_text (x = 12, y = 3, show.legend = F,
             label = paste ("y = 10.08 - 0.26 x \n R-squared= 0.94, p < 2.2e-16"), color="black") -> A

##Predicted genome copies
### Using MODEL 8 to predict using different levels of the factor cycler
data.std.lm$predicted<- 10^predict(lm.SCCyc)
data.std.lm$residuals<- 10^residuals(lm.SCCyc)

data.std.lm%>%
  dplyr::select(predicted, Cycler)%>%
  group_by(Cycler) %>%
  get_summary_stats(predicted, type = "mean_sd")

data.std.lm%>%
  dplyr::select(predicted, Cycler)%>%
  anova_test(predicted ~ Cycler)-> cycler.aov 

data.std.lm%>%
  dplyr::select(predicted, Cycler)%>%
  pairwise_t_test(predicted ~ Cycler, p.adjust.method = "bonferroni")-> cycler.pwc

# Show adjusted p-values
cycler.pwc%>%
  add_xy_position(x = "Cycler")%>%
  mutate(y.position= log10(y.position))-> cycler.pwc

data.std.lm%>%
  dplyr::select(predicted, Cycler)%>%
  ggboxplot(x = "Cycler", y = "predicted", color = "black", 
            fill = "Cycler", palette =c("#00BA38", "#F8766D", "#619CFF"), ylab = "log10 Predicted Eimeria Genome copies") +
  yscale("log10")+
  stat_pvalue_manual(cycler.pwc, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(cycler.aov, detailed = TRUE),
       caption = get_pwc_label(cycler.pwc), tag = "B)")+
  theme_bw()+
  theme(text = element_text(size=20), legend.position= "top")+
  font("caption", size = 14)+
  font("subtitle", size = 14)-> tmp.fig.3

##Linear model Genome copies modeled by Oocyst count 
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Genome_copies)) +
  geom_smooth(method = "lm", se = F, color= "black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "b", x= "log 10 *Eimeria* Oocysts Count (Flotation)", y= "log 10 *Eimeria* genome copies (qPCR)")+
  theme_bw() +
  theme(text = element_text(size=20), legend.position= "top",
        axis.title.x = element_markdown(), axis.title.y = element_markdown())+
  annotation_logticks(sides = "bl")-> B

##Model 11: Genome copies modeled by Oocyst count and cycle 
lm.SCOoc<- lm(log10(Genome_copies_ngDNA)~log10(Oocyst_count)+Cycler, data.std.lm)

## ### Figure 1 Final Standard curves 
#pdf(file = "fig/Figure_1.pdf", width = 8, height = 10)
#grid.arrange(A, B)
#dev.off()
rm(A,B, lm.SCOoc, cycler.aov, cycler.pwc)
## If it is necessary some of the previous figures could be included as supplementary  

##Mean comparison standards against NTC (Supplementary 1)
set.seed(2020)
data.std%>%
  dplyr::select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Date)%>%
  dplyr::filter(Task%in%c("Standard", "NTC"))%>%
  ggplot(aes(x = Sample.Name, y = Ct)) +
  scale_x_discrete(name = "Standard",
                   labels= c("Eimeria_10_0"= "Oocysts 10^0", "Eimeria_10_1"= "Oocysts 10^1",
                             "Eimeria_10_2"= "Oocysts 10^2", "Eimeria_10_3"= "Oocysts 10^3",
                             "Eimeria_10_4"= "Oocysts 10^4", "Eimeria_10_5"= "Oocysts 10^5",
                             "Eimeria_10_6"= "Oocysts 10^6", "H2O"= "NTC")) +
  scale_y_continuous(name = "Ct")+ 
  geom_jitter(shape=21, position=position_jitter(0.2), color= "black", size= 5, alpha= 0.5,
              aes(fill= Cycler))+
  theme_bw() +
  theme(text = element_text(size=16),legend.position = "top")+
  theme(axis.text.x = element_text(angle=-90))+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange",
               shape=16, size=0.5, color="black")+
  labs(tag = "a")+
  geom_hline(yintercept = 30, linetype = 2)+
  stat_compare_means(method = "anova",
                     aes(label = paste0(..method.., ",\n","p=",..p.format..)),
                     label.y= 33, label.x = 7)+
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "H2O", 
                     label.y = c(40, 39, 33, 30.5, 24, 20, 16, 0)) -> Supp_1

###Determine that 10^0 and 10^1 measurements are basically like NTC when all the information is taken into account

### Tm as complementary reference for Negative samples 
set.seed(2020)
data.std%>%
  dplyr::select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Date)%>%
  dplyr::filter(Task%in%c("Standard", "NTC"))%>%
  ggplot(aes(x = Sample.Name, y = Tm)) +
  scale_x_discrete(name = "Standard",
                   labels= c("Eimeria_10_0"= "Oocysts 10^0", "Eimeria_10_1"= "Oocysts 10^1",
                             "Eimeria_10_2"= "Oocysts 10^2", "Eimeria_10_3"= "Oocysts 10^3",
                             "Eimeria_10_4"= "Oocysts 10^4", "Eimeria_10_5"= "Oocysts 10^5",
                             "Eimeria_10_6"= "Oocysts 10^6", "H2O"= "NTC")) +
  scale_y_continuous(name = "Tm")+ 
  geom_jitter(shape=21, position=position_jitter(0.2), color= "black", size= 5, alpha= 0.5,
              aes(fill= Cycler))+
  theme_bw() +
  theme(text = element_text(size=16),legend.position = "none")+
  theme(axis.text.x = element_text(angle=-90))+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange",
               shape=16, size=0.5, color="black")+
  labs(tag = "b")+
  geom_hline(yintercept = 75, linetype = 2)+
  stat_compare_means(method = "anova",
                     aes(label = paste0(..method.., ", ","p=",..p.format..)),
                     label.y= 76, label.x = 6)-> Supp_2

### Estimate mean Eimeria Tm
data.std%>%
  dplyr::select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Date)%>%
  dplyr::filter(Task%in%c("Standard"))%>%
  dplyr::filter(Cycler%in%c("BioRad"))%>%
  dplyr::select(Sample.Name,Task,Std_series,Ct,Oocyst_count,Tm)%>%
  dplyr::filter(complete.cases(.))%>%
  dplyr::summarise(mean = mean(Tm), sd= sd(Tm), n = n())%>%
  dplyr::mutate(se = sd / sqrt(n),
                lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
                upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)%>%
  dplyr::mutate(upper.ran = mean + (2*sd), 
                lower.ran = mean - (2*sd))

#pdf(file = "fig/Supplementary_1.pdf", width = 10, height = 15)
#grid.arrange(Supp_1, Supp_2)
#dev.off()
rm(Supp_1, Supp_2)

######### Infection experiment data############
## Define real positive and negatives based on Tm 

### Estimate mean Eimeria Tm for positive controls
data.inf.exp%>%
  dplyr::select(Task,labels,Tm)%>%
  dplyr::filter(Task%in%c("Pos_Ctrl"))%>%
  dplyr::filter(complete.cases(.))%>%
  dplyr::mutate(Tm = as.numeric(Tm))%>%
  #dplyr::slice(1,4,7,11,13)%>%
  dplyr::summarise(mean = mean(Tm), sd= sd(Tm), n = n())%>%
  dplyr::mutate(se = sd / sqrt(n),
                lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
                upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)%>%
  dplyr::mutate(upper.ran = mean + (2*sd), 
                lower.ran = mean - (2*sd))

##Positive controls have a Tm at:
#mean        sd n  se lower.ci upper.ci upper.ran lower.ran
#74.1 0.8944272 5 0.4 72.98942 75.21058  75.88885  72.31115

#there is a second Tm higher than 80°C from an unspecific product

##The mean Tm is at 74.1°C and sd of 0.89
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
data.inf.exp %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                      (Tm >= 76 | Tm <= 72.2)  ~ "Negative",
                                      (Tm >=72.3 | Tm <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm == 3 ~ "Positive",
                                      Count.Tm != 3 ~ "Negative"))%>%
  dplyr::left_join(data.inf.exp, by= "labels")-> data.inf.exp 

##Estimate number of genome copies with qPCR Ct value (Model 8)
data.inf.exp$Genome_copies<- 10^predict(lm.SCCyc, data.inf.exp)

##Summarize genome copies by sample  
data.inf.exp %>%
  dplyr::filter(Infection=="Positive")%>% ## Select true positives
  dplyr::select(Genome_copies,labels)%>% # select variables to summarize
  na.omit()%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise_each(funs(Genome_copies_min = min, Genome_copies_q25 = quantile(., 0.25),
                             Genome_copies_median = median, Genome_copies_q75 = quantile(., 0.75), 
                             Genome_copies_max = max, Genome_copies_mean = mean, Genome_copies_sd = sd)) -> Sum.inf

##Summarize Tm by sample 
data.inf.exp %>%
  dplyr::filter(Infection=="Positive")%>%
  dplyr::mutate(Tm = as.numeric(Tm))%>%## Select true positives
  dplyr::select(Tm,labels)%>% # select variables to summarize
  na.omit()%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise_each(funs(Tm_mean = mean, Tm_sd = sd))%>%
  dplyr::left_join(Sum.inf, by= "labels")-> Sum.inf

##Join summarized data
data.inf.exp<- left_join(data.inf.exp, Sum.inf, by= "labels")

##Eliminate an unprocessed sample and controls
data.inf.exp%>%
  dplyr::mutate(Genome_copies = case_when((Infection=="Negative")~ 0,
                                          TRUE~ Genome_copies))%>% ## Make Negative zero
  dplyr::select(labels, Genome_copies_mean, Tm_mean, Infection)%>%
  dplyr::filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl"))%>% ## Replace NAs in real negative samples to 0 
  dplyr::mutate(Genome_copies_mean= replace_na(Genome_copies_mean, 0))%>%
  dplyr::mutate(Tm_mean= replace_na(Tm_mean, 0))%>%
  ##Get unique labels from qPCR data
  dplyr::distinct(labels, .keep_all = TRUE)-> data.inf.exp

##Merging Infection experiment oocyst and weight loss data with qPCR data
##Check differences between two dataframes
setdiff(sample.data$labels, data.inf.exp$labels)

###Join all the data in the same dataframe
sdt<- left_join(sample.data, data.inf.exp, by="labels") ## Add qPCR data
###Tiny adjustment  
sdt$dpi<- as.factor(sdt$dpi)

sdt%>%
  dplyr::mutate(Genome_copies_ngDNA= Genome_copies_mean/50, ## copies by ng of fecal DNA considering 1uL from 50 ng/uL DNA
                DNA_sample= Conc_DNA*30, ## Estimate total gDNA of sample considering 30uL of elution buffer
                DNA_g_feces= DNA_sample/fecweight_DNA,
                ## Transform it to ng fecal DNA by g of faeces
                Genome_copies_gFaeces= Genome_copies_ngDNA*DNA_g_feces) -> sdt ## Estimate genome copies by g of faeces

##Transform to zero OPGs for DPI 1 and 2 (Dile and 3)
sdt$OPG[sdt$dpi==1] <- 0
sdt$OPG[sdt$dpi==2] <- 0  
sdt$OPG[sdt$dpi==3] <- 0

##Remove dataframes with data not related to the infection experiment data that won't be used in the following scripts
rm(data.std, data.std.lm, Sum.inf, 
   tmp.fig.1, tmp.fig.2, tmp.fig.3)












