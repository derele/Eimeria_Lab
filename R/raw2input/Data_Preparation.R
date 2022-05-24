## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Lab
## Script adapted from Victor's 1_Data_preparation

library(tidyverse)
library(dplyr)

## Quant2_E57 includes both genotype and flotation data merged unlike in Victor's case
##Hence, Victor's R script is adapted accordingly 

##Load data
#Quant2_E57 data has flotation and genotype data
#Nanadrop_2.0 has the DNA concentration data of our extractions
##E88_fecweight has the weight of fecal matter used in the DNA extraction
sample.data <- read.csv("data/Experiment_results/Quant_Eimeria/Quant2_E57_copy.csv")
Nanodrop <- read.csv ("data/Experiment_results/Quant_Eimeria/Fecal_DNA_Extractions/Nanodrop_2.0.csv")
fecweight <- read.csv ("data/Experiment_results/Quant_Eimeria/Fecal_DNA_Extractions/E88_fecweight.csv", header=T, na.strings=c("","NA"))

## Filter out the E64 (E. ferrisi) mice in Quant2_E57 file
sample.data <- sample.data %>% 
  dplyr::filter(!primary_infection == "E64") %>%
  dplyr::select(-X)

##rename the entries under the column 'Sample_ID' in file Nanodrop_2.0 to 
##E57aXXX to allow merge with file Quant2_E57
Nanodrop$labels = paste('E57a',Nanodrop$Sample_ID, sep = '')

##rename column 'DNA_conc' of Nanadrop_2.0 to Conc_Data to allow cohesion with Victor's data
Nanodrop <- Nanodrop %>% 
  rename(Conc_DNA = DNA_Conc)

########CHECK WITH PROFESSOR
##rename column 'feces_weight' of Quant2_E57 to 'fecweight_flot' to allow cohesion with Victor's data
sample.data <- sample.data %>% 
  rename(fecweight_flot = feces_weight)

##merge files Quant2_E57, Nanadrop_2.0 and fecweight
sample.data = left_join(sample.data, Nanodrop) %>%
  left_join(., fecweight)

##remove column 'Sample_ID' from merged file
sample.data = subset(sample.data, select = -c(Sample_ID) )

########CHECK WITH PROFESSOR
##R identifies the entries under column decweight_DNA as 'NULL'
#Hence, class changed to numeric
sample.data$fecweight_DNA<-sapply(sample.data$fecweight_DNA, as.numeric)

  ###Use function from Alice to calculate OPG
calculateOPG <- function(sample.data){
  sample.data$mean_Neubauer <- 
    (sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4) / 4
                                           # NB! Limit of detection = 1 oocysts
  sample.data$mean_Neubauer[sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4 == 1] <- 0
  sample.data$oocysts.per.tube <- sample.data$mean_Neubauer * 10000 * sample.data$dilution
  sample.data$OPG <- sample.data$oocysts.per.tube / sample.data$fecweight_flot
  ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
  sample.data$oocysts.per.tube[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
  sample.data$OPG[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
  return(sample.data)
}

sample.data <- calculateOPG(sample.data = sample.data)

##Check for spaces
sample.data$EH_ID <- gsub(pattern = " ", replacement = "", x = sample.data$EH_ID)

###Estimate microbiota density (MD) from Contijoch et al. 2019
### MD = total DNA per sample (µg)/ mg of fresh feces
##considering 30µL of elution volume 
sample.data %>%
  mutate(Total_DNA = (sample.data$Conc_DNA*30)*0.001) -> sample.data ### Add a new variable that will contain total DNA extracted per sample in µg

sample.data %>%
  mutate(Microbial_density = sample.data$Total_DNA/(sample.data$fecweight_DNA*1000)) -> sample.data ### Total DNA extracted per sample in µg by feces weight in mg

sample.data$labels<- as.vector(sample.data$labels)
rownames(sample.data) <- make.unique(sample.data$labels)

##Generate table 1
sample.data %>%
  dplyr::group_by(EH_ID)%>%
  dplyr:: select(EH_ID, weight_dpi0,ageAtdpi0expe1a, Sex, Locality, Genome, mouse_strain)%>%
  dplyr::distinct()%>%
  dplyr::rename(Weight_at_DPI_0= weight_dpi0, Age_at_DPI_0 = ageAtdpi0expe1a, Strain = mouse_strain)-> tmp1

write.csv(tmp1, "data/Experiment_results/Quant_Eimeria/Tables/tmp1.csv")
rm(tmp1)



