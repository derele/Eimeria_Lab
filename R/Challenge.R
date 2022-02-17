### This could be an example of how to access particular subsets of
### the data
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tidyr)

#select columns: 
basics <- c("EH_ID", "mouse_strain", "experiment", "primary_infection", 
            "challenge_infection", "labels", "dpi", "infection", "infection_history")

weight_loss <- c("weight", "weight_dpi0", "relative_weight")

oocysts_counts <- c("feces_weight", "oocyst_sq1", "oocyst_sq2", "oocyst_sq3",
                    "oocyst_sq4", "dilution", "OOC", "OO4sq")

intensity_qPCR <- c("Eim_MC", "delta")

cewe_elisa <- "IFNy_CEWE"

mes_elisa <- "IFNy_MES"

#reading the overview table. In each row there is a link to the raw data for each experiment
OV <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv")

## ## Only the challenge experiments
#a list of the names of each experiments
# you can use it later, to select the experiments from ov 
ChallengeEx  <- c("E57", "E10", "E11", "P4", "P3")

## ## download and append the weigth tables
#lapply: applies a function to every element of the list
#we select the challenge experimetns and the weight columns 
#we apply to every element of the list the function read.csv
W <- lapply(OV[OV$Experiment%in%ChallengeEx, "weight"], read.csv)

#reduce: works on 1st and 2nd element, produces result and then uses the result 
#with the 3rd element and so on
#we apply in this case the function rbind
Weight <- Reduce(rbind, W)

## ## Same for shedding
O <- lapply(OV[OV$Experiment%in%ChallengeEx, "shedding"], read.csv)

Oocysts <- Reduce(rbind, O)

## ## some don't agree
table(Oocysts$labels%in%Weight$labels)
table(Weight$labels%in%Oocysts$labels)

## there are more in in the weights table which are not found in the
## oocysts though
Weight[!Weight$labels%in%Oocysts$labels, ]

## But if the weight is NA the mouse was dead
table(Weight[!Weight$labels%in%Oocysts$labels &
                 is.na(Weight$weight), "dpi"])
### confirmed by the late dpi of these!

Results <- merge(Weight, Oocysts, all=TRUE)

## IDs sometimes with "_" sometimes without
Results$EH_ID <- gsub("LM_", "LM", Results$EH_ID)

## For now we hav to exclude those E11aBMI which have no mouse IDs
Results <- Results[!is.na(Results$EH_ID), ]

## Same for design
D <- lapply(OV[OV$Experiment%in%ChallengeEx, "design"], read.csv)

Des.cols <- Reduce(intersect, lapply(D, colnames))

Design <- Reduce(rbind, lapply(D, "[", Des.cols))

## remove all whitespaces
Design %>%
    mutate(across(where(is.character), str_trim)) ->
    Design

## IDs sometimes with "_" sometimes without
Design$EH_ID <- gsub("LM_", "LM", Design$EH_ID)

#### NOW ALSO REMOVE EXPERIMENT column?
### Not for now
ALL <- merge(Design, Results, all=TRUE, by="EH_ID")
table(ALL$experiment.x, ALL$experiment.y)

## let's call experiment 5 and 7 57, like in all tables

Design$experiment[Design$experiment%in%c("E5", "E7")] <- "E57"
Results$experiment[Results$experiment%in%c("E5", "E7")] <- "E57"
ALL <- merge(Design, Results, all=TRUE)

### some mice don't have an infection history 
forgotten.mice <- unique(ALL[is.na(ALL$primary_infection), "EH_ID"])
## one mouse "LM0295" got lost?


Design$EH_ID[!Design$EH_ID%in%Results$EH_ID]
### Two samples with no result records
### [1] "LM0205" "LM0290"


## all the data has a mouse ID
table(is.na(ALL$EH_ID))

## but feces weights are ZERO T ( when the oocyst counts are NA ->
## the mouse was dead)!
table(ALL$feces_weight==0)

##  only one prolbem remaining after asking Alice for better data...
ALL[which(ALL$feces_weight==0 &
              !is.na(ALL$oocyst_sq2)), ]

ALL %>% filter(!is.na(challenge_infection)) %>%
    rowwise() %>% mutate(OO4sq = rowSums(across(starts_with("oocyst_")))) %>%
    ## 0.1µl per square -> *10.000 to scale up to ml
    mutate(OOC=(OO4sq/4*10000)/dilution) %>%
    ## we have ZEROS in feces weight (also when we counted oocysts) so we better don't
    ## calculate OPG for now but just max (see below)
    ## mutate(OPG=OOC/feces_weight) %>%
    ## also re-calculate relative weight, as this seems to have errors
    ## from a spreadsheet program (wtf!)
    mutate(relative_weight= weight/weight_dpi0*100) %>%
    ## also look at this for OPG above (by uncommenting)
    ## select(feces_weight, starts_with("oocyst_"), OO4sq, OOC) %>%
    ## look at this for controlling the weight calculation
    ## select(EH_ID, dpi, infection, weight, relative_weight) %>%
    ## print(n=40)
    ## the E88 innoculum in E57 challenge infection was "not
    ## working", these mice are basically unifected controls
    mutate(challenge_infection=ifelse(!experiment%in%"E57",
                                      challenge_infection,
                                      ifelse(challenge_infection%in%"E88", "UNI",
                                             challenge_infection))) %>%
    ## then correct the infection history
    mutate(infection_history=paste0(primary_infection, "_",
                                    challenge_infection)) ->
    ALL

#here is a function to join the following merged data frames to ALL
join_to_ALL <- function(x) {
    left_join(ALL, unique(x), by = c(intersect(colnames(x), colnames(ALL)))) 
}


#download and append the infection intensity tables (qPCR)
I <- lapply(OV[OV$Experiment%in%ChallengeEx, "infection_intensity"], read.csv)

#now combine the infection intensity tables
Intensity <- Reduce(rbind, I)

## Corrrect wrong IDs
Intensity$EH_ID <- gsub("LM_", "LM", Intensity$EH_ID)

## Melting curve sometimes positve sometimes TRUE, and the opposite
Intensity$Eim_MC <- gsub("pos", TRUE, Intensity$Eim_MC)
Intensity$Eim_MC <- gsub("neg", FALSE, Intensity$Eim_MC)

#Eim_MC type is character, change to logical
Intensity$Eim_MC <- as.logical(Intensity$Eim_MC)


#Combined all of the qPCR data for the challenge infections in "Intensity"
#Questions on how to go on: 
#to join this data to the "ALL" file I need to match for EH_ID, Experiment and dpi

##summarize for the maximum dpi for the challenge infections
ALL_sum_max_dpi <- ALL %>%
    group_by(EH_ID, infection, experiment) %>%
    summarize(max(dpi)) 

ALL_sum_max_dpi <- unique(ALL_sum_max_dpi)


##turn the table into wide format,  to see if both max dpi exists for challenge and primary
ALL_sum_max_dpi <- ALL_sum_max_dpi %>%
    pivot_wider(names_from = 'max(dpi)', values_from = infection) %>%
    rename(dpi_8 = '8', dpi_11 = '11')

# Make a death column to show when the mice for sacrificed
#if it is chal_8, then the mouse was sacrificed during the secondary challenge infections
# on day 8
#if it is chal_11, then the mouse was sacrificed during the primary infections
# on day 11
ALL_sum_max_dpi <- ALL_sum_max_dpi %>%
    mutate(death = 
               case_when(
                   !is.na(dpi_8) ~ "chal_8",
                   TRUE ~ "prim_11"
               ))

#now join the Intensity file to the summarized ALL_sum_max_dpi
Intensity <- ALL_sum_max_dpi %>%
    left_join(Intensity, by = c(intersect(colnames(Intensity), colnames(ALL_sum_max_dpi)))) %>%
    mutate(infection =
               case_when(
                   death == "chal_8" ~ "challenge",
                   TRUE ~ "primary"
               )) 

Intensity <- Intensity %>%
    select(EH_ID, experiment, death, Eim_MC, delta, infection)

#Now join the Intensity to the ALL file, while taking account of the death variable
ALL <- join_to_ALL(Intensity)

### Step: Join the CEWE_ELISA to our challenge infections file
## ## download and append the CEWE_ELISA tables
#lapply: applies a function to every element of the list
C <- OV[OV$Experiment%in%ChallengeEx, "CEWE_ELISA"] 

#I have to apply the read.csv to vector elemnts wich contain the raw data, therefore
#I have to first select from the OV file the lines with actual links to the raw files
C <- lapply(C[c(1,2,5)], read.csv)

CEWE_ELISA <- Reduce(rbind, C)

#next step clean the Mouse ID
## IDs sometimes with "_" sometimes without
CEWE_ELISA$EH_ID <- gsub("LM_", "LM", CEWE_ELISA$EH_ID)

#change IFNy to IFNy_cewe to show origin of the measurement
CEWE_ELISA <- CEWE_ELISA %>% rename(IFNy_CEWE = IFNy) 

#merge with ALL
ALL <- join_to_ALL(CEWE_ELISA)

##Joining data on ELISA from the mesenterial lymphnodes
#download and append the data from the ELISA's - Mesentrial Lymphnodes
#step by step, as there are empty spaces and this won't work
#C <- OV[OV$Experiment%in%ChallengeEx, "CEWE_ELISA"] 
M <- OV[OV$Experiment%in%ChallengeEx, "MES_ELISA"]

#I have to apply the read.csv to vector elemnts wich contain the raw data, therefore
#I have to first select from the OV file the lines with actual links to the raw files
#as there is only one row line with row data, lapply is not required 
MES_ELISA <- read.csv(M[[1]]) 

MES_ELISA <- MES_ELISA %>% select(EH_ID, IFNy)

## Corrrect wrong IDs
MES_ELISA$EH_ID <- gsub("LM_", "LM", MES_ELISA$EH_ID)

#change IFNy to IFNy_cewe to show origin of the measurement
<<<<<<< HEAD
MES_ELISA <- MES_ELISA %>% 
=======
MES_ELISA <- MES_ELISA %>% rename(IFNy_MES = IFNy) %>% 
>>>>>>> 45e31f5569e102c14159f0807e54fdbd5a8362cd
    mutate(experiment = "P4")

write.csv(ALL, "data/Experiment_results/P4_082020_Eim_MES_ELISA.csv", row.names=FALSE)

#change IFNy to IFNy_cewe to show origin of the measurement
MES_ELISA <- MES_ELISA %>% rename(IFNy_MES = IFNy)
#Now join the MES_ELISA to the ALL file
ALL <- join_to_ALL(MES_ELISA)

## Moving on to the gene expression data
#download and append the gene expression data
G <- OV[OV$Experiment%in%ChallengeEx, "gene_expression"]

G <- lapply(G[c(2,5)], read.csv)

#now combine the infection intensity tables
Gene_Expression <- Reduce(rbind, G) %>%
    select(EH_ID, CXCR3, IRG6, IL.12) #remove unecessary column X

## Corrrect wrong IDs
Gene_Expression$EH_ID <- gsub("LM_", "LM", Gene_Expression$EH_ID)

#join to the ALL file
ALL <- join_to_ALL(Gene_Expression)

#What's next? 
#FACS!

#download and append the FACS data
F <- OV[OV$Experiment%in%ChallengeEx, "FACS"]

F <- lapply(F[c(1,2,4)], read.csv)

F[[1]]



#go back and add an experiment tag for the cewe_elisa files and mes_elisa!

write.csv(ALL, "data_products/Challenge_infections.csv", row.names=FALSE)


