### This could be an example of how to access particular subsets of
### the data
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)

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


## For an analysis of immune protection we want the following
## categories in one column
ALL$infection_type <- NA
## primary_E88
ALL$infection_type[ALL$infection_history%in%"UNI_E88" &
                   ALL$infection%in%"challenge"] <-  "primary_E88"
ALL$infection_type[ALL$primary_infection%in%"E88" &
                   ALL$infection%in%"primary"] <-  "primary_E88"
## homologous_E88
ALL$infection_type[ALL$infection_history%in%"E88_E88" &
                   ALL$infection%in%"challenge"] <-  "homologous_E88"
## heterologous_E88 
ALL$infection_type[ALL$infection_history%in%"E64_E88" &
                   ALL$infection%in%"challenge"] <-  "heterologous_E88"
## primary_E64  ("UNI_E64" and challenge in infection) E64_* and primary in infection
ALL$infection_type[ALL$infection_history%in%"UNI_E64" &
                   ALL$infection%in%"challenge"] <-  "primary_E64"
ALL$infection_type[ALL$primary_infection%in%"E64" &
                   ALL$infection%in%"primary"] <-  "primary_E64"
## homologous_E64
ALL$infection_type[ALL$infection_history%in%"E64_E64" &
                   ALL$infection%in%"challenge"] <-  "homologous_E64"
## heterologous_E64
ALL$infection_type[ALL$infection_history%in%"E88_E64" &
                   ALL$infection%in%"challenge"] <-  "heterologous_E64"
## the remaining should be UNI?!
ALL$infection_type[is.na(ALL$infection_type)] <- "UNI"

## ## download and append the infection_intensity tables
#lapply: applies a function to every element of the list
#we select the challenge experimetns and the qpcr columns 
#we apply to every element of the list the function read.csv
I <- lapply(OV[OV$Experiment%in%ChallengeEx, "infection_intensity"], read.csv)

#write a function to left join 
join_my_tables <- function (x, y) {
    x %>%
        left_join(y, by = c(names_common(x, y)))
}

I[[2]] <- I[[2]] %>%
    select(colnames(I[[1]]))

write.csv(I[[2]], "data/Experiment_results/E1_012017_Eim_CEWE_qPCR.csv", row.names = FALSE)


#reduce: works on 1st and 2nd element, produces result and then uses the result 
#with the 3rd element and so on
#we apply in this case the function rbind
Intensity <- Reduce(rbind,I)

ALL <- join_my_tables(ALL, I)
write.csv(ALL, "data_products/Challenge_infections.csv", row.names=FALSE)

