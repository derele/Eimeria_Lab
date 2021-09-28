### This could be an example of how to access particular subsets of
### the data
library(RCurl)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)

OV <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv"))

loadFromGH <- function(URL){
    if(url.exists(URL)){
        U <- getURL(URL)
        read.csv(text = U)
    } else {
        message("URL \"", URL, "\" does not exist")
    }
}

## ## Only the challenge experiments
ChallengeEx  <- c("E57", "E10", "E11")

## ## download and append the weigth tables
W <- lapply(OV[OV$Experiment%in%ChallengeEx, "weight"], loadFromGH)
Weight <- Reduce(rbind, W)

## ## Same for shedding
O <- lapply(OV[OV$Experiment%in%ChallengeEx, "shedding"], loadFromGH)
Oocysts <- Reduce(rbind, O)

## ## some don't agree
table(Oocysts$labels%in%Weight$labels)
table(Weight$labels%in%Oocysts$labels)

### FIXME!!!  It's experiment E11!!!
Oocysts[!Oocysts$labels%in%Weight$labels, ]
## FIXED by Franzi, just two labels missing: E11aBMI and E11aJOY

## there are more in in the weights table which are not found in the
## oocysts though
Weight[!Weight$labels%in%Oocysts$labels, ]
## Let's export this for Franzi to have a look:
write.csv(Weight[!Weight$labels%in%Oocysts$labels, ],
          "data/Experiment_results/FIXME_E11_Weight.csv", row.names=FALSE)

Results <- merge(Weight, Oocysts, all=TRUE)

## IDs sometimes with "_" sometimes without
Results$EH_ID <- gsub("LM_", "LM", Results$EH_ID)

## For now we hav to exclude those E11aBMI and E11aJOY which have no
## mouse IDs
Results <- Results[!is.na(Results$EH_ID), ]

## Same for design
D <- lapply(OV[OV$Experiment%in%ChallengeEx, "design"], loadFromGH)

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


write.csv(ALL, "data_products/Challenge_infections.csv", row.names=FALSE)

