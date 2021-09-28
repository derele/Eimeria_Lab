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
Weight[!Weight$labels%in%Oocysts$labels, ]

Results <- merge(Weight, Oocysts, all=TRUE)

## IDs sometimes with "_" sometimes without
Results$EH_ID <- gsub("LM_", "LM", Results$EH_ID)

## BUT we now have some NAs in the mouse IDs here:
table(is.na(Results$EH_ID))
Results[is.na(Results$EH_ID), ]
## These are labels for which in E11 (both first and challenge
## infection) no weight was in the table and thus tube labels have no
## mouse EH_ID association. Again: FIX ME!!!

## For now we hav to exclude those
Results <- Results[!is.na(Results$EH_ID), ]

## Same for design
D <- lapply(OV[OV$Experiment%in%ChallengeEx, "design"], loadFromGH)

Des.cols <- Reduce(intersect, lapply(D, colnames))

Design <- Reduce(rbind, lapply(D, "[", Des.cols))

## remove all whitespaces
Design %>%
    mutate(across(where(is.character), str_trim)) ->
    Design

Design[Design$experiment%in%c("E5", "E7"),]

## IDs sometimes with "_" sometimes without
Design$EH_ID <- gsub("LM_", "LM", Design$EH_ID)

#### NOW ALSO REMOVE EXPERIMENT AND ADD E7 Design to Design!

ALL <- merge(Design, Results, all=TRUE, by="EH_ID")

### some mice don't have an infection history 
forgotten.mice <- unique(ALL[is.na(ALL$primary_infection), "EH_ID"])

length(intersect(Design$EH_ID, Results$EH_ID))

length(union(Design$EH_ID, Results$EH_ID))

Design$EH_ID[!Design$EH_ID%in%Results$EH_ID]
### Two samples with now result records
### [1] "LM0205" "LM0290"

table(ALL$experiment.x, ALL$experiment.y)


## write.csv(ALL, "data_products/Challenge_infections.csv")


## oldE11 <- read.csv("data/Experiment_results/E11_Oocyst_CSV.csv")

## newE11 <- read.delim("data/Experiment_results/E11_Oocyst_updateCSV.csv", sep=";")

## newE11$labels <- paste0("E11", newE11$batch, newE11$labels)
## newE11$dilution <- 1

write.csv(newE11[,colnames(oldE11)],
          "data/Experiment_results/E11_Oocyst_updateCSV.csv", row.names=FALSE)

