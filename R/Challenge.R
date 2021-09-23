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

## Only the challenge experiments
ChallengeEx  <- c("E57", "E10", "E11")

## download and append the weigth tables
W <- lapply(OV[OV$Experiment%in%ChallengeEx, "weight"], loadFromGH)
Weight <- Reduce(rbind, W)

## Same for shedding
O <- lapply(OV[OV$Experiment%in%ChallengeEx, "shedding"], loadFromGH)
Oocysts <- Reduce(rbind, O)


Results <- merge(Weight, Oocysts, all=TRUE)
## IDs sometimes with "_" sometimes without
Results$EH_ID <- gsub("LM_", "LM", Results$EH_ID)

## We have some NAs in the mouse IDs here: 
table(is.na(Results$EH_ID))
Results[is.na(Results$EH_ID), ]

## Labels for which in E11 (both first and challenge infection) no
## weight was in the table and thus tube labels have no mouse EH_ID
## association. FIX ME!!!

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

## IDs sometimes with "_" sometimes without
Design$EH_ID <- gsub("LM_", "LM", Design$EH_ID)

ALL <- merge(Design, Results, all=TRUE)

write.csv(ALL, "data_products/Challenge_infections.csv")
