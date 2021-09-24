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
## W <- lapply(OV[OV$Experiment%in%ChallengeEx, "weight"], loadFromGH)
## Weight <- Reduce(rbind, W)

## ## messed up some lables for E7 challenge
## Weight[grepl("aE7", Weight$labels),]


## ## Same for shedding
O <- lapply(OV[OV$Experiment%in%ChallengeEx, "shedding"], loadFromGH)
Oocysts <- Reduce(rbind, O)

## ## The same for oocysts
## head(Oocysts[grepl("aE7", Oocysts$labels), ])

## ## some don't agree
## table(Oocysts$labels%in%Weight$labels)
## table(Weight$labels%in%Oocysts$labels)

## ## But that's not the cause for the bigger problems as result are
## ## still okay, as both the labels are messed up in the same way:

## Results <- merge(Weight, Oocysts, all=TRUE)
## ## IDs sometimes with "_" sometimes without
## Results$EH_ID <- gsub("LM_", "LM", Results$EH_ID)

E57OO <- read.csv("data/Experiment_results/E5_062018_Eim_oocyst.csv")
E57OO <- E57OO[, !colnames(E57OO)%in%"X"]

E57OO[nchar(E57OO$labels)==3, "labels"] <-
    paste0("E57a", E57OO[nchar(E57OO$labels)==3, "labels"])

E57OO$labels <- gsub("E7a", "E57bx", E57OO$labels)

E57OO$labels <- gsub("E7b", "E57by", E57OO$labels)

E57OO$dilution <- as.numeric(gsub(",", ".", E57OO$dilution))

write.csv(E57OO, "data/Experiment_results/E57_xxxxx_Eim_oocyst.csv")


## BUT we now have some NAs in the mouse IDs here:
table(is.na(Results$EH_ID))
Results[is.na(Results$EH_ID), ]
## These are labels for which in E11 (both first and challenge
## infection) no weight was in the table and thus tube labels have no
## mouse EH_ID association. FIX ME!!!

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


write.csv(ALL, "data_products/Challenge_infections.csv")


