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

Oocysts[!Oocysts$labels%in%Weight$labels, ]
## FIXED by Franzi, just one labels missing: E11aBMI

## there are more in in the weights table which are not found in the
## oocysts though
Weight[!Weight$labels%in%Oocysts$labels, ]

## But if the weight is NA the mouse was dead
table(Weight[!Weight$labels%in%Oocysts$labels &
             is.na(Weight$weight), "dpi"])
### confirmed by the late dpi of these!


## we have this prolem remaining:
OOcystProblem <- Weight[!Weight$labels%in%Oocysts$labels &
                        !is.na(Weight$weight),]

## Let's export this for Franzi to have a look:
write.csv(OOcystProblem,
          "data/Experiment_results/FIXME_E11_OOcysts.csv", row.names=FALSE)

table(OOcystProblem[,"dpi"])
## seems like oocyst from dip 4 in primary infection have been omitted
## Franzi Says one box has been missing for counting

Results <- merge(Weight, Oocysts, all=TRUE)

## IDs sometimes with "_" sometimes without
Results$EH_ID <- gsub("LM_", "LM", Results$EH_ID)

## For now we hav to exclude those E11aBMI which have no mouse IDs
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

## all the data has a mouse ID
table(is.na(ALL$EH_ID))

## but feces weights are ZERO T
table(ALL$feces_weight==0)

table(ALL[which(ALL$feces_weight==0), "experiment"])
table(ALL[which(ALL$feces_weight==0), "dpi"])
table(ALL[which(ALL$feces_weight==0), "EH_ID"])

### need to this data to correct the E57 data 

## ## Alice
AE5 <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_RelatedParasitesResTol/master/data/ExpeDF_005_Alice.csv")

## ## make agreeing labels
## AE5$labelsNew <- paste0("E57a", AE5$labels)

## ## Missing feces weight
## WE57 <- read.csv("data/Experiment_results/E57_xxxxx_Eim_record.csv")

## replace <- WE57$labels[WE57$labels%in%AE5$labelsNew]

## ## so we replace all the primary infection
## table(nchar(replace))
## table(WE57[WE57$labels%in%replace, "infection"])

## NEW57.1 <- WE57[!WE57$labels%in%replace, ]

## ## and leave only challenge infection in place
## table(NEW57.1$infection)

## NEW57.2<- AE5[AE5$labelsNew%in%replace,
##               c("EH_ID", "labelsNew", "weight", "weight_dpi0",
##                 "relativeWeight", "fecweight",
##                 "dpi")]

## NEW57.2$experiment <- "E5"
## NEW57.2$infection <- "primary"

## NEW57.2 <- NEW57.2[, c("EH_ID", "labelsNew", "experiment", "weight",
##                        "weight_dpi0", "relativeWeight", "fecweight",
##                        "dpi", "infection")]

## colnames(NEW57.2) <- colnames(NEW57.1)

## NEW57 <- rbind(NEW57.2, NEW57.1)

## dupes <- NEW57$labels[duplicated(NEW57$labels)]

## foo <- NEW57[NEW57$labels%in%dupes, ]

## foo$EH_ID <- NULL

## dim(foo)
## dim(unique(foo))

## NEW57 <- NEW57[!duplicated(NEW57$labels), ]

## write.csv(NEW57, "data/Experiment_results/E57_xxxxx_Eim_record.csv",
##           row.names=FALSE)


## for this just one feces weight is ZERO
ALL[ALL$EH_ID%in%"LM0236",]

## for this just four feces weights are ZERO
ALL[ALL$EH_ID%in%"LM0274",]

## they are all from primary, let's ask Alice
table(ALL[which(ALL$feces_weight==0), "infection"])
