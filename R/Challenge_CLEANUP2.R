### This could be an example of how to access particular subsets of
### the data
library(RCurl)
library(dplyr)
library(magrittr)
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

ChallengeEx  <- c("E57", "E10", "E11")
ExpNames <- OV$Experiment[OV$Experiment%in%ChallengeEx]

## create a list of dataframes for the weight data subsetting for only
## the experiments with the E139 used
W <- lapply(OV[OV$Experiment%in%ChallengeEx, "weight"], loadFromGH)
names(W) <- ExpNames

## Same for shedding
O <- lapply(OV[OV$Experiment%in%ChallengeEx, "shedding"], loadFromGH)
names(O) <- ExpNames

## Same for design
D <- lapply(OV[OV$Experiment%in%ChallengeEx, "design"], loadFromGH)
names(D) <- ExpNames


W[["E57"]]$infection <- ifelse(W[["E57"]]$experiment%in%"E5", "primary", "challenge")

W[["E11"]]$infection <- ifelse(W[["E11"]]$batch%in%"a", "primary", "challenge")
W[["E10"]]$infection <- ifelse(W[["E10"]]$batch%in%"a", "primary", "challenge")


W.cols <- c("EH_ID", "labels", "experiment", "weight", "weight_dpi0", "relative_weight",
            "feces_weight", "dpi", "infection")


write.csv(W[["E57"]], "data/Experiment_results/E57_xxxxx_Eim_record.csv")

write.csv(W[["E10"]], "data/Experiment_results/E10_112020_Eim_record.csv")

write.csv(W[["E11"]], "data/Experiment_results/E11_012021_Eim_record.csv")
