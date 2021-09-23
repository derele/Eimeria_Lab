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

## Only the challenge experiments
ChallengeEx  <- c("E57", "E10", "E11")
ExpNames <- OV$Experiment[OV$Experiment%in%ChallengeEx]

W <- lapply(OV[OV$Experiment%in%ChallengeEx, "weight"], loadFromGH)
names(W) <- ExpNames


## Same for shedding
O <- lapply(OV[OV$Experiment%in%ChallengeEx, "shedding"], loadFromGH)
names(O) <- ExpNames

O[["E10"]] <- O[["E10"]][, !colnames(O[["E10"]])%in%"X"]

O[["E11"]]$labels <- paste0("E11", O[["E11"]]$batch, O[["E11"]]$labels)
O[["E11"]]$dilution <- 1

OO.cols <- Reduce(intersect, lapply(O, colnames))

O[["E10"]] <- O[["E10"]][, OO.cols]
O[["E11"]] <- O[["E11"]][, OO.cols]
O[["E57"]] <- O[["E57"]][, OO.cols]

write.csv(O[["E57"]], "data/Experiment_results/E57_xxxxx_Eim_oocyst.csv",
          row.names=FALSE)

write.csv(O[["E10"]], "data/Experiment_results/E10_112020_Eim_oocyst.csv",
          row.names=FALSE)

write.csv(O[["E11"]], "data/Experiment_results/E11_Oocyst_CSV.csv",
          row.names=FALSE)       

## Same for design
D <- lapply(OV[OV$Experiment%in%ChallengeEx, "design"], loadFromGH)
names(D) <- ExpNames

