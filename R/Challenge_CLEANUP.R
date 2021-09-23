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

ChallengeEx  <- c("E5", "E7", "E10", "E11")
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


head(W[["E11"]])
tail(W[["E11"]])

head(O[["E11"]])

head(O[["E10"]])

tail(O[["E5"]])


O[["E11"]]$labels <- paste0("E11", O[["E11"]]$batch, O[["E11"]]$labels)

## here we have a complication: E7 a and b batches are really batches
## not primary and challenge infection. We'll code this in x/y to
## avoid confusion with a and b primary/challenge (wrongly called
## "batches" in E10 and E11. 

O[["E5"]]$labels <- paste0("E57", "a", O[["E5"]]$labels)

O[["E7"]]$labels <- gsub("E7a", "E7x", O[["E7"]]$labels)
O[["E7"]]$labels <- gsub("E7b", "E7y", O[["E7"]]$labels)

O[["E7"]]$labels <- gsub("E7", "E57b", O[["E7"]]$labels)

oocyst.Cols <- intersect(colnames(O[["E5"]]), colnames(O[["E7"]]))

## Now we can concatenate those
OE57 <- rbind(O[["E5"]][,oocyst.Cols], O[["E7"]][,oocyst.Cols])
OE57$dilution <- as.numeric(gsub(",", ".", OE57$dilution))


## Now the weight data
W[["E5"]]$labels <- paste0("E57a", W[["E5"]]$labels)


W[["E7"]]$labels <- gsub("E7a", "E7x", W[["E7"]]$labels)
W[["E7"]]$labels <- gsub("E7b", "E7y", W[["E7"]]$labels)
W[["E7"]]$labels <- gsub("E7", "E57b", W[["E7"]]$labels)

weight.Cols <- intersect(colnames(W[["E5"]]), colnames(W[["E7"]]))

## Now we can concatenate those
WO57 <- rbind(W[["E5"]][, weight.Cols], W[["E7"]][, weight.Cols])


DE57  <- merge(D[["E5"]], D[["E7"]], all.x=TRUE)

write.csv(DE57, "data/Experimental_design/E57_xxxxx_Eim_DESIGN.csv", row.names = FALSE)
write.csv(WO57, "data/Experiment_results/E57_xxxxx_Eim_record.csv", row.names = FALSE)
write.csv(OE57, "data/Experiment_results/E57_xxxxx_Eim_oocyst.csv", row.names = FALSE)


## c("Date", "E64", "E88", "E139", "Eflab", "UNI", "no_of_mice", 
## "CLS", "WDS", "Experiment", "weight", "shedding", "infection_intensity", 
## "CEWE_ELISA", "MES_ELISA", "gene_expression", "FACS", "design"
## )

