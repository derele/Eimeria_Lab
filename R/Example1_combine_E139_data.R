library(RCurl)
library(dplyr)
library(magrittr)

OV <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv"))

loadFromGH <- function(URL){
    if(url.exists(URL)){
        U <- getURL(URL)
        read.csv(text = U)
    } else {
        message("URL \"", URL, "\" does not exist")
    }
}

## create a list of dataframes for the weight data subsetting for only
## the experiments with the E139 used
E139W <- lapply(OV[OV$E139, "weight"], loadFromGH)

## URL for an empty cell in the table does not exist, that's fine, all
## others are read. Remove the empty element in the list
E139W <- E139W[!unlist(lapply(E139W, is.null))]

## Same for shedding
E139Shed <- lapply(OV[OV$E139, "shedding"], loadFromGH)
E139Shed <- E139Shed[!unlist(lapply(E139Shed, is.null))]

Reduce(intersect, lapply(E139W, colnames))
## "labels" is the only common column name!  FIX this!!

Reduce(intersect, lapply(E139Shed, colnames))
lapply(E139Shed, colnames)


## Once the files all have the same (number of) columns (and names),
## you can do something like:
#Reduce(rbind, E139W)

# until separate standardised files are made, this should work

keep <- c("EH_ID", "weight", "faces_weight", "relative_weight", "weight_dpi0")
weight <- lapply(E139W, function(x) subset(x, select = intersect(keep, colnames(x))))
