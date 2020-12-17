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

## create a list of dataframes for the weight data subsetting for only
## the experiments with the E139 used
W <- lapply(OV[, "weight"], loadFromGH)
W <- W[!unlist(lapply(W, is.null))]
### TODO: weight table missing for 4 experiments!!!(URL for an empty
### cell in the table does not exist)

## Same for shedding
O <- lapply(OV[, "shedding"], loadFromGH)
O <- O[!unlist(lapply(O, is.null))]
### TODO: shedding table missing for 4 experiments!!!(URL for an empty
### cell in the table does not exist)

## Same for design
D <- lapply(OV[, "design"], loadFromGH)
D <- D[!unlist(lapply(D, is.null))]
### TODO: design table missing for 3 experiments!!!(URL for an empty
### cell or NA in the table does not exist)

Wcolnames <- Reduce(intersect, lapply(W, colnames))
Wcolnames
## TODO: only EH_ID is common to all weigth tables as a colname

Ocolnames <- Reduce(intersect, lapply(O, colnames))
## TODO: NOTHIN is common to all shedding tables as a colname

Dcolnames <- Reduce(intersect, lapply(D, colnames))
## TODO: only "EH_ID", "primary_infection" and "mouse_strain" are
## common for design tables


### WITHOUT common colnames this can't work so generally for now it
### could work to soem extend if we were not that greeding wanting to
### combine all datsets already!!! (see Example1_combine_E139_data.R
### for working only with datasets having this strain). 

DR <- lapply(D, "[", Dcolnames)
DA <- Reduce(rbind, DR)

### TODO: The data has a HUGE OMISSION that needs to be fixed asap: we
### need a column to make tell us whether the weight (and shedding) is
### for primary or secondary infection!!

### ONE WAY TO EASIELY do that: put a "challenge_infection" in every
### design table and set it NA in every design table. Then you can set
### the controls (no challenge) in those tables to NA too! Or
### challenge to "unifected/Unif". Just to be not confusing: also add
### a column measured_infection (with values "challenge" or "primary")


WR <- lapply(W, "[", Wcolnames)
WA <- Reduce(rbind, WR)

OR <- lapply(O, "[", Ocolnames)
OA <- Reduce(rbind, OR)

### TODO, well... MAKE THE ABOVE work by sanitizing the data!!
