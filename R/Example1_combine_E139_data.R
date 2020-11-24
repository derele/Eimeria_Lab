library(RCurl)

OV <- read.csv("./Eimeria_Lab_overview.csv")

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


Reduce(intersect, lapply(E139W, colnames))
## weight is the only common column name!  FIX this!!


Reduce(intersect, lapply(E139W[1:2], colnames))
### four column names common in the first two files

## Once the files all have the same (number of) columns (and names),
## you can do something like:

## Reduce(rbind, E139W)

