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
E139W <- lapply(OV[OV$E139, "weight"], loadFromGH)

## URL for an empty cell in the table does not exist, that's fine, all
## others are read. Remove the empty element in the list
E139W <- E139W[!unlist(lapply(E139W, is.null))]

## Same for shedding
E139Shed <- lapply(OV[OV$E139, "shedding"], loadFromGH)
E139Shed <- E139Shed[!unlist(lapply(E139Shed, is.null))]

W139colnames <- Reduce(intersect, lapply(E139W, colnames))
S139colnames <- Reduce(intersect, lapply(E139Shed, colnames))

E139W <- lapply(E139W, "[", W139colnames)
E139Shed <- lapply(E139Shed, "[", S139colnames)

W139 <- Reduce(rbind, E139W)
S139 <- Reduce(rbind, E139Shed)

### calculating max weight loss for each mouse
as_tibble(W139) %>%
    group_by(EH_ID) %>%
    slice_min(n=1, order_by=weight) %>%
    left_join(OV[, c("Experiment", "Date")],
              by= c("experiment"="Experiment")) %>%
    group_by(experiment) %>%
    mutate(meanMaxWLat = mean(dpi)) %>%
    mutate(meanMaxWL = mean(1-(weight/weight_dpi0))*100) %>%
    select(experiment, Date, meanMaxWLat, meanMaxWL) %>%
    slice_head(n=1) %>%
    mutate(Date=as.Date(paste0("01/", Date), format="%d/%m/%Y")) ->
    maxWL

pdf("example_fig/E139_evol_maxWL.pdf", width=6, height=5)
ggplot(maxWL, aes(Date, meanMaxWLat, size= meanMaxWL)) +
    geom_point(color="red") +
    scale_y_continuous("maximal weightloss at dpi (average)") +
    scale_x_date("date of the experiment") +
    scale_size_continuous("mean maximal\nweightloss in %")
dev.off()

