## Our challenge infection data product had the problem (scope, that
## is) to only include mice that were challenge infected. Our quant2
## project also needs the mice from experiment E5 (primary infection
## in E57, in which E7 was the challenge infection) which were not
## challenge infected but sacrificed at the end of the primary. 

## reading the overview table
OV <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv")

## identify all the files we want
files <- OV[OV$Experiment%in%"E57", c("weight", "shedding", "design")]

## download them in a lapply loop
allTab <- lapply(files, read.csv)

## ## to check the tables from these
## lapply(allTab, head)

## ## fix a problem with the "experiment entries". 
allTab[[1]]$experiment <- gsub("E5|E7", "E57", allTab[[1]]$experiment)
allTab[[3]]$experiment <- gsub("E5|E7", "E57", allTab[[3]]$experiment)

## ## merge those files (the lazy way)
AT <- Reduce(merge, allTab)

lapply(allTab, nrow)
nrow(AT)
## -> Yes the merge loses exactly 300 samples!

### Which samples are those? My guess is they have different labels
### (or mouse IDs) in the shedding and weight data. 
setdiff(allTab[[1]]$labels, allTab[[2]]$labels)
## -> no this is not the case!

setdiff(allTab[[1]]$EH_ID, allTab[[3]]$EH_ID)
## oh, wow 25 mouse IDs don't overlap between weight data and
## experimental design

table(allTab[[1]]$EH_ID%in%allTab[[3]]$EH_ID)
## 300 samples (a number exactly matching the missing samples) have a
## mouse ID that's not in the design

## What are those?
unique(allTab[[1]]$EH_ID[!allTab[[1]]$EH_ID%in%allTab[[3]]$EH_ID])

table(allTab[[3]]$EH_ID%in%allTab[[1]]$EH_ID)
## the design has 27 mouse IDs that are not in the weight

## What are those?
allTab[[3]]$EH_ID[!allTab[[3]]$EH_ID%in%allTab[[1]]$EH_ID]

## hah, they just have some whitspces around them!!!
## We fix it and then do the merge again!!
allTab[[3]]$EH_ID <- gsub(" *", "", allTab[[3]]$EH_ID)

AT <- Reduce(merge, allTab)

lapply(allTab, nrow)
nrow(AT)
## fixed, sorry for missing this! I guess this created a lot of
## additional work...

## limit to ony primary infection
AT <- AT[AT$infection%in%"primary", ]

## order for better overview
AT <- AT[order(AT$EH_ID, AT$dpi), ]

## those are the E. falciformis infected the mice we (Dile) will be
## working with (first)
unique(AT[AT$primary_infection%in%"E88", "EH_ID"])

## here is a function to find a particular sample

gimmeSample <- function (x) { 
    AT[grepl(x, AT$labels), c("EH_ID", "labels", "dpi", "infection")]
}

gimmeSample("ACU")

## and we write a new data product
write.csv(AT, "data_products/Quant2_E57.csv")
