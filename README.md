# Eimeria_Lab repository

### This repository is for storage of clean data from experiments conducted at AG Heitlinger. 

# Structure:
## [data](https://github.com/derele/Eimeria_Lab/tree/master/data) = contains all cleaned up data generated during our experiments and templates for tables
### [Experimental_design](https://github.com/derele/Eimeria_Lab/tree/master/data/Experimental_design) = mouse information sheets containing attributes such as: sex, strain, date of birth, EH_ID (unique identifier) and infection strain.
### [Experiment_results](https://github.com/derele/Eimeria_Lab/tree/master/data/Experiment_results) = clean tables of results obtained from a given experiment and assay/observation
### [Templates](https://github.com/derele/Eimeria_Lab/tree/master/data/templates) = examples of what corresponding tables should look like
#### [mouse_paperwork_mandatory](https://github.com/derele/Eimeria_Lab/tree/master/data/templates/mouse_paperwork_mandatory) = This folder contains all the necessary files for setting up an infection experiment in our mouse facilities. There is a complete protocol for animal handling, sampling, euqipment and facility handling. Accompanied by a cage placement template to keep a full track of experimental setup and a !Score Sheet! for each mouse to keep track of animal health. This sheet is a legal requirement that we must fulfill as per our animal handling license.

## [data_access_code](https://github.com/derele/Eimeria_Lab/tree/master/data_access_code) = contains examples of R code related to accessing and reading information from the [data](https://github.com/derele/Eimeria_Lab/tree/master/data) folder.

## [data_creation_code](https://github.com/derele/Eimeria_Lab/tree/master/data_creation_code) = contains examples of R code related to processing raw data and making it suitable to be in the [data](https://github.com/derele/Eimeria_Lab/tree/master/data) folder.

# 1. Accessing data:
### All clean and standardised tables are kept in an [overview table](https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv), which can be used to check information availability and retreive said tables. 


## 1.1. General description:

All data in this repository has been processed and saved as a clean table according to the corresponding [template](https://github.com/derele/Eimeria_Lab/tree/master/data/templates)

All file names contain information to distinguish the number of an experiment, date of infection (MMYYYY), infection agent (e.g. Eimeria, Crypto, etc.)and a format.
In addition, elements of the name may designate: design (outline of mouse information, infection and label), [tissue type](https://github.com/derele/Eimeria_Lab/blob/master/Tissue_labels.csv), oocyst, record (record of mouse weight over dpi) or assay type (ELISA, qPCR, FACS, etc.). 

Exmaple: E7_112018_Eim_CEWE_ELISA.csv
This means the table contains information generated from Experiment 7, infections started in November 2018, mice were infected with Eimeria, the tissue used in the assay was Caecum, the essay was ELISA and the table is in a .csv format.

## 1.2. Examples:
### 1.2.1. Example 1 getting mean maximal wightloss for E139 isolate over multiple experiments
```r
# load libraries
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

```
### 1.2.1. Example 2 

# 2. Adding data:
## 2.1. General description:
Each file should be named according to the tamplate of:
ExperimentNumber_MonthYear_Pathogen_Neccessary_Information_Not_Too_Long.format
E.g.: E1_012017_Eim_Caecum_cDNA.csv
E = Experiment, P = Passaging

Raw data should be stored here, processed using code saved here as well and both should be subsequently deleted once a clean table exists. The raw data and code should be both commited and pushed to git to keep track of events. Commit messages should contain information on what files are being handled.

General rule is:
1. upload raw data table
2. upload code to process raw data table
3. upload clean data table
4. delete raw data table and code

The column names are preordained by the [create_design_table](https://raw.githubusercontent.com/derele/Eimeria_Lab/master/R/create_design_table.R), [create_oocyst_table](https://raw.githubusercontent.com/derele/Eimeria_Lab/master/R/create_oocyst_table.R) and [create_record_table](https://raw.githubusercontent.com/derele/Eimeria_Lab/master/R/create_record_table.R) scripts.
## These should contain the following columns for "design":
### EH_ID (unique mouse identifier), mouse_strain (NMRI, SWISS, PWD, BUSNA, etc.), primary_infection (Eimeria strain) and experiment (unique experiment identifier). + any other available information about the mice challenge_infection (if reinfected), infection_history (if reinfected)

## For "record":
### experiment(unique experiment identifier), EH_ID (unique mouse identifier), labels (unique mouse-timepoint identifier), weight (g), weight_dpi0 (weight on day of infection), relative_weight (percentage change in weight from weight_dpi0), feces_weight, dpi (days post infection) and dpi_dissection (dpi at which mouse was dissected). 

## For "oocyst":
### labels (unique mouse-timepoint identifier), experiment (unique experiment identifier), oocyst_sq1, oocys_sq2, oocyst_sq3, oocyst_sq4 (4 Neubauer counting chambers) and dilution.

## For infection intensity qPCRs:
### labels, EH_ID, Target, Ct (cycle threshold), Ct_SD (standard deviation within sample duplicates/triplicates), dissection_dpi, Eim_MC (melting curve (positive or negative)), Amp (is amplification good?)

## 2.2. Examples:
### 2.2.1. Creating a design table for experiment E10
```r

# create design table
# Load information table
library(RCurl)
# load in initial dataset from GitHub (must be raw.)
infoTable = read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_creation_code/E10_112020_Eim_INFO.csv")
# check last experiment and get highest EH_ID
lastEH_ID <- "LM0399"
# divide dataset into groups as desired
## e.g.: 100 mice would divide into c(rep("isolate1", 25), rep("isolate2", 25), 
##                                  rep("isolate3", 25), rep("uninfected", 25))
# for challenge infection plans this has to be spread further
primary_infection <- c(rep("E64", 11), rep("E88", 11), rep("UNI", 10))
challenge_infection <- c(rep("E64", 11), rep("E88", 11), rep("UNI", 10))
# number of mice
Nmice = nrow(infoTable)
#Give EH_IDs
num = as.numeric(sub("LM", "", lastEH_ID))
num = num + (1:(Nmice))
EH_ID = paste0("LM", sprintf("%04d", num))
#Assign infection isolate
designTable <- data.frame(primary_infection = primary_infection,
                          challenge_infection = challenge_infection,
                          EH_ID= EH_ID)
# Spread names and infections randomly among mice (restrict by total amount of mice)
infoTable$EH_ID <- sample(EH_ID, size = 32)
infoTable$primary_infection <- sample(primary_infection, size = 32)
infoTable$challenge_infection <- sample(challenge_infection, size = 32)
# merge
finaldesignTable <- merge(infoTable, designTable, all.x = T)
####### check necessary columns at https://github.com/derele/Eimeria_Lab
# add experiment column
finaldesignTable$experiment <- "E10"
# rename columns to match other design tables as stated in the repo
names(finaldesignTable)[names(finaldesignTable) == "Strain"] <- "mouse_strain"
finaldesignTable$infection_history <- paste(finaldesignTable$primary_infection, finaldesignTable$challenge_infection, sep = ":")


# write out
write.csv(finaldesignTable,
          "~/GitHub/Eimeria_Lab/data_creation_code/design_table_creation_example.csv",
          row.names = F, quote = F )
```

### 2.2.2. Adding qPCR data


DO NOT edit the data in these folders. You can clone the repos, work with the data and generate new tables as you please. The envisioned structure is that each manuscript has it's own repository. To avoid hardcoding, it is recommended to load tables via GitHub raw links using the "RCurl" package. See 1.2.1 Example
