Eimeria_Lab repository

This repository is for storage of clean data from experiments conducted at AG Heitlinger. 
# Structure:
## [data](https://github.com/derele/Eimeria_Lab/tree/master/data) = contains all cleaned up data generated during our experiments and templates for tables
### [Experimental_design](https://github.com/derele/Eimeria_Lab/tree/master/data/Experimental_design) = mouse information sheets containing attributes such as: sex, strain, date of birth, EH_ID, InfectionStrain.
### [Experiment_results](https://github.com/derele/Eimeria_Lab/tree/master/data/Experiment_results) = clean tables of results obtained from a given experiment and assay/observation
### [Templates](https://github.com/derele/Eimeria_Lab/tree/master/data/templates) = examples of what corresponding tables should look like
#### [mouse_paperwork_mandatory](https://github.com/derele/Eimeria_Lab/tree/master/templates/mouse_paperwork_mandatory) = This folder contains all the necessary files for setting up an infection experiment in our mouse facilities. There is a complete protocol for animal handling, sampling, euqipment and facility handling. Accompanied by a cage placement template to keep a full track of experimental setup and a !Score Sheet! for each mouse to keep track of animal health. This sheet is a legal requirement that we must fulfill as per our animal handling license.

## [data_access_code](https://github.com/derele/Eimeria_Lab/tree/master/data_access_code) = contains examples of R code related to accessing and reading information from the [data](https://github.com/derele/Eimeria_Lab/tree/master/data) folder.

## [data_creation_code](https://github.com/derele/Eimeria_Lab/tree/master/data_creation_code) = contains examples of R code related to processing raw data and making it suitable to be in the [data](https://github.com/derele/Eimeria_Lab/tree/master/data) folder.

# 1. Accessing data:
## 1.1. General description:

All data in this repository has been processed and saved as a clean table according to the corresponding [template](https://github.com/derele/Eimeria_Lab/tree/master/data/templates)

All file names contain information to distinguish the number of an experiment, date of infection (MMYYYY), infection agent (e.g. Eimeria, Crypto, etc.)and a format.
In addition, elements of the name may designate: design (outline of mouse information, infection and label), [tissue type](https://github.com/derele/Eimeria_Lab/blob/master/Tissue_labels.csv), oocyst, record (record of mouse weight over dpi) or assay type (ELISA, qPCR, FACS, etc.). 

Exmaple: E7_112018_Eim_CEWE_ELISA.csv
This means the table contains information generated from Experiment 7, infections started in November 2018, mice were infected with Eimeria, the tissue used in the assay was Caecum, the essay was ELISA and the table is in a .csv format.

## 1.2. Examples:
### 1.2.1. Example 1 (extracting maximum weightloss and maximum oocysts shed for the E139 and E64 strains across all experiments where these strains are present)
```r
# load libraries
library(RCurl)

# load 
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

The column names are preordained by the "makeDesignTable.R" and "makeRecordTable.R".
These should contain the following columns for "design":
### EH_ID (unique mouse identifier), mouse_strain (NMRI, SWISS, PWD, BUSNA, etc.) and experiment (unique experiment identifier). + any other available information about the mice
primary_infection (Eimeria strain), challenge_infection (if reinfected), infection_history (if reinfected)

For "record":
### EH_ID (unique mouse identifier), labels (unique timepoint identifier), weight (g), weight_dpi0 (weight on day of infection), weightloss, relative_weight (percentage change in weight from weight_dpi0), feces_weight and dpi (days post infection). 

For "oocyst":
### label (unique timepoint identifier), experiment (unique experiment identifier), oocyst_sq1, oocys_sq2, oocyst_sq3, oocyst_sq4, oocyst_mean, OPG, dilution.

and these for infection intensity qPCRs:
### label, EH_ID, delta, dpi, Eim_MC (melting curve (positive or negative)), Amp (is amplification good?)

## 2.2. Examples:
### 2.2.1. Adding genotype data

### 2.2.2. Adding qPCR data


DO NOT edit the data in these folders. You can clone the repos, work with the data and generate new tables as you 
please. The envisioned structure is that each manuscript has it's own repository. To avoid hardcoding, it is recommended to load tables via GitHub raw links using the "RCurl" package. See 1.2.1 Example

### 2. Experiment design
Fill the experiment design with the help of the information table,
with the help of R function "makeDesignTable.R"

ex:

```r

myDesignTable2 <- makeDesignTable(myseed = 1234 ,
                                  pathToInfoTable = "../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv",
                                  firstEH_Id = "LM0145")

# Separate equally between Mouse_strains
library(experiment)

expe <- randomize(data = myDesignTable2, group = c("E64", "E139"),
                   indx = myDesignTable2$EH_id, block = myDesignTable2$Strain)

trt <- data.frame(infection_isolate = expe$treatment)
trt$EH_id <- rownames(trt)
rownames(trt) <- NULL

# Create design table
designTable <- merge(trt, myDesignTable2)
print(table(designTable$Sex, designTable$Strain, designTable$infection_isolate))

write.csv(designTable,     "../data/2_designTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_DESIGN.csv", row.names = F)
```



**makeDesignTable.R** is a general function taking as input an *INFO.csv* data.frame,
and creating a DESIGN.csv one

**selectBasedOnAge.R** was used in April 2018 to split our groups in 2. Hard coded.
