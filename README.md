Eimeria_Lab repository

This repository is for storage of clean data from experiments conducted at AG Heitlinger. 
# Structure:
## data :
### Experimental_design = mouse information sheets containing attributes such as: sex, strain, date of birth, EH_ID, InfectionStrain.
### Experiment_results = clean tables of results obtained from a given experiment and assay/observation
### Templates = examples of what corresponsing tables should look like
#### mouse_paperwork_mandatory = This folder contains all the necessary files for setting up an infection experiment in our
mouse facilities. There is a complete protocol for animal handling, sampling, euqipment and facility handling. Accompanied
by a cage placement template to keep a full track of experimental setup and a !Score Sheet! for each mouse to keep track of
animal health. This sheet is a legal requirement that we must fulfill as per our animal handling license.



# 1. Accessing data:
## 1.1. General description:

All data in this repository has been processed and saved as a clean table according to the corresponding template in 
https://github.com/derele/Eimeria_Lab/tree/master/data/Templates/

All file names contain information to distinguish the number of an experiment, date of infection (MMYYYY), infection agent (e.g. Eimeria, Crypto,
etc.), format.
In addition, elements of the name may designate: design (outline of mouse information, infection and label), info (mouse information sheet),
tissue type (see table ......), oocyst, record (record of mouse weight over dpi) or assay type (ELISA, qPCR, FACS, etc.). 

Exmaple: E7_112018_Eim_CEWE_ELISA.csv
This means the table contains information generated from Experiment 7, infections started in November 2018, mice were infected with Eimeria, the
tissue used in the assay was Caecum, the essay was ELISA and tthe table is in a .csv format.

## 1.2. Examples
### 1.2.1. Example 1 (tabulate number of mice per experiment)

### 1.2.1. Example 2 

# 2. Adding data:
## 2.1. General description:
Each file should be named according to the tamplate of:
ExperimentNumber_MonthYear_Pathogen_Neccessary_Information_Not_Too_Long.format
E.g.: E1_012017_Eim_Caecum_cDNA.csv
E = Experiment, P = Passaging

Raw data should be stored here, processed using code saved here as well and both should be subsequently deleted once a clean table exists.
The raw data and code should be both commited and pushed to git to keep track of events. Commit messages should contain information on what
files are being handled.

General rule is:
1. upload raw data table
2. upload code to process raw data table
3. upload clean data table
4. delete raw data table and code

## 2.2. Examples
### 2.2.1. Adding genotype data

### 2.2.2. Adding qPCR data








R = This folder contains data processing scripts for our experiments. Here you can find all the necessary code 
for how previous work was processed and help yourself to useful functions to make your analysis easier and !replicable!. 
In that spirit, as all scripts and data from the entire Eimeria_Lab are on GitHub, please write your script 
in a way to load the data from there as well. And same as before, name your scripts by the template + 
whatever they are made for.

	e.g.: 
	#load libraries for loading raw github files
	library(httr)
	library(RCurl)
	#read in cell counts (FACS) data
	cell.countsURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_cell_counts_processed.csv"
	cell.counts <- read.csv(text = getURL(cell.countsURL)) 
	................................................................
	#after processing and cleaning, write to the appropriate repo and Git it
	write.csv(E7, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACS_clean.csv", quote = FALSE)

.git = A folder for GitHub use when initializing new repo or cloning an existing one. I wouldn't touch it unless you know
what you're doing.

DO NOT edit the raw data in these folders. You can clone the repos, work with the data and generate new tables as you 
please. The envisioned structure is that each manuscript has it's own folder/repo, the raw data is loaded from the raw.
GitHub files and saved in the manuscript folder. Then all can be edited and analysed there.

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

### 3. Retrieving data



## Historic

### Exp001. Parental and F1 wild mice Ploen (x2) cross infection (Eflab, E64) (Francisca)

### Pass001 Nov 2017, passaging 4 isolates (Eflab, E88, E139, E64) in NMRI. 2 mice per cage

### Exp002. NMRI 4 strains passaging, March 2018

### Exp003 & Exp004 Parental wild mice (x4) cross infection (E64, E139) (Vivian) 

#### Exp003 First batch: 02/05/2018

#### Exp004 Second batch:tba

## Diverse R codes

**makeDesignTable.R** is a general function taking as input an *INFO.csv* data.frame,
and creating a DESIGN.csv one

**selectBasedOnAge.R** was used in April 2018 to split our groups in 2. Hard coded.
