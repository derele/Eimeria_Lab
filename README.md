Record of infection experiments performed in AG Heitlinger

## Protocol

### 1.Information
Eimeria_Lab repository

This repository is for sorage of all !relevant material! on handling mice in our facilities and experiments conducted there. 
Each file should be named according to the tamplate of: ExperimentNumber_MonthYear_Pathogen_Neccessary_Information_Not_Too_Long.format
E.g.: E1_012017_Eim_Caecum_cDNA.csv
E = Experiment, P = Passaging

Folders:
data = Main folder containing 1_information Tables, 2_design Tables, 3_recording Tables
	
	1_infromationTables = Usually large tables containing infromation on mice when they come in. Must contain "Strain", "Birth       date", "Sex", "HybridStatus". 

	2_designTables = Tables made based on 1_informationTables, showing the plan for infection experiments 
	"InfectionStrain" and has to have added "EH_ID". Must at minimum contain columns from 1_informationTables specified above.

	3_recordingTables = A large folder containing record tables of any measurements taken during experiments or passaging. 
	Only raw data should be contained here (work in progress) and named consistently. If many repeated measurements are taken 
	(e.g. raw qPCR outputs), please make a folder dedicated to the raw data using the naming template above. However, 
	each file in the folder should have a unique name.

figures = A folder containing mostly R generated graphics of data processing, showing trends and distributions. 
Please label each file using the same template as everywhere else and don't forget to include Title, legend and axis names. 
Check redundancy of your own figures and delete/replace appropriately. Basic descriptive graphics are ok 
but a more focused one should be kept for the manuscript repos.

mouse_paperwork_mandatory = This folder contains all the necessary files for setting up an infection experiment in our mouse facilities. There is a complete protocol for animal handling, sampling, euqipment and facility handling. Accompanied by a cage placement template to keep a full track of experimental setup and a !Score Sheet! for each mouse to keep track of animal health. This sheet is a legal requirement that we must fulfill as per our animal handling license.

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

.git = A folder for GitHub use when initializing new repo or cloning an existing one. I wouldn't touch it unless you know what you're doing.

DO NOT edit the raw data in these folders. You can clone the repos, work with the data and generate new tables as you please. The envisioned structure is that each manuscript has it's own folder/repo, the raw data is loaded from the raw. GitHub files and saved in the manuscript folder. Then all can be edited and analysed there.

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
