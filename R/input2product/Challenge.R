### This could be an example of how to access particular subsets of
### the data
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tidyr)
library(plyr)
library(data.table)
library(psych)
library(tibble)

#select columns: 
basics <- c("EH_ID", "mouse_strain", "experiment", "primary_infection", 
            "challenge_infection", "labels", "dpi", "infection", "infection_history", "sex", "death")

weight_loss <- c("weight", "weight_dpi0", "relative_weight")

oocysts_counts <- c("feces_weight", "oocyst_sq1", "oocyst_sq2", "oocyst_sq3",
                    "oocyst_sq4", "dilution", "OOC", "OO4sq", "OPG_O")

intensity_qPCR <- c("Eim_MC", "delta")

cewe_elisa <- "IFNy_CEWE"

mes_elisa <- "IFNy_MES"

gene_expr <- c("CXCR3", "IRG6", "IL.12")

CellCount.cols <- c("Position", "CD4", "Treg", "Div_Treg", "Treg17", "Th1", "Div_Th1", "Th17", 
                    "Div_Th17", "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8",  
                    "Treg_prop", "IL17A_CD4", "batch")  

IFNy_faeces <- c("IFNy_FEC", "Caecum")

Gene.Exp.cols   <- c("IFNy", "IL.12", "IRG6", "CXCR3", "IL.6", "GBP2",
                     "IL.10", "IL.13", "IL.10", "IL.13", "IL1RN",
                     "CXCR3", "CASP1", "CXCL9", 
                     "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
                     "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

#reading the overview table. In each row there is a link to the raw data for each experiment
OV <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv")

## ## Only the challenge experiments
#a list of the names of each experiments
# you can use it later, to select the experiments from ov 
ChallengeEx  <- c("E57", "E10", "E11", "P4", "P3")

## ## download and append the weigth tables
#lapply: applies a function to every element of the list
#we select the challenge experimetns and the weight columns 
#we apply to every element of the list the function read.csv
W <- lapply(OV[OV$Experiment%in%ChallengeEx, "weight"], read.csv)

#reduce: works on 1st and 2nd element, produces result and then uses the result 
#with the 3rd element and so on
#we apply in this case the function rbind
Weight <- Reduce(rbind, W)
rm(W)

## ## Same for shedding
O <- lapply(OV[OV$Experiment%in%ChallengeEx, "shedding"], read.csv)

Oocysts <- Reduce(rbind, O)
rm(O)

## ## some don't agree
table(Oocysts$labels%in%Weight$labels)
table(Weight$labels%in%Oocysts$labels)

## there are more in in the weights table which are not found in the
## oocysts though
Weight[!Weight$labels%in%Oocysts$labels, ]

## But if the weight is NA the mouse was dead
table(Weight[!Weight$labels%in%Oocysts$labels &
                 is.na(Weight$weight), "dpi"])
### confirmed by the late dpi of these!

Results <- merge(Weight, Oocysts, all=TRUE)

## IDs sometimes with "_" sometimes without
Results$EH_ID <- gsub("LM_", "LM", Results$EH_ID)

## For now we hav to exclude those E11aBMI which have no mouse IDs
Results <- Results[!is.na(Results$EH_ID), ]

## Same for design
D <- lapply(OV[OV$Experiment%in%ChallengeEx, "design"], read.csv)

Des.cols <- Reduce(intersect, lapply(D, colnames))

Design <- Reduce(rbind, lapply(D, "[", Des.cols))

rm(D)

## remove all whitespaces
Design %>%
    mutate(across(where(is.character), str_trim)) ->
    Design

## IDs sometimes with "_" sometimes without
Design$EH_ID <- gsub("LM_", "LM", Design$EH_ID)

#### NOW ALSO REMOVE EXPERIMENT column?
### Not for now
ALL <- merge(Design, Results, all=TRUE, by="EH_ID")
table(ALL$experiment.x, ALL$experiment.y)

## let's call experiment 5 and 7 57, like in all tables

Design$experiment[Design$experiment%in%c("E5", "E7")] <- "E57"
Results$experiment[Results$experiment%in%c("E5", "E7")] <- "E57"
ALL <- merge(Design, Results, all=TRUE)

### some mice don't have an infection history 
forgotten.mice <- unique(ALL[is.na(ALL$primary_infection), "EH_ID"])
## one mouse "LM0295" got lost?


Design$EH_ID[!Design$EH_ID%in%Results$EH_ID]
### Two samples with no result records
### [1] "LM0205" "LM0290"


## all the data has a mouse ID
table(is.na(ALL$EH_ID))

## but feces weights are ZERO T ( when the oocyst counts are NA ->
## the mouse was dead)!
table(ALL$feces_weight == 0)

##  only one prolbem remaining after asking Alice for better data...
ALL[which(ALL$feces_weight == 0 &
              !is.na(ALL$oocyst_sq2)), ]

ALL %>% filter(!is.na(challenge_infection)) %>%
    rowwise() %>% 
    dplyr::mutate(OO4sq = rowSums(across(starts_with("oocyst_")))) %>%
    ## 0.1µl per square -> *10.000 to scale up to ml
    mutate(OOC = (OO4sq/4*10000)/dilution) %>%
    ## we have ZEROS in feces weight (also when we counted oocysts) so we better don't
    ## calculate OPG for now but just max (see below)
    ## mutate(OPG=OOC/feces_weight) %>%
    ## also re-calculate relative weight, as this seems to have errors
    ## from a spreadsheet program (wtf!)
    mutate(relative_weight= weight/weight_dpi0*100) %>%
    ## also look at this for OPG above (by uncommenting)
    ## select(feces_weight, starts_with("oocyst_"), OO4sq, OOC) %>%
    ## look at this for controlling the weight calculation
    ## select(EH_ID, dpi, infection, weight, relative_weight) %>%
    ## print(n=40)
    ## the E88 innoculum in E57 challenge infection was "not
    ## working", these mice are basically unifected controls
    mutate(challenge_infection=ifelse(!experiment%in%"E57",
                                      challenge_infection,
                                      ifelse(challenge_infection%in%"E88", "UNI",
                                             challenge_infection))) %>%
    ## then correct the infection history
    mutate(infection_history=paste0(primary_infection, "_",
                                    challenge_infection)) ->
    ALL

#here is a function to join the following merged data frames to ALL
join_to_ALL <- function(x) {
    left_join(ALL, unique(x), by = c(intersect(colnames(x), colnames(ALL)))) 
}

#download and append the infection intensity tables (qPCR)
I <- lapply(OV[OV$Experiment%in%ChallengeEx, "infection_intensity"], read.csv)

#now combine the infection intensity tables
Intensity <- Reduce(rbind, I)
rm(I)

## Corrrect wrong IDs
Intensity$EH_ID <- gsub("LM_", "LM", Intensity$EH_ID)

## Melting curve sometimes positve sometimes TRUE, and the opposite
Intensity$Eim_MC <- gsub("pos", TRUE, Intensity$Eim_MC)
Intensity$Eim_MC <- gsub("neg", FALSE, Intensity$Eim_MC)

#Eim_MC type is character, change to logical
Intensity$Eim_MC <- as.logical(Intensity$Eim_MC)

ALL <- join_to_ALL(Intensity)



### Step: Join the CEWE_ELISA to our challenge infections file
## ## download and append the CEWE_ELISA tables
#lapply: applies a function to every element of the list
C <- OV[OV$Experiment %in% ChallengeEx, "CEWE_ELISA"] 

#I have to apply the read.csv to vector elements which contain the raw data, therefore
#I have to first select from the OV file the lines with actual link to the raw files
C <- lapply(C[c(1,2,5)], read.csv)

CEWE_ELISA <- Reduce(rbind, C)

rm(C)

#Next step clean the Mouse ID
## IDs sometimes with "_" sometimes without
CEWE_ELISA$EH_ID <- gsub("LM_", "LM", CEWE_ELISA$EH_ID)

#Change IFNy to IFNy_cewe to show origin of the measurement
CEWE_ELISA <- CEWE_ELISA %>% dplyr::rename(IFNy_CEWE = IFNy) 

#merge with ALL
ALL <- join_to_ALL(CEWE_ELISA)

##Joining data on ELISA from the mesenterial lymphnodes
#download and append the data from the ELISA's - Mesentrial Lymphnodes
#step by step, as there are empty spaces and this won't work
M <- OV[OV$Experiment%in%ChallengeEx, "MES_ELISA"]

#I have to apply the read.csv to vector elements which contain the raw data, therefore
#I have to first select from the OV file the lines with actual links to the raw files
#as there is only one row line with row data, lapply is not required here
MES_ELISA <- read.csv(M[[1]])
rm(M)

## Corrrect wrong IDs
MES_ELISA$EH_ID <- gsub("LM_", "LM", MES_ELISA$EH_ID)

#change IFNy to IFNy_cewe to show origin of the measurement
MES_ELISA <- MES_ELISA %>% dplyr::rename(IFNy_MES = IFNy)

#Now join the MES_ELISA to the ALL file
ALL <- join_to_ALL(MES_ELISA)

## Moving on to the gene expression data
#download and append the gene expression data
G <- OV[OV$Experiment%in%ChallengeEx, "gene_expression"]

G <- lapply(G[c(2,5)], read.csv)

#now combine the infection intensity tables
Gene_Expression <- Reduce(rbind, G)

rm(G)

## Corrrect wrong IDs
Gene_Expression$EH_ID <- gsub("LM_", "LM", Gene_Expression$EH_ID)

# Merge to ALL
ALL <- join_to_ALL(Gene_Expression)

#remove the duplicates from the "ALL" file
ALL <- unique(ALL)

# Add the gene expression data
## using part of Finn's code
#### 2016-2019 (Luke)
IFC1 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC1.csv")
IFC2 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC2.csv")
IFC3 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC3.csv")
IFC4 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC4.csv")
IFC5 <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experiment_results/IFC5.csv")

IFC  <- bind_rows(IFC1, IFC2, IFC3, IFC4, IFC5)

rm(IFC1, IFC2, IFC3, IFC4, IFC5)

# filter out AA
IFC <- IFC %>% filter(str_starts(EH_ID, "LM"))

## remove unsuccessful amplifications 
## == 999, equivalent to bad quality, should not be used
## Luke's Version:
IFC <- subset(IFC, IFC$Call != "Flag")

IFC <- IFC %>% dplyr::group_by(EH_ID, Target) %>% 
  dplyr::summarise(Ct = mean(Value)) 
 IFC <- distinct(IFC)

## separate data using the pivot_wider()
## turns Gene Expression Markers into individual columns (from Target),
## values taken (from Ct)
IFC <- pivot_wider(IFC, names_from = "Target", values_from = "Ct")

## We need to add Variance and n() again - if I add variance before pivoting, 
# it will lead to numerous lines, because we have numerous variances for the different markers, 
# same for the n()... either we create a variance variable for EACH marker individually 
# or we don't include it for the pivoting... 


## rename columns that have incorrect column names
colnames(IFC)[colnames(IFC)%in%"IL6"]   <- "IL.6"
colnames(IFC)[colnames(IFC)%in%"IL10"]  <- "IL.10"
colnames(IFC)[colnames(IFC)%in%"IL12A"]  <- "IL.12A"
colnames(IFC)[colnames(IFC)%in%"IL13"]  <- "IL.13"
colnames(IFC)[colnames(IFC)%in%"IL17A"]  <- "IL.17A"
colnames(IFC)[colnames(IFC)%in%"IFNG"]  <- "IFNy"

sapply(IFC, function(x) sum(is.na(x)))

IFC <- IFC %>%
  dplyr::rename(CXCR3_bio = CXCR3)

ALL <- join_to_ALL(IFC)

rm(IFC)

#What's next? 
#FACS!

#download and append the FACS data
F <- OV[OV$Experiment %in% ChallengeEx, "FACS"]

F <- lapply(F[c(1,2,4)], read.csv)
F <- list(F[[1]], F[[2]], F[[3]])
#F[[1]] <- F[[1]] %>% mutate(infection = "challenge") #P4
#F[[2]] <- F[[2]] %>% mutate(infection = "challenge") #E57
F[[3]] <- F[[3]] %>% dplyr::select(-c("X", "weight", "weight_dpi0", "relative_weight", "feces_weight", "infection", "labels")) #E11
F[[3]] <- unique(F[[3]])

#merge the different facs files together
FACS <- Reduce(bind_rows, F)

rm(F)

#rename the opg counts to opg_o = opg counts old (can compare it alter to see if they are the same with the all File)
FACS <- FACS %>% dplyr::rename(OPG_O = OPG)

## Corrrect wrong IDs
FACS$EH_ID <- gsub("LM_", "LM", FACS$EH_ID)
names(FACS) <- gsub(" ", "_", names(FACS))

#keep only the unique mice
FACS <- unique(FACS)

#remove unecessary columns from the raw data. All the unecessary columns with some NAs will mess 
#my next merges
FACS <- FACS %>% dplyr::select(-c("labels", "delta", "IFNy_CEWE", "CXCR3",
                           "IRG6", "IL.12")) 

#I assume all the mice in the facs data are in the challenge infection! 
#THis could potentially solve my problems with merging
FACS <- FACS %>% dplyr::mutate(infection = "challenge")

#how many are duplicated?
length(unique(FACS$EH_ID)) #85
sum(duplicated(FACS$EH_ID)) #33

#How many observations do I expect? Amount of mice + one time the duplicated ones
2951 + 33 # = 2984 observations

#let's make some summaries of our data
my_summary_data <- FACS %>%
  dplyr::group_by(EH_ID) %>%
  dplyr::summarise(Count = n())

my_summary_data <- my_summary_data %>%
  dplyr::group_by(Count) %>%
  dplyr::summarise(Count2 = n())

#indeed we have 33 duplicates of the EH_ID. This is to be expected, as for some mice 
#we had replicates, as we measured from two positions the Immune measures


#how many columns do I expect?
length(colnames(ALL)) #ALL has 77 columns
length(colnames(FACS)) #FACS has 24 columns
length(intersect(colnames(ALL), colnames(FACS))) #they have 4 common columns

77 + 24 - 4 #I expect 97 columns in my new merged file

ALL <- ALL %>% left_join(unique(FACS), by = intersect(colnames(ALL), colnames(FACS)), copy = FALSE)

length(colnames(ALL)) #indeed 97!


#check the primary and challenge infection columns
colnames(ALL)

#check the unique values of each column
unique(ALL$primary_infection)
#"E88"   "E64"   "Eflab" "E139"  "UNI"  
#check if Eflab ise E88 - checked!

#In column position - replace nas with Mln, as the nas are coming from the 
#mesenterial lymphnodes
ALL$Position[is.na(ALL$Position)] <- "mLN"

#Combined all of the qPCR data for the challenge infections in "Intensity"

#We want to create the variable death. The variable death will show when the mice died 
#To be expected: If the mouse was sacrificed during the primary infection 
#the tag will be prim_11 (dpi = 11)
#If in the other case the mouse was sacrificed on day 8 as planned in the challenge infections
#it will get the tag chal_8
#In this way we will know to which observation the experimental data should match
#according to the day the mice was sacrificed

# death = challenge / primary (did the mouse die during the primary or 
# during the challenge infection

# summarize data by mouse and infection 
# and then creat the variables:
#max_dpi = maximum dpi that the mouse that the mouse reached for each infection
# challenge or primara
# maximum oocysts for each infection type 
# maximum weight loss for each infection type 
mouse_infection <- as_tibble(ALL) %>%
  dplyr::group_by(EH_ID, infection) %>%
  dplyr::summarise(max_dpi = max(dpi, na.rm = TRUE),
                   max_OOC = max(OOC, na.rm=TRUE),
                   max_WL = min(relative_weight, na.rm=TRUE))


#death variable
# select mice that only appear once
# these mice have died in the primary infection
A <- as_tibble(mouse_infection) %>%
  filter(EH_ID %in% mouse_infection$EH_ID
         [which(!(duplicated(mouse_infection$EH_ID) |
                    duplicated(mouse_infection$EH_ID, fromLast = TRUE)))]) %>%
  #now create a new column that shows that the mice died in the primary infection
  mutate(death = "primary")

B <- mouse_infection[!(mouse_infection$infection == "primary"),] %>%
  mutate(death = "challenge")

di <- rbind(A, B) %>%
  select(-infection)


# join to ALL
ALL <- ALL %>% left_join(di, by = "EH_ID")

rm(A, B, mouse_infection, di)
#Remove column OPG_O with not checked old oocyst counts

#check the missing sex in the design

### Categorizing mouse strains by hybrid status
# Mouse_strains: (hybrid status)

# BUSNA_STRA: outbred hybrids
# STRA_BUSNA: outbred hybrids
# SCHUNT_SCHUNT: parental strains
# PWD_SCHUNT: outbred hybrids
# STRA_STRA: parental strains
# STRA_SCHUNT: outbred hybrids
# PWD_BUSNA: outbred hybrids
# BUSNA_PWD: outbred hybrids
# PWD_PWD: parental strains
# SCHUNT_PWD: outbred hybrids
# SCHUNT_STRA: outbred hybrids
# STRA_SCHUNT: outbred hybrids


ALL <- ALL %>% 
  group_by(EH_ID) %>%
  dplyr::mutate(hybrid_status = case_when(
    mouse_strain == "BUSNA_STRA" ~ "F1 hybrid",
    mouse_strain == "STRA_BUSNA" ~ "F1 hybrid",
    mouse_strain == "SCHUNT_SCHUNT" ~ "F0 M. m. domesticus",
    mouse_strain == "PWD_SCHUNT" ~ "F1 hybrid",
    mouse_strain == "STRA_STRA" ~ "F0 M. m. domesticus",
    mouse_strain == "STRA_SCHUNT" ~ "F1 M. m. domesticus",
    mouse_strain == "PWD_BUSNA" ~ "F1 M. m. musculus",
    mouse_strain == "BUSNA_PWD" ~ "F1 M. m. musculus",
    mouse_strain == "PWD_PWD" ~ "F0 M. m. musculus",
    mouse_strain == "SCHUNT_PWD" ~ "F1 hybrid",
    mouse_strain == "SCHUNT_STRA" ~ "F1 M. m. domesticus",
    mouse_strain == "STRA_SCHUNT" ~ "F1 M. m. domesticus",
    mouse_strain == "BUSNA_BUSNA" ~ "F0 M. m. musculus",
    mouse_strain == "NMRI" ~ "other",
    mouse_strain == "BUSNA_STRA" ~ "F1 hybrid",
    TRUE ~ ""))



write.csv(ALL,"data_products/Challenge_infections.csv", row.names=FALSE)

