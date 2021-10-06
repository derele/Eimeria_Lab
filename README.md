# Eimeria_Lab repository

### This repository is for storage of clean data from experiments conducted at AG Heitlinger. 

# Structure:

## [data_products](https://github.com/derele/Eimeria_Lab/tree/master/data_products)

This older contains clened data sets (more or less, discuss, e.g. in
issues!) ready for use data sets for analyses. 

## [data_products/Challenge_infections.csv](https://github.com/derele/Eimeria_Lab/tree/master/data_products/Challenge_infections.csv)

Contains data for challenge (repeated) infections performed between
2017 and 2019. The data product is structure in the folowing columns:

- EH_ID
- experiment
- mouse_strain
- primary_infection
- challenge_infection
- infection_history
- labels
- weight 
- weight_dpi0
- relative_weight
- feces_weight 
- dpi 
- infection 
- oocyst_sq1, oocyst_sq2, oocyst_sq3, oocyst_sq4
- dilution
- OO4sq
- OOC
- infection_type


```{r }
as_tibble(CI) %>%
    group_by(EH_ID, infection) %>%
    summarize(max_OOC = max(OOC, na.rm=TRUE),
              max_WL = min(relative_weight, na.rm=TRUE),
              experiment = unique(experiment),
              mouse_strain= unique(mouse_strain),
              primary_infection=unique(primary_infection),
              challenge_infection=unique(challenge_infection),
              infection_history=unique(infection_history),
              infection_type=unique(infection_type),
              experiment=unique(experiment)) ->
    CIMouse
```

## [data](https://github.com/derele/Eimeria_Lab/tree/master/data) = contains all cleaned up data generated during our experiments and templates for tables

### [Experimental_design](https://github.com/derele/Eimeria_Lab/tree/master/data/Experimental_design) = mouse information sheets containing attributes such as: sex, strain, date of birth, EH_ID (unique identifier) and infection strain.
### [Experiment_results](https://github.com/derele/Eimeria_Lab/tree/master/data/Experiment_results) = clean tables of results obtained from a given experiment and assay/observation

# 1. Accessing data (and compilling it into data_products):
### Pointers to clean and standardised raw data are given in the [overview table](https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv). This  can be used to compile the latest raw data into data_products. 


