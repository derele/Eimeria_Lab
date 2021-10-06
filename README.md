# Eimeria_Lab repository

This repository is for the storage of data from experiments conducted
in the research group of Emanuel Heitlinger.

## [data_products](https://github.com/derele/Eimeria_Lab/tree/master/data_products)

This older contains clened data sets (more or less, discuss, e.g. in
issues!) ready for use data sets for analyses.

### [data_products/Challenge_infections.csv](https://github.com/derele/Eimeria_Lab/tree/master/data_products/Challenge_infections.csv)

Contains data for challenge (repeated) infections performed between
2017 and 2019. The data product is structure in the folowing columns:

- EH_ID: the unique identifier of the mouse
- experiment: the experiment as numbered in the overview table
- mouse_strain: the strain (inbred or outbred) of the mouse
- primary_infection: The Eimeria strain used for the primary infection
- challenge_infection: The Eimeria strain used for the challenge infection
- infection_history: The resulting infection history
- labels: the unique label of the fecal sample at a particular dpi
- weight: the weight of the mouse at this dpi
- weight_dpi0: the weight at the day of infection
- relative_weight: the weight of the mouse at this dpi relative to the
  weight at dpi0
- feces_weight: the weight of the feces collected at this dpi
- dpi: days post infection at which samples and data in this row were taken
- infection: the infection (primary or challenge) this row/dpi
  corresponds to
- oocyst_sq1, oocyst_sq2, oocyst_sq3, oocyst_sq4: the raw values for
  squared during oocyst counting
- dilution: the amount of PBS the feces (with it's relative weight)
  was dissolved in
- OO4sq: the sum of oocysts in the four counting squares
- OOC: the overall number of oocysts in in the feces (of a particular
  weight) at this dpi
- infection_type: what kind of infection are we looking at (challenge
  or primary, homologous or heterologous immunization). This is
  differently coded to infection, as here UNI:E88 (first uninfected,
  then infected with E88) would count as a "primaryE88" infection


In order to summarize data by mouse and infection you will use code
similar to the below in your analysis.

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


# Accessing raw data (and compilling it into data_products)

Pointers to clean and standardised raw data are given in the [overview
table](https://raw.githubusercontent.com/derele/Eimeria_Lab/master/Eimeria_Lab_overview.csv).
This can be used to compile the latest raw data into data_products. An
example of how this can be done is available in
[Challenge.R](https://raw.githubusercontent.com/derele/Eimeria_Lab/master/R/Challenge.R)

### [data](https://github.com/derele/Eimeria_Lab/tree/master/data) 
Contains somewhat cleaned (but still "raw") data generated during our
experiments. The most relevant of these raw data sets are indicated in
the overview tables as explained above for further compilation into
data products.

#### [Experimental_design](https://github.com/derele/Eimeria_Lab/tree/master/data/Experimental_design) 
Mouse information sheets containing attributes such as: sex, strain,
date of birth, EH_ID (unique identifier) and infection strain.

#### [Experiment_results](https://github.com/derele/Eimeria_Lab/tree/master/data/Experiment_results) 
Clean tables of results obtained from a given experiment and
assay/observation



