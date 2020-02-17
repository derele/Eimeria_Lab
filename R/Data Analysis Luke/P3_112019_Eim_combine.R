# script for combining all clean datasets into one
library(httr)
library(Rmisc)
library(RCurl)

# load in weight and oocysts
P3_weight <- "https://raw.githubusercontent.com/LubomirBednar/Manuscript_1/master/clean_data/P3_112019_Eim_weight_complete.csv"
P3_weight <- read.csv(text = getURL(P3_weight))
# load in ELISAs FEC
P3_FEC_ELISA <- "https://raw.githubusercontent.com/LubomirBednar/Manuscript_1/master/clean_data/P3_112019_Eim_feces_ELISA1_complete.csv"
P3_FEC_ELISA <- read.csv(text = getURL(P3_fecal_ELISA))
#load in ELISAs CEWE
P3_CEWE_ELISA <- ""
P3_CEWE_ELISA <- read.csv(text = getURL(P3_CEWE_ELISA))
# load in RT-qPCRs
P3_RTqPCR_CEWE <- "https://raw.githubusercontent.com/LubomirBednar/Manuscript_1/master/clean_data/P3_112019_Eim_RT-qPCR_complete.csv"
P3_RTqPCR_CEWE <- read.csv(text = getURL(P3_RTqPCR_CEWE))

