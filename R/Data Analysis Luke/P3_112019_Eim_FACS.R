# P3 FACS
library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(ggplot2)
library(flowCore)
library(flowCytBioc)
library(ggcyto)

Afs <- read.flowSet(path = "~/Eimeria_Lab/data/3_recordingTables/P3b_112019_Eim_FACS/",
                    package = "flowCore", pattern = "a(LN)") 
markernames(Afs)
sampleNames(Afs)

Pfs <- read.flowSet(path = "~/Eimeria_Lab/data/3_recordingTables/P3b_112019_Eim_FACS/",
                    package = "flowCore", pattern = "p(LN)")
markernames(Pfs)
sampleNames(Pfs)

Sfs <- read.flowSet(path = "~/Eimeria_Lab/data/3_recordingTables/P3b_112019_Eim_FACS/",
                    package = "flowCore", pattern = "spleen")
markernames(Sfs)
sampleNames(Sfs)


autoplot(Afs, "CD4")