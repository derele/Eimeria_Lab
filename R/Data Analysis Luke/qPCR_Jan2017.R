#install.packages("RDML")
#install.packages("xml2")
end_point <- read.csv("../lubomir/Documents/Eimeria_Lab/Jan2017/Caecum_RT-qPCR/18.04.2019/FAM=SYBR/Luke_2019-04-18 15-46-21_CC011384 -  End Point Results_FAM.csv", sep = ";")
melt_deriv <- read.csv("../lubomir/Documents/Eimeria_Lab/Jan2017/Caecum_RT-qPCR/18.04.2019/FAM=SYBR/Luke_2019-04-18 15-46-21_CC011384 -  Melt Curve Derivative Results_FAM.csv", sep = ";")
melt_RFU <- read.csv("../lubomir/Documents/Eimeria_Lab/Jan2017/Caecum_RT-qPCR/18.04.2019/FAM=SYBR/Luke_2019-04-18 15-46-21_CC011384 -  Melt Curve RFU Results_FAM.csv", sep = ";")
quant_amp <- read.csv("../lubomir/Documents/Eimeria_Lab/Jan2017/Caecum_RT-qPCR/18.04.2019/FAM=SYBR/Luke_2019-04-18 15-46-21_CC011384 -  Quantification Amplification Results_FAM.csv", sep = ";")
#check structure
str(end_point)
str(melt_deriv)
str(melt_RFU)
str(quant_amp)
#convert factors to numeric, the "," from original .csv is the problem
end_point <- data.frame(lapply(end_point, function(x) {gsub(",", ".", x)}))
str(end_point)
end_point[,'End.RFU'] <- as.numeric(as.character(end_point[,'End.RFU']))

melt_deriv <- data.frame(lapply(melt_deriv, function(x) {gsub(",", ".", x)}))
str(melt_deriv)
melt_deriv[] <- lapply(melt_deriv, function(y) {
  if(is.factor(y)) as.numeric(as.character(y)) else y
})
sapply(melt_deriv, class)
str(melt_deriv)

melt_RFU <- data.frame(lapply(melt_RFU, function(x) {gsub(",", ".", x)}))
str(melt_RFU)
melt_RFU[] <- lapply(melt_RFU, function(y) {
  if(is.factor(y)) as.numeric(as.character(y)) else y
})
sapply(melt_RFU, class)
str(melt_RFU)

quant_amp <- data.frame(lapply(quant_amp, function(x) {gsub(",", ".", x)}))
str(quant_amp)
melt_deriv[] <- lapply(melt_deriv, function(y) {
  if(is.factor(y)) as.numeric(as.character(y)) else y
})
sapply(melt_deriv, class)
str(melt_deriv)

library(ggplot2)

