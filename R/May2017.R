#########################
## Info about experiment:
## 25 mice (LM_0104 was deleted as she died at dpi1)
## 272 mice.days data points
#########################
library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

######### I Experiment design:
# Info on the mouse genotype AND on the infection strain:
ExpePlanDF <- read.csv("https://raw.githubusercontent.com/alicebalard/Eimeria_Lab/master/data_raw/Experiment_Table_raw_Ploen_May2017.csv")

# Print the expe design:
as.data.frame.matrix(table(ExpePlanDF$strain, ExpePlanDF$Inf_strain))

# Load data:
mydata_oocysts <- read.csv("https://raw.githubusercontent.com/alicebalard/Eimeria_Lab/master/data_raw/Oocyst_coding_PloenEx_May2017.csv")

# Remove LM_0104:
mydata_oocysts <- mydata_oocysts[which(mydata_oocysts$X != "LM0104"),]

# Check:
table(is.na(mydata_oocysts[grep("raw.count", names(mydata_oocysts))]))
table(is.na(mydata_oocysts[grep("fec.weight", names(mydata_oocysts))]))
## Error corrected (same amount of fecal weight and oocysts counts)
#######################################################################################!!!!!!!!!!!!!!!!!!!!!!!!!

mydata_weight <- read.csv("https://raw.githubusercontent.com/alicebalard/Eimeria_Lab/master/data_raw/mouse%20weight.csv")

# Weight of the feces separated for PCR:
subsampleDF <- read.csv("https://raw.githubusercontent.com/alicebalard/Eimeria_Lab/master/data_raw/Feces_kept_separated_for_PCR.csv")
subsampleDF <- na.omit(subsampleDF)

# Calculate how many animals at start:
summary_table_at_dpi <- function(n){
  A = data.frame(mouse_number = ExpePlanDF$EH_id, ExpePlanDF$strain, ExpePlanDF$Inf_strain)
  B = data.frame(mouse_number = mydata_weight$mouse_no., mydata_weight[paste0("mouse_weight.g._dpi", n)])
  C = merge(A, B, by = "mouse_number")
  C = na.omit(C)
  table(C$ExpePlanDF.strain, C$ExpePlanDF.Inf_strain)
}

summary_table_at_dpi(11)
 
######### II Oocysts shedding:

# Initialise Newdata
Newdata <- data.frame(NULL)

# Run on all the dpi in 1 go!
for (DPI in 0:11){
  mydata_oocyststemp <- cbind(mydata_oocysts$X,mydata_oocysts[, grep(paste0("dpi.", paste(DPI, "", sep="\\.")), colnames(mydata_oocysts))])
  names(mydata_oocyststemp)<- c("mouse.ID","label","fec.weight","dil", "raw.count")
  # add column
  mydata_oocyststemp$dpi<- DPI
  Newdata <- rbind(Newdata, mydata_oocyststemp)
}

# Remove useless objects:
rm(mydata_oocyststemp, DPI)

# Rm the dil column, useless:
Newdata <- Newdata[!names(Newdata) %in% "dil"]

# Rm empty lines (when the mice died, no more counting) :
Newdata <- na.omit(Newdata)

# Merge the 2 dataframe by labels:
names(subsampleDF)[names(subsampleDF) %in% "Label"] <- "label"

Oo_Df <- na.omit(merge(Newdata, subsampleDF, by = "label"))

## Calculate the oocysts in the subsampling:
Oo_Df$Total_feces <- Oo_Df$subsample_g + Oo_Df$fec.weight

# Calculate the oocysts number in 1 mL:
Oo_Df$oocysts.per.tube <- Oo_Df$raw.count * 10000

# Calculate the oocyst per gram in the processed part of the feces::
Oo_Df$oocysts.per.g <- Oo_Df$oocysts.per.tube / Oo_Df$fec.weight

#Put my dataframe at the good format:
Oo_Df <- Oo_Df[c("label", "mouse.ID", "dpi", "oocysts.per.g", "oocysts.per.tube")]

# Add the info on the mouse genotype AND on the infection strain:
names(Oo_Df)[names(Oo_Df) %in% "mouse.ID"] <- "EH_id"

Oo_Df <- merge(Oo_Df, ExpePlanDF, by = "EH_id")

# Visual check-up:
table(Oo_Df$EH_id, Oo_Df$dpi)

######### III. Weight of the mice:
mydata_weight_long <- colnames(mydata_weight) <- gsub("mouse_weight.g.", "weight", colnames(mydata_weight))
mydata_weight_long <- melt(mydata_weight, id=c("mouse_no.", "Infection_date"))  # convert to long format

colnames(mydata_weight_long)[colnames(mydata_weight_long)%in%c("mouse_no.", "variable", "value")] <-
  c("LM_id", "dpi", "weight")

mydata_weight_long$dpi <- gsub("weight_", "", mydata_weight_long$dpi)

# Now let's do it relative to the weight of D1:
mydata_weight_relative <- data.frame(LM_id=mydata_weight[,1], mydata_weight[ ,3:14]/mydata_weight[,4]*100)
colnames(mydata_weight_relative) <- gsub("mouse_weight.g.", "rel_w", colnames(mydata_weight_relative))
mydata_weight_relative_long <- melt(mydata_weight_relative, id=c("LM_id"))  # convert to long format

colnames(mydata_weight_relative_long)[colnames(mydata_weight_relative_long)%in%c("variable", "value")] <-
  c("dpi", "rel.weight")

mydata_weight_relative_long$dpi <- gsub("weight_", "", mydata_weight_relative_long$dpi)

Alldata <- merge(mydata_weight_long, mydata_weight_relative_long, all=TRUE)

Alldata <- merge(Alldata, ExpePlanDF, by.x="LM_id", by.y="EH_id", all=TRUE)

Alldata$dpi <- reorder(Alldata$dpi, as.numeric(gsub("dpi", "", Alldata$dpi)))

# Rename same mouse labels:
names(Alldata)[names(Alldata) %in% "LM_id"] <- "EH_id"

# Rename correct dpi:
Alldata$dpi <- gsub(pattern = "dpi", replacement = "", x = Alldata$dpi)

Total_Franci <- merge(Alldata, Oo_Df)

# write.csv(x = Total_Franci, file = "../data_clean/May2017_crossing_infection.csv", row.names = F)


## Extra for DNA extraction : tubes with maximum oocysts
library(dplyr)

tabhighoo <- function(inf){
  al <- inf[which((inf$oocysts.per.tube) %in% tail(sort(inf$oocysts.per.tube), 20)), c(1,2,10,15,20,21,22), ]
  arrange(al, -oocysts.per.tube)
}

tabhighoo(Total_Franci[which(Total_Franci$Inf_strain == "EI64"),])
tabhighoo(Eflab <- Total_Franci[which(Total_Franci$Inf_strain == "Eflab"),])
