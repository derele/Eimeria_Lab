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
ExpePlanDF <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/1_informationTables/Exp001_May2017_crossing_infection_INFO.csv")

# Print the expe design:
as.data.frame.matrix(table(ExpePlanDF$strain, ExpePlanDF$Inf_strain))

# Load data:
mydata_oocysts <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/raw/Exp001_RAW_May2017_WildPloen_crossInfection_RECORDoocysts.csv")

# Remove LM_0104 (died before)
mydata_oocysts <- mydata_oocysts[which(mydata_oocysts$X != "LM0104"),]

# Check:
table(is.na(mydata_oocysts[grep("raw.count", names(mydata_oocysts))]))
table(is.na(mydata_oocysts[grep("fec.weight", names(mydata_oocysts))]))
## Error corrected (same amount of fecal weight and oocysts counts)
#######################################################################################!!!!!!!!!!!!!!!!!!!!!!!!!

mydata_weight <- read.csv("../../../data/3_recordingTables/raw/Expe001_RAW_mouse weight.csv")

# Calculate how many animals at start:
summary_table_at_dpi <- function(n){
  A = data.frame(mouse_number = ExpePlanDF$EH_id, ExpePlanDF$strain, ExpePlanDF$Inf_strain)
  B = data.frame(mouse_number = mydata_weight$mouse_no., mydata_weight[paste0("mouse_weight.g._dpi", n)])
  C = merge(A, B, by = "mouse_number")
  C = na.omit(C)
  table(C$ExpePlanDF.strain, C$ExpePlanDF.Inf_strain)
}
summary_table_at_dpi(1)
summary_table_at_dpi(11)
 
######### II Oocysts shedding:

# Initialise Newdata
Oo_Df <- data.frame(NULL)

# Run on all the dpi in 1 go!
for (DPI in 0:11){
  mydata_oocyststemp <- cbind(mydata_oocysts$X,mydata_oocysts[, grep(paste0("dpi.", paste(DPI, "", sep="\\.")), colnames(mydata_oocysts))])
  names(mydata_oocyststemp)<- c("mouse.ID","label","fec.weight","dil", "raw.count")
  # add column
  mydata_oocyststemp$dpi<- DPI
  Oo_Df <- rbind(Oo_Df, mydata_oocyststemp)
}

# Remove useless objects:
rm(mydata_oocyststemp, DPI)

# Rm the dil column, useless:
Oo_Df <- Oo_Df[!names(Oo_Df) %in% "dil"]

# Rm empty lines (when the mice died, no more counting) :
Oo_Df <- na.omit(Oo_Df)

# Correct error
Oo_Df[Oo_Df$fec.weight > 10, "fec.weight"] <- Oo_Df[Oo_Df$fec.weight > 10, "fec.weight"] / 1000

# Calculate the oocysts number in 1 mL:
Oo_Df$oocysts.per.tube <- Oo_Df$raw.count * 10000

# Calculate the oocyst per gram in the processed part of the feces::
Oo_Df$oocysts.per.g <- Oo_Df$oocysts.per.tube / Oo_Df$fec.weight

#Put my dataframe at the good format:
Oo_Df <- Oo_Df[c("label", "mouse.ID", "dpi", "oocysts.per.g", "oocysts.per.tube", "fec.weight")]

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

## HERE careful change names if rerun
write.csv(x = Total_Franci,
          file = "../../../data/3_recordingTables/Exp001_May2017_crossing_infection.csv",
          row.names = F)
