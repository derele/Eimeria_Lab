# E7 ELISA FEC cleanup
library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(drc)
library(data.table)

############################################## first plate
# standards 
# LABELS are missing
E1_std <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA1_std.csv"
E1_std <- read.csv(text = getURL(E1_std))
# samples
E1_samples <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA1_samples.csv"
E1_samples <- read.csv(text = getURL(E1_samples))

###### use drc to construct standard curve and pinpoint protein content

model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E1_std)
plot(model1)

E1<-ED(model1, E1_samples$OD, type="absolute", display=F)
row.names(E1) <- E1_samples$label

points(y=E1_samples$OD,x=E1[,1],col="lightblue",pch=19,cex=2)
text(y =E1_samples$OD, x = E1[,1], labels=E1_samples$label, data=E1, cex=0.9, font=2)

E1 <- data.frame(E1)
colnames(E1)[1] <- "IFNy"
E1 <- dplyr::select(E1, IFNy)
setDT(E1, keep.rownames = TRUE)[]
colnames(E1)[1] <- "labels"
# the 0 replacement is just temporary
E1[E1=="NaN"]<- 0
# write out
write.csv(E1, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA1_clean.csv")
write.csv(E1, "./Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA1_clean.csv")
########################################################### second plate
# standards
E2_std <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA2_std.csv"
E2_std <- read.csv(text = getURL(E2_std))
# samples
E2_samples <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA2_samples.csv"
E2_samples <- read.csv(text = getURL(E2_samples))

###### use drc to construct standard curve and pinpoint protein content

model2<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E2_std)
plot(model2)

E2<-ED(model2, E2_samples$OD, type="absolute", display=F)
row.names(E2) <- E2_samples$label

points(y=E2_samples$OD,x=E2[,1],col="lightblue",pch=19,cex=2)
text(y =E2_samples$OD, x = E2[,1], labels=E2_samples$label, data=E2, cex=0.9, font=2)

E2 <- data.frame(E2)
colnames(E2)[1] <- "IFNy"
E2 <- dplyr::select(E2, IFNy)
setDT(E2, keep.rownames = TRUE)[]
colnames(E2)[1] <- "labels"
# the 0 replacement is just temporary
E2[E2=="NaN"]<- 0

write.csv(E2, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA2_clean.csv")
write.csv(E2, "./Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISA2_clean.csv")

############################ merge datasets
ELISA_complete <- rbind(E1, E2)
#check for duplicates
ELISA_complete <- filter(ELISA_complete)
# write out
write.csv(ELISA_complete, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISAs_complete.csv")
write.csv(ELISA_complete, "./Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FEC_ELISAs/E7_112018_Eim_FEC_ELISAs_complete.csv")
