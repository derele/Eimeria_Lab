# E7 FEC ELISA scleanup script
library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(drc)
library(data.table)

##### add clean tables
E1_std <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_P3_Eim_CEWE_ELISA1_std.csv"
E1_std <- read.csv(text = getURL(E1_std))

E1_samples <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_CEWE_ELISAs/E7_112018_Eim_CEWE_ELISA1_samples.csv"
E1_samples <- read.csv(text = getURL(E1_samples))

###### use drc to construct standard curve and pinpointprotein content

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
colnames(E1)[1] <- "label"
E1 <- merge(E1, E1_samples)
E1$OD <- NULL

write.csv(E1, "./Eimeria_Lab/data/3_recordingTables/E7_112019_Eim_CEWE_ELISAs/EP_112019_Eim_CEWE_ELISA1_complete.csv")
write.csv(E1, "C:/Users/Luke Bednar/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CEWE_ELISAs/E7_112018_Eim_CEWE_ELISA1_complete.csv")
