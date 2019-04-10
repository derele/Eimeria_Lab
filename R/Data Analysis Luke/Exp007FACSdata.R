ANT <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
POS <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_posteriorMLN.csv")

?read.csv
#cleanup artifact columns
ANT$X.6 <- NULL
ANT$X.10 <- NULL
POS$X <- NULL
POS$X.6 <- NULL
POS$X.9 <- NULL
#cleanup artifact rows#
ANT <- ANT[-c(1, 2, 3, 31, 32, 33), ]
POS <- POS[-c(1, 2, 30, 31), ]
#name columns properly
names(ANT) <- c("Sample","ThCD4p","TcCD8p","Th1IFN-gp_in_CD4p","Th17IL-17Ap_in_CD4p","Tc1IFN-gp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_ K-67p_in_Foxp3p","RORgtp_in_Foxp3p","Th1T-betp_in_CD4pFoxp3n","Dividing_Ki-67p_in_T-betp","Th17RORgp_in_CD4pFoxp3n","DividingKi-67p_in_RORgtp")
names(POS) <- c("Sample","ThCD4p","TcCD8p","Th1IFN-gp_in_CD4p","Th17IL-17Ap_in_CD4p","Tc1IFN-gp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_ K-67p_in_Foxp3p","RORgtp_in_Foxp3p","Th1T-betp_in_CD4pFoxp3n","Dividing_Ki-67p_in_T-betp","Th17RORgp_in_CD4pFoxp3n","DividingKi-67p_in_RORgtp")

#set columns to numbers
library("dplyr")
library("magrittr")


str(ANT)

#plot
boxplot(split(ANT$Sample, ANT$`Th%CD4+`, ANT$`Tc%CD8+`))
library(ggplot2)
