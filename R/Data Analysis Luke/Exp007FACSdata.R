#personal path
#ANT <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")

#IZW path
#ANT <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_posteriorMLN.csv")

#HU path
ANT <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
POS <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_posteriorMLN.csv")
#cleanup artifact columns (check before using, csv reads different now and then)#
ANT$X.6 <- NULL
ANT$X.10 <- NULL
POS$X.6 <- NULL
POS$X.9 <- NULL
POS$X <- NULL
#cleanup artifact rows#
ANT <- ANT[-c(1, 2, 3, 31, 32, 33), ]
POS <- POS[-c(1, 2, 30, 31), ]
#name columns properly (check before using, csv reads different now and then#
names(ANT) <- c("Sample","ThCD4p","TcCD8p","Th1IFN-gp_in_CD4p","Th17IL-17Ap_in_CD4p",
                "Tc1IFN-gp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_ K-67p_in_Foxp3p",
                "RORgtp_in_Foxp3p","Th1T-betp_in_CD4pFoxp3n","Dividing_Ki-67p_in_T-betp",
                "Th17RORgp_in_CD4pFoxp3n","DividingKi-67p_in_RORgtp")
names(POS) <- c("Sample","ThCD4p","TcCD8p","Th1IFN-gp_in_CD4p","Th17IL-17Ap_in_CD4p",
                "Tc1IFN-gp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_ K-67p_in_Foxp3p",
                "RORgtp_in_Foxp3p","Th1T-betp_in_CD4pFoxp3n","Dividing_Ki-67p_in_T-betp",
                "Th17RORgp_in_CD4pFoxp3n","DividingKi-67p_in_RORgtp")

#set columns to numbers
library("dplyr")
library("magrittr")
str(ANT)
#ANT
ANT$ThCD4p <- as.numeric(as.character(ANT$ThCD4p))
ANT$TcCD8p <- as.numeric(as.character(ANT$TcCD8p))
ANT$`Th1IFN-gp_in_CD4p` <- as.numeric(as.character(ANT$`Th1IFN-gp_in_CD4p`))
ANT$`Th17IL-17Ap_in_CD4p` <- as.numeric(as.character(ANT$`Th17IL-17Ap_in_CD4p`))
ANT$`Tc1IFN-gp_in_CD8p` <- as.numeric(as.character(ANT$`Tc1IFN-gp_in_CD8p`))
ANT$Treg_Foxp3_in_CD4p <- as.numeric(as.character(ANT$Treg_Foxp3_in_CD4p))
ANT$`Dividing_ K-67p_in_Foxp3p` <- as.numeric(as.character(ANT$`Dividing_ K-67p_in_Foxp3p`))
ANT$RORgtp_in_Foxp3p <- as.numeric(as.character(ANT$RORgtp_in_Foxp3p))
ANT$`Th1T-betp_in_CD4pFoxp3n` <- as.numeric(as.character(ANT$`Th1T-betp_in_CD4pFoxp3n`))
ANT$`Dividing_Ki-67p_in_T-betp` <- as.numeric(as.character(ANT$`Dividing_Ki-67p_in_T-betp`))
ANT$Th17RORgp_in_CD4pFoxp3n <- as.numeric(as.character(ANT$Th17RORgp_in_CD4pFoxp3n))
ANT$`DividingKi-67p_in_RORgtp` <- as.numeric(as.character(ANT$`DividingKi-67p_in_RORgtp`))
#POS
POS$ThCD4p <- as.numeric(as.character(POS$ThCD4p))
POS$TcCD8p <- as.numeric(as.character(POS$TcCD8p))
POS$`Th1IFN-gp_in_CD4p` <- as.numeric(as.character(POS$`Th1IFN-gp_in_CD4p`))
POS$`Th17IL-17Ap_in_CD4p` <- as.numeric(as.character(POS$`Th17IL-17Ap_in_CD4p`))
POS$`Tc1IFN-gp_in_CD8p` <- as.numeric(as.character(POS$`Tc1IFN-gp_in_CD8p`))
POS$Treg_Foxp3_in_CD4p <- as.numeric(as.character(POS$Treg_Foxp3_in_CD4p))
POS$`Dividing_ K-67p_in_Foxp3p` <- as.numeric(as.character(POS$`Dividing_ K-67p_in_Foxp3p`))
POS$RORgtp_in_Foxp3p <- as.numeric(as.character(POS$RORgtp_in_Foxp3p))
POS$`Th1T-betp_in_CD4pFoxp3n` <- as.numeric(as.character(POS$`Th1T-betp_in_CD4pFoxp3n`))
POS$`Dividing_Ki-67p_in_T-betp` <- as.numeric(as.character(POS$`Dividing_Ki-67p_in_T-betp`))
POS$Th17RORgp_in_CD4pFoxp3n <- as.numeric(as.character(POS$Th17RORgp_in_CD4pFoxp3n))
POS$`DividingKi-67p_in_RORgtp` <- as.numeric(as.character(POS$`DividingKi-67p_in_RORgtp`))
#check structure
str(ANT)
str(POS)
#mutate "Sample" to character
i <- sapply(ANT, is.factor)
ANT[i] <- lapply(ANT[i], as.character)

i <- sapply(POS, is.factor)
POS[i] <- lapply(POS[i], as.character)

#reshape df
library(reshape2)
cANT <- melt(ANT, id=c("Sample"), measure.vars = c("ThCD4p","TcCD8p","Th1IFN-gp_in_CD4p","Th17IL-17Ap_in_CD4p","Tc1IFN-gp_in_CD8p",
                                                   "Treg_Foxp3_in_CD4p","Dividing_ K-67p_in_Foxp3p","RORgtp_in_Foxp3p","Th1T-betp_in_CD4pFoxp3n",
                                                   "Dividing_Ki-67p_in_T-betp","Th17RORgp_in_CD4pFoxp3n","DividingKi-67p_in_RORgtp") )
cPOS <- melt(ANT, id=c("Sample"), measure.vars = c("ThCD4p","TcCD8p","Th1IFN-gp_in_CD4p","Th17IL-17Ap_in_CD4p","Tc1IFN-gp_in_CD8p",
                                                   "Treg_Foxp3_in_CD4p","Dividing_ K-67p_in_Foxp3p","RORgtp_in_Foxp3p","Th1T-betp_in_CD4pFoxp3n",
                                                   "Dividing_Ki-67p_in_T-betp","Th17RORgp_in_CD4pFoxp3n","DividingKi-67p_in_RORgtp") )

#reduce labels
library(stringr)
cANT$Sample <- gsub(cANT$Sample, pattern="Anterior", replacement='')
cANT$Sample <- gsub(cANT$Sample, pattern=".fcs", replacement='')
cPOS$Sample <- gsub(cPOS$Sample, pattern="Anterior", replacement='')
cPOS$Sample <- gsub(cPOS$Sample, pattern=".fcs", replacement='')

#plot
library(ggplot2)

ggplot(cANT, aes(Sample, value)) +   
  geom_bar(aes(fill = variable), position = "fill", stat="identity") +
  labs (title = "Anterior MLN", x = "Samples", y = "% of population") +
  coord_flip()

ggplot(cPOS, aes(Sample, value)) +   
  geom_bar(aes(fill = variable), position = "fill", stat="identity") +
  labs (title = "Posterior MLN", x = "Samples", y = "% of population") + 
  coord_flip()
