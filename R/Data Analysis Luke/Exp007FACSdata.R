#LOAD, CLEAN UP AND PROCESS DATA#
#PC path
ANT <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
POS <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_posteriorMLN.csv")

#laptop path
#ANT <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")

#IZW path
#ANT <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_posteriorMLN.csv")

#HU path
#ANT <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_posteriorMLN.csv")

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
names(ANT) <- c("Sample","ThCD4p","TcCD8p","Th1IFNgp_in_CD4p","Th17IL17Ap_in_CD4p",
                "Tc1IFNgp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_Ki67p_in_Foxp3p",
                "RORgtp_in_Foxp3p","Th1Tbetp_in_CD4pFoxp3n","Dividing_Ki67p_in_Tbetp",
                "Th17RORgp_in_CD4pFoxp3n","Dividing_Ki67p_in_RORgtp")
names(POS) <- c("Sample","ThCD4p","TcCD8p","Th1IFNgp_in_CD4p","Th17IL17Ap_in_CD4p",
                "Tc1IFNgp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_Ki67p_in_Foxp3p",
                "RORgtp_in_Foxp3p","Th1Tbetp_in_CD4pFoxp3n","Dividing_Ki67p_in_Tbetp",
                "Th17RORgp_in_CD4pFoxp3n","Dividing_Ki67p_in_RORgtp")

#set columns to numbers
library("dplyr")
library("magrittr")
str(ANT)
#ANT
ANT$ThCD4p <- as.numeric(as.character(ANT$ThCD4p))
ANT$TcCD8p <- as.numeric(as.character(ANT$TcCD8p))
ANT$Th1IFNgp_in_CD4p <- as.numeric(as.character(ANT$Th1IFNgp_in_CD4p))
ANT$Th17IL17Ap_in_CD4p <- as.numeric(as.character(ANT$Th17IL17Ap_in_CD4p))
ANT$Tc1IFNgp_in_CD8p <- as.numeric(as.character(ANT$Tc1IFNgp_in_CD8p))
ANT$Treg_Foxp3_in_CD4p <- as.numeric(as.character(ANT$Treg_Foxp3_in_CD4p))
ANT$Dividing_Ki67p_in_Foxp3p <- as.numeric(as.character(ANT$Dividing_Ki67p_in_Foxp3p))
ANT$RORgtp_in_Foxp3p <- as.numeric(as.character(ANT$RORgtp_in_Foxp3p))
ANT$Th1Tbetp_in_CD4pFoxp3n <- as.numeric(as.character(ANT$Th1Tbetp_in_CD4pFoxp3n))
ANT$Dividing_Ki67p_in_Tbetp <- as.numeric(as.character(ANT$Dividing_Ki67p_in_Tbetp))
ANT$Th17RORgp_in_CD4pFoxp3n <- as.numeric(as.character(ANT$Th17RORgp_in_CD4pFoxp3n))
ANT$Dividing_Ki67p_in_RORgtp <- as.numeric(as.character(ANT$Dividing_Ki67p_in_RORgtp))

#POS
ANT$ThCD4p <- as.numeric(as.character(ANT$ThCD4p))
ANT$TcCD8p <- as.numeric(as.character(ANT$TcCD8p))
ANT$Th1IFNgp_in_CD4p <- as.numeric(as.character(ANT$Th1IFNgp_in_CD4p))
ANT$Th17IL17Ap_in_CD4p <- as.numeric(as.character(ANT$Th17IL17Ap_in_CD4p))
ANT$Tc1IFNgp_in_CD8p <- as.numeric(as.character(ANT$Tc1IFNgp_in_CD8p))
ANT$Treg_Foxp3_in_CD4p <- as.numeric(as.character(ANT$Treg_Foxp3_in_CD4p))
ANT$Dividing_Ki67p_in_Foxp3p <- as.numeric(as.character(ANT$Dividing_Ki67p_in_Foxp3p))
ANT$RORgtp_in_Foxp3p <- as.numeric(as.character(ANT$RORgtp_in_Foxp3p))
ANT$Th1Tbetp_in_CD4pFoxp3n <- as.numeric(as.character(ANT$Th1Tbetp_in_CD4pFoxp3n))
ANT$Dividing_Ki67p_in_Tbetp <- as.numeric(as.character(ANT$Dividing_Ki67p_in_Tbetp))
ANT$Th17RORgp_in_CD4pFoxp3n <- as.numeric(as.character(ANT$Th17RORgp_in_CD4pFoxp3n))
ANT$Dividing_Ki67p_in_RORgtp <- as.numeric(as.character(ANT$Dividing_Ki67p_in_RORgtp))

#check structure
str(ANT)
str(POS)

#mutate "Sample" to character
i <- sapply(ANT, is.factor)
ANT[i] <- lapply(ANT[i], as.character)
str(ANT)

i <- sapply(POS, is.factor)
POS[i] <- lapply(POS[i], as.character)

#reshape df
library(reshape2)
library(meltt)
library(MASS)

cANT <- melt(ANT, id = c("Sample"))
cPOS <- melt(POS, id = c("Sample"))

#reduce labels
library(stringr)
cANT$Sample <- gsub(cANT$Sample, pattern="Anterior", replacement='')
cANT$Sample <- gsub(cANT$Sample, pattern=".fcs", replacement='')

cPOS$Sample <- gsub(cPOS$Sample, pattern="Posterior", replacement='')
cPOS$Sample <- gsub(cPOS$Sample, pattern=".fcs", replacement='')


#orientation plot
library(ggplot2)

ggplot(cANT, aes(Sample, value)) +   
  geom_bar(aes(fill = variable), position = "fill", stat="identity") +
  labs (title = "Anterior MLN", x = "Samples", y = "% of population") +
  coord_flip()

ggplot(cPOS, aes(Sample, value)) +   
  geom_bar(aes(fill = variable), position = "fill", stat="identity") +
  labs (title = "Posterior MLN", x = "Samples", y = "% of population") + 
  coord_flip()

#####################################################################################################################
#select T-cell subsets from ANT
CD4 <- cANT %>% filter(variable == "ThCD4p")
CD8 <- cANT %>% filter(variable == "TcCD8p")

#
Th1_IFNgp_in_CD4 <- cANT %>% filter(variable == "Th1IFNgp_in_CD4p")
Th1_IL17Ap_in_CD4 <- cANT %>% filter(variable == "Th17IL17Ap_in_CD4p")
Treg_Foxp3_in_CD4 <- cANT %>% filter(variable == "Treg_Foxp3_in_CD4p")

Tc1_IFNgp_in_CD8 <- cANT %>% filter(variable == "Tc1IFNgp_in_CD8p")

Th1_Tbet_in_CD4Foxp3n <- cANT %>% filter(variable == "Th1Tbetp_in_CD4pFoxp3n")
Th17RORgp_in_CD4Foxp3n <- cANT %>% filter(variable == "Th17RORgp_in_CD4pFoxp3n")

RORgtp_in_Foxp3p <- cANT %>% filter(variable == "RORgtp_in_Foxp3p")

Dividing_K67p_in_Foxp3p <- cANT %>% filter(variable == "Dividing_Ki67p_in_Foxp3p")

Dividing_Ki67p_in_Tbetp <- cANT %>% filter(variable == "Dividing_Ki67p_in_Tbetp")

Dividing_Ki67p_in_RORgtp <- cANT %>% filter(variable == "Dividing_Ki67p_in_RORgtp")

#create populations with subsets
CDANT <- rbind(CD4, CD8)
CD4T <- rbind(Th1_IFNgp_in_CD4, Th1_IL17Ap_in_CD4, Treg_Foxp3_in_CD4)

ggplot(CD4T, aes(Sample, value)) + 
  geom_bar(aes(fill = variable), position = "stack", stat = "identity") +
  coord_flip()





#basic bar
ggplot(CDANT, aes(Sample, value)) + 
  geom_bar(aes(fill = variable), position = "dodge", stat = "identity") +
  coord_flip()
#by T-cells
ggplot(CDANT, aes(variable, value)) + 
  geom_bar(aes(fill = Sample), position = "dodge", stat = "identity") +
  coord_flip()

#####################################################################################################################

#select T-cell subsets from POS
CD4 <- cANT %>% filter(variable == "ThCD4p")
CD8 <- cANT %>% filter(variable == "TcCD8p")

#
Th1_IFNgp_in_CD4 <- cANT %>% filter(variable == "Th1IFNgp_in_CD4p")
Th1_IL17Ap_in_CD4 <- cANT %>% filter(variable == "Th17IL17Ap_in_CD4p")
Treg_Foxp3_in_CD4 <- cANT %>% filter(variable == "Treg_Foxp3_in_CD4p")

Tc1_IFNgp_in_CD8 <- cANT %>% filter(variable == "Tc1IFNgp_in_CD8p")

Th1_Tbet_in_CD4Foxp3n <- cANT %>% filter(variable == "Th1Tbetp_in_CD4pFoxp3n")
Th17RORgp_in_CD4Foxp3n <- cANT %>% filter(variable == "Th17RORgp_in_CD4pFoxp3n")

RORgtp_in_Foxp3p <- cANT %>% filter(variable == "RORgtp_in_Foxp3p")

Dividing_K67p_in_Foxp3p <- cANT %>% filter(variable == "Dividing_Ki67p_in_Foxp3p")

Dividing_Ki67p_in_Tbetp <- cANT %>% filter(variable == "Dividing_Ki67p_in_Tbetp")

Dividing_Ki67p_in_RORgtp <- cANT %>% filter(variable == "Dividing_Ki67p_in_RORgtp")

#create populations with subsets
CDANT <- rbind(CD4, CD8)
CD4T <- rbind(Th1_IFNgp_in_CD4, Th1_IL17Ap_in_CD4, Treg_Foxp3_in_CD4)

ggplot(CD4T, aes(Sample, value)) + 
  geom_bar(aes(fill = variable), position = "stack", stat = "identity") +
  coord_flip()





#basic bar
ggplot(CDANT, aes(Sample, value)) + 
  geom_bar(aes(fill = variable), position = "dodge", stat = "identity") +
  coord_flip()
#by T-cells
ggplot(CDANT, aes(variable, value)) + 
  geom_bar(aes(fill = Sample), position = "dodge", stat = "identity") +
  coord_flip()