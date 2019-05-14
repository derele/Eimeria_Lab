#LOAD, CLEAN UP AND PROCESS DATA#
##########################################################################################################################
#PC path
ANT <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_anteriorMLN.csv")
POS <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_posteriorMLN.csv")

#laptop path
#ANT <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")

#IZW path
#ANT <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_posteriorMLN.csv")

#HU path
#ANT <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_posteriorMLN.csv")

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
str(POS)
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
POS$ThCD4p <- as.numeric(as.character(POS$ThCD4p))
POS$TcCD8p <- as.numeric(as.character(POS$TcCD8p))
POS$Th1IFNgp_in_CD4p <- as.numeric(as.character(POS$Th1IFNgp_in_CD4p))
POS$Th17IL17Ap_in_CD4p <- as.numeric(as.character(POS$Th17IL17Ap_in_CD4p))
POS$Tc1IFNgp_in_CD8p <- as.numeric(as.character(POS$Tc1IFNgp_in_CD8p))
POS$Treg_Foxp3_in_CD4p <- as.numeric(as.character(POS$Treg_Foxp3_in_CD4p))
POS$Dividing_Ki67p_in_Foxp3p <- as.numeric(as.character(POS$Dividing_Ki67p_in_Foxp3p))
POS$RORgtp_in_Foxp3p <- as.numeric(as.character(POS$RORgtp_in_Foxp3p))
POS$Th1Tbetp_in_CD4pFoxp3n <- as.numeric(as.character(POS$Th1Tbetp_in_CD4pFoxp3n))
POS$Dividing_Ki67p_in_Tbetp <- as.numeric(as.character(POS$Dividing_Ki67p_in_Tbetp))
POS$Th17RORgp_in_CD4pFoxp3n <- as.numeric(as.character(POS$Th17RORgp_in_CD4pFoxp3n))
POS$Dividing_Ki67p_in_RORgtp <- as.numeric(as.character(POS$Dividing_Ki67p_in_RORgtp))

#check structure
str(ANT)
str(POS)

#mutate "Sample" to character
i <- sapply(ANT, is.factor)
ANT[i] <- lapply(ANT[i], as.character)
str(ANT)

i <- sapply(POS, is.factor)
POS[i] <- lapply(POS[i], as.character)
str(POS)

#reduce labels
library(stringr)
#ANT$Sample <- gsub(ANT$Sample, pattern="Anterior", replacement='')
#ANT$Sample <- gsub(ANT$Sample, pattern=".fcs", replacement='')
#
#POS$Sample <- gsub(POS$Sample, pattern="Posterior", replacement='')
#POS$Sample <- gsub(POS$Sample, pattern=".fcs", replacement='')

#extract Mouse_ID from that mess and paste in "LM02" to standardize with our data structure
x = ANT$Sample
AIDs <- data.frame(Sample = x, EH_ID = sapply(strsplit(x,"_"), function(f)f[2]), Position = sapply(strsplit(x, " "), function(f)f[2]))
AIDs$EH_ID <- paste0("LM02", AIDs$EH_ID)

y = POS$Sample
PIDs <- data.frame(Sample = y, EH_ID = sapply(strsplit(y,"_"), function(f)f[2]), Position = sapply(strsplit(y, " "), function(f)f[2]))
PIDs$EH_ID <- paste0("LM02", PIDs$EH_ID)

#merge
ANT <- merge(AIDs, ANT)
POS <- merge(PIDs, POS)
MLNs <- rbind(ANT, POS)

#####################################################################################################################################
#introduce parasitological data
# HU: 
#Exp007 <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge.csv")
# Home PC: 
Exp007 <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge.csv")

Exp007_FACS <- merge(Exp007, MLNs)
Exp007_FACS$X <- NULL
str(Exp007_FACS)
#reshape for model
library(reshape)
library(reshape2)
library(meltt)
library(MASS)

wide <- dcast(Exp007_FACS, value.var = "dpi", "ThCD4p", "TcCD8p", "Th1IFNgp_in_CD4p", "Th17IL17Ap_in_CD4p", "Tc1IFNgp_in_CD8p", "Treg_Foxp3_in_CD4p", "Dividing_Ki67p_in_Foxp3p", "RORgtp_in_Foxp3p", "Th1Tbetp_in_CD4pFoxp3n", "Dividing_Ki67p_in_Tbetp", "Th17RORgp_in_CD4pFoxp3n", "Dividing_Ki67p_in_RORgtp")

ModelFACS <- recast(Exp007_FACS, id.var = "weight_dpi0", "fecweight", "ThCD4p")







#reshape df


cANT <- melt(ANT, id = c("Sample"))
cPOS <- melt(POS, id = c("Sample"))


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
A_CD4 <- cANT %>% filter(variable == "ThCD4p")
A_CD8 <- cANT %>% filter(variable == "TcCD8p")

#
A_Th1_IFNgp_in_CD4 <- cANT %>% filter(variable == "Th1IFNgp_in_CD4p")
A_Th1_IL17Ap_in_CD4 <- cANT %>% filter(variable == "Th17IL17Ap_in_CD4p")
A_Treg_Foxp3_in_CD4 <- cANT %>% filter(variable == "Treg_Foxp3_in_CD4p")

A_Tc1_IFNgp_in_CD8 <- cANT %>% filter(variable == "Tc1IFNgp_in_CD8p")

A_Th1_Tbet_in_CD4Foxp3n <- cANT %>% filter(variable == "Th1Tbetp_in_CD4pFoxp3n")
A_Th17RORgp_in_CD4Foxp3n <- cANT %>% filter(variable == "Th17RORgp_in_CD4pFoxp3n")

A_RORgtp_in_Foxp3p <- cANT %>% filter(variable == "RORgtp_in_Foxp3p")

A_Dividing_K67p_in_Foxp3p <- cANT %>% filter(variable == "Dividing_Ki67p_in_Foxp3p")

A_Dividing_Ki67p_in_Tbetp <- cANT %>% filter(variable == "Dividing_Ki67p_in_Tbetp")

A_Dividing_Ki67p_in_RORgtp <- cANT %>% filter(variable == "Dividing_Ki67p_in_RORgtp")

#create populations with subsets
A_CDANT <- rbind(CD4, CD8)
A_CD4T <- rbind(Th1_IFNgp_in_CD4, Th1_IL17Ap_in_CD4, Treg_Foxp3_in_CD4)

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
P_CD4 <- cANT %>% filter(variable == "ThCD4p")
P_CD8 <- cANT %>% filter(variable == "TcCD8p")

#
P_Th1_IFNgp_in_CD4 <- cANT %>% filter(variable == "Th1IFNgp_in_CD4p")
P_Th1_IL17Ap_in_CD4 <- cANT %>% filter(variable == "Th17IL17Ap_in_CD4p")
P_Treg_Foxp3_in_CD4 <- cANT %>% filter(variable == "Treg_Foxp3_in_CD4p")

P_Tc1_IFNgp_in_CD8 <- cANT %>% filter(variable == "Tc1IFNgp_in_CD8p")

P_Th1_Tbet_in_CD4Foxp3n <- cANT %>% filter(variable == "Th1Tbetp_in_CD4pFoxp3n")
P_Th17RORgp_in_CD4Foxp3n <- cANT %>% filter(variable == "Th17RORgp_in_CD4pFoxp3n")

P_RORgtp_in_Foxp3p <- cANT %>% filter(variable == "RORgtp_in_Foxp3p")

P_Dividing_K67p_in_Foxp3p <- cANT %>% filter(variable == "Dividing_Ki67p_in_Foxp3p")

P_Dividing_Ki67p_in_Tbetp <- cANT %>% filter(variable == "Dividing_Ki67p_in_Tbetp")

P_Dividing_Ki67p_in_RORgtp <- cANT %>% filter(variable == "Dividing_Ki67p_in_RORgtp")

#create populations with subsets
P_CDANT <- rbind(CD4, CD8)
P_CD4T <- rbind(Th1_IFNgp_in_CD4, Th1_IL17Ap_in_CD4, Treg_Foxp3_in_CD4)

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
