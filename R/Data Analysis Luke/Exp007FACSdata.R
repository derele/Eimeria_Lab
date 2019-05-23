#LOAD, CLEAN UP AND PROCESS DATA#
##########################################################################################################################
#PC path
#ANT <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_posteriorMLN.csv")

#laptop path
#ANT <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")

#IZW path
#ANT <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_posteriorMLN.csv")

#HU path
ANT <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_anteriorMLN.csv")
POS <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_posteriorMLN.csv")

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
ANT_num.convert <- function(data) {
  
  char.columns <- colnames(data)[1]
  cols <- colnames(data)[2:ncol(data)]
  
  #make empty matrix/data.frame
  df <- data.frame()
  
  for(a in 1:length(cols))
  {
    df[1:nrow(data),a] <- as.numeric(as.character(data[,cols[a]]))
  }
  colnames(df) <- cols
  new.df <- data.frame(Sample = data[,char.columns], df)
  
  return(new.df)
}

ANT <- ANT_num.convert(data = ANT)
POS <- ANT_num.convert(data = POS)

#function defined
# ANT$ThCD4p <- as.numeric(as.character(ANT$ThCD4p))
# ANT$TcCD8p <- as.numeric(as.character(ANT$TcCD8p))
# ANT$Th1IFNgp_in_CD4p <- as.numeric(as.character(ANT$Th1IFNgp_in_CD4p))
# ANT$Th17IL17Ap_in_CD4p <- as.numeric(as.character(ANT$Th17IL17Ap_in_CD4p))
# ANT$Tc1IFNgp_in_CD8p <- as.numeric(as.character(ANT$Tc1IFNgp_in_CD8p))
# ANT$Treg_Foxp3_in_CD4p <- as.numeric(as.character(ANT$Treg_Foxp3_in_CD4p))
# ANT$Dividing_Ki67p_in_Foxp3p <- as.numeric(as.character(ANT$Dividing_Ki67p_in_Foxp3p))
# ANT$RORgtp_in_Foxp3p <- as.numeric(as.character(ANT$RORgtp_in_Foxp3p))
# ANT$Th1Tbetp_in_CD4pFoxp3n <- as.numeric(as.character(ANT$Th1Tbetp_in_CD4pFoxp3n))
# ANT$Dividing_Ki67p_in_Tbetp <- as.numeric(as.character(ANT$Dividing_Ki67p_in_Tbetp))
# ANT$Th17RORgp_in_CD4pFoxp3n <- as.numeric(as.character(ANT$Th17RORgp_in_CD4pFoxp3n))
# ANT$Dividing_Ki67p_in_RORgtp <- as.numeric(as.character(ANT$Dividing_Ki67p_in_RORgtp))
# 
# #POS
# POS$ThCD4p <- as.numeric(as.character(POS$ThCD4p))
# POS$TcCD8p <- as.numeric(as.character(POS$TcCD8p))
# POS$Th1IFNgp_in_CD4p <- as.numeric(as.character(POS$Th1IFNgp_in_CD4p))
# POS$Th17IL17Ap_in_CD4p <- as.numeric(as.character(POS$Th17IL17Ap_in_CD4p))
# POS$Tc1IFNgp_in_CD8p <- as.numeric(as.character(POS$Tc1IFNgp_in_CD8p))
# POS$Treg_Foxp3_in_CD4p <- as.numeric(as.character(POS$Treg_Foxp3_in_CD4p))
# POS$Dividing_Ki67p_in_Foxp3p <- as.numeric(as.character(POS$Dividing_Ki67p_in_Foxp3p))
# POS$RORgtp_in_Foxp3p <- as.numeric(as.character(POS$RORgtp_in_Foxp3p))
# POS$Th1Tbetp_in_CD4pFoxp3n <- as.numeric(as.character(POS$Th1Tbetp_in_CD4pFoxp3n))
# POS$Dividing_Ki67p_in_Tbetp <- as.numeric(as.character(POS$Dividing_Ki67p_in_Tbetp))
# POS$Th17RORgp_in_CD4pFoxp3n <- as.numeric(as.character(POS$Th17RORgp_in_CD4pFoxp3n))
# POS$Dividing_Ki67p_in_RORgtp <- as.numeric(as.character(POS$Dividing_Ki67p_in_RORgtp))

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
Exp007 <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge.csv")

#IZW path
#Exp007 <- read.csv("../luke/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge")
#Home PC: 
#Exp007 <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge.csv")

Exp007_FACS <- merge(Exp007, MLNs)
Exp007_FACS$X <- NULL
str(Exp007_FACS)
#reshape for model
library(reshape)
library(reshape2)
library(meltt)
library(MASS)

#reshape df
A_CD4 <- tibble(ANT$EH_ID, ANT$ThCD4p)
A_CD8 <- tibble(ANT$EH_ID, ANT$TcCD8p)
A_tpop <- merge(A_CD4, A_CD8)

#select cell population names
facs.measure.cols <- colnames(Exp007_FACS)[14:ncol(Exp007_FACS)]

## apply summary over all interesting columns
apply(Exp007_FACS[,facs.measure.cols], 2, summary)

## for one interesing column only for Posterior and dpi8 example selection
summary(Exp007_FACS[Exp007_FACS$dpi%in%8 & Exp007_FACS$Position%in%"Posterior","Tc1IFNgp_in_CD8p"])

#select dpi8 from FACS
dpi8POS <- Exp007_FACS[Exp007_FACS$dpi%in%8 & Exp007_FACS$Position%in%"Posterior",]
dpi8ANT <- Exp007_FACS[Exp007_FACS$dpi%in%8 & Exp007_FACS$Position%in%"Anterior",]

#infection strain effect on CD8 INF positive cell populations
tapply(dpi8POS$Tc1IFNgp_in_CD8p, dpi8POS$InfectionStrain, summary)
tapply(dpi8ANT$Tc1IFNgp_in_CD8p, dpi8ANT$InfectionStrain, summary)

#create list of cell populations summaries infection strains
cell.sumariesANT <- lapply(facs.measure.cols, function (x){
         tapply(dpi8ANT[, x], dpi8ANT$InfectionStrain, summary)
})

cell.sumariesPOS <- lapply(facs.measure.cols, function (x){
  tapply(dpi8POS[, x], dpi8POS$InfectionStrain, summary)
})
#set identical names
names(cell.sumariesANT) <- facs.measure.cols
names(cell.sumariesPOS) <- facs.measure.cols

#medians comparisons between infection strains
cell.mediansANT <- lapply(facs.measure.cols, function (x){
  tapply(dpi8ANT[, x], dpi8ANT$InfectionStrain, median)
})

cell.mediansPOS <- lapply(facs.measure.cols, function (x){
  tapply(dpi8POS[, x], dpi8POS$InfectionStrain, median)
})
#set identical names
names(cell.mediansANT) <- facs.measure.cols
names(cell.mediansPOS) <- facs.measure.cols

#non-parametric Whitney Mann tests of tcell pops
cell.testsANT <- lapply(facs.measure.cols, function (x){
  wilcox.test(dpi8ANT[dpi8ANT$InfectionStrain%in%"E88", x], dpi8ANT[dpi8ANT$InfectionStrain%in%"E64", x])
})

cell.testsPOS <- lapply(facs.measure.cols, function (x){
  wilcox.test(dpi8POS[dpi8POS$InfectionStrain%in%"E88", x], dpi8POS[dpi8POS$InfectionStrain%in%"E64", x])
})
#set identical names
names(cell.testsANT) <- facs.measure.cols
names(cell.testsPOS) <- facs.measure.cols

#store intra strain wilcox results 
wilcox_medians_ANT <- cbind(unlist(cell.mediansANT), rep(unlist(lapply(cell.testsANT, "[", "p.value")), each=2))
colnames(wilcox_medians_ANT)
colnames(wilcox_medians_ANT) <- c("population_percentages", "p_value")

wilcox_medians_POS <- cbind(unlist(cell.mediansPOS), rep(unlist(lapply(cell.testsPOS, "[", "p.value")), each=2))
colnames(wilcox_medians_POS)
colnames(wilcox_medians_POS) <- c("population_percentages", "p_value")

### combine two factors
as.factor(letters):as.factor(rev(letters))


#Melt?
#orientation plot

library(ggplot2)
