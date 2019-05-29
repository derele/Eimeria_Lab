#LOAD, CLEAN UP AND PROCESS DATA#
##########################################################################################################################
#PC path
#ANT <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/FACS/CD40L_assays_Exp007_posteriorMLN.csv")

#laptop path
#ANT <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")

#IZW path
ANT <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_anteriorMLN.csv")
POS <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_posteriorMLN.csv")

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

#ANT
ANT_num.convert <- function(data) {
  
  char.columns <- colnames(data)[1:3]
  cols <- colnames(data)[4:ncol(data)]
  
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


#mutate "Sample" to character
i <- sapply(ANT, is.factor)
ANT[i] <- lapply(ANT[i], as.character)


i <- sapply(POS, is.factor)
POS[i] <- lapply(POS[i], as.character)

CELLS <- rbind(ANT, POS)


#extract Mouse_ID from that mess and paste in "LM02" to standardize with our data structure
CELLS$EH_ID <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM02\\2", CELLS$Sample.Sample)
CELLS$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", CELLS$Sample.Sample)


#####################################################################################################################################
#introduce parasitological data
# HU: 
#Exp007 <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7")

#IZW path
Exp007 <- read.csv("../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv")
#Home PC: 
#Exp007 <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge.csv")

Exp007_FACS <- merge(Exp007, CELLS)
Exp007_FACS$X <- NULL
str(Exp007_FACS)
#reshape for model
library(reshape)
library(reshape2)
library(meltt)
library(MASS)


#select cell population names
facs.measure.cols <- c("ThCD4p", "TcCD8p", "Th1IFNgp_in_CD4p", "Th17IL17Ap_in_CD4p", 
                       "Tc1IFNgp_in_CD8p", "Treg_Foxp3_in_CD4p", "Dividing_Ki67p_in_Foxp3p", 
                       "RORgtp_in_Foxp3p", "Th1Tbetp_in_CD4pFoxp3n", "Dividing_Ki67p_in_Tbetp", 
                       "Th17RORgp_in_CD4pFoxp3n", "Dividing_Ki67p_in_RORgtp")


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

cell.mediansANT_HIST <- lapply(facs.measure.cols, function (x){
  tapply(dpi8ANT[, x], dpi8ANT$primary:dpi8ANT$challenge, median)
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


round(wilcox_medians_ANT, 3)


### combine two factors
as.factor(letters):as.factor(rev(letters))


#Melt?
#orientation plot

library(ggplot2)
