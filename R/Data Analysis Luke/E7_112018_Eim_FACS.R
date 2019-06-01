#LOAD, CLEAN UP AND PROCESS DATA#
##########################################################################################################################
#PC path
ANT <- read.csv("./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_anteriorMLN.csv")
POS <- read.csv("./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_posteriorMLN.csv")

#laptop path
#ANT <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")
#POS <- read.csv(file = "../lubomir/Documents/Eimeria_Lab/data/3_recordingTables/Exp007/CD40L_assays_Exp007_anteriorMLN.csv")

#IZW path
#ANT <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_anteriorMLN.csv")
#POS <- read.csv(file = "../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_posteriorMLN.csv")

#HU path
#ANT <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_anteriorMLN.csv")
#POS <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_CD40L_assays_posteriorMLN.csv")

#remove mouse 293 as it was mixed both posterior and anterior
ANT <- ANT[-c(52),]

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

#convert 
# num.convert <- function(data) {
#   
#   char.columns <- colnames(data)[1:3]
#   cols <- colnames(data)[4:ncol(data)]
#   
#   #make empty matrix/data.frame
#   df <- data.frame()
#   
#   for(a in 1:length(cols))
#   {
#     df[1:nrow(data),a] <- as.numeric(as.character(data[,cols[a]]))
#   }
#   colnames(df) <- cols
#   new.df <- data.frame(Sample = data[,char.columns], df)
#   
#   return(new.df)
# }
# 
# ANT <- num.convert(data = ANT)
# POS <- num.convert(data = POS)


#mutate "Sample" to character
i <- sapply(ANT, is.factor)
ANT[i] <- lapply(ANT[i], as.character)


i <- sapply(POS, is.factor)
POS[i] <- lapply(POS[i], as.character)

CELLS <- rbind(ANT, POS)


#extract Mouse_ID from that mess and paste in "LM02" to standardize with our data structure
CELLS$EH_ID <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM02\\2", CELLS$Sample)
CELLS$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", CELLS$Sample)


#####################################################################################################################################
#introduce parasitological data
# HU: 
#E7 <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv")

#IZW path
#Exp007 <- read.csv("../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv")

#Home PC: 
E7 <- read.csv("./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv")

E7_FACS <- merge(E7, CELLS)
E7_FACS$X <- NULL

#select cell population names
facs.measure.cols <- c("ThCD4p", "TcCD8p", "Th1IFNgp_in_CD4p", "Th17IL17Ap_in_CD4p", 
                       "Tc1IFNgp_in_CD8p", "Treg_Foxp3_in_CD4p", "Dividing_Ki67p_in_Foxp3p", 
                       "RORgtp_in_Foxp3p", "Th1Tbetp_in_CD4pFoxp3n", "Dividing_Ki67p_in_Tbetp", 
                       "Th17RORgp_in_CD4pFoxp3n", "Dividing_Ki67p_in_RORgtp")


#select dpi8 from FACS
dpi8POS <- E7_FACS[E7_FACS$dpi%in%8 & E7_FACS$Position%in%"Posterior",]
dpi8ANT <- E7_FACS[E7_FACS$dpi%in%8 & E7_FACS$Position%in%"Anterior",]

#primary infection strain effect on CD8 INF positive cell populations
CD8_INF_POS_prim <- tapply(dpi8POS$Tc1IFNgp_in_CD8p, dpi8POS$primary, summary)
CD8_INF_ANT_prim <- tapply(dpi8ANT$Tc1IFNgp_in_CD8p, dpi8ANT$primary, summary)

#challenge infection strain effect on CD8 INF positive cell populations
CD8_INF_POS_cha <- tapply(dpi8POS$Tc1IFNgp_in_CD8p, dpi8POS$challenge, summary)
CD8_INF_ANT_cha <- tapply(dpi8ANT$Tc1IFNgp_in_CD8p, dpi8ANT$challenge, summary)

######################################## wilcox of medians (infection strains comparison primary)

#create list of cell populations summaries infection strains
cell.sumariesANT.prim <- lapply(facs.measure.cols, function (x){
  tapply(dpi8ANT[, x], dpi8ANT$primary, summary)
})

cell.sumariesPOS.prim <- lapply(facs.measure.cols, function (x){
  tapply(dpi8POS[, x], dpi8POS$primary, summary)
})


#set identical names
names(cell.sumariesANT.prim) <- facs.measure.cols
names(cell.sumariesPOS.prim) <- facs.measure.cols

#medians comparisons between infection strains
cell.mediansANT.prim <- lapply(facs.measure.cols, function (x){
  tapply(dpi8ANT[, x], dpi8ANT$primary, median)
})

cell.mediansPOS.prim <- lapply(facs.measure.cols, function (x){
  tapply(dpi8POS[, x], dpi8POS$primary, median)
})

#set identical names
names(cell.mediansANT.prim) <- facs.measure.cols
names(cell.mediansPOS.prim) <- facs.measure.cols

#non-parametric Whitney Mann tests of tcell pops
cell.testsANT.prim <- lapply(facs.measure.cols, function (x){
  wilcox.test(dpi8ANT[dpi8ANT$primary%in%"E88", x], dpi8ANT[dpi8ANT$primary%in%"E64", x])
})


cell.testsPOS.prim <- lapply(facs.measure.cols, function (x){
  wilcox.test(dpi8POS[dpi8POS$primary%in%"E88", x], dpi8POS[dpi8POS$primary%in%"E64", x])
})


#set identical names
names(cell.testsANT.prim) <- facs.measure.cols
names(cell.testsPOS.prim) <- facs.measure.cols


#store intra strain wilcox results 
wilcox_medians_ANT.prim <- cbind(unlist(cell.mediansANT.prim), rep(unlist(lapply(cell.testsANT.prim, "[", "p.value")), each=2))
colnames(wilcox_medians_ANT.prim)
colnames(wilcox_medians_ANT.prim) <- c("population_percentages", "p_value")

wilcox_medians_POS.prim <- cbind(unlist(cell.mediansPOS.prim), rep(unlist(lapply(cell.testsPOS.prim, "[", "p.value")), each=2))
colnames(wilcox_medians_POS.prim)
colnames(wilcox_medians_POS.prim) <- c("population_percentages", "p_value")

##### compare between results of wilcox (ANT vs POS)

##### row differences
library(matrixTests)
wilcox_medians_ANT.prim <- as.data.frame(wilcox_medians_ANT.prim)
row_kruskalwallis(t(wilcox_medians_ANT.prim), g = wilcox_medians_ANT.prim$population_percentages)

#sapply(intersect(rownames(wilcox_medians_ANT.prim), rownames(wilcox_medians_POS.prim)), 
 #      function(x) Kruskal(wilcox_medians_ANT.prim[x,], wilcox_medians_POS.prim[x,]))

######################################## wilcox of medians (infection strains comparison challenge)

#create list of cell populations summaries infection strains
cell.sumariesANT.cha <- lapply(facs.measure.cols, function (x){
  tapply(dpi8ANT[, x], dpi8ANT$challenge, summary)
})

cell.sumariesPOS.cha <- lapply(facs.measure.cols, function (x){
  tapply(dpi8POS[, x], dpi8POS$challenge, summary)
})

#set identical names
names(cell.sumariesANT.cha) <- facs.measure.cols
names(cell.sumariesPOS.cha) <- facs.measure.cols

#medians comparisons between infection strains
cell.mediansANT.cha <- lapply(facs.measure.cols, function (x){
  tapply(dpi8ANT[, x], dpi8ANT$challenge, median)
})

cell.mediansPOS.cha <- lapply(facs.measure.cols, function (x){
  tapply(dpi8POS[, x], dpi8POS$challenge, median)
})

#set identical names
names(cell.mediansANT.cha) <- facs.measure.cols
names(cell.mediansPOS.cha) <- facs.measure.cols

#non-parametric Whitney Mann tests of tcell pops
cell.testsANT.cha <- lapply(facs.measure.cols, function (x){
  wilcox.test(dpi8ANT[dpi8ANT$challenge%in%"E88", x], dpi8ANT[dpi8ANT$challenge%in%"E64", x])
})

cell.testsPOS.cha <- lapply(facs.measure.cols, function (x){
  wilcox.test(dpi8POS[dpi8POS$challenge%in%"E88", x], dpi8POS[dpi8POS$challenge%in%"E64", x])
})

#set identical names
names(cell.testsANT.cha) <- facs.measure.cols
names(cell.testsPOS.cha) <- facs.measure.cols

#store intra strain wilcox results 
wilcox_medians_ANT.cha <- cbind(unlist(cell.mediansANT.cha), rep(unlist(lapply(cell.testsANT.cha, "[", "p.value")), each=2))
colnames(wilcox_medians_ANT.cha)
colnames(wilcox_medians_ANT.cha) <- c("population_percentages", "p_value")

wilcox_medians_POS.cha <- cbind(unlist(cell.mediansPOS.cha), rep(unlist(lapply(cell.testsPOS.cha, "[", "p.value")), each=2))
colnames(wilcox_medians_POS.cha)
colnames(wilcox_medians_POS.cha) <- c("population_percentages", "p_value")

W_Med_ANT.cha <- round(wilcox_medians_ANT.cha, 3)
W_Med_POS.cha <- round(wilcox_medians_POS.cha, 3)

###################################################cell medians vs infection history
names(cell.mediansANT_HIST) <- facs.measure.cols
names(cell.mediansPOS_HIST) <- facs.measure.cols

cell.mediansANT_HIST <- lapply(facs.measure.cols, function (x){
  tapply(dpi8ANT[, x], dpi8ANT$primary:dpi8ANT$challenge, median)
})

cell.mediansPOS_HIST <- lapply(facs.measure.cols, function (x){
  tapply(dpi8POS[, x], dpi8POS$primary:dpi8ANT$challenge, median)
})

#######test with man whitney
library(magrittr)
names(cell.mediansANT_HIST) <- facs.measure.cols

cell.mediansANT.history <- data.frame(matrix(unlist(cell.mediansANT_HIST), nrow=length(cell.mediansANT_HIST), byrow=T))
cell.mediansANT.history <- set_rownames(cell.mediansANT.history, facs.measure.cols)
cell.mediansANT.history <- set_colnames(cell.mediansANT.history, c("E64:E64", "E64:E88", "E88:E64", "E88:E88"))
strains <- c(names(cell.mediansANT.history))

cell.mediansPOS.history <- data.frame(matrix(unlist(cell.mediansPOS_HIST), nrow=length(cell.mediansPOS_HIST), byrow=T))
cell.mediansPOS.history <- set_rownames(cell.mediansPOS.history, facs.measure.cols)
cell.mediansPOS.history <- set_colnames(cell.mediansPOS.history, c("E64:E64", "E64:E88", "E88:E64", "E88:E88"))
strains <- c(names(cell.mediansPOS.history))


#figure this out (test cell populations between infection histories)
 #cell.tests.history <- lapply(strains, function (x){
  # wilcox.test(cell.mediansANT.history[x, ], cell.mediansPOS.history[x, ])
 #})

 
 
 

# cell.testsANT.history <- lapply(cell.mediansANT.history, function (x){
#   wilcox.test(x$`E64:E64`, x$`E64:E88`, x$`E88:E64`,  x$`E88:E88`)
# })
# 
# #then rename
# names(cell.testsANT.cha) <- facs.measure.cols

