#LOAD, CLEAN UP AND PROCESS DATA#
##########################################################################################################################

library(httr)
library(RCurl)

ANTfileUrl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_CD40L_assays_anteriorMLN.csv"

ANT <- read.csv(text=getURL(ANTfileUrl))

POSfileUrl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_CD40L_assays_posteriorMLN.csv"

POS <- read.csv(text=getURL(POSfileUrl))

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
library(dplyr)
library(magrittr)

CELLS <- rbind(ANT, POS)

#extract Mouse_ID from that mess and paste in "LM02" to standardize with our data structure
CELLS$EH_ID <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM02\\2", CELLS$Sample)
CELLS$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", CELLS$Sample)

#remove mouse 293 as it was mixed both posterior and anterior
CELLS <- CELLS[!CELLS$EH_ID%in%"LM0293",]

#####################################################################################################################################
#introduce parasitological data
paraURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_complete.csv"

E7 <- read.csv(text=getURL(paraURL))
E7 <- merge(E7[E7$dpi%in%8,], CELLS)

E7$infHistory <- E7$primary:E7$challenge

##select cell population names
facs.measure.cols <- c("ThCD4p", "TcCD8p", "Th1IFNgp_in_CD4p", "Th17IL17Ap_in_CD4p", 
                       "Tc1IFNgp_in_CD8p", "Treg_Foxp3_in_CD4p", "Dividing_Ki67p_in_Foxp3p", 
                       "RORgtp_in_Foxp3p", "Th1Tbetp_in_CD4pFoxp3n", "Dividing_Ki67p_in_Tbetp", 
                       "Th17RORgp_in_CD4pFoxp3n", "Dividing_Ki67p_in_RORgtp")

## tabulate  medians for different infection histories and anterior vs posterior
## create list of cell populations summaries infection strains
cell.medians <- lapply(facs.measure.cols, function (x){
  tapply(E7[, x], list(E7$infHistory, as.factor(E7$Position)), median)
})

## test differences between infection histories (by positon anterior/posterior)
cell.testsPOS <- lapply(facs.measure.cols, function (x){
  kruskal.test(E7[E7$Position%in%"Posterior", x], E7[E7$Position%in%"Posterior", "infHistory"])
})

## set cell type names
names(cell.testsPOS) <- facs.measure.cols
cell.testsPOS


## test differences between infection histories (by positon anterior/posterior)
cell.testsANT <- lapply(facs.measure.cols, function (x){
  kruskal.test(E7[E7$Position%in%"Anterior", x], E7[E7$Position%in%"Anterior", "infHistory"])
})

## #check distribution
library(ggplot2)



plotCells <- function (col){
    ggplot(E7, aes(infHistory, get(col))) +
        geom_boxplot() +
        geom_jitter(width=0.2) +
        facet_wrap(~Position) +
        ggtitle(col)
}

facs_boxplots <- lapply(facs.measure.cols, plotCells)
names(facs_boxplots) <-  facs.measure.cols

for(i in seq_along(facs_boxplots)){
    pdf(paste0(names(facs_boxplots)[[i]], ".pdf"))
    plot(facs_boxplots[[i]])
    dev.off()
}



## ## remember to check the huge outlier value in Treg_Foxp3_in_CD4p
## ## this would exclude it graphically...
## ggplot(E7, aes(infHistory, Treg_Foxp3_in_CD4p)) +
##     geom_boxplot() +
##     geom_jitter(width=0.2) +
##     facet_wrap(~Position) +
##     ggtitle("Treg_Foxp3_in_CD4p") +
##     ylim(0, 14)
