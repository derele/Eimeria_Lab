#load up libraries and raw csv of cell counts
library(httr)
library(RCurl)

CELLSfileUrl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_cell_counts_raw.csv"
CELLS <- read.csv(text=getURL(CELLSfileUrl))

# include percentage counts in the df
ANTfileUrl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_CD40L_assays_anteriorMLN.csv"
ANT <- read.csv(text=getURL(ANTfileUrl))

POSfileUrl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_CD40L_assays_posteriorMLN.csv"
POS <- read.csv(text=getURL(POSfileUrl))

#rename columns sensibly
names(ANT) <- c("Sample","ThCD4p","TcCD8p","Th1IFNgp_in_CD4p","Th17IL17Ap_in_CD4p",
                "Tc1IFNgp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_Ki67p_in_Foxp3p",
                "RORgtp_in_Foxp3p","Th1Tbetp_in_CD4pFoxp3n","Dividing_Ki67p_in_Tbetp",
                "Th17RORgp_in_CD4pFoxp3n","Dividing_Ki67p_in_RORgtp")
names(POS) <- c("Sample","ThCD4p","TcCD8p","Th1IFNgp_in_CD4p","Th17IL17Ap_in_CD4p",
                "Tc1IFNgp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_Ki67p_in_Foxp3p",
                "RORgtp_in_Foxp3p","Th1Tbetp_in_CD4pFoxp3n","Dividing_Ki67p_in_Tbetp",
                "Th17RORgp_in_CD4pFoxp3n","Dividing_Ki67p_in_RORgtp")

#bind ANT and POS together + include actual cell counts
PERCENTAGES <- rbind(ANT, POS)
cell.counts <- cbind(CELLS, PERCENTAGES)

#extract Mouse_ID from that mess and paste in "LM02" to standardize with our data structure
cell.counts$EH_ID <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "LM02\\2", cell.counts$Sample)
cell.counts$Position <- gsub("\\d+: (Anterior|Posterior) LN_(\\d{2})_\\d{3}.fcs", "\\1", cell.counts$Sample)

#need to make caculations for each row per cell subset from the originating one
#total cell counts subset to thCD4p and TcCD8p
cell.counts$ThCD4p.cells <- with(cell.counts, (cell_counts/100)*ThCD4p)
cell.counts$TcCD8p.cells <- with(cell.counts, (cell_counts/100)*TcCD8p)

#ThCD4p subsets to Th1IFNgp_in_CD4p, Th17IL17Ap_in_CD4p, Treg_Foxp3_in_CD4p
cell.counts$Th1IFNgp_in_CD4p.cells <- with(cell.counts, (ThCD4p.cells/100)*Th1IFNgp_in_CD4p)
cell.counts$Th17IL17Ap_in_CD4p.cells <- with(cell.counts, (ThCD4p.cells/100)*Th17IL17Ap_in_CD4p)
cell.counts$Treg_Foxp3_in_CD4p.cells <- with(cell.counts, (ThCD4p.cells/100)*Treg_Foxp3_in_CD4p)

#TcCD8p subsets to Tc1IFNgp_in_CD8p
cell.counts$Tc1IFNgp_in_CD8p.cells <- with(cell.counts, (TcCD8p.cells/100)*Tc1IFNgp_in_CD8p)

#Treg_Foxp3_in_CD4p subsets to RORgtp_in_Foxp3p
cell.counts$RORgtp_in_Foxp3p.cells <- with(cell.counts, (Treg_Foxp3_in_CD4p.cells/100)*RORgtp_in_Foxp3p)

#RORgtp_in_Foxp3p subsets to Dividing_Ki67p_in_Foxp3p and Dividing_Ki67p_in_RORgtp
cell.counts$Dividing_Ki67p_in_Foxp3p.cells <- with(cell.counts, (RORgtp_in_Foxp3p.cells/100)*Dividing_Ki67p_in_Foxp3p)
cell.counts$Dividing_Ki67p_in_RORgtp.cells <- with(cell.counts, (RORgtp_in_Foxp3p.cells/100)*Dividing_Ki67p_in_RORgtp)

#ThCD4p_Foxp3n (doesn't exist here) subsets to Th1Tbetp_in_CD4pFoxp3n, Th17RORgp_in_CD4pFoxp3n
cell.counts$ThCD4p_Foxp3n <- with(cell.counts, (100-(Th17IL17Ap_in_CD4p + Th1IFNgp_in_CD4p + Treg_Foxp3_in_CD4p)))
cell.counts$ThCD4p_Foxp3n.cells <- with(cell.counts, (ThCD4p.cells - Treg_Foxp3_in_CD4p))
cell.counts$Th1Tbetp_in_CD4pFoxp3n.cells <- with(cell.counts, (ThCD4p_Foxp3n.cells/100)*Th1Tbetp_in_CD4pFoxp3n)
cell.counts$Th17RORgp_in_CD4pFoxp3n.cells <- with(cell.counts, (ThCD4p_Foxp3n.cells/100)*Th17RORgp_in_CD4pFoxp3n)

#Th1Tbetp_in_CD4pFoxp3n subsets to Dividing_Ki67p_in_Tbetp
cell.counts$Dividing_Ki67p_in_Tbetp.cells <- with(cell.counts, (Th1Tbetp_in_CD4pFoxp3n.cells/100)*Dividing_Ki67p_in_Tbetp)

#clean up table before saving
#remove mouse column X (artifact)
cell.counts$X = NULL
cell.counts[,4] = NULL # only because of identical column names

#write the .csv (rewrite Hongweis raw (all data included)) 
write.csv(cell.counts, "./Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACS_cell_counts_processed.csv", quote = FALSE)
