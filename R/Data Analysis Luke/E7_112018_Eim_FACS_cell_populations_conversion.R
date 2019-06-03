#load up libraries and raw csv of cell counts
library(httr)
library(RCurl)

CELLSfileUrl <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_cell_counts.csv"
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
#ThCD4p subsets to Th1IFNgp_in_CD4p, Th17IL17Ap_in_CD4p, Treg_Foxp3_in_CD4p
#TcCD8p subsets to Tc1IFNgp_in_CD8p
#Treg_Foxp3_in_CD4p subsets to RORgtp_in_Foxp3p
#RORgtp_in_Foxp3p subsets to Dividing_Ki67p_in_Foxp3p and Dividing_Ki67p_in_RORgtp
#ThCD4p_Foxp3n (doesn't exist here) subsets to Th1Tbetp_in_CD4pFoxp3n, Th17RORgp_in_CD4pFoxp3n
#Th1Tbetp_in_CD4pFoxp3n subsets to Dividing_Ki67p_in_Tbetp

