#personal path
#DF <- read.csv(file = "../luke/Documents/Mm_norm_counts.csv")

#IZW path
#DF <- read.csv(file = "../luke/Documents/Mm_norm_counts.csv")

#HU path
DF <- read.csv(file = "../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/Mm_norm_counts.csv")

IFNg <- DF["ENSMUSG00000055170",]
TGFb <- DF["ENSMUSG00000002603",]
STAT6 <- DF["ENSMUSG00000002147",]
CXCL9 <- DF["ENSMUSG00000029417",]
IL10 <- DF["ENSMUSG00000016529",]
TNFa <- DF["ENSMUSG00000024401",]
CXCR3 <- DF["ENSMUSG00000050232",]
IRG6A <- DF["ENSMUSG00000054072",]
GBP2 <- DF["ENSMUSG00000028270",]
MyD88 <- DF["ENSMUSG00000032508",]
Ticam1 <- DF["ENSMUSG00000047123",]

CYT <- rbind(IFNg, TGFb, STAT6, CXCL9, IL10, TNFa, CXCR3, IRG6A, GBP2, MyD88, Ticam1)

#install.packages("data.table")
library(data.table)
library(reshape2)

setDT(CYT, keep.rownames = TRUE)
row.names(CYT) <- c("IFN-y", "TGF-b", "STAT6", "CXCL9", "IL-10", "TNFa", "CXCR3", "IRG6A", "GBP2", "MyD88", "Ticam1")

cCYT <- melt(CYT, id.vars = c("rn"))
cCYT[cCYT=="ENSMUSG00000055170"]<-"IFNg"
cCYT[cCYT=="ENSMUSG00000002603"]<-"TGFb"
cCYT[cCYT=="ENSMUSG00000002147"]<-"STAT6"
cCYT[cCYT=="ENSMUSG00000029417"]<-"CXCL9"
cCYT[cCYT=="ENSMUSG00000016529"]<-"IL10"
cCYT[cCYT=="ENSMUSG00000024401"]<-"TNFa"
cCYT[cCYT=="ENSMUSG00000050232"]<-"CXCR3"
cCYT[cCYT=="ENSMUSG00000054072"]<-"IRG6A"
cCYT[cCYT=="ENSMUSG00000028270"]<-"GBP2"
cCYT[cCYT=="ENSMUSG00000032508"]<-"MyD88"
cCYT[cCYT=="ENSMUSG00000047123"]<-"Ticam1"

x <- c("C57BL6_1stInf_0dpi_rep1", "C57BL6_1stInf_0dpi_rep2", "C57BL6_1stInf_5dpi_rep1", "C57BL6_1stInf_5dpi_rep2", 
             "NMRI_1stInf_3dpi_rep1", "NMRI_1stInf_3dpi_rep2", "NMRI_1stInf_5dpi_rep1", "NMRI_1stInf_5dpi_rep2", 
             "NMRI_1stInf_5dpi_rep3", "NMRI_1stInf_7dpi_rep1", "NMRI_1stInf_7dpi_rep2", "Rag_1stInf_0dpi_rep1", 
             "Rag_1stInf_0dpi_rep2", "Rag_1stInf_5dpi_rep1", "Rag_1stInf_5dpi_rep2", "C57BL6_2ndInf_5dpi_rep1", 
             "NMRI_2ndInf_0dpi_rep1", "NMRI_2ndInf_0dpi_rep2", "NMRI_2ndInf_3dpi_rep2", "NMRI_2ndInf_5dpi_rep1", 
             "NMRI_2ndInf_7dpi_rep1", "Rag_2ndInf_5dpi_rep1", "NMRI_2ndInf_7dpi_rep2")

spl <- sapply(strsplit(x,"_"), function(x)x[1])
lames <- data.frame(variable = x, strain = sapply(strsplit(x,"_"), function(x)x[1]), infection = sapply(strsplit(x,"_"), function(x)x[2]), dpi = sapply(strsplit(x,"_"), function(x)x[3]), replicate = sapply(strsplit(x,"_"), function(x)x[4]))

scCYT <- merge(cCYT, lames, by = "variable")
install.packages("hexbin")

library(ggplot2)
#start plot#
ggplot(scCYT, mapping = aes(x = rn , y = value, color = variable)) + 
  geom_point() +
  coord_flip()
#mouse strains#
ggplot(scCYT, mapping = aes(x = rn , y = value, color = strain)) + 
  geom_boxplot()
#infection type#
ggplot(scCYT, mapping = aes(x = rn , y = value, color = infection)) + 
  geom_boxplot()
#dpi#

