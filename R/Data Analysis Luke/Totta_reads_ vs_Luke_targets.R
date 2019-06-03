#libraries
library(httr)
library(RCurl)
library(data.table)
library(reshape2)
library(ggplot2)
#read in csv
targets <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_Mm_norm_counts.csv"
targets <- read.csv(text=getURL(targets))

#pick target codes
IFNg <- targets["ENSMUSG00000055170",]
TGFb <- targets["ENSMUSG00000002603",]
STAT6 <- targets["ENSMUSG00000002147",]
CXCL9 <- targets["ENSMUSG00000029417",]
IL10 <- targets["ENSMUSG00000016529",]
TNFa <- targets["ENSMUSG00000024401",]
CXCR3 <- targets["ENSMUSG00000050232",]
IRG6A <- targets["ENSMUSG00000054072",]
GBP2 <- targets["ENSMUSG00000028270",]
MyD88 <- targets["ENSMUSG00000032508",]
Ticam1 <- targets["ENSMUSG00000047123",]

CYT <- rbind(IFNg, TGFb, STAT6, CXCL9, IL10, TNFa, CXCR3, IRG6A, GBP2, MyD88, Ticam1)
setDT(CYT, keep.rownames = TRUE)
row.names(CYT) <- c("IFN-y", "TGF-b", "STAT6", "CXCL9", "IL-10", "TNFa", "CXCR3", "IRG6A", "GBP2", "MyD88", "Ticam1")

#make long data and rename rows
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

#use colnames insteda of writing them out#
x = colnames(targets)

#split names and exract into columns
lames <- data.frame(variable = x, strain = sapply(strsplit(x,"_"), function(f)f[1]), infection = sapply(strsplit(x,"_"), function(f)f[2]),
                    dpi = sapply(strsplit(x,"_"), function(f)f[3]), replicate = sapply(strsplit(x,"_"), function(f)f[4]))
scCYT <- merge(cCYT, lames, by = "variable")

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

