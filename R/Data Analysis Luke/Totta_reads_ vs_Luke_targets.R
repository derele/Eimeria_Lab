#personal path
#DF <- read.csv(file = "../luke/Documents/Mm_norm_counts.csv")

#IZW path
DF <- read.csv(file = "../luke/Documents/Mm_norm_counts.csv")

#HU path
#DF <- read.csv(file = "../luke/Documents/Mm_norm_counts.csv")

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

setDT(CYT, keep.rownames = TRUE)
                
#install.packages("data.table")
library(data.table)
library(reshape2)

row.names(CYT) <- c("IFN-y", "TGF-b", "STAT6", "CXCL9", "IL-10", "TNFa", "CXCR3", "IRG6A", "GBP2", "MyD88", "Ticam1")

cCYT <- melt(CYT, id.vars = c("rn"))


library(ggplot2)
ggplot(cCYT, aes(cCYT$variable, cCYT$value, color =cCYT$rn)) +
  geom_point() +
  coord_flip()
