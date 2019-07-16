## Analysis of experiment performed on NMRI mice in January 2017
## Emanuel Heitlinger
## Alice Balard
## Lubomir Bednar
library(ggplot2)
library(lme4)
library(lmerTest)
library(lsmeans) # "mvtnorm" dependency is not available through base repos,
#use install.packages("https://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-8.tar.gz", repos=NULL)
library(devtools) # devtools::install_github("strengejacke/strengejacke")
library(strengejacke)
library(plyr)
library(coin)
library(ggeffects)
library(httr)
library(RCurl)
library(Rmisc)

#### ### get the data --------------------------------------------

## for control: general experimental setup -------------------
#stab <- read.csv("./Experiment_Table_raw_NMRI_Jan2017.csv")
stabURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/E1_012017_Eim_Experiment_Table_raw_NMRI.csv"
stab <- read.csv(text=getURL(stabURL))
stab <- subset(stab, !stab$inf.strain%in%"EI70", drop=TRUE)
stab$inf.strain <- as.factor(as.character(stab$inf.strain))

## weight --------------------------------------------------------
weightURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_record.csv"
weight <- read.csv(text=getURL(weightURL))
## weight was only obtained for mice dissected at or after 7dpi
## correcting the names
names(weight)[names(weight)%in%"Mouse.ID"] <- "EH_ID"
names(stab)[names(stab)%in%"EH_ID"] <- "EH_ID"

## EI70 infections were not followed up as likely double infections
weight <- weight[!weight$inf.strain%in%"EI70", ] 

## remove non meaningful columns
weight <- weight[, !colnames(weight)%in%c("Mouse.number", "Date.of.Birth")]

#set weight as percentages
p.weight <- apply(weight[, grepl("^Day.*g", colnames(weight))], 2,
                  function (x) {
                      (x/weight$Day.1_.g)*100
                  })
#append the weight percentages
colnames(p.weight) <- gsub(".g", ".p", colnames(p.weight))

weight <- cbind(weight, p.weight[, 2:ncol(p.weight)])

#create long format data
weight.long <- reshape(weight,
                       direction = "long",
                       idvar = "EH_ID", ids = EH_ID,
                       varying = list(grep(".p$", colnames(weight), value=TRUE)),
                       timevar="dpi_of_perc",
                       v.names = "perc_of_dpi1", 
                       times = grep(".p$", colnames(weight), value=TRUE))
#extract the  weight-dpi information from names
weight.long$dpi_of_perc <- as.numeric(gsub("Day\\.?(\\d+)_\\.p", "\\1",
                                           weight.long$dpi_of_perc))

weight.long <- weight.long[, !grepl("^Day.", names(weight.long)) ]

### oocysts ---------------------------------------------------------
#oocysts <- read.csv("./Clean_oocyst_data.csv")
oocystsURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_Clean_oocyst_data.csv"
oocysts <- read.csv(text=getURL(oocystsURL))

#add a column for total oocyst numbers
oocysts$Total.oocysts.g <- ((oocysts$Count..8.Neubauer.squares. / 8)*
                            10000 * 2) / oocysts$used.in.flotation..g.

## correcting the names
names(oocysts)[names(oocysts)%in%c("Sample.ID", "dpi")] <- c("EH_ID", "dpi_count")
oocysts <- merge(oocysts, stab, all.y=TRUE)

all.data <- merge(oocysts[, c("EH_ID", "dpi_count", "Total.oocysts.g",
                              "dpi.diss", "inf.strain")],
                  weight.long,
                  by.x=c("EH_ID", "dpi_count", "dpi.diss", "inf.strain"),
                  by.y=c("EH_ID", "dpi_of_perc", "dpi.diss", "inf.strain"),
                  all=TRUE)

## we can observe some pattern:
## we have any measurements only for mice killed at or after 7dpi
## we have oocyst counts only for mice kileed at or after 7dpi

## in other words: 
## weight is not reported (measured) for mice killed at 3 and 5dpi
## oocyst counts are not repored (measured) for mice killed at 3, 5 and 7 dpi


############# --------------- qPCR for DNA -------------------

#Rtissue <- read.csv("./Eimeria-NMRI_Relative quantification_clean.csv")
RtissueURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_NMRI_Relative_quantification_clean.csv"
Rtissue <- read.csv(text=getURL(RtissueURL))

## only the means
RtMeans <- Rtissue[seq(1, nrow(Rtissue), by=2),
                   c("Sample", "Cq.Mean", "Cq.Mean.1")]

RtMeans <- Rtissue[seq(1, nrow(Rtissue), by=2),
                   c("Sample", "Cq.Mean", "Cq.Mean.1")]
names(RtMeans) <- c("EH_ID", "Mouse_gDNA", "Eimeria_mDNA")

RtMeans$EH_ID <- toupper(RtMeans$EH_ID)
## LM0065 was measured twice with the same outcome
RtMeans <- RtMeans[!duplicated(RtMeans$EH_ID),]
RtMeans <- merge(RtMeans, stab, all=TRUE)
RtMeans$dpi.diss <- gsub("dip$", "dpi", RtMeans$dpi.diss)
RtMeans$dpi_count<- as.numeric(gsub("dpi", "", RtMeans$dpi.diss))

all.data <- merge(RtMeans, all.data, all=TRUE)

## remove the mice with no data whatsoever
all.data <- all.data[!is.na(all.data$dpi_count), ]

all.data$PH.delta <- all.data$Mouse_gDNA - all.data$Eimeria_mDNA

max.neg <- max(all.data[all.data$inf.strain%in%"Uninf", "PH.delta"],
               na.rm=TRUE)
min.neg <- min(all.data$PH.delta, na.rm=TRUE)

LLD <- mean(all.data[all.data$inf.strain%in%"Uninf", "PH.delta"], na.rm=TRUE)+
    2*(sd(all.data[all.data$inf.strain%in%"Uninf", "PH.delta"], na.rm=TRUE))

no.log.lld <- 2^LLD
## per 100 transcripts (Dreisatz ;-))
(no.log.lld*1/(1-no.log.lld))*100

maxPHd <- max(all.data$PH.delta, na.rm=TRUE)
no.log.max.PHd <- 2^maxPHd
(no.log.max.PHd*1/(1-no.log.max.PHd))*100

max(all.data[all.data$inf.strain%in%"EI64", "PH.delta"], na.rm=TRUE)

## generate summary data for plotting ---------------------------
w.means <- ddply(all.data,
                 c("dpi_count", "inf.strain"),
                 summarize,
                 N    =  sum(!is.na(perc_of_dpi1)),
                 sd    =  sd(perc_of_dpi1, na.rm=TRUE),
                 mean = mean(perc_of_dpi1, na.rm=TRUE))


o.means <- ddply(all.data[!all.data$inf.strain%in%"Uninf", ],
                 c("dpi_count", "inf.strain"),
                 summarize,
                 N    =  sum(!is.na(Total.oocysts.g)),
                 sd    =  sd(Total.oocysts.g, na.rm=TRUE),
                 mean = mean(Total.oocysts.g, na.rm=TRUE))


### plot summary data -------------------------------------


pdf("figures/weight_retained.pdf", width=8, height=4)
ggplot(w.means, aes(dpi_count, mean, group=inf.strain,
                    color=inf.strain)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2,
                  position = "dodge") +
    geom_text(aes(x=dpi_count,
                  y=50+as.numeric(w.means$inf.strain)*2.5, label=N))+
    scale_y_continuous("mean weight retained as percent of 1 dpi")+
    scale_x_continuous("days post infection (dpi)", breaks=1:11,
                       labels=1:11, limits=c(2, 11.5))+
    scale_colour_brewer(palette = "Dark2")+
    theme_bw()
dev.off()


pdf("figures/oocysts_shed.pdf", width=8, height=4)
ggplot(o.means, aes(dpi_count, mean, group=inf.strain,
                    color=inf.strain)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2,
                  position = "dodge") +
    geom_text(aes(x=dpi_count,
                  y=-1.2*10^5+
                      (as.numeric(as.factor(o.means$inf.strain))*10^5)*0.25,
                  label=N))+
    scale_y_continuous("number of oocysts shed per gramm feces")+
    scale_x_continuous("days post infection (dpi)", breaks=1:11,
                       labels=1:11, limits=c(2, 11.5))+
    scale_colour_brewer(palette = "Dark2")+
    theme_bw()
dev.off()


pdf("figures/Genomic_pPCR.pdf", width=8, height=4)
ggplot(all.data, aes(dpi_count, PH.delta, color=inf.strain)) +
    geom_point() +
    geom_smooth(se=FALSE) +
    scale_y_continuous("delta parasite-host DNA")+
    scale_x_continuous("days post infection (dpi)",
                       breaks=c(3, 5, 7, 9, 11),
                       labels=c("3dpi", "5dpi", "7dpi", "9dpi", "11dpi")) +
    scale_colour_brewer(palette = "Dark2")+
    geom_ribbon(aes(ymax = LLD, ymin = min.neg-0.2, color=NA),
                data = , alpha=.2)+
    geom_hline(yintercept = LLD, linetype=2)+
    theme_bw()
dev.off()


## ------------------------- statistical testing of weight and oocyst data 

## testing weight loss ----- ----------------------------------
test.day.for.strain <- function(day, strain){
    foo <- subset(all.data,
                  inf.strain%in%strain & dpi_count%in%day)
    test <- wilcox_test(perc_of_dpi1~as.factor(inf.strain),
                   data=foo,
                   distribution = "exact")
    list(test, length(foo$perc_of_dpi1[!is.na(foo$perc_of_dpi1)]))
}

test.day.for.strain(day="4", strain=c("EI64" ,"Uninf"))
test.day.for.strain("5", c("EI64" ,"Uninf"))

test.day.for.strain("8", c("Efwild" ,"Uninf"))
test.day.for.strain("9", c("Efwild" ,"Uninf"))

test.day.for.strain("8", c("Eflab" ,"Uninf"))
test.day.for.strain("9", c("Eflab" ,"Uninf"))

foo <- subset(all.data,
              inf.strain%in%"EI64" & dpi_count%in%"5"|
              inf.strain%in%"Efwild" & dpi_count%in%"9")

length(foo$perc_of_dpi1[!is.na(foo$perc_of_dpi1)])

wilcox_test(perc_of_dpi1~as.factor(inf.strain),
            data=foo,
            distribution = "exact")

foo <- subset(all.data,
              inf.strain%in%"EI64" & dpi_count%in%"5"|
              inf.strain%in%"Eflab" & dpi_count%in%"9")

length(foo$perc_of_dpi1[!is.na(foo$perc_of_dpi1)])

wilcox_test(perc_of_dpi1~as.factor(inf.strain),
            data=foo,
            distribution = "exact")

### oocyst shedding ----------------------------------------------

## ferW reduction from dpi 6 to 7
ferWsub67 <- subset(all.data,
                    inf.strain%in%"EI64"&dpi_count%in%c(6,7))

nrow(ferWsub67)

wilcox_test(Total.oocysts.g~as.factor(dpi_count),
            data=ferWsub67,
            distribution = "exact")


## peak day difference for falciformis
falMaxWL <- subset(all.data,
                   inf.strain%in%"Eflab"&dpi_count%in%8|
                   inf.strain%in%"Efwild"&dpi_count%in%9)

nrow(falMaxWL)

wilcox_test(Total.oocysts.g~as.factor(inf.strain),
            data=falMaxWL,
            distribution = "exact")

## peak day differnce between fal and fer
falferMax <- subset(all.data,
##                    inf.strain%in%"Eflab"&dpi_count%in%8|
                    inf.strain%in%"Efwild"&dpi_count%in%9|
                    inf.strain%in%"EI64"&dpi_count%in%6)

wilcox_test(Total.oocysts.g~as.factor(inf.strain),
            data=falferMax,
            distribution = "exact")

## peak day differnce between fal and fer
falferMax <- subset(all.data,
                    inf.strain%in%"Eflab"&dpi_count%in%8|
##                    inf.strain%in%"Efwild"&dpi_count%in%9|
                    inf.strain%in%"EI64"&dpi_count%in%6)

wilcox_test(Total.oocysts.g~as.factor(inf.strain),
            data=falferMax,
            distribution = "exact")


## ## get the sample size
## by(stab, stab$inf.strain, function (x)
##     table(is.na(x$Eimeria_mDNA-x$Mouse_gDNA)))

## Efal <- subset(stab, inf.strain%in%c("Efwild", "Eflab") &
##                      dpi==9)

## ## Trying to predict weight loss with parasite DNA abundance...
## d7 <- lm(X7DPI_count~PH.delta, data=Efal)
## summary(d7)

## d7 <- lm(Day7_.p~PH.delta, data=Efal)
## summary(d7)
## ## .. not like this but maybe arrange the data and do it in a more
## ## sophisticated way (predict day x influencing e.g. day x+1)

############ - Gene expression data -----------------------

GE.files <- list.files(pattern="^GE_")

GeMeans.l <- lapply(GE.files, function (file) {
    data <- read.csv(file)
    data <- data[seq(1, nrow(Rtissue), by=2), c("Sample", "Gene", "NE")]
    data[!is.na(data$Gene), ]
})

GeMeans <- Reduce(rbind, GeMeans.l)
GeMeans$Sample <- toupper(GeMeans$Sample)

### Correction
## removing an empty row
GeMeans <- GeMeans[!GeMeans$Gene%in%"",]
GeMeans$Gene <- toupper(GeMeans$Gene)
## standard naming
names(GeMeans)[names(GeMeans)%in%"Sample"] <- "EH_ID"

## wide dateset for merging in overall table
GeMeans.wide <- reshape(GeMeans, timevar = "Gene", idvar = "EH_ID", direction = "wide")

M <- merge(GeMeans, stab, all=TRUE)
M.wide <- merge(GeMeans.wide, stab, all=TRUE)

pdf("figures/Cytokines.pdf", width=12, height=4)
ggplot(subset(M, nchar(M$Gene)>2), aes(dpi, NE, color=inf.strain)) +
    geom_jitter(width=0.2) +
    geom_smooth(se=FALSE) +
    scale_x_continuous(breaks=c(3, 5, 7, 9, 11), labels=c("3dpi", "5dpi", "7dip", "9dpi", "11dpi")) +
    facet_wrap(~Gene, scales="free_y", nrow=2)+
    scale_colour_brewer(palette = "Dark2") +
    theme_bw()
dev.off()


pdf("figures/Cytokines_vs_Inf.pdf", width=12, height=4)
ggplot(subset(M, nchar(M$Gene)>2), aes((Mouse_gDNA-Eimeria_mDNA), NE, color=inf.strain)) +
    geom_jitter(width=0.2) +
##    geom_smooth(se=FALSE) +
##    scale_x_continuous(limits=c(-5, 9)) +
    facet_wrap(~Gene, scales="free_y", nrow=2)
dev.off()


## Contrasting against Eflab
modCXCL9 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"CXCL9"))
summary(modCXCL9)

modIL10 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL10"))
summary(modIL10)

modIL12 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL12"))

modIL6 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL6"))
summary(modIL6)

modINFG <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"INFG"))
summary(modINFG)

modSTAT6 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"STAT6"))
summary(modSTAT6)

modTGFB <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"TGFB"))
summary(modTGFB)

modTNFA <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"TNFA"))
summary(modTNFA)

sjt.lmer(modCXCL9, modIL10, modIL12, modIL6,
         modINFG, modSTAT6, modTGFB, modTNFA, file="table_VS_Eflab(itercept).html",
         depvar.labels=c("CXCL9", "IL10", "IL12", "IL6",
                         "INFG", "STAT6", "TGFB", "TNFA"))

## Now contrasting against negative control
M$inf.strain = factor(M$inf.strain, levels(M$inf.strain)[c(4,1:3)])

modCXCL9 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"CXCL9"))
summary(modCXCL9)

modIL10 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL10"))
summary(modIL10)

modIL12 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL12"))
summary(modIL12)

modIL6 <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"IL6"))
summary(modIL6)

modINFG <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"INFG"))
summary(modINFG)

modSTAT6 <- lmer(NE~inf.strain  +(1|dpi.diss), data=subset(M, M$Gene%in%"STAT6"))
summary(modSTAT6)

modTGFB <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"TGFB"))
summary(modTGFB)

modTNFA <- lmer(NE~inf.strain + (1|dpi.diss), data=subset(M, M$Gene%in%"TNFA"))
summary(modTNFA)

sjt.lmer(modCXCL9, modIL10, modIL12, modIL6,
         modINFG, modSTAT6, modTGFB, modTNFA, file="table_VS_non_infected(itercept).html",
         depvar.labels=c("CXCL9", "IL10", "IL12", "IL6",
                         "INFG", "STAT6", "TGFB", "TNFA"))

hist <- read.csv("./Histo.csv")

hist$mMLS <- rowMeans(hist[, grepl("^MLS", colnames(hist))])
hist$sMLS <- rowSums(hist[, grepl("^MLS", colnames(hist))])

names(hist)[names(hist)%in%"Sample.ID"] <- "EH_ID"
stab <- merge(stab, hist, all=TRUE)

## no need to subset because Enas does not report MLS for uninfected
## subset(stab, inf.strain%in%c("Efwild", "Eflab", "EI64"))
modS <- glm(sMLS~PH.delta*inf.strain, 
            data=stab, family=poisson)
summary(modS)

sjt.glm(modS, file="table_glm_PHvsMLS.html")


## An analysis involving dpi is for now left out (not included).
##

## modSRAND <- glmer(sMLS~PH.delta+inf.strain +
##                       (1+PH.delta|dpi) +
##                      (1+PH.delta|inf.strain:dpi),
##                  data=stab,
##                family=poisson)

## summary(modSRAND)
## sjt.glmer(modSRAND, file="table_glmm_PHvsMLS.html")
## sjp.glmer(modSRAND)

## modSFIX <- glm(sMLS~PH.delta+(inf.strain*dpi.diss),
##                data=stab,
##                family=poisson)
## summary(modSFIX)

pdf("figures/PHvsMLS.pdf", width=8, height=5)
foo <- ggpredict(modS, terms=c("PH.delta", "inf.strain"))
ggplot() +
    geom_line(data=foo, aes(x, predicted, group=group, color=group)) +
    geom_ribbon(data=foo, aes(x, ymin=conf.low, ymax=conf.high,
                              fill=group), alpha=0.4) +
    geom_point(data=stab, aes(x=PH.delta, y=sMLS, color=inf.strain,
                              shape=dpi.diss), size=3) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_continuous("parasite-host DNA log-ratio",
                       limits=c(-10, max(foo$x)))+
    scale_y_continuous("lesions observed in histological sections")+
    theme_bw()
dev.off()

## Infiltration scores at different dpi 
tapply(stab$Score1, stab$inf.strain:stab$dpi.diss, print)

FlamMod <- lmer(Score1~inf.strain + (1|dpi.diss), data=stab)
summary(FlamMod)

difflsmeans(FlamMod, test.effs = "inf.strain")

FlamModmls <- lmer(Score1~sMLS*inf.strain + (1|dpi.diss), data=stab)
summary(FlamModmls)

## which peak is first
diffthat <- by(all.data, all.data$EH_ID, function (x){
    minW <- min(x[, "perc_of_dpi1"], na.rm=TRUE)
    WminW <- which(x[, "perc_of_dpi1"] == minW)
    dayminW <- mean(x[WminW, "dpi_count"])
    maxO <- max(x[, "Total.oocysts.g"], na.rm=TRUE)
    WmaxO <- which(x[, "Total.oocysts.g"] == maxO)
    daymaxO <- mean(x[WmaxO, "dpi_count"])
    cbind(OO=daymaxO, WL=dayminW)
})

diffthis <- do.call(rbind, diffthat)
rownames(diffthis) <- unique(all.data$EH_ID)
diffthis <- diffthis[rowSums(is.na(diffthis))==0, ]
diffthis <- merge(stab, diffthis, by.x="EH_ID", by.y=0)

diffthis$EH_ID <- as.factor(as.character(diffthis$EH_ID))

library(reshape)
difflong <- melt(diffthis)

pdf("figures/peaks.pdf", width=8, height=4)
ggplot(difflong, aes(y=EH_ID, x=value, color=variable)) +
    geom_point(size=4, alpha=0.5) +
    facet_wrap(~inf.strain) +
    scale_x_continuous("days post infection", breaks=3:11) +
    scale_color_manual(name="Measured peak in",
                         labels=c("oocyst shedding", "weight loss"),
                         values=c("red", "blue")) +
    scale_y_discrete("mouse ID") +
    theme_bw()
dev.off()


pdf("figures/peaks_alt.pdf", width=6, height=6)
pd <- position_dodge(0.15)
ggplot(data=difflong, aes(x=value, y=variable, group=EH_ID)) +
    geom_line(color="red", position=pd)+
    geom_point(position=pd) +
    facet_wrap(~inf.strain) +
    scale_x_continuous("days post infection", breaks=3:11) +
    scale_y_discrete(label=c("peak oocyst\nshedding", "peak weight\nloss"))+
    coord_flip() +
    theme_bw()
dev.off()


