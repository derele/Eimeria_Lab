### Code to analysis
## 1) Correlation among DNA quantification methods
## 2) Course of infection by DNA quantification method

library("ggplot2")
#library("dplyr")
library("vegan")
library("gridExtra")
#library("tscount")
#library(gplots)
library(lme4)
library(MASS)
#library("stargazer")
library("plm")
library("AER")
library("nparLD")

source("R/1_Data_preparation.R")
source("R/2_qPCR_data_preparation.R")

###Get Panel A, B and C from script 3
if(!exists(c("A","B", "C"))){
  source("R/3_Data_analysis_qPCR_flotation.R")
}

str(sdt)

###### Question 1: Can DNA and oocyst predict weightloss & and is there an interaction between the 2?
# 1.1 Weightloss: as a rate
# 1.2 weightloss as the peak --> ALICE

#### classical glmm with ID as random effect
weightglmm <- lmer(weightloss~Genome_copies_gFaeces*OPG + (1|EH_ID), data= sdt)
summary(weightglmm)

#### linear regression, does not include time or ID.
myweight_lm <- lm(weightloss~Genome_copies_gFaeces*OPG, data=sdt)
coeftest(myweight_lm, vcov=vcovHC, type="HC1")

#### Linear regression with ID as fixed effect #######################
# eliminates the risk of biases due to omitted factors that vary across animals but not over time
# linear regression with ID as a fixed effect: reports dummy coefficients for each ID. which is annoying

myweight_lmf <- lm(weightloss~Genome_copies_gFaeces*OPG + EH_ID - 1, data=sdt)
#coeftest(myweight_lm, vcov=vcovHC, type="HC1")
summary(myweight_lmf)

#OLS to the demeaned data
##obtain subject demeaned data using non-zero dataframe 
sdt_Sdemeaned <- with(sdt.nozero, data.frame(weightloss=weightloss- ave(weightloss, EH_ID),
                                             Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                                             OPG=OPG-ave(OPG, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                                             dpi=dpi,
                                             EH_ID=EH_ID))
# estimate the regression
summary(lm(weightloss~Genome_copies_gFaeces*OPG, data=sdt_Sdemeaned))

# time demeaned data using non-zero dataframe
sdt_Tdemeaned <- with(sdt.nozero, data.frame(weightloss=weightloss- ave(weightloss, dpi),
                                             Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, dpi, FUN=function(x) mean(x, na.rm=T)),
                                             OPG=OPG-ave(OPG, dpi, FUN=function(x) mean(x, na.rm=T))))
# estimate the regression
summary(lm(weightloss~Genome_copies_gFaeces*OPG, data=sdt_Tdemeaned))

#### Regression with time and ID fixed effects###################################
# Controlling for variables that are constant across entities but vary over time
# can be done by including time fixed effects.
# The combined model (time and ID fixed effects) allows to eliminate bias from
#unobservables that change over time but are constant over entities and it controls
#for factors that differ across entities but are constant over time. Such models
#can be estimated using the OLS algorithm

sdt_STdemeaned <- with(sdt_Sdemeaned, data.frame(weightloss=weightloss- ave(weightloss, dpi),
                                                 Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, dpi, FUN=function(x) mean(x, na.rm=T)),
                                                 OPG=OPG-ave(OPG, dpi, FUN=function(x) mean(x, na.rm=T)),
                                                 dpi=dpi,
                                                 EH_ID=EH_ID))

summary(sdt_STdemeaned)

cor.test(sdt_STdemeaned$Genome_copies_gFaeces, sdt_STdemeaned$OPG, method="spearman") #--> Reported in the plot (Fig. 3 panel D)

cor.test (sdt.nozero$Genome_copies_gFaeces, sdt.nozero$OPG, method="spearman") #--> Reported in the plot (Fig. 3 panel A)

# plotting correlations
sdt_STdemeaned%>%
  ggplot(aes(y=(Genome_copies_gFaeces+ max(na.omit(Genome_copies_gFaeces))), x=(OPG+max(na.omit(OPG)))))+
  geom_point(shape=21, size=5, alpha=0.75, aes(fill= dpi))+
  scale_fill_manual(values = colores, guide= "none")+
  scale_x_log10(name = "log10 (Oocysts/g Faeces + max) \n (Flotation)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces + max) \n (qPCR)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  labs(tag= "d")+
  theme_bw()+
  stat_cor(label.y = 10,  label.x = 6.22, method = "spearman",
           aes(label= paste("rho","'='", ..r.., ..p.label.., sep= "~` `~")))+
  theme(text = element_text(size=16), legend.position = "none")+
  annotation_logticks()-> D

##To visualize it externally 
#ggsave(filename = "../../Rplots.pdf", D)

AB<- ggarrange(A, B, common.legend = TRUE, ncol = 2, nrow = 1)
CD<- ggarrange(C, D, ncol = 2, nrow = 1)

##Figure 3# Genome copies predicted by OPG overall and by dpi
#ggarrange(AB, CD, ncol = 1, nrow = 2)-> tmp.fig
#ggsave(file = "fig/Figure_3.pdf", tmp.fig, width = 13.5, height = 10.5, dpi = 600)


########Stopped script apllication here 
---------------------------------------------------------------------------------------------------

##obtain subject demeaned data using the complete dataframe 
sdt_Sdemeaned <- with(sdt, data.frame(weightloss=weightloss- ave(weightloss, EH_ID),
                                      Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                                      OPG=OPG-ave(OPG, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                                      dpi=dpi,
                                      EH_ID=EH_ID))

##obtain subject and time demeaned data using the complete dataframe 
sdt_STdemeaned <- with(sdt_Sdemeaned, data.frame(weightloss=weightloss- ave(weightloss, dpi),
                                                 Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, dpi, FUN=function(x) mean(x, na.rm=T)),
                                                 OPG=OPG-ave(OPG, dpi, FUN=function(x) mean(x, na.rm=T)),
                                                 dpi=dpi,
                                                 EH_ID=EH_ID))

# estimate the regression
sdtST=na.omit(sdt_STdemeaned[,c("Genome_copies_gFaeces", "OPG", "weightloss", "dpi")])
ST.lm=lm(weightloss~Genome_copies_gFaeces*OPG, data=sdtST)
int.lm=lm(weightloss~1, data=sdtST)

OPG.lm=lm(weightloss~Genome_copies_gFaeces, data=sdtST)
GC.lm=lm(weightloss~OPG, data=sdtST)
I.lm=lm(weightloss~Genome_copies_gFaeces+OPG, data=sdtST)

anova(ST.lm, GC.lm, test="LRT")
anova(ST.lm, OPG.lm, test="LRT")
anova(ST.lm, I.lm, test="LRT")


lrtest(GC.lm, ST.lm) #--> goes to table Predicted models of host health by Eimeria DNA and oocysts in faeces (Genome copies)
lrtest(OPG.lm, ST.lm) #--> goes to table Predicted models of host health by Eimeria DNA and oocysts in faeces (OPG)
lrtest(I.lm, ST.lm) #--> goes to table Predicted models of host health by Eimeria DNA and oocysts in faeces (Interaction)
summary(ST.lm) #--> goes to table Predicted models of host health by Eimeria DNA and oocysts in faeces (linear model) 

lrtest(ST.lm, int.lm)


library(relaimpo)
calc.relimp(I.lm)
calc.relimp(I.lm, rela=TRUE)

## Preliminary plots
ggplot(sdtST, aes(x=log(1+Genome_copies_gFaeces), y=weightloss))+
  geom_point(size=2, alpha=0.8)+
  annotate("text", x=10, y=15, label="F=25.1, p<0.001", hjust = "left")+
  labs(x="Genome copies (log+1)", y="Weight loss")+
  theme_classic()-> tmp.fig.1

ggplot(sdtST, aes(x=log(1+OPG), y=weightloss))+
  geom_point(size=2, alpha=0.8)+
  annotate("text", x=10, y=15, label="F=3.9, p=0.05", hjust = "left")+
  labs(x="OPG (log 1+)", y="Weight loss")+
  theme_classic()-> tmp.fig.2


# plot residuals
d <- sdtST[c("OPG", "Genome_copies_gFaeces", "weightloss", "dpi")]

d$predicted <- predict(ST.lm)   # Save the predicted values
d$residuals <- residuals(ST.lm) # Save the residual values

#https://drsimonj.svbtle.com/visualising-residuals

# change: color by day of max weight loss
d <- d %>% 
  gather(key = "iv", value = "x", -weightloss, -predicted, -residuals, 
         -dpi) # Get data into shape
d$dpi <- as.factor(d$dpi)

ResDemSusana_1 <- d[d$iv %in% "Genome_copies_gFaeces",] %>% 
  ggplot(aes(x = x, y = weightloss)) +  # Note use of `x` here and next line
  geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
  geom_point(aes(fill = dpi, alpha = abs(residuals)), size = 2.5, shape=21, col=1) +
  scale_alpha(range = c(0.1, 1), guide = F) +
  geom_point(aes(y = predicted), shape = 1) +
  xlab("Genome copies per gram of faeces")+
  ylab("Weight loss relative to DPI 0 (%)") +
  theme_bw() + 
  labs(tag = "c", fill= "DPI")+
  theme(text = element_text(size=16), legend.position = "top")+ 
  guides(fill = guide_legend(nrow = 1))

##Extract legend
legend <- cowplot::get_legend(ResDemSusana_1)

##Remove legend
ResDemSusana_1 <- d[d$iv %in% "Genome_copies_gFaeces",] %>% 
  ggplot(aes(x = x, y = weightloss)) +  # Note use of `x` here and next line
  geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
  geom_point(aes(fill = dpi, alpha = abs(residuals)), size = 2.5, shape=21, col=1) +
  scale_alpha(range = c(0.1, 1), guide = F) +
  geom_point(aes(y = predicted), shape = 1) +
  xlab("Genome copies per gram of faeces")+
  ylab("Weight loss relative to DPI 0 (%)") +
  theme_bw() + 
  labs(tag = "c")+
  theme(text = element_text(size=16), legend.position = "none")

ResDemSusana_2 <- d[d$iv %in% "OPG",] %>% 
  ggplot(aes(x = x, y = weightloss)) +  # Note use of `x` here and next line
  geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
  geom_point(aes(fill = dpi, alpha = abs(residuals)), size = 2.5, shape=21, col=1) +
  scale_alpha(range = c(0.1, 1), guide = F) +
  geom_point(aes(y = predicted), shape = 1) +
  xlab("Oocysts per gram of faeces")+
  ylab("Weight loss relative to DPI 0 (%)") +
  theme_bw() +
  labs(tag = "d")+
  theme(text = element_text(size=16), legend.position = "none")

ResDemSusana <- ggarrange(ResDemSusana_1, ResDemSusana_2, ncol = 2, nrow = 1)

saveRDS(ResDemSusana, file="fig/ResDemSusana.rds")

saveRDS(legend, file="fig/legend.rds")

ggplot(sdtST, aes(x = logGC, y = weightloss)) +  # Set up canvas with outcome variable on y-axis
  geom_segment(aes(xend = logGC, yend = predicted), alpha = .2) +  # alpha to fade lines
  geom_point() +
  geom_point(aes(y = predicted), shape = 1) +
  theme_classic()-> tmp.fig.3  # Add theme for cleaner look

# interestingly when we include time fixed effects (controls for effects that change
#over time and not because of ID) we don't get a significant effect for OPG and the
#interaction between Genome_copies and OPG.

### non parametric analysis of longitudinal data in factorial experiments
### Brunner et al. 2002
# LD-F1 design fers to the experimental design with one sub-plot factor
#(longitudinal data forone homogeneous group of subjects).
#1 hypothesis: no time effect

ex.f1 <- ld.f1(y=mytab$weightloss, time=mytab$dpi, subject=mytab$EH_ID, time.order=c(0,1,2,3,4,5,6,7,8,9,10))

summary(ex.f1)
plot(ex.f1)
print(ex.f1)

########################################
#######################################
########### Question 2: do DNA and OOcyst predict infection stage?

mytab$inf <- NA

for (i in 1:nrow(mytab)){
  if (mytab$dpi[i]== 1|| mytab$dpi[i]==2||mytab$dpi[i]==3||mytab$dpi[i]==4) {
    mytab$inf[i] <- "early"
  } else if (mytab$dpi[i]==0){
    mytab$inf[i] <- "NI"
  } else if (mytab$dpi[i]==5||mytab$dpi[i]==6||mytab$dpi[i]==7) {
    mytab$inf[i] <- "peak"
  } else if (mytab$dpi[i]==8||mytab$dpi[i]==9||mytab$dpi[i]==10) {
    mytab$inf[i] <- "late"
  } else {mytab$inf[i] <- NA}
}


##########################################

library(ggpubr)

# granger causality
# hypothesis: DNA does not cause OPG for all individuals

pgrangertest(OPG~Genome_copies_mean, data=mytab, index=c("EH_ID", "dpi"))

pgrangertest(Genome_copies_mean~OPG, data=mytab, index=c("EH_ID", "dpi"))


# plm does not report any dummy variable which is cool
myweigh_plm <- plm(weightloss~Genome_copies_mean * OPG,
                   data=mytab,
                   index=c("EH_ID", "dpi"),
                   model="within",
                   effect = "twoways")


## testing with time series: create a mean for each individualö
meantab <- sdt[,c("Genome_copies_mean", "OPG", "dpi", "weightloss")]
meantab <- unique(meantab)

head(meantab)

meantab$dpi <- as.factor(meantab$dpi)

meantab  <- aggregate(meantab[, c(1, 2, 4)], list(meantab$dpi), mean, na.action=na.exclude)

tab <- meantab %>% group_by(dpi)

mytab

mtab <- with(mytab, data.frame(weightloss=ave(weightloss, dpi),
                               Genome_copies_mean=ave(Genome_copies_mean, dpi, FUN=function(x) mean(x, na.rm=T)),
                               OPG=ave(OPG, dpi, FUN=function(x) mean(x, na.rm=T))))

mtab <- unique(mtab)
mtab$dpi <- seq(0,10,1)

# let's just pretend that OPG at day 1 and 2 are 0
mtab[2, 3] <- 0
mtab[3, 3] <- 0


acf(na.omit(mtab$Genome_copies_mean))

acf(na.omit(mtab$OPG))

acf(na.omit(mtab$weightloss))

##### to do:
## stage of infection
## residual plots for LM's for supplementary material
