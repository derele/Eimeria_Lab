Challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")


# Chapter 2. T-tests and Wilcoxon

# Load package libraries

library(dunn.test)             # R package for Dunn's test post hoc comparisons of Kruskal-Wallis data
library(ggplot2)               # R package library for better graphic procedures
library(ggpubr)                # R package for ggarrange() command, producing multi-panel graphics
library(gridExtra)             # pre-requisite for ggpubr library
library(lme4)                  # R package for linear mixed models
library(multcomp)              # R package to make multiple comparison adjustments for ANOVA-type problems
library(PMCMR)                 # R Package for post hoc tests of Kruskal-Wallis and Friedman's test

# Fig 2a. Student's T-test

PCRneg = rnorm(n=18,mean=15,sd=2.2)                    # Generate random normal data for PCRneg and PCRpos groups
PCRpos = rnorm(n=18,mean=23,sd=2.4)
TotIgG = c(PCRneg,PCRpos)                              # Combine PCRneg and PCRpos into single column of data
Infect = c(rep("PCR neg",18),rep("PCR pos",18))        # Create new variable "Infect" to label TotIgG as PCRneg or PCRpos
dat2a  = data.frame(TotIgG,Infect)                     # Create a small data table for Fig 2a

t.test(x=PCRneg,y=PCRpos)                              # Compute T-test at the "command line"

result = t.test(formula=TotIgG ~ Infect,data=dat2a)    # Store T-test in "result" variable
print(result)                                          # Display T-test result
print(result$statistic)                                # Recall only the T-statistic from the results

# Create Fig 2a

fig2a  = ggplot(dat2a, aes(x=Infect,y=TotIgG,group=Infect)) + theme_classic() + theme(legend.position="none") + 
  labs(title="A",x="",y="Total IgG") + geom_dotplot(binaxis="y",binwidth=0.5,stackdir="center",fill="grey") + 
  geom_text(x=1,y=25,label=paste("p = ",formatC(result$p.value,format="e",digits=2)),colour="black") + 
  stat_summary(fun.y=mean,fun.ymin=mean,fun.ymax=mean,geom="crossbar",colour="black") + 
  stat_summary(fun.ymin=function(y) mean(y) - sd(y),fun.ymax=function(y) mean(y) + sd(y),geom="errorbar",colour="black",width=0.3)

# Fig 2b. Mann-Whitney Test

WT        = abs(rnorm(n=50,mean=0.7,sd=3.4))          # Generate random data and re-organize as done in last example
KO        = abs(rnorm(n=80,mean=0.1,sd=0.8))          # Low means and abs() function make the data non-normal
PctTBXpos = c(WT,KO)                                  # Combine responses for WT and KO strains
Strain    = c(rep("WT",50),rep("KO",80))              # Define strain label data
dat2b     = data.frame(PctTBXpos,Strain)              # Create data frame

wilcox.test(x=WT,y=KO)                                # Compute Mann-Whitney test at command line

result = wilcox.test(formula=PctTBXpos ~ Strain,data=dat2b)    # Store and recall MW-test as shown previously
print(result)
print(result$p.value)

# Create Fig 2b

fig2b  = ggplot(dat2b, aes(x=Strain,y=PctTBXpos,group=Strain)) + theme_classic() + theme(legend.position="none") + 
  labs(title="B",x="",y="%TBX+") + expand_limits(y=c(0,10)) + geom_dotplot(binaxis="y",binwidth=0.15,stackdir="center",fill="grey") + 
  geom_text(x=1,y=8,label=paste("p = ",formatC(result$p.value,format="e",digits=2)),colour="black") + 
  geom_boxplot(fill="white",col="black",alpha=0) 

