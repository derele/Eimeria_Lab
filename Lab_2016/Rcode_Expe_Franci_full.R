#########################
## Info about experiment:
## 25 mice (LM_0104 was deleted as she died at dpi1)
## 272 mice.days data points
#########################
library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

setwd("/home/alice/Schreibtisch/Franci+Phoungs_Projects/Inf_expe_Franci_April_2017/")

######### I Experiment deign:
# Info on the mouse genotype AND on the infection strain:
ExpePlanDF <- read.csv("Experiment_Table_raw_Ploen_May2017.csv")
# Print the expe design:
as.data.frame.matrix(table(ExpePlanDF$strain, ExpePlanDF$Inf_strain))

# Load data:
mydata_oocysts <- read.csv("Oocyst_coding_PloenEx_May2017.csv")

# Remove LM_0104:
mydata_oocysts <- mydata_oocysts[which(mydata_oocysts$X != "LM0104"),]

# Check:
table(is.na(mydata_oocysts[grep("raw.count", names(mydata_oocysts))]))
table(is.na(mydata_oocysts[grep("fec.weight", names(mydata_oocysts))]))
## Error corrected (same amount of fecal weight and oocysts counts)
#######################################################################################!!!!!!!!!!!!!!!!!!!!!!!!!

mydata_weight <- read.csv("mouse weight.csv")

# Weight of the feces separated for PCR:
subsampleDF <- read.csv("Feces_kept_separated_for_PCR.csv")
subsampleDF <- na.omit(subsampleDF)

# Calculate how many animals at start:
summary_table_at_dpi <- function(n){
  A = data.frame(mouse_number = ExpePlanDF$EH_id, ExpePlanDF$strain, ExpePlanDF$Inf_strain)
  B = data.frame(mouse_number = mydata_weight$mouse_no., mydata_weight[paste0("mouse_weight.g._dpi", n)])
  C = merge(A, B, by = "mouse_number")
  C = na.omit(C)
  table(C$ExpePlanDF.strain, C$ExpePlanDF.Inf_strain)
}

summary_table_at_dpi(11)
 
######### II Oocysts shedding:

# Initialise Newdata
Newdata <- data.frame(NULL)

# Run on all the dpi in 1 go!
for (DPI in 0:11){
  mydata_oocyststemp <- cbind(mydata_oocysts$X,mydata_oocysts[, grep(paste0("dpi.", paste(DPI, "", sep="\\.")), colnames(mydata_oocysts))])
  names(mydata_oocyststemp)<- c("mouse.ID","label","fec.weight","dil", "raw.count")
  # add column
  mydata_oocyststemp$dpi<- DPI
  Newdata <- rbind(Newdata, mydata_oocyststemp)
}

# Remove useless objects:
rm(mydata_oocyststemp, DPI)

# Rm the dil column, useless:
Newdata <- Newdata[!names(Newdata) %in% "dil"]

# Rm empty lines (when the mice died, no more counting) :
Newdata <- na.omit(Newdata)

# Merge the 2 dataframe by labels:
names(subsampleDF)[1] <- "label"

Oo_Df <- na.omit(merge(Newdata, subsampleDF, by = "label"))

## Calculate the oocysts in the subsampling:
Oo_Df$Total_feces <- Oo_Df$subsample_g + Oo_Df$fec.weight

# Calculate the oocysts number in 1 mL:
Oo_Df$abscount <- Oo_Df$raw.count * 10000

# Calculate the oocyst per gram in the processed part of the feces::
Oo_Df$oocysts.per.g <- Oo_Df$abscount / Oo_Df$fec.weight

#Put my dataframe at the good format:
Oo_Df <- Oo_Df[c("mouse.ID", "dpi", "oocysts.per.g")]

# Add the info on the mouse genotype AND on the infection strain:
names(Oo_Df)[1] <- "EH_id"

Oo_Df <- merge(Oo_Df, ExpePlanDF, by = "EH_id")

Oo_Df <- data.frame(EH_id = Oo_Df$EH_id, 
                    dpi = Oo_Df$dpi, 
                    oocysts.per.g = Oo_Df$oocysts.per.g, 
                    Inf_strain = Oo_Df$Inf_strain,
                    strain = Oo_Df$strain)

# Visual check-up:
table(Oo_Df$EH_id, Oo_Df$dpi)

# To keep all plot with the same theme (for Alice s TAC):
theme_alice <- theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        text = element_text(size = 20))

## PLOT mice strains:
ggplot(Oo_Df, aes(x=dpi, y=oocysts.per.g, group = strain, col = strain))+
  geom_smooth(aes(fill = strain), alpha = 0.2)+
  ggtitle("Oocyst count along the experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "Loess smoothing + 95% CI")+
  scale_x_continuous(breaks = 0:11) +
  facet_wrap(~Inf_strain)+
  geom_jitter(width=0.1, size=5, pch = 21, color = "black", aes(fill = strain), alpha = 0.78) +
  labs(x = "Oocyst per gram", y = "Day post infection") +
  scale_fill_discrete(labels = c("Eastern mice", "Hybrids", "Western mice")) +
  scale_y_continuous(labels = scientific) +
  theme_alice

# Violin plots of the total sum of oocysts collected during 11 days: 
sum.oocysts <- do.call("rbind", by(Oo_Df, Oo_Df$EH_id, function (x){
  x$sum.oo <- sum(x$oocysts.per.g, na.rm=TRUE)
  x
}))

sum.oocysts <- sum.oocysts[!duplicated(sum.oocysts$EH_id),]
# NB: some mice died before!!

ggplot(sum.oocysts, aes(strain, sum.oo)) +
  ggtitle("Sum of oocysts shed during the experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "PWD : HMHZ-eastern-like mice; WSB : HMHZ-west-like mice; WP : HMHZ-hybrids-like mice \n EI64 : wild 'eastern' Eimeria, Eflab : lab 'western-like' Eimeria")+
  geom_violin(color = "black")+
  facet_wrap(~Inf_strain) +
  geom_jitter(width=0.1, size=7, alpha = 0.8, pch = 21, aes(fill = strain)) +
  scale_y_continuous(labels = scientific) +
  theme_alice

######### III. Weight of the mice:
mydata_weight_long <- colnames(mydata_weight) <- gsub("mouse_weight.g.", "weight", colnames(mydata_weight))
mydata_weight_long <- melt(mydata_weight, id=c("mouse_no.", "Infection_date"))  # convert to long format

colnames(mydata_weight_long)[colnames(mydata_weight_long)%in%c("mouse_no.", "variable", "value")] <-
  c("LM_id", "dpi", "weight")

mydata_weight_long$dpi <- gsub("weight_", "", mydata_weight_long$dpi)

# Now let's do it relative to the weight of D1:
mydata_weight_relative <- data.frame(LM_id=mydata_weight[,1], mydata_weight[ ,3:14]/mydata_weight[,4]*100)
colnames(mydata_weight_relative) <- gsub("mouse_weight.g.", "rel_w", colnames(mydata_weight_relative))
mydata_weight_relative_long <- melt(mydata_weight_relative, id=c("LM_id"))  # convert to long format

colnames(mydata_weight_relative_long)[colnames(mydata_weight_relative_long)%in%c("variable", "value")] <-
  c("dpi", "rel.weight")

mydata_weight_relative_long$dpi <- gsub("weight_", "", mydata_weight_relative_long$dpi)

Alldata <- merge(mydata_weight_long, mydata_weight_relative_long, all=TRUE)

Alldata <- merge(Alldata, ExpePlanDF, by.x="LM_id", by.y="EH_id", all=TRUE)

Alldata$dpi <- reorder(Alldata$dpi, as.numeric(gsub("dpi", "", Alldata$dpi)))

# Plot to follow the weight : if < 80% weight, mouse has to be sacrificed
ggplot(data=Alldata,
       aes(x=dpi, y=rel.weight, group=LM_id, color=LM_id)) +
  geom_line()+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=2)+
  geom_text(data=subset(Alldata, dpi == "dpi9" | dpi == "dpi11" ),
            aes(label=LM_id))+
  theme(legend.position="none")

# maximum weight lost before death
max.loss <- do.call("rbind", by(Alldata, Alldata$LM_id, function (x){
  m.loss <- which(x$rel.weight==min(x$rel.weight, na.rm=TRUE))
  x[m.loss,]
}))

table(max.loss$dpi, max.loss$strain)

table(max.loss$dpi, max.loss$Inf_strain)

table(max.loss$dpi, max.loss$Inf_strain, max.loss$strain)

tapply(max.loss$rel.weight, max.loss$Inf_strain:max.loss$strain, mean)

ggplot(max.loss, aes(strain, rel.weight, color=Inf_strain)) +
  ggtitle("Relative weight lost during 11 days of experiment, according to the mice genotype, for both Eimeria infection strains", 
          subtitle = "PWD : HMHZ-eastern-like mice; WSB : HMHZ-west-like mice; WP : HMHZ-hybrids-like mice \n EI64 : wild 'eastern' Eimeria, Eflab : lab 'western-like' Eimeria")+
  geom_violin(color = "black")+
  facet_wrap(~Inf_strain) +
  geom_jitter(width=0.1, size=7, pch = 21, color = "black", aes(fill = strain), alpha = 0.8) +
  theme_alice

summary(glm(rel.weight~strain + Inf_strain, data=max.loss))

summary(glm(rel.weight~strain, data = max.loss[max.loss$Inf_strain == "EI64",]))

## using offsets and real weight instead of relative weight for modeling

###################
## Plots weight and oocyst counts for E64:

names(sum.oocysts)[1] <- "LM_id"
E64DF <- merge(max.loss, sum.oocysts, by = "LM_id")
E64DF <- E64DF[E64DF$Inf_strain.y == "EI64", ]
E64DF$weight.loss.max <- 100 - E64DF$rel.weight

p1 <- ggplot(E64DF, aes(strain.y, sum.oo)) +
  ggtitle("Sum of oocysts shed during the experiment") +
  geom_violin(color = "black") + 
  geom_jitter(width=0.1, size=7, alpha = 0.8, pch = 21, aes(fill = strain.y)) +
  scale_y_continuous(labels = scientific) +
  labs(x = "", y = "") + 
  theme_alice +
  theme(legend.position="none")

p2 <- ggplot(E64DF, aes(strain.y, weight.loss.max)) +
  ggtitle("Max weight loss reached during the experiment (%)") +
  geom_violin(color = "black") + 
  geom_jitter(width=0.1, size=7, alpha = 0.8, pch = 21, aes(fill = strain.y)) +
  labs(x = "", y = "") + 
  theme_alice

grid.arrange(p1, p2, ncol = 2)





