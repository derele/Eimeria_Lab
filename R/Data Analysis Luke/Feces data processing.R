# setwd("Repositories/Eimeria_Lab/R/Data Analysis Luke/")

library(Rmisc)

#load csv files at HU#
#Exp007a_design <- read.csv("../luke/Repositories/Eimeria_Lab/data/2_designTables/Exp007a_design.csv")
#Exp007b_design <- read.csv("../luke/Repositories/Eimeria_Lab/data/2_designTables/Exp007b_design.csv")
#E7aF <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/Exp007a/Exp_007a_feces.csv", row.names = NULL)
#E7bF <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/Exp007b/Exp_007b_feces.csv", row.names = NULL)


#load csv files at home (Win)
Exp007a_design <- read.csv("./Eimeria_Lab/data/2_designTables/Exp007a_design.csv")
Exp007b_design <- read.csv("./Eimeria_Lab/data/2_designTables/Exp007b_design.csv")
E7aF <- read.csv("Eimeria_Lab/data/3_recordingTables/Exp007/Exp007a/Exp_007a_feces.csv")
E7bF <- read.csv("Eimeria_Lab/data/3_recordingTables/Exp007/Exp007b/Exp_007b_feces.csv")

#load csv at home (win)
#Exp007a_design <- read.csv("./Eimeria_Lab/data/2_designTables/Exp007a_design.csv")
#Exp007b_design <- read.csv("./Eimeria_Lab/data/2_designTables/Exp007b_design.csv")
#E7aF <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/Exp007a/Exp_007a_feces.csv")
#E7bF <- read.csv("./Eimeria_Lab/data/3_recordingTables/Exp007/Exp007b/Exp_007b_feces.csv") 

#the columns we want to keep
col2keep <- c("Strain", "HybridStatus", "InfectionStrain", "EH_ID")

Exp007a_design <- Exp007a_design[col2keep]
Exp007b_design <- Exp007b_design[col2keep]

# rename EH_id to EH_ID#
names(Exp007a_design)[names(Exp007a_design) == "EH_id"] <- "EH_ID"
names(Exp007b_design)[names(Exp007b_design) == "EH_id"] <- "EH_ID"

# let's make one big fat Expe007 design table (E88 = 31 entries, E64 = 38 entries)
Exp007_design <- rbind(Exp007a_design, Exp007b_design)

# remove shit columns
#E7aF <- E7aF[-grep(pattern = "X", x = names(E7aF))]
#E7bF <- E7bF[-grep(pattern = "X", x = names(E7bF))]

# keep the batch information
E7aF$batch <- "october2018"
E7bF$batch <- "december2018"

# Make one big fat table Expe 7
Exp007_record <- rbind(E7aF, E7bF)

# Merge all, #
Exp007 <- merge(Exp007_design, Exp007_record)

# Split in 2 infection batches (Efalciformis and Eferrisi are studied separetely)
Exp007_E88 <- Exp007[Exp007$InfectionStrain %in% "E88",]
Exp007_E64 <- Exp007[Exp007$InfectionStrain %in% "E64",]

#convert Wchange to numeric (291 and 295 are NAs)#
str(Exp007_E64)
Exp007_E64[,9] <- sapply(Exp007_E64[,9], as.numeric)
Exp007_E64[,7] <- sapply(Exp007_E64[,7], as.numeric)
str(Exp007_E64)
str(Exp007_E88)
Exp007_E88[,9] <- sapply(Exp007_E88[,9], as.numeric)
Exp007_E88[,7] <- sapply(Exp007_E88[,7], as.numeric)
str(Exp007_E88)

#export HU
#write.csv(Exp007, "../luke/Repositories/Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge", quote = FALSE)

#export home (Win)
write.csv(Exp007, "./Eimeria_Lab/data/3_recordingTables/Exp007/Exp_007_Merge", quote = FALSE)

Exp007

# calculate summary statistics on weight loss NEEDS FIXING #
#WLoss_007_E64 <- summarise(Exp007_E64, measurevar = "Wchange",
#                                   groupvars=c("HybridStatus", "dpi", "batch"), na.rm = T)
#WLoss_007_E64$ci[is.na(WLoss_007_E64$ci)] <- 0
#
#WLoss_007_E88 <- summarySE(Exp007_E88, measurevar = "Wchange",
#                                   groupvars=c("HybridStatus", "dpi", "batch"), na.rm = T)
#WLoss_007_E88$ci[is.na(WLoss_007_E88$ci)] <- 0

#convert Exp007 Wchange to numeric, replace dpi0 0% with 100%#
#Exp007[,9] <- sapply(Exp007[,9], as.numeric)
#Exp007$Wchange[Exp007$Wchange == 0] <- 100

#Summary and confidence intervals (feces only)#
#library(dplyr)
#DF <- Exp007 %>%
 # group_by(dpi, InfectionStrain, HybridStatus) %>%
  #dplyr::summarise(mean.weight = mean(Wchange, na.rm = T),
   #                sd.weight = sd(Wchange, na.rm = TRUE),
    #               n.weight = n()) %>%
  #dplyr::mutate(se.weight = sd.weight / sqrt(n.weight),
     #           lower.ci.weight = mean.weight - qt(1 - (0.05 / 2), 
      #                                             n.weight - 1) * se.weight,
       #         upper.ci.weight = mean.weight + qt(1 - (0.05 / 2), 
        #                                           n.weight - 1) * se.weight)

#Plot the Exp007#
library(ggplot2)

ggplot(DF, aes(x = dpi, y = mean.weight))+
  geom_errorbar(aes(ymin = lower.ci.weight, ymax = upper.ci.weight, 
                    col = HybridStatus),
                alpha = 0.5) +
  geom_line(aes(group = HybridStatus, col = HybridStatus), size = 2) +
  geom_point(aes(fill = HybridStatus), size=4, pch = 21, color = "black") +
  facet_grid(. ~ InfectionStrain) +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "relative weight")


