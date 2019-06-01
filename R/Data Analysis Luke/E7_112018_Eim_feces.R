# setwd("Repositories/Eimeria_Lab/R/Data Analysis Luke/")

library(Rmisc)

#load csv files at HU#
#E7a_design <- read.csv("../luke/Repositories/Eimeria_Lab/data/2_designTables/E7a_112018_Eim_design.csv")
#E7b_design <- read.csv("../luke/Repositories/Eimeria_Lab/data/2_designTables/E7b_112018_Eim_design.csv")

#load csv files at home (Win)
E7a_design <- read.csv("../Documents/Eimeria_Lab/data/2_designTables/E7a_112018_Eim_design.csv")
E7b_design <- read.csv("../Documents/Eimeria_Lab/data/2_designTables/E7b_112018_Eim_design.csv")
E7aF <- read.csv("../Documents/Eimeria_Lab/data/3_recordingTables/E7a_112018_Eim_feces.csv")
E7bF <- read.csv("../Documents/Eimeria_Lab/data/3_recordingTables/E7b_112018_Eim_feces.csv") 

#load csv at IZW
#E7a_design <- read.csv("../luke/Documents/Eimeria_Lab/data/2_designTables/E7a_112018_Eim_design.csv")
# E7b_design <- read.csv("../luke/Documents/Eimeria_Lab/data/2_designTables/E7b_112018_Eim_design.csv")
# E7aF <- read.csv("../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7a_112018_Eim_feces.csv")
# E7bF <- read.csv("../luke/Documents/Eimeria_Lab/data/3_recordingTables/E7b_112018_Eim_feces.csv") 

#the columns we want to keep
col2keep <- c("Strain", "HybridStatus", "EH_ID")

E7a_design <- E7a_design[col2keep]
E7b_design <- E7b_design[col2keep]

# rename EH_id to EH_ID#
names(E7a_design)[names(E7a_design) == "EH_id"] <- "EH_ID"
names(E7b_design)[names(E7b_design) == "EH_id"] <- "EH_ID"

#add infection history IZW
# history_a <- read.csv("../luke/Documents/Eimeria_Lab/data/2_designTables/E7a_112018_Eim_infection.history.csv")
# history_b <- read.csv("../luke/Documents/Eimeria_Lab/data/2_designTables/E7b_112018_Eim_infection.history.csv")
# history <- rbind(history_a, history_b)

#add infection history HU
#history_a <- read.csv("../luke/Repositories/Eimeria_Lab/data/2_designTables/E7a_112018_Eim_infection.history.csv")
#history_b <- read.csv("../luke/Repositories/Eimeria_Lab/data/2_designTables/E7b_112018_Eim_infection.history.csv")
#history <- rbind(history_a, history_b)

#add infection history home Win
history_a <- read.csv("../Documents/Eimeria_Lab/data/2_designTables/E7a_112018_Eim_infection.history.csv")
history_b <- read.csv("../Documents/Eimeria_Lab/data/2_designTables/E7b_112018_Eim_infection.history.csv")
history <- rbind(history_a, history_b)


# let's make one big fat Expe007 design table (E88 = 31 entries, E64 = 38 entries)
E7_design <- rbind(E7a_design, E7b_design)

E7_design <- merge(E7_design, history, by = "EH_ID")

# remove shit columns
#E7aF <- E7aF[-grep(pattern = "X", x = names(E7aF))]
#E7bF <- E7bF[-grep(pattern = "X", x = names(E7bF))]

#load feces data HU
#E7aF <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7a_112018_Eim_feces.csv", row.names = NULL)
#E7bF <- read.csv("../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7b_112018_Eim_feces.csv", row.names = NULL)

#load feces data home Win
E7aF <- read.csv("../Documents/Eimeria_Lab/data/3_recordingTables/E7a_112018_Eim_feces.csv", row.names = NULL)
E7bF <- read.csv("../Documents/Eimeria_Lab/data/3_recordingTables/E7b_112018_Eim_feces.csv", row.names = NULL)

# keep the batch information
E7aF$batch <- "october2018"
E7bF$batch <- "december2018"

# Make one big fat table Expe 7
E7_record <- rbind(E7aF, E7bF)

# Merge all, #
E7 <- merge(E7_design, E7_record)

#export HU
#write.csv(E7, "../luke/Repositories/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv", quote = FALSE)

#export IZW
#write.csv(Exp007, "../luke/Documents//Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv", quote = FALSE", quote = FALSE)

#export home (Win)
write.csv(Exp007, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_complete.csv", quote = FALSEe)

# Split in 2 infection batches (Efalciformis and Eferrisi are studied separetely)
E7_E88 <- E7[E7$primary %in% "E88",]
E7_E64 <- E7[E7$primary %in% "E64",]

#convert Wchange to numeric (291 and 295 are NAs)#
# str(E7_E64)
# Exp007_E64[,9] <- sapply(Exp007_E64[,9], as.numeric)
# Exp007_E64[,7] <- sapply(Exp007_E64[,7], as.numeric)
# str(Exp007_E64)
# str(Exp007_E88)
# Exp007_E88[,9] <- sapply(Exp007_E88[,9], as.numeric)
# Exp007_E88[,7] <- sapply(Exp007_E88[,7], as.numeric)
# str(Exp007_E88)




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


