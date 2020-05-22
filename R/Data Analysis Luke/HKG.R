library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

E7 <- read.csv("C:/Users/Luke Bednar/Eimeria_Lab/data/3_recordingTables/HKG_E7.csv")
HZ16 <- read.csv("C:/Users/Luke Bednar/Eimeria_Lab/data/3_recordingTables/HKG_HZ16-17.csv")
HZ18 <- read.csv("C:/Users/Luke Bednar/Eimeria_Lab/data/3_recordingTables/HKG_HZ18.csv")

E7 <- select(E7, Mouse_ID, Target, RT.Ct, EXP)
HZ16 <- select(HZ16, Mouse_ID, Target, RT.Ct, EXP)
HZ18 <- select(HZ18, Mouse_ID, Target, RT.Ct, EXP)

HZ18$Target <- as.character(HZ18$Target)
HZ18$Target[HZ18$Target == "beta-Actin"] <- "B-actin"

HKG <- rbind(E7, HZ16)
HKG <- rbind(HKG, HZ18)
HKG <- distinct(HKG)



HKGL <- spread(HKG, Target, RT.Ct)
names(HKGL)[names(HKGL) == "B-actin"] <- "B.actin"

ggplot(HKGL, aes(x = B.actin, y = GAPDH, color = EXP)) +
  geom_point() + 
  geom_smooth(method='lm')

ggscatter(HKGL, x = "B.actin", y = "GAPDH", color = "EXP", add = "reg.line") +
  facet_wrap(~EXP)+
  stat_cor(label.x = 3, label.y = 34) +
  stat_regline_equation(label.x = 3, label.y = 32)
