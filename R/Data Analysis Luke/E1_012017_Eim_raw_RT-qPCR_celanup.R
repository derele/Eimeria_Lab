# raw RT-qPCR cleanup for E1_012017, additional data from CEC
library(RCurl)
library(httr)
library(Rmisc)

#load in data from GitHub
Samples_36_37_38_39url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/20_05_2019_Samples_36_37_38_39/Luke_2019_05_21_Enas2017samples_36_37_38_39%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_36_37_38_39 <- read.csv(text = getURL(Samples_36_37_38_39url), sep = ";")

Samples_23_24_25_26url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/03_07_2019_Samples_23_24_25_26/admin_2019-07-03%2013-01-15_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_23_24_25_26 <- read.csv(text = getURL(Samples_23_24_25_26url), sep = ";")

Samples_C_27_74_95url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/20_05_2019_Samples_C_27_74_95/admin_2019-05-20%2011-56-24_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_C_27_74_95 <- read.csv(text = getURL(Samples_C_27_74_95url), sep = ";")

Samples_30_29_54_92url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/08_07_2019_Samples_30_29_54_92/admin_2019-07-08%2010-09-26_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_30_29_54_92 <- read.csv(text = getURL(Samples_30_29_54_92url), sep = ";")

Samples_40_41_43_46url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/27_05_2019_Samples_40_41_43_46/admin_2019-05-27%2011-30-32_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_40_41_43_46 <- read.csv(text = getURL(Samples_40_41_43_46url), sep = ";")

Samples_42_44_47_50url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/28_05_2019_Samples_42_44_47_50/admin_2019-05-28%2011-37-44_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_42_44_47_50 <- read.csv(text = getURL(Samples_42_44_47_50url), sep = ";")

Samples_49_32_28_34url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/10_07_2019_Samples_49_32_28_34/admin_2019-07-10%2015-30-23_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_49_32_28_34 <- read.csv(text = getURL(Samples_49_32_28_34url), sep = ";")

Samples_51_55_56_58url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/28_05_2019_Samples_51_55_56_58/admin_2019-05-28%2014-06-31_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_51_55_56_58 <- read.csv(text = getURL(Samples_51_55_56_58url), sep = ";")

Samples_59_60_61_62url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/29_05_2019_Samples_59_60_61_62/admin_2019-05-29%2010-26-36_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_59_60_61_62 <- read.csv(text = getURL(Samples_59_60_61_62url), sep = ";")

Samples_57_94_70_76url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/08_07_2019_Samples_57_94_70_76/admin_2019-07-08%2011-57-32_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_57_94_70_76 <- read.csv(text = getURL(Samples_57_94_70_76url), sep = ";")

Samples_63_64_65_66url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/29_05_2019_Samples_63_64_65_66/admin_2019-05-29%2014-29-53_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_63_64_65_66 <- read.csv(text = getURL(Samples_63_64_65_66url), sep = ";")

Samples_67_68_69_73url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/01_07_2019_Samples_67_68_69_73/admin_2019-07-01%2011-00-22_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_67_68_69_73 <- read.csv(text = getURL(Samples_67_68_69_73url), sep = ";")

Samples_45_71_72_78url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/11_07_2019_Samples_45_71_72_78/admin_2019-07-11%2010-39-22_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_45_71_72_78 <- read.csv(text = getURL(Samples_45_71_72_78url), sep = ";")

Samples_75_77_79_80url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/01_07_2019_Samples_75_77_79_80/admin_2019-07-01%2011-00-22_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_75_77_79_80 <- read.csv(text = getURL(Samples_75_77_79_80url), sep = ";")

Samples_81_82_83_84url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/02_07_2019_Samples_81_82_83_84/admin_2019-07-02%2011-00-22_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_81_82_83_84 <- read.csv(text = getURL(Samples_81_82_83_84url), sep = ";")

Samples_85_86_88_89url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/02_07_2019_Samples_85_86_88_89/admin_2019-07-02%2011-00-22_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_85_86_88_89 <- read.csv(text = getURL(Samples_85_86_88_89url), sep = ";")

Samples_90_91_93_22url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/03_07_2019_Samples_90_91_93_22/admin_2019-07-03%2010-38-19_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_90_91_93_22 <- read.csv(text = getURL(Samples_90_91_93_22url), sep = ";")

Samples_87_73_53_31url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/10_07_2019_%20Samples_87_73_53_31/admin_2019-07-10%2013-29-23_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_87_73_53_31 <- read.csv(text = getURL(Samples_87_73_53_31url), sep = ";")

Samples_48url<- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E1_012017_Eim_RT-Caecum_qPCR/11_07_2019_Samples_48/admin_2019-07-11%2014-19-57_CC011384%20-%20%20Quantification%20Cq%20Results_0.csv"
Samples_48 <- read.csv(text = getURL(Samples_48url), sep = ";")

#meger it all
All <- rbind(Samples_48, Samples_87_73_53_31, Samples_90_91_93_22, Samples_85_86_88_89, Samples_81_82_83_84, Samples_75_77_79_80, 
             Samples_45_71_72_78, Samples_67_68_69_73, Samples_63_64_65_66, Samples_57_94_70_76, Samples_59_60_61_62, Samples_51_55_56_58,
             Samples_49_32_28_34, Samples_42_44_47_50, Samples_40_41_43_46, Samples_30_29_54_92, Samples_C_27_74_95, Samples_23_24_25_26,
             Samples_36_37_38_39)
#Delete rubbish columns
All$X <- NULL
All$Biological.Set.Name <- NULL
All$Cq.Std..Dev <- NULL
All$Starting.Quantity..SQ. <- NULL
All$Log.Starting.Quantity <- NULL
All$Well.Note <- NULL
All$SQ.Mean <- NULL
All$Cq <- NULL
All$SQ.Std..Dev <- NULL
#convert "n.def" and 0 to NAs, Unknowns to Cecum, Empty Target to Nothing and Ppia to Pos Ctrl
All[All=="n. def."] <- NA
All[All=="0"] <- NA
All$Content <- as.character(All$Content)
All$Content[All$Content == "Unkn"] <- "Cecum"
All$Target <- as.character(All$Target)
All$Target[All$Target == ""] <- "Nothing"
#sets Ppia accidentally to NA, doen't matter, just convert NAs to "Pos Ctrl" ## Definitely fix later
All <- within(All, Content[Target=="Ppia"] <- ("Pos Ctrl"[Target=="Ppia"]))
All$Content[is.na(x = All$Content)] <- "Pos Ctrl"
#Omit NAs
All <- na.omit(All)
#add 00 to LM names
All$Sample <- sub("LM", "LM00", All$Sample )
#replace commas with dots for decimal values
All$Cq.Mean <- gsub(",", '.', All$Cq.Mean, fixed = T)
#convert number values from factor to numeric
All$Cq.Mean <- as.numeric(All$Cq.Mean)
#write, omit row names for compatible raw read
write.csv(All, file = "~/Eimeria_Lab/data/3_recordingTables/E1_012017_Eim_RT-qPCR_clean.csv", row.names = FALSE)
