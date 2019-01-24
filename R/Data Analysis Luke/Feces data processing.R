setwd("Repositories/Eimeria_Lab/R/Data Analysis Luke/")
#load csv files#
Inf2a <- read.csv("Inf2a_Exp005.DESIGN.csv")
Inf2b <- read.csv("Inf2b_Exp005.DESIGN.csv")
E7aF <- read.csv("Exp_007a feces.csv")
E7bF <- read.csv("Exp_007b feces.csv")
#remove uneccessary columns#
E7aFC <- E7aF[, -c(8:18)]
E7bFC <- E7bF[, -c(8:18)]
Inf2aC <- Inf2a[, -c(1:2)]
Inf2aC <- Inf2aC[, -c(2:25)]
Inf2aC <- Inf2aC [, -c(3:4)]
Inf2bC <- Inf2b[, -c(1:2)]
Inf2bC <- Inf2bC[, -c(2:25)]
Inf2bC <- Inf2bC [, -c(3:4)]
#rename EH_id to EH_ID#
names(Inf2aC)[names(Inf2aC) == "EH_id"] <- "EH_ID"
names(Inf2bC)[names(Inf2bC) == "EH_id"] <- "EH_ID"
#merge feces and infection dfs#
rbind.data.frame(Inf2aC)
