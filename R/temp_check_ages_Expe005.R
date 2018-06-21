mydata <- read.csv("../data/1_informationTables/rawData/Exp005_WDS_21.6.18.csv")

# mydata <- mydata[-grep("REZERVA", mydata$NOTE),]

mydata$miceGenome <- as.factor(paste0(mydata$Strain, "_", mydata$Sir))

# first date : "04/20/2018"
mydata$Born <- as.Date(mydata$Born, "%m/%d/%Y")

# dates of each infection experiment: 
# dpi0expe1a <- as.Date("2018-07-09")
# dpi0expe1b <- as.Date("2018-07-10")
# dpi0expe2a <- as.Date("2018-07-23")
dpi0expe2b <- as.Date("2018-07-24")

# mydata$ageAtdpi0expe1a <- as.numeric(dpi0expe1a - mydata$Born) / 7
# mydata$ageAtdpi0expe1b <- (dpi0expe1b - mydata$Born) / 7
# mydata$ageAtdpi0expe2a <- (dpi0expe2a - mydata$Born) / 7
mydata$ageAtdpi0expe2b <- as.numeric(dpi0expe2b - mydata$Born) / 7

library(dplyr)
aggData <- mydata %>% 
  dplyr::group_by(miceGenome, ageAtdpi0expe2b) %>%
  dplyr::summarise(n = n()) %>%
  data.frame()

library(ggplot2)
ggplot(aggData,
       aes(x = ageAtdpi0expe2b, y = miceGenome)) + 
  geom_point(aes(fill = miceGenome),
             pch = 21, alpha = .5, size = 10) +
  geom_text(aes(label = n), size=5) + 
    scale_x_continuous(breaks = 0 : 100) +
  geom_vline(xintercept = 12, linetype = 2) +
  geom_vline(xintercept = 14, linetype = 3) +
  theme_bw()
