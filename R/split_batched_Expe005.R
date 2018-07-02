mydata <- read.csv("../data/1_informationTables/Exp005_WDS_21.6.18.csv")

# first date : "04/20/2018"
mydata$Born <- as.Date(mydata$Born, "%m/%d/%Y")

# dates of each infection experiment: 
dpi0expe1a <- as.Date("2018-07-09")

mydata$ageAtdpi0expe1a <- as.numeric(dpi0expe1a - mydata$Born) / 7

# order by age at infection first expe
mydata = mydata[with(mydata, order(-ageAtdpi0expe1a)), ]

# Split by mice strain
mydataList = split(mydata, mydata$Strain)

## We want 48 mice for the first group (2 infection strains) and the rest after
## Parental strains : 6 and rest
## Hybrids : 3 and rest
# 4 * 6 + 8 * 3

library(dplyr)

sampleGroup <- function(df, howmany){
  ss <- sample(1:nrow(df), size = howmany, replace = FALSE)
  df[ss, "InfectionStrain"] <- 'E64'
  df[-ss, "InfectionStrain"] <- 'E88'
  return(df)
}

selectedMice <- lapply(
  mydataList, 
  function(data){
    howManyRows <- nrow(data)
    takeHead <- 3
    if (data$HybridStatus == "parental strains") {
      takeHead <- 6
    }
    takeTail <- howManyRows - takeHead
    howManyToSampleHead <- round(runif(1)) * takeHead %% 2 + as.integer(takeHead / 2)
    howManyToSampleTail <- round(runif(1)) * takeTail %%2 + as.integer(takeTail / 2)
    data[1:takeHead, "Batch"] <- 1
    data[seq(takeHead+1, howManyRows), "Batch"] <- 2
    dataHead <- sampleGroup(head(data, takeHead), howManyToSampleHead)
    dataTail <- sampleGroup(tail(data, takeTail), howManyToSampleTail)
    return(rbind(dataHead, dataTail))
  }
) %>% bind_rows()

# plop <- sample(1:nrow(df), size = 3, replace = FALSE)


library(dplyr)
aggData <- selectedMice %>% 
  dplyr::group_by(Strain, ageAtdpi0expe1a, Batch) %>%
  dplyr::summarise(n = n()) %>%
  data.frame()

library(ggplot2)
ggplot(selectedMice,
       aes(x = ageAtdpi0expe1a, y = Strain)) + 
  geom_jitter(aes(fill = as.factor(Batch)),
             size = 3, pch =21) +
  facet_grid(.~InfectionStrain) +
  scale_x_continuous(breaks = 0 : 100) +
  theme_bw()

write.csv(selectedMice, "../data/1_informationTables/Exp005_WDS_21.6.18.csv", row.names = F)
