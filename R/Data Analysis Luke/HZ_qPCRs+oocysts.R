#### process HZ infection intensities, gene expressions and add oocysts









# load in wild data for comparison
HZ18 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ18_complete.csv"
HZ18 <- read.csv(text = getURL(HZ18))
names(HZ18)[names(HZ18) == "inf"] <- "Caecum"
HZ18$Caecum <- as.character(HZ18$Caecum)
HZ18$Caecum[HZ18$Caecum == "TRUE"] <- "pos"
HZ18$Caecum[HZ18$Caecum == "FALSE"] <- "neg"
names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"
# graph to compare NE between infected and non-infected
ggplot(HZ18, aes(x = HI, y = NE, color = Caecum)) +
  geom_point() +
  facet_wrap("Target") +
  geom_smooth(method = "lm")

ggplot(HZ18, aes(x = HI, y = NE, color = Eimeria.subspecies)) +
  geom_point() +
  facet_wrap("Target") +
  geom_smooth(method = "lm")

ggplot(HZ18, aes(x = delta, y = NE, color = Eimeria.subspecies)) +
  geom_point() +
  facet_wrap("Target") +
  geom_smooth(method = "lm")

###################################################################
# HZ18 <- HZ18[!(HZ18$Caecum == "neg"),]

#process to graph together
HZ18 <- HZ18[!(HZ18$Target=="GBP2"),]
HZ18 <- HZ18[!(HZ18$Target=="IL-6"),]
i <- sapply(HZ18, is.factor)
HZ18[i] <- lapply(HZ18[i], as.character)
HZ18$Target[HZ18$Target == "IL-12b"] <- "IL-12"
HZ18$inf <- NULL

names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"

ggplot(HZ18, aes(x = HI, y = NE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap("Target")

ggplot(HZ18, aes(x = HI, y = delta)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap("Target")

#remove negs
E7 <- E7[!(E7$Caecum == "neg"),]
names(HZ18)[names(HZ18) == "deltaCtMmE_tissue"] <- "delta"
E7 %>% mutate_if(is.factor, as.character) -> E7
E7$Target[E7$Target == "IL.12"] <- "IL-12"


# combine in one table and distinguish as batches (rename Mouse_ID to EH_ID for sake of merging)
HZ18$batch <- "HZ18"
E7$batch <- "E7"
complete$batch <- "E7"
names(HZ18)[names(HZ18) == "Mouse_ID"] <- "EH_ID"
All <- bind_rows(HZ18, E7)
#remove dodgy outliers
# All <- All[-c(127, 129), ]

ggplot(All, aes(x = NE, y = delta, color = batch)) +
  geom_point() +
  facet_wrap("Target") +
  coord_flip() +
  geom_smooth(method = "lm")

# just complete minus outliers
complete <-complete[!(complete$NE > 0),]
All <-All[!(All$NE > 0),]