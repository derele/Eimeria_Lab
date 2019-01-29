# NMRI passaging 4 strains January 2019
# Alice Balard

library(directlabels)
source("functions4InfectionExperiments.R")

# Import data 
Expe008_record <- read.csv("../data/3_recordingTables/Exp008_Jan2019_NMRI_4isolatespassaging.csv") 
Expe008_design <- read.csv("../data/2_designTables/Exp008_NMRI_DESIGN_jan2019.csv")

Expe008 <- merge(Expe008_design, Expe008_record)

# Calculate weight loss
Expe008 <- calculateWeightLoss(Expe008, startingDay = 0)

Expe008$EimeriaSpecies[Expe008$Isolate %in% c("E64", "E139")] <- "E. ferrisi"
Expe008$EimeriaSpecies[Expe008$Isolate %in% c("E88", "Eflab")] <- "E. falciformis"

# Plot follow weight
ggplot(Expe008, aes(x = dpi, y = relativeWeight)) +
  geom_line(aes(group = EH_ID, col = Isolate)) +
  geom_point(size=4, pch = 21, aes(fill = Isolate))+
  mytheme +
  scale_x_continuous(breaks = 0:11, name = "Day post infection (dpi)") +
  scale_y_continuous(name = "Weight (g)") +
  theme(strip.text.y = element_text(size = 15)) +
  ggtitle("Weight along experiment per individual") +
  facet_grid(.~EimeriaSpecies)
