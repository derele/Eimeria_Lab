## Experiment 007 Reinfection experiment Nov-Dec 2018
## Eimeria isolates : E64 (E. ferrisi), E88 (E. falciformis)

# Import our data
Expe007a_oocysts <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp_007a_oocyst_counts.csv")
Expe007a_weight <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/Exp_007a_feces.csv")
# add b tables here... TBC
Expe007a_design <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/Exp007a_design.csv")
Expe007b_design <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/2_designTables/Exp007b_design.csv")

# First, merge the tables
Expe007a <- merge(Expe007a_oocysts, Expe007a_weight)
Expe007a <- merge(Expe007a, Expe007a_design)

# Then prepare the oocysts counts

Expe007a$totalOocysts <- ((Expe007a$oocyst_1 + 
                             Expe007a$oocyst_2 + 
                             Expe007a$oocyst_3 +
                             Expe007a$oocyst_4) / 4) * 
  10000 * # because volume chamber
  Expe007a$volume_PBS_mL

Expe007a$OPG <- Expe007a$totalOocysts / Expe007a$fecweight 

# keep in mind that some of the mice were at MORE than 2 infections !!!

# Make plot cause we can
# install.packages("ggplot2")
library(ggplot2) # call the package

ggplot(Expe007a, aes(x = dpi, y = OPG, group = EH_ID, col = EH_ID)) +
  geom_point() +
  geom_line() +
  theme_bw() # pretty


