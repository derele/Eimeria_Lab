mytab <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/1_informationTables/April2018_wildmice_Eferrisi_INFO.csv")

# check for duplicates
duplicated(mytab$original_mouse_id)

# table to organise by factor
table(mytab$sex)

table(mytab$Mouse_strain)

table(mytab$sex, mytab$Mouse_strain)
# Perfect :)

# Check the ages
birthdays <- as.Date(mytab$born)
birthdays <- format(birthdays, format = "%d-%m-%Y")

birthdays <- sub("18-", replacement = "2018-", birthdays)
birthdays <- sub("17-", replacement = "2017-", birthdays)

# calculate age at infection
infectionDate <- "2018-05-02"

# Then we got the age
age <- difftime(infectionDate, birthdays, units = "weeks")
age <- as.numeric(age)

# Histogrammes
hist(age, breaks = 50)

# How many we have above 11 weeks?
mytab$ageAtInfection <- age

table(mytab$ageAtInfection >= 12)

# to select the 25 good mice
G1Tab <- mytab[mytab$ageAtInfection >= 12,]
table(G1Tab$sex, G1Tab$Mouse_strain)

## Too many STRA, let's select the oldest M and F only
x <- tail(sort(G1Tab$ageAtInfection[G1Tab$Mouse_strain %in% "STRA"]), 6)

keepStra <- G1Tab$original_mouse_id[G1Tab$ageAtInfection %in% x & G1Tab$Mouse_strain %in% "STRA"]

G1Tab <- G1Tab[G1Tab$original_mouse_id %in% keepStra | !G1Tab$Mouse_strain %in% "STRA",]
table(G1Tab$sex, G1Tab$Mouse_strain)

# Some plot
library(ggplot2)
library(ggrepel)

ggplot(mytab, aes(x = Mouse_strain, y = ageAtInfection, col = sex,
                  label = original_mouse_id)) +
  geom_boxplot() +
  geom_label_repel() +
  theme_bw()

# Which STRA do we want?
tail(sort(mytab$ageAtInfection[mytab$Mouse_strain %in% "STRA"]), 6)


