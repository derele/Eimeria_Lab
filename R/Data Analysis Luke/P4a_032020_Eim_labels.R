# make design tables
P4a <- read.csv("~/Eimeria_Lab/data/2_designTables/P4a_012020_Eim_design.csv")

P4a$labels <- sample(combn(LETTERS, 3, paste, collapse = ""), nrow(P4a))

write.csv(P4a, "~/Eimeria_Lab/data/2_designTables/P4a_012020_Eim_design.csv")
