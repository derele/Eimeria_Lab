# Eimeria_Lab
A. Balard

Record of infection experiments performed in AG Heitlinger

## Protocol

### 1.Information
Fill the information table with known informations about the mice

### 2. Experiment design
Fill the experiment design with the help of the information table,
with the help of R function "makeDesignTable.R"

ex:

```r

myDesignTable2 <- makeDesignTable(myseed = 1234 ,
                                  pathToInfoTable = "../data/1_informationTables/Exp004_May2018_wildmice_Eferrisi_secondbatch_INFO.csv",
                                  firstEH_Id = "LM0145")

# Separate equally between Mouse_strains
library(experiment)

expe <- randomize(data = myDesignTable2, group = c("E64", "E139"),
                   indx = myDesignTable2$EH_id, block = myDesignTable2$Strain)

trt <- data.frame(infection_isolate = expe$treatment)
trt$EH_id <- rownames(trt)
rownames(trt) <- NULL

# Create design table
designTable <- merge(trt, myDesignTable2)
print(table(designTable$Sex, designTable$Strain, designTable$infection_isolate))

write.csv(designTable,     "../data/2_designTables/Exp004_June2018_wildmice_Eferrisi_Secondbatch_DESIGN.csv", row.names = F)
```

## Historic

### Exp001. Parental and F1 wild mice Ploen (x2) cross infection (Eflab, E64) (Francisca)

### Pass001 Nov 2017, passaging 4 isolates (Eflab, E88, E139, E64) in NMRI. 2 mice per cage

### Exp002. NMRI 4 strains passaging, March 2018

### Exp003 & Exp004 Parental wild mice (x4) cross infection (E64, E139) (Vivian) 

#### Exp003 First batch: 02/05/2018

#### Exp004 Second batch:tba

## Diverse R codes

**makeDesignTable.R** is a general function taking as input an *INFO.csv* data.frame,
and creating a DESIGN.csv one

**selectBasedOnAge.R** was used in April 2018 to split our groups in 2. Hard coded.
