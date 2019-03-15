# code Vivian & Franci & Alice's experiments
# March 2019
## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85
source("myFunctions.R")
#### Load data ####
source("loadExpe001toExpe005.R")

### Part 1: Franci, preliminary experiment
plotsExpe1 <- makeIntermPlots(ExpeDF_001)
plotsExpe1[[1]]
plotsExpe1[[2]]

## Passaging in NMRI
plotsExpe2 <- makeIntermPlots(ExpeDF_002)
plotsExpe2[[1]]
plotsExpe2[[2]]

## Parental inbred
plotsExpe3_4 <- makeIntermPlots(ExpeDF_003_4)
plotsExpe3_4[[1]]
plotsExpe3_4[[2]]

## Full design (F0, F1 outbred, F1 hybrids)
plotsExpe5 <- makeIntermPlots(ExpeDF_005[ExpeDF_005$Batch == 1,])
plotsExpe5[[1]]
plotsExpe5[[2]]

## ALL TOGETHER
ALL_Expe <- merge(ExpeDF_001, ExpeDF_002, all = T)
ALL_Expe <- merge(ALL_Expe, ExpeDF_003_4, all = T)
ALL_Expe <- merge(ALL_Expe, ExpeDF_005[ExpeDF_005$Batch == 1,], all = T)

plotsALL <- makeIntermPlots(ALL_Expe)
plotsALL[[1]] + coord_cartesian(ylim = c(-10,20))
makeIntermPlots



tolerance_001 <- prepFinalTolComp(tolerance_001)

## Test of hybrid effect
tolerance_001$HybridStatus <- "inbred"
tolerance_001[tolerance_001$Mouse_strain == "WP", "HybridStatus"] <- "hybrids"

### 1.1 On weight loss
myboxPlot(tolerance = tolerance_001,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "HybridStatus", groupName = "Hybrid status") +
  scale_fill_manual(values = c("darkorchid", "yellow", "pink"))

myMod(mydata = tolerance_001,
      response = tolerance_001$maxweightloss, 
      group = tolerance_001$HybridStatus)

### 1.2 On max of oocysts
myboxPlot(tolerance = tolerance_001,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "HybridStatus", groupName = "Hybrid status") +
  scale_fill_manual(values = c("darkorchid", "yellow", "pink"))

myMod(mydata = tolerance_001,
      response = tolerance_001$maxOPG_inmillion, 
      group = tolerance_001$HybridStatus)

### 1.2.3 On tolerance
myboxPlot(tolerance = tolerance_001,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "HybridStatus", groupName = "Hybrid status") 

myMod(mydata = tolerance_001,
      response = tolerance_001$toleranceFactor, 
      group = tolerance_001$HybridStatus)

## Part 2. Strain effect

tolerance_001_E64 <- tolerance_001[tolerance_001$infection_isolate %in% "E.ferrisi (E64)",] 

### 1.1 On weight loss
pdf("fig1.pdf", width = 4, height = 5)
myboxPlot(tolerance = tolerance_001_E64,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "Mouse_genotype", groupName = "Mouse genotype")
dev.off()
# Mmm-Mmd_Hybrid (WP)-MMd_F0 (Ws-Ws) 3.753162e-05 3.753162e-05
# MMm_F0 (Pw-Pw)-MMd_F0 (Ws-Ws) 2.900146e-04 2.900146e-04

myMod(mydata = tolerance_001,
      response = tolerance_001$maxweightloss, 
      group = tolerance_001$Mouse_genotype)

### 1.2 On sum of oocysts shed along expe
pdf("fig2.pdf", width = 4, height = 5)
myboxPlot(tolerance = tolerance_001_E64,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "Mouse_genotype", groupName = "Mouse genotype")
dev.off()
# MMm_F0 (Pw-Pw)-Mmm-Mmd_Hybrid (WP) 0.04848581 0.04848581

myMod(mydata = tolerance_001,
      response = tolerance_001$maxOPG_inmillion, 
      group = tolerance_001$Mouse_genotype)

### 1.2.3 On tolerance
pdf("fig3.pdf", width = 4, height = 5)
myboxPlot(tolerance = tolerance_001_E64,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "Mouse_genotype", groupName = "Mouse genotype") 
dev.off()
#nothing significant
myMod(mydata = tolerance_001,
      response = tolerance_001$toleranceFactor, 
      group = tolerance_001$Mouse_genotype)

### Part 2: Vivian
### comparison tolerance for Vivian
tolerance_comparison <- merge(tolerance_003_4, tolerance_005, all = T)
# Remove individuals that DIDN'T shed oocysts at all! (no infection)
tolerance_comparison <- tolerance_comparison[tolerance_comparison$sum.opg != 0,]
# Let's remove double infections (Exp_005_2a and Exp_005_2b)
tolerance_comparison <- tolerance_comparison[
  !tolerance_comparison$Exp_ID %in% c("Exp_005_2a", "Exp_005_2b"),]
# clean table
tolerance_comparison <- prepFinalTolComp(tolerance_comparison)

## Check first strain
ggplot(tolerance_comparison, aes(maxOPG_inmillion, maxweightloss)) +
  geom_point(size = 5, shape = 21, alpha = 0.6,
             aes(fill = Mouse_genotype)) +
  facet_grid(.~infection_isolate, scales = "free_x") +
  scale_x_log10() # x axis in log 10 to see the fold differences

#######################
# Statistical modelling
#######################

## Preliminary : can we merge E139 and E64 as E.ferrisi?? --> YES
# tolerance_comparisonFERRISI <- tolerance_comparison[tolerance_comparison$infection_isolate != "E.falciformis (E88)",]
# myMod(mydata = tolerance_comparisonFERRISI,
#       response = tolerance_comparisonFERRISI$maxweightloss, 
#       group = tolerance_comparisonFERRISI$Mouse_genotype)
# myMod(mydata = tolerance_comparisonFERRISI,
#       response = tolerance_comparisonFERRISI$maxOPG, 
#       group = tolerance_comparisonFERRISI$Mouse_genotype)
# myMod(mydata = tolerance_comparisonFERRISI,
#       response = tolerance_comparisonFERRISI$toleranceFactor, 
#       group = tolerance_comparisonFERRISI$Mouse_genotype)
## -> no difference between E139 and E64 at the genotype level for all 3 measurements

tolerance_comparison$infection_isolate <- 
  as.character(tolerance_comparison$infection_isolate)

tolerance_comparison[
  tolerance_comparison$infection_isolate %in% 
    c("E.ferrisi (E64)", "E.ferrisi (E139)"), "infection_isolate"] <- "E.ferrisi (E64|E139)"

table(tolerance_comparison$HybridStatus,
      tolerance_comparison$infection_isolate)

## Part 1. Hybrid effect

### 1.1 On weight loss
myboxPlot(tolerance = tolerance_comparison,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "HybridStatus", groupName = "Hybrid status") +
  scale_fill_manual(values = c("darkorchid", "yellow", "pink"))

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxweightloss, 
           group = tolerance_comparison$HybridStatus)

### 1.2 On max of oocysts
myboxPlot(tolerance = tolerance_comparison,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "HybridStatus", groupName = "Hybrid status")

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxOPG_inmillion, 
           group = tolerance_comparison$HybridStatus)

### 1.2.3 On tolerance
myboxPlot(tolerance = tolerance_comparison,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "HybridStatus", groupName = "Hybrid status") 

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$toleranceFactor, 
           group = tolerance_comparison$HybridStatus)

## Part 2. Strain effect

table(tolerance_comparison$Mouse_genotype,
      tolerance_comparison$infection_isolate)

### 1.1 On weight loss
myboxPlot(tolerance = tolerance_comparison,
          response = "maxweightloss", respName = "Maximum relative weight loss (%)",
          group = "Mouse_genotype", groupName = "Mouse genotype")

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxweightloss, 
      group = tolerance_comparison$Mouse_genotype)

### 1.2 On sum of oocysts shed along expe
myboxPlot(tolerance = tolerance_comparison,
          response = "maxOPG_inmillion", respName = "Max OPG (million)",
          group = "Mouse_genotype", groupName = "Mouse genotype")

myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$maxOPG_inmillion, 
      group = tolerance_comparison$Mouse_genotype)

### 1.2.3 On tolerance
myboxPlot(tolerance = tolerance_comparison,
          response = "toleranceFactor", respName = "tolerance factor",
          group = "Mouse_genotype", groupName = "Mouse genotype") 
myMod(mydata = tolerance_comparison,
      response = tolerance_comparison$toleranceFactor, 
           group = tolerance_comparison$Mouse_genotype)



