## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Lab
## Script adapted from Victor's 3_Data_analysis_qPCR_flotation


### Code to analyses
## 1) Correlation among qPCR and oocyst flotation quantification
### library(ggsci)

library(dplyr)

if(!exists("sample.data")){
  source("R/raw2input/Data_Preparation.R")
}

if(!exists("sdt")){
    source("R/raw2input/qPCR_Data_Preparation.R")
}

##select the completed 22 mice 
sdt <- sdt[sdt$EH_ID %in% c("LM0180", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229","LM0238", "LM0191",
                     "LM0240","LM0246","LM0247","LM0194","LM0190","LM0199","LM0197","LM0198","LM0188","LM0254","LM0192"),]
sample.data <- sample.data[sample.data$EH_ID %in% c("LM0180", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229","LM0238", "LM0191",
                            "LM0240","LM0246","LM0247","LM0194","LM0190","LM0199","LM0197","LM0198","LM0188","LM0254","LM0192"),]


###Let's start plotting and analysing the data!
### 1) Course of infection 
##Genome copies/g of faeces
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  filter(!is.na(Genome_copies_gFaeces))%>%
  wilcox_test(Genome_copies_gFaeces ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "data/Experiment_results/Quant_Eimeria/Tables/Genome_copies_gFaeces_DPI_Comparison.csv")
stats.test%>%
  dplyr::mutate(y.position = log10(y.position))%>%
  dplyr::mutate(dpi = c("1","2","3","4", "5","6", "7", "8", "9", "10","11"))-> stats.test

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces, Infection)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  filter(!is.na(Genome_copies_gFaeces))%>% 
  #filter(Genome_copies_gFaeces!=0)%>% 
  ggplot(aes(x= dpi, y= Genome_copies_gFaeces+1))+
  scale_y_log10("log10 (Genome copies/g Faeces + 1) (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot()+
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Infection, fill= dpi), color= "black")+
  scale_shape_manual(values = c(21, 24))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  labs(tag= "a", shape= "qPCR status (Melting curve)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  guides(fill=FALSE)+
  stat_pvalue_manual(stats.test, label= "p.adj.signif", x= "dpi", y.position = 100000000000)-> A

##Significant mean difference from DPI 2!!!!! not all samples but some sowed signal!!!!

## Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
# sdt%>%
#   filter(dpi%in%c("0","1", "2", "3","4", "5","6", "7", "8", "9", "10"))%>%
#   dplyr::select(EH_ID, dpi,OPG)%>%
#   dplyr::arrange(EH_ID)%>%
#   dplyr::arrange(dpi)%>% ##for comparison 
#   dplyr::mutate(OPG= OPG+1)%>% ##To check
#   filter(!is.na(OPG))->#%>% 
  #wilcox_test(OPG ~ dpi, alternative = "two.sided", ref.group = "0")%>%
#  adjust_pvalue(method = "bonferroni") %>%
# add_significance()%>%
#  add_xy_position(x = "dpi")#-> stats.test

##Save statistical analysis
#x <- stats.test
#x$groups<- NULL
#write.csv(x, "Tables/OPG_DPI_Comparison.csv")

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= OPG+1))+
  scale_y_log10("log10 (Oocysts/g Faeces + 1) (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  scale_color_brewer(palette = "Paired")+
  labs(tag= "b")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "none")+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> B

##Significant mean difference from day 4 and on... Basically DPI 0 to 3 No OPG and DPI 4 equal to 10!

####calculate weightloss using function
weightloss <- function(sdt){
  sdt$weightloss <-
    (sdt$weight_dpi0 - sdt$weight)
  return(sdt)
}
sdt <- weightloss(sdt)

##Weight loss 
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi, weightloss)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  wilcox_test(weightloss ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "data/Experiment_results/Quant_Eimeria/Tables/Weightloss_DPI_Comparison.csv")

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi,weightloss)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= weightloss))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  xlab("Day post infection")+
  ylab("Weight loss (%)")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  scale_color_brewer(palette = "Paired")+
  labs(tag= "c", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> C

##Figure 2: Course of Eimeria Infection in genome copies, OPG, and weight loss
pdf(file = "data/Experiment_results/Quant_Eimeria/Figures/Figure_2.pdf", width = 10, height = 15)
grid.arrange(A,B,C)
dev.off()
rm(A,B,C, x, stats.test)


#####STOPPED SCRIPT APPLICATION AT THIS POINT################ -Dilé

-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------





### 2) Correlation among Eimeria quantification methods
###Aiming to state a quantitative measure of the relationship between both measurements 
###Generate a data set without samples with zero OPG counts and/or zero genome copies 
sdt%>%
  filter(!(OPG== 0))%>%
  filter(!(Genome_copies_gFaeces== 0))-> sdt.nozero

##Model 1: Genome copies/g faeces modeled by OPG
DNAbyOPG <- lm(log10(Genome_copies_gFaeces)~log10(OPG),
               data = sdt.nozero, na.action = na.exclude)
summary(DNAbyOPG)

### Plot and extract estimates
require(sjPlot)

sdt.nozero$predictedM1 <- predict(DNAbyOPG)   # Save the predicted values
sdt.nozero$residualsM1 <- residuals(DNAbyOPG) # Save the residual values

##Plot model
##Assign the colors for dpi and keep consistency with previous plots
colores<- c("4"="#00BD5C", "5"= "#00C1A7", "6"= "#00BADE", "7"= "#00A6FF", 
            "8" = "#B385FF", "9"= "#EF67EB", "10" = "#FF63B6")

####Genome copies modeled by OPGs 
sdt.nozero%>%
  ggplot(aes(OPG, Genome_copies_gFaeces))+
  geom_smooth(method = lm, col= "black")+
  scale_x_log10(name = "log10 (Oocysts/g Faeces) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces)  \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  scale_fill_manual(values = colores, guide= "none")+
  labs(tag= "a", fill= "DPI")+
  theme_bw()+
  #stat_cor(label.y = 5,  label.x = 6, method = "spearman",
  #         aes(label= paste("rho","'='", ..r.., ..p.label.., sep= "~` `~")))+
  theme(text = element_text(size=16), legend.position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  annotation_logticks()-> A

##To visualize it 
#ggsave(filename = "Rplots.pdf", A)

##Plot residuals
##Mean residuals for plot 
sdt.nozero%>%
  group_by(dpi) %>% 
  summarise(residualsM1_mean = mean(na.omit(residualsM1)))%>%
  inner_join(sdt.nozero, by= "dpi")%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(dpi, residualsM1, residualsM1_mean)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  mutate(residualim= 0)%>%
  ggplot(aes(x= dpi, y= residualsM1))+
  geom_jitter(width = 0.5, shape=21, size=2.5, aes(fill= dpi), alpha= 0.75, color= "black")+
  scale_fill_manual(values = colores, guide= "none")+
  geom_segment(aes(y= residualsM1_mean, yend= residualim, xend= dpi), color= "black", size= 1) +
  geom_point(aes(x = dpi, y = residualsM1_mean), size=4)+
  #geom_rect(aes(xmin=-0.1,xmax=4.5,ymin=-Inf,ymax=Inf),alpha=0.01,fill="grey")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  xlab("Day post infection")+
  scale_y_continuous(name = "Residuals\n (Genome copies/g Faeces)")+
  labs(tag= "b")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")-> B

##To visualize it externally 
#ggsave(filename = "Rplots.pdf", B)

##Model 2: Genome copies/g faeces modeled by OPG without DPI interaction
DNAbyOPG_dpi <- lm(log10(Genome_copies_gFaeces)~log10(OPG)+dpi,
                   data = sdt.nozero, na.action = na.exclude)
summary(DNAbyOPG_dpi)

##Model 3: Genome copies/g faeces modeled by OPG with DPI interaction
DNAbyOPGxdpi <- lm(log10(Genome_copies_gFaeces)~log10(OPG)*dpi,
                   data = sdt.nozero, na.action = na.exclude)
summary(DNAbyOPGxdpi)

### Plot and extract estimates
sjPlot:: tab_model(DNAbyOPGxdpi, 
                   terms = c("log10(OPG):dpi5", "log10(OPG):dpi6", "log10(OPG):dpi7", 
                             "log10(OPG):dpi8", "log10(OPG):dpi9", "log10(OPG):dpi10"))

plot_model(DNAbyOPGxdpi, terms = c("log10(OPG):dpi5", "log10(OPG):dpi6", "log10(OPG):dpi7", 
                                   "log10(OPG):dpi8", "log10(OPG):dpi9", "log10(OPG):dpi10"))-> tmp.fig.1

#ggsave(filename = "Rplots.pdf", tmp.fig.1)

##Comparison of models
# test difference LRT or anova
lrtest(DNAbyOPG, DNAbyOPG_dpi) #--> Report this table in the results 
lrtest(DNAbyOPG, DNAbyOPGxdpi) #--> Report this table in the results 
lrtest(DNAbyOPG_dpi, DNAbyOPGxdpi) #--> Report this table in the results 

###Model 4: GLMM with dpi as random factor
require(lme4)

DNAbyOPG_dpi_glmm <- lmer(log10(Genome_copies_gFaeces)~log10(OPG) + (1|dpi),
                          data = sdt.nozero, na.action = na.exclude, REML=TRUE)

summary(DNAbyOPG_dpi_glmm)
sjPlot:: tab_model(DNAbyOPG_dpi_glmm)

### Plot estimates
## Random effect estimates
plot_model(DNAbyOPG_dpi_glmm, type = "re", show.values = TRUE)-> tmp.fig.2
#ggsave(filename = "Rplots.pdf", tmp.fig.2)

##Plot model by DPI
sdt.nozero%>%
  mutate(dpi = fct_relevel(dpi, "0","1", "2", "3", "4", "5", 
                                   "6", "7", "8", "9", "10"))%>%
  ggplot(aes(OPG, Genome_copies_gFaeces, fill=dpi))+
  geom_point(shape=21, size=5) +
  geom_smooth(method = lm, se=FALSE, aes(OPG, Genome_copies_gFaeces, color=dpi))+
  scale_color_manual(values = colores, guide= "none")+
  scale_fill_manual(values = colores, guide= "none")+
  scale_x_log10(name = "log10 (Oocysts/g Faeces) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces) \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  labs(tag= "c")+
  theme(text = element_text(size=16), legend.position = "none")+
  annotation_logticks()-> C

##To visualize it externally 
#ggsave(filename = "Rplots.pdf", C)

##Extraction of coefficients by dpi
sdt.nozero%>% 
  nest(-dpi)%>% 
  mutate(cor=map(data,~lm(log10(Genome_copies_gFaeces) ~ log10(OPG), data = .))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()-> x

x$data<- NULL
x$cor<- NULL
x%>%
  arrange(dpi)%>%
  filter(term != "(Intercept)")-> fitted_models_dpi

#write.csv(fitted_models_dpi, "Tables/Q1_OPG_DNA_estimates_DPI.csv",  row.names = F)

sdt.nozero%>% 
  nest(-dpi)%>% 
  mutate(cor=map(data,~cor.test(log10(.x$Genome_copies_gFaeces), log10(.x$OPG), method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()-> x

x$data<- NULL
x$cor<- NULL
x%>%
  arrange(dpi)-> corOPGbyDNA_DPI

#write.csv(corOPGbyDNA_DPI, "Tables/Q1_OPG_DNA_Correlation_DPI.csv",  row.names = F)
##Non significant correlation between measurements by DPI

rm(x, tmp.fig.1, tmp.fig.2, DNAbyOPG, DNAbyOPG_dpi, DNAbyOPG_dpi_glmm, DNAbyOPGxdpi)