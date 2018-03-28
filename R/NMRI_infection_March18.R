source("../R/summarySE.R")
library(ggplot2)

# import file
mydata <- read.csv(file = "../data_raw/NMRIMarch2018/NMRI_Mäuse.csv")

# Calculate the weight loss compared to dpi0 automatically
for(i in 0:11){
  mydata[[paste0("dpi.",i,".weightloss")]] <- 
    (1 - (mydata[[paste0("dpi.",i,".weight")]] / mydata$dpi.0.weight)) * 100
}

# Change to long format
Newdata <- data.frame(NULL)

for (DPI in 0:11){
  temp <- cbind(Mouse_ID = mydata$mouse.ID, 
                mydata[, grep(paste0("dpi.", paste(DPI, "", sep="\\.")), colnames(mydata))])
  names(temp) <- gsub(paste0("dpi.", paste(DPI, "", sep="\\.")), 
                      replacement = "", 
                      x = names(temp)) 
  temp$dpi <- DPI 
  Newdata <- rbind(Newdata, temp)
}

# Remove useless objects:
rm(temp, DPI, i)

# Add groups and sex
Newdata$infection_isolate[Newdata$Mouse_ID %in% 1:6] <- "EfLab"
Newdata$infection_isolate[Newdata$Mouse_ID %in% 7:12] <- "E64"
Newdata$infection_isolate[Newdata$Mouse_ID %in% 13:18] <- "E88"
Newdata$infection_isolate[Newdata$Mouse_ID %in% 19:23] <- "E139"

# And sex

# Plot all mice
ggplot(Newdata, aes(x = dpi, y = weightloss, color = infection_isolate, group = Mouse_ID )) +
  geom_point() +
  geom_line() +
  theme_bw()

## REMOVING NON INFECTED (to do when the oocysts are counted)
Newdata <- Newdata[Newdata$Mouse_ID != 14,]

# Summary table & plot
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
summary_data <- summarySE(Newdata, measurevar="weightloss", groupvars=c("infection_isolate","dpi"))
summary_data

ggplot(summary_data, aes(x = dpi, y = weightloss, colour = infection_isolate, group = infection_isolate)) + 
  geom_errorbar(aes(ymin = weightloss-ci, ymax = weightloss+ci), width=1, position = position_dodge(0.2), size = 2) +
  geom_line(position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2), size=3) +
  theme_bw() +
  theme(legend.position=c(.9, .2), legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +
  scale_color_manual(values=c("#336600", "#66CC00", "#FF3300", "#660000"), 
                     name="Infection\nstrains")+
  scale_x_continuous(breaks = 0:11, name = "Day post infection" )+ 
  scale_y_continuous(name = "Weight loss relative to infection day")
