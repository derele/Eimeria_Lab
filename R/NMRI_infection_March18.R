source("../R/summarySE.R")
library(ggplot2)

# import file
mydata <- read.csv(file = "../data_raw/NMRIMarch2018/NMRI_Mäuse.csv")

# Calculate the weight loss compared to dpi0 automatically
for(i in 0:11){
  mydata[[paste0("dpi.",i,".weightRelativeToInfection(%)")]] <- 
    (mydata[[paste0("dpi.",i,".weight")]] / mydata$dpi.0.weight) * 100
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
Newdata$sex <- "male"
Newdata$sex[Newdata$Mouse_ID %in% c(1,2,10,11,12,16,17,18,21,22,23)] <- "female"

# Plot all mice
ggplot(Newdata, aes(x = dpi, y = `weightRelativeToInfection(%)`, 
                    color = infection_isolate, fill = infection_isolate, group = Mouse_ID )) +
  geom_line() +
  geom_point(pch = 21, size = 4, col = "white") +
  theme_bw()

# Summary table & plot
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
summary_data <- summarySE(Newdata, measurevar="weightRelativeToInfection(%)", groupvars=c("infection_isolate","dpi"))
summary_data

ggplot(summary_data, aes(x = dpi, y = `weightRelativeToInfection(%)`, 
                         colour = infection_isolate, fill = infection_isolate, group = infection_isolate)) + 
  geom_errorbar(aes(ymin = `weightRelativeToInfection(%)`-ci, 
                    ymax = `weightRelativeToInfection(%)`+ci), 
                width=0, position = position_dodge(0.2), size = 1, alpha = .5) +
  geom_line(position = position_dodge(0.2), size = 2) +
  geom_point(pch = 21, position = position_dodge(0.2), size=3, col = "white") +
  theme_bw() +
  theme(legend.position=c(.9, .2), legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +
  scale_color_manual(values=c("#336600", "#66CC00", "#FF3300", "#660000"), 
                     name="Infection\nstrains")+
  scale_fill_manual(values=c("#336600", "#66CC00", "#FF3300", "#660000"), 
                     name="Infection\nstrains")+
  scale_x_continuous(breaks = 0:11, name = "Day post infection" )+ 
  scale_y_continuous(name = "Weight loss relative to infection day")

## Oocysts shedding
Newdata$Neubauer1 <- NA
Newdata$Neubauer2 <- NA
Newdata$Neubauer3 <- NA
Newdata$Neubauer4 <- NA
Newdata$mean <- NA
Newdata$dilution_ml <- NA
Newdata$oocysts_number <- NA
Newdata$OPG <- NA

write.csv(x = Newdata, file = "../data_raw/NMRIMarch2018/followUp.csv", row.names = F)
