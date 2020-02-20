E7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_complete.csv"))
E7$X <- NULL

P3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_COMPLETE.csv"))
colnames(P3)[2] <- "labels"
P3$X <- NULL
P3$Caecum <- NA
E7$comment <- NULL

ggplot(subset(x = E7, dpi == 8), aes(x = IFNy_CEWE, y = IFNy_FEC, color = challenge)) +
  geom_point() +
  facet_wrap("primary", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

ggplot(P3, aes(y = IFNy_CEWE, x = delta, color = challenge)) +
  geom_point() +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")


ggplot(P3, aes(y = IFNy_CEWE, x = delta, color = N.oocyst)) +
  geom_point() +
  facet_wrap("challenge") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

ggplot(E7, aes(y = IFNy_FEC, x = dpi , color = EH_ID)) +
  geom_point() +
  geom_line() +
  facet_wrap("challenge") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")

ggplot(P3, aes(wloss, IFNy_CEWE)) +
  geom_point()

P3_88 <- P3[P3$challenge == "E88",]
P3_64 <- P3[P3$challenge == "E64",]
P3_UNI <- P3[P3$challenge == "UNI",]

mod88 <- lm(wloss ~ IFNy_CEWE, data = P3_88)
mod64 <- lm(wloss ~ IFNy_CEWE, data = P3_64)
modUNI <- lm(wloss ~ IFNy_CEWE, data = P3_UNI)

ggplot(P3_88, aes(OPG, IFNy_CEWE)) +
  geom_abline(data = mod88) + 
  geom_point() 

ggplot(P3_64, aes(OPG, IFNy_CEWE)) +
  geom_abline(data = mod64) +
  geom_point()

ggplot(P3_UNI, aes(OPG, IFNy_CEWE)) +
  geom_abline(data = modUNI) +
  geom_point()

ggplot(P3, aes(y = IFNy_CEWE, x = OPG)) +
  geom_text_repel(aes(label=delta)) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")
