library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(mgcv)
library(lubridate)
library(patchwork)
library(gratia)
library(corrplot)


extraction_data<-extraction_data|>
  mutate(depth_cm=as.numeric(str_extract(V2, "(?<=_)\\d+\\.?\\d*(?=CM)")))#extract depth from last two characters before 'cm'
filtered<-extraction_data|>
  filter(!is.na(depth_cm))#filter out na
filtered$totalDNA<-as.numeric(filtered$totalDNA)#turn DNA into numeric
# filtered<-filtered|>
#   rename(PostWash=V13)
# filtered<-filtered|>
#   rename(WashMin=V29)
# filtered<-filtered|>
#   rename(ExDate=V8) 
filtered$ExDate<-as.character(filtered$ExDate)
filtered$ExDate<-dmy(filtered$ExDate)#turn date into day/month/year
filtered$ExDate_numeric <- as.numeric(filtered$ExDate - min(filtered$ExDate, na.rm = TRUE))
filtered$WashMin<-as.numeric(filtered$WashMin)#calculate wash time
filtered$PostWash<-as.numeric(as.character(filtered$PostWash))#convert post wash mass

# ggplot(data=filtered, aes(x=depth_cm,y=totalDNA, color=V1))+
#   #already renamed V25 to totalDNA
#   geom_point()+
#   geom_smooth(method = "lm", se = FALSE, color = "blue")+
#   scale_y_continuous(breaks = seq(100, 1000, by = 100), limits = c(0, 1000))+
#   theme_minimal()+
#   labs(x="Core Depth(cm)", y="Total DNA(ng)")

model <- gam(totalDNA ~ depth_cm, data = filtered)#basic linear model
summary(model)

mod_sm<-gam(totalDNA ~ s(depth_cm)+ExDate_numeric+WashMin+s(PostWash), data = filtered, na.action = na.exclude)#s for smooth,not enough distinct values for washmin and date extracted
summary(mod_sm) #GAM for basic before any filtering

filtered$residuals<-residuals(mod_sm)
threshold<-2*sd(filtered$residuals, na.rm=TRUE)#calculate residuals to filter out outliers
filtered_clean<-filtered|>
  filter(abs(residuals)<threshold)
gam_clean <- gam(totalDNA~depth_cm+factor(ExDate)+WashMin+s(PostWash), data = filtercleanwo2)#use factor() for date
summary(gam_clean)#0.639 r squared

filtercleanwo2<-filtered_clean|>
  filter(V1!="SR2422_SBB1_04GC_SECT2")#filter out sect2

p1<-ggplot(data=filtercleanwo2, aes(x=depth_cm,y=totalDNA))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue")+
  scale_y_continuous(breaks = seq(100, 1000, by = 100), limits = c(0, 1000))+
  theme_minimal()+
  labs(x="Core Depth(cm)", y="Total DNA(ng)")
p2<-ggplot(data=filtercleanwo2, aes(x=factor(ExDate),y=totalDNA))+
  geom_violin()+
  scale_y_continuous(breaks = seq(100, 1000, by = 100), limits = c(0, 1000))+
  theme_minimal()+
  labs(x="Date Extracted", y="Total DNA(ng)")
p3<-ggplot(data=filtercleanwo2, aes(x=factor(WashMin), y=totalDNA))+
  geom_violin()+
  scale_y_continuous(breaks = seq(100, 1000, by = 100), limits = c(0, 1000))+
  theme_minimal()+
  labs(x="PBS Wash Length(min)", y="Total DNA(ng)")
p1+p2+p3

filtercleanwo2$gam_pred <- predict(gam_clean)
plot(gam_clean, pages = 3, residuals = TRUE, shade = TRUE, seWithMean = TRUE)#check GAM effect plot
draw(gam_clean)

ggplot(filtercleanwo2, aes(x=depth_cm, y=totalDNA))+
  geom_point()+
  geom_line(aes(y=gam_pred),color="blue")+
  theme_minimal()

#correlation between pre and post wash mass
ggplot(filtercleanwo2, aes(x=V12, y=PostWash, color=factor(WashMin)))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Pre Wash Weight(g)", y="Post Wash Weight(g)")

#correlation between pre/post wash mass g and total DNA
filtercleanwo2$PreWash<-as.numeric(filtercleanwo2$V12)
ggplot(filtercleanwo2, aes(x=PreWash, y=totalDNA, color=PostWash))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Pre Wash Weight(g)", y="Total DNA extracted(ng)")


#add correlation plots, pearson's, 
corrmat<-cor(filtercleanwo2[, c("totalDNA", "depth_cm","WashMin")], use = "complete.obs")
corrplot(corrmat, method = "circle", type = "upper", tl.col = "black", addCoef.col = "black")
#correlation between pre/post wash weight and total DNA
corrWeight<-cor(filtercleanwo2[,c("totalDNA", "PreWash", "PostWash")], use="complete.obs")
corrplot(corrWeight,method="circle", type="upper", tl.col="black", addCoef.col="black")

#ggplot(data=filtered, aes(x=PostWash,y=totalDNA))+
  #geom_point()
  
