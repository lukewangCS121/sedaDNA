library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(mgcv)
library(lubridate)
# install.packages("gratia")
# library(gratia)


extraction_data<-extraction_data|>
  mutate(depth_cm=as.numeric(str_extract(V2, "(?<=_)\\d+\\.?\\d*(?=CM)")))
filtered<-extraction_data|>
  filter(!is.na(depth_cm))
filtered$totalDNA<-as.numeric(filtered$totalDNA)
# filtered<-filtered|>
#   rename(PostWash=V13)
# filtered<-filtered|>
#   rename(WashMin=V29)
# filtered<-filtered|>
#   rename(ExDate=V8) 
filtered$ExDate<-as.character(filtered$ExDate)
filtered$ExDate<-dmy(filtered$ExDate)
filtered$ExDate_numeric <- as.numeric(filtered$ExDate - min(filtered$ExDate, na.rm = TRUE))
filtered$WashMin<-as.numeric(filtered$WashMin)
filtered$PostWash<-as.numeric(as.character(filtered$PostWash))

# ggplot(data=filtered, aes(x=depth_cm,y=totalDNA, color=V1))+
#   #already renamed V25 to totalDNA
#   geom_point()+
#   geom_smooth(method = "lm", se = FALSE, color = "blue")+
#   scale_y_continuous(breaks = seq(100, 1000, by = 100), limits = c(0, 1000))+
#   theme_minimal()+
#   labs(x="Core Depth(cm)", y="Total DNA(ng)")

model <- gam(totalDNA ~ depth_cm, data = filtered)
summary(model)
mod_sm<-gam(totalDNA ~ s(depth_cm)+ExDate_numeric+WashMin+s(PostWash), data = filtered, na.action = na.exclude)#s for smooth,not enough distinct values for washmin and date extracted
summary(mod_sm)
filtered$residuals<-residuals(mod_sm)
threshold<-2*sd(filtered$residuals, na.rm=TRUE)
filtered_clean<-filtered|>
  filter(abs(residuals)<threshold)
gam_clean <- gam(totalDNA ~s(depth_cm)+ExDate_numeric+WashMin+s(PostWash), data = filtered_clean)
summary(gam_clean)#0.573 r squared

ggplot(data=filtered_clean, aes(x=depth_cm,y=totalDNA, color=V1))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue")+
  scale_y_continuous(breaks = seq(100, 1000, by = 100), limits = c(0, 1000))+
  theme_minimal()+
  labs(x="Core Depth(cm)", y="Total DNA(ng)")

filtered_clean$gam_pred <- predict(gam_clean)
plot(gam_clean, pages = 3, residuals = TRUE, shade = TRUE, seWithMean = TRUE)

#ggplot(data=filtered, aes(x=PostWash,y=totalDNA))+
  #geom_point()
  
