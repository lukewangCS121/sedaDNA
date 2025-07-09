library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(mgcv)
library(lubridate)
library(patchwork)
library(gratia)
library(corrplot)
hasDepth<-LukeREU|>
  filter(!is.na(Core.Depth.cm.)&!is.na(total.DNA..ng.))
# hasDepth$total.DNA..ng.<-as.numeric(hasDepth$total.DNA..ng.)
# hasDepth$Core.Depth.cm.<-as.numeric(hasDepth$Core.Depth.cm.) 
# hasDepth$PostWashMass<-as.numeric(hasDepth$PostWashMass)
#DNA over total depth
ggplot(data=hasDepth, aes(x=Core.Depth.cm., y=total.DNA..ng., color=factor(Site)))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue")+
  labs(x="Core Depth(cm)", y="Total DNA(ng)")
dnaDepthLM<-lm(total.DNA..ng.~Core.Depth.cm., data=hasDepth)
# hasDepth<-hasDepth|>
#   rename(PostWashMass=Sample.mass.post.PBS.wash..added.by.Max.)
summary(dnaDepthLM)#R squared of 0.1029

#TODO:high oxygen vs low oxygen vs medium oxygen multiplot

#gamplot:
hasDepth$gam_pred<-predict(dnaDepthGAM)
plot(dnaDepthGAM, pages = 3, residuals = TRUE, shade = TRUE, seWithMean = TRUE)
draw(dnaDepthGAM)

#elemental
hasDepth$RbNormal=hasDepth$Rb.Ka.Area/sum(hasDepth$Rb.Ka.Area,na.rm=TRUE)
ggplot(data=hasDepth,aes(x=RbNormal,y=total.DNA..ng.,na.rm=TRUE))+
  geom_point()#supposed to be well
hasDepth$SrNormal=hasDepth$Sr.Ka.Area/sum(hasDepth$Sr.Ka.Area,na.rm=TRUE)
ggplot(data=hasDepth,aes(x=SrNormal,y=total.DNA..ng.,na.rm=TRUE))+
  geom_point()
hasDepth$YNormal=hasDepth$Y.Ka.Area/sum(hasDepth$Y.Ka.Area,na.rm=TRUE)
ggplot(data=hasDepth,aes(x=YNormal,y=total.DNA..ng.,na.rm=TRUE))+
  geom_point()
hasDepth$ZrNormal=hasDepth$Zr.Ka.Area/sum(hasDepth$Zr.Ka.Area,na.rm=TRUE)
ggplot(data=hasDepth,aes(x=ZrNormal,y=total.DNA..ng.,na.rm=TRUE))+
  geom_point()

#linear models for elements:
RbLM<-lm(total.DNA..ng.~Core.Depth.cm.+RbNormal+SrNormal+YNormal+ZrNormal, data=hasDepth)
summary(RbLM)
# EleGAM<-gam(total.DNA..ng.~Core.Depth.cm.+s(RbNormal)+s(SrNormal)+s(YNormal)+s(ZrNormal),data=hasDepth)
# summary(EleGAM)
# plot(EleGAM, pages = 3, residuals = TRUE, shade = TRUE, seWithMean = TRUE)


hasDNA<-LukeREU|>
  filter(!is.na(total.DNA..ng.))
hasDNA$total.DNA..ng.<-as.numeric(LukeREU$total.DNA..ng.)
hasDNA<-hasDNA|>
  filter((total.DNA..ng.)!=0)
#Which date had the most successful extractions
ggplot(data=hasDNA, aes(x=Date.of.extraction))+
   geom_bar()
#Which sites have the most DNA
ggplot(data=hasDNA, aes(x=factor(Site)))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
#average DNA per site
hasDNA|>
  group_by(Site)|>
  summarize(MeanDNABySite=mean(total.DNA..ng.))