library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(mgcv)
library(lubridate)
library(patchwork)
library(gratia)
library(corrplot)
library(effects)

MC10KV<-read_csv("/Users/lukewang/Downloads/processed_sbb1_xrf_data_6_27_2025/SBB1_MC1_10kv.csv")

MC10KV <- MC10KV%>%
  rename("real_depth_cm" = "Real Depth(cm)")

MC10KVpivot <- MC10KV %>%
  pivot_longer(cols= ends_with("Area"), 
               names_to = "element",
               values_to = "element_counts")
ggplot(data=MC10KVpivot,aes(x=`real_depth_cm`,y=log(`element_counts`), color = element))+
  geom_line()

MC10kv_summary<-MC10KVpivot|>
  mutate(depth_bin=floor(real_depth_cm))|>
  group_by(depth_bin,element)|>
  summarize(avg_element_counts=mean(element_counts,na.rm=TRUE),.groups = "drop")
#0 cm means 0.1,0.2,0.3... basically 0 is 0-1, 1 is 1-2 on and on

ggplot(data=filter(MC10kv_summary,element!="Rh-La-Inc Area"),aes(x=depth_bin,y=avg_element_counts,color=element))+
  geom_line()

MC10KV_percent<-MC10kv_summary|>
  group_by(element)|>
  mutate(total_mass=sum(avg_element_counts),
         percent=(avg_element_counts/total_mass)*100)
ggplot(data=filter(MC10KV_percent, element!="Rh-La-Inc Area" & element!='Ar-Ka Area'&element!="Cr-Ka Area"), aes(x=depth_bin, y=percent, color=element))+
  geom_line()#excluded elements with negative values
MC10KV_percent<-MC10KV_percent|>
  filter(element!="Rh-La-Inc Area" & element!='Ar-Ka Area'&element!="Cr-Ka Area")
onlySBB1<-hasDepth|>
  filter(Site=="Santa Barbara Basin Site 1")

# MC_5cm<-MC10KV_percent|>
#   mutate(depth_5cm = floor(depth_bin / 5) * 5)
# MC_5cm_sum<-MC_5cm|>
#   group_by(depth_5cm, element)|>
#   summarize(abundance_percent_5cm=sum(percent))
# MC_5cm_sum<-MC_5cm_sum|>
#   rename("Core_Depth"="depth_5cm")

MaxSample<-read_csv("/Users/lukewang/Downloads/Luke Wang Copy of aces_cruise_data_SR2422 - extraction_data (2).csv")
MaxSample<-MaxSample|>
  rename("depth_bin"="CoreDepth(cm)")
Dna_joined<-left_join(MC10KV_percent,MaxSample,by="depth_bin")

Dna_joined<-Dna_joined|>
  rename("Total_DNA"="total DNA (ng)")
Dna_joined<-Dna_joined|>
  filter(!is.na("Sample ID")&Total_DNA!=0)

dna_wide<-Dna_joined|>
  group_by(depth_bin, element) %>%
  summarise(percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from=element, values_from=percent)
wide_joined<-left_join(dna_wide,MaxSample)
wide_joined<-wide_joined|>
  rename("Total_DNA"='total DNA (ng)')
wide_joined$Total_DNA<-as.numeric(wide_joined$Total_DNA)
ggplot(wide_joined, aes(x=`V -Ka Area`, y=Total_DNA))+
  geom_point()

ggplot(Dna_joined, aes(x=percent,y=as.numeric(Total_DNA), color=element))+
  geom_line()

#TODO: column for each element, each row is a cm, entry is average. build GLM of DNA~depth+add every element that matters
na_lm<-lm(Total_DNA~`Ca-Ka Area`+`Mg-Ka Area`+`Fe-Ka Area`+`P -Ka Area`+`S -Ka Area`+`Cl-Ka Area`+`Al-Ka Area`+`K -Ka Area`+`Mn-Ka Area`+`Si-Ka Area`+`Ti-Ka Area`+depth_bin, 
          data=wide_joined)#no Ba, Na, Rh,V
summary(na_lm)

clean_df <- wide_joined %>%
  rename_with(~ gsub("[- ]", "_", .x))

element_gam<-gam(Total_DNA~depth_bin+s(Ca_Ka_Area)+Mg_Ka_Area+Fe_Ka_Area+P__Ka_Area+S__Ka_Area+Cl_Ka_Area+Al_Ka_Area+Mn_Ka_Area+s(Si_Ka_Area)+Ba_La_Area+s(Rh_La_Area)+s(V__Ka_Area),
                 data=clean_df)
summary(element_gam)
plot(element_gam, pages = 1, residuals = TRUE, se = TRUE)
draw(element_gam)

visreg(element_gam,overlay=FALSE)


