library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
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
  summarize(avg_element_counts=mean(element_counts,na.rm=TRUE),.groups = "drop")#0 cm means 0.1,0.2,0.3... basically 0 is 0-1, 1 is 1-2 on and on

ggplot(data=filter(MC10kv_summary,element!="Rh-La-Inc Area"),aes(x=depth_bin,y=avg_element_counts,color=element))+
  geom_line()

