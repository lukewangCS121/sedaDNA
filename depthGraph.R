library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
extraction_data<-extraction_data|>
  mutate(depth_cm=as.numeric(str_extract(V2, "\\d{2}(?=CM)")))
filtered<-extraction_data|>
  filter(!is.na(depth_cm))
filtered$totalDNA<-as.numeric(filtered$totalDNA)
ggplot(data=filtered, aes(x=depth_cm,y=totalDNA, color=V1))+
  #already renamed V25 to totalDNA
  geom_point()+
  scale_y_continuous(breaks = seq(100, 1000, by = 100), limits = c(0, 1000))+
  theme_minimal()
