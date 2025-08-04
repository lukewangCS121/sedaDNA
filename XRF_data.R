library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(lme4)
library(stringr)
library(broom)
library(multcomp)


MC10KV<-read_csv("/Users/lukewang/Downloads/processed_sbb1_xrf_data_6_27_2025/SBB1_MC1_10kv.csv")

MC30KV<-read_csv("/Users/lukewang/Downloads/processed_sbb1_xrf_data_6_27_2025/SBB1_MC1_30kv.csv")

MC50KV<-read_csv("/Users/lukewang/Downloads/processed_sbb1_xrf_data_6_27_2025/SBB1_MC1_50kv.csv")

MC10K_avg<-MC10KV|>
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

MC30K_avg<-MC30KV|>
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

MC50K_avg<-MC50KV|>
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

MC10K_avg$source <- "MC10K"
MC30K_avg$source <- "MC30K"
MC50K_avg$source <- "MC50K"

# Combine all into one long-format table
combined_avg <- bind_rows(MC10K_avg, MC30K_avg, MC50K_avg)

# Step 2: Pivot to long format
long_avg <- combined_avg %>%
  pivot_longer(
    cols = -source,
    names_to = "element",
    values_to = "avg_value"
  )

# Step 3: For each element, keep the row with the highest average
highest_avg_elements <- long_avg %>%
  group_by(element) %>%
  slice_max(avg_value, n = 1) %>%
  ungroup()
highest_avg_elements<-highest_avg_elements|>
  filter(str_detect(element, "Area$"))


MC10KV <- MC10KV%>%
  rename("real_depth_cm" = "Real Depth(cm)")

MC30KV <- MC30KV%>%
  rename("real_depth_cm" = "Real Depth(cm)")

MC50KV <- MC50KV%>%
  rename("real_depth_cm" = "Real Depth(cm)")

MC10KVpivot <- MC10KV %>%
  pivot_longer(cols= ends_with("Area"), 
               names_to = "element",
               values_to = "element_counts")

MC30KVpivot<- MC30KV %>%
  pivot_longer(cols= ends_with("Area"), 
               names_to = "element",
               values_to = "element_counts")

MC50KVpivot<- MC50KV %>%
  pivot_longer(cols= ends_with("Area"), 
               names_to = "element",
               values_to = "element_counts")

ggplot(data=MC10KVpivot,aes(x=`real_depth_cm`,y=log(`element_counts`), color = element))+
  geom_line()
#average into 1cm depths
MC10kv_summary<-MC10KVpivot|>
  mutate(depth_bin=floor(real_depth_cm))|>
  group_by(depth_bin,element)|>
  summarize(avg_element_counts=mean(element_counts,na.rm=TRUE),.groups = "drop")

MC30kv_summary<-MC30KVpivot|>
  mutate(depth_bin=floor(real_depth_cm))|>
  group_by(depth_bin,element)|>
  summarize(avg_element_counts=mean(element_counts,na.rm=TRUE),.groups = "drop")

MC50kv_summary<-MC50KVpivot|>
  mutate(depth_bin=floor(real_depth_cm))|>
  group_by(depth_bin,element)|>
  summarize(avg_element_counts=mean(element_counts,na.rm=TRUE),.groups = "drop")
#0 cm means 0.1,0.2,0.3... basically 0 is 0-1, 1 is 1-2 on and on

#average element by percent abundance of total
MC10KV_percent<-MC10kv_summary|>
  group_by(element)|>
  mutate(total_mass=sum(avg_element_counts),
         percent=(avg_element_counts/total_mass)*100)

MC30KV_percent<-MC30kv_summary|>
  group_by(element)|>
  mutate(total_mass=sum(avg_element_counts),
         percent=(avg_element_counts/total_mass)*100)

MC50KV_percent<-MC50kv_summary|>
  group_by(element)|>
  mutate(total_mass=sum(avg_element_counts),
         percent=(avg_element_counts/total_mass)*100)

MC10KV_percent <- MC10KV_percent %>% mutate(source = "MC10K")
MC30KV_percent <- MC30KV_percent %>% mutate(source = "MC30K")
MC50KV_percent <- MC50KV_percent %>% mutate(source = "MC50K")

# Combine all three percent tables
all_percent <- bind_rows(MC10KV_percent, MC30KV_percent, MC50KV_percent)

# Join to get the percent for each element from the filtered average table
final_element_df <- highest_avg_elements %>%
  left_join(all_percent, by = c("source", "element"))
final_element_df<-final_element_df|>
  filter(avg_value>0&element!="Ag-Ka Area"&element!="Rh-Ka Inc Area"&element!="Rh-Ka Coh Area"&element!="Cr-Ka Area")

final_wide<-final_element_df%>%
group_by(depth_bin, element) %>%
  summarise(percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from=element, values_from=percent)

#exclude out elemental with negative abundances
MC10KV_percent<-MC10KV_percent|>
  filter(element!="Rh-La-Inc Area" & element!='Ar-Ka Area'&element!="Cr-Ka Area")


# MC_5cm<-MC10KV_percent|>
#   mutate(depth_5cm = floor(depth_bin / 5) * 5)
# MC_5cm_sum<-MC_5cm|>
#   group_by(depth_5cm, element)|>
#   summarize(abundance_percent_5cm=sum(percent))
# MC_5cm_sum<-MC_5cm_sum|>
#   rename("Core_Depth"="depth_5cm")

#combine abundance with DNA yields
MaxSample<-read_csv("/Users/lukewang/Downloads/Luke Wang Copy of aces_cruise_data_SR2422 - extraction_data (2).csv")
MaxSample<-MaxSample|>
  rename("depth_bin"="CoreDepth(cm)")
Dna_joined<-left_join(MC10KV_percent,MaxSample,by="depth_bin")
#make 30KV into each row for each depth and columns for elements
MC30KV_row_depth<-MC30KV_percent|>
  group_by(depth_bin, element)|>
  summarise(percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from=element, values_from=percent)

MC50kv_row_depth<-MC50KV_percent|>
  group_by(depth_bin, element)|>
  summarise(percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from=element, values_from=percent)

Dna_joined<-Dna_joined|>
  rename("Total_DNA"="total DNA (ng)")
Dna_joined<-Dna_joined|>
  filter(!is.na("Sample ID")&Total_DNA!=0)

alex_fixed_dna_joined<-left_join(final_wide,MaxSample,by="depth_bin")
alex_fixed_dna_joined<-alex_fixed_dna_joined|>
  rename("Total_DNA"='total DNA (ng)')
alex_fixed_dna_joined$Total_DNA<-as.numeric(alex_fixed_dna_joined$Total_DNA)
alex_fixed_dna_joined<- alex_fixed_dna_joined %>%
  rename_with(~ gsub("[- ]", "_", .x))
alex_fixed_dna_joined<-alex_fixed_dna_joined%>%
  select(-contains("Rh"))


#combine table into each row being for a depth
# dna_wide<-Dna_joined|>
#   group_by(depth_bin, element) %>%
#   summarise(percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
#   pivot_wider(names_from=element, values_from=percent)
# wide_joined<-left_join(dna_wide,MaxSample)
# wide_joined<-wide_joined|>
#   rename("Total_DNA"='total DNA (ng)')
# wide_joined$Total_DNA<-as.numeric(wide_joined$Total_DNA)
# 
# all_energy_levels<-left_join(MC50kv_row_depth,new_wide,by="depth_bin")
# all_energy_levels<- all_energy_levels %>%
#   rename_with(~ gsub("[- ]", "_", .x))



# 
# new_wide<-left_join(MC30KV_row_depth,wide_joined,by="depth_bin")
# new_wide<-new_wide|>
#   select(-`Fe-Ka Area.x`,-`Mn-Ka Area.x`)
# new_wide<-new_wide|>
#   filter(!is.na(`Al-Ka Area`))|>
#   rename("Total_DNA"="total DNA (ng)")
# new_wide$Total_DNA<-as.numeric(new_wide$Total_DNA)
# new_clean_df <- new_wide %>%
#   rename_with(~ gsub("[- ]", "_", .x))


#LM of all elements and depth of 10kv
element_pred<-alex_fixed_dna_joined|>
  select(contains("_Area"))|>
  colnames()
all_pred<-c("depth_bin",element_pred)
formula <- as.formula(paste("Total_DNA ~", paste(all_pred, collapse = " + ")))
#making a formula for all elements across all levels
all_ele_lm<-lm(formula,data=alex_fixed_dna_joined)
summary(all_ele_lm)


all_30kv_lm<-lm(formula, data=new_clean_df)
summary(all_30kv_lm)

all_lm<-lm(Total_DNA~`Ca-Ka Area`+`Mg-Ka Area`+`Fe-Ka Area`+`P -Ka Area`+`S -Ka Area`+`Cl-Ka Area`+`Al-Ka Area`+`K -Ka Area`+`Mn-Ka Area`+`Si-Ka Area`+`Ti-Ka Area`+`Ba-La Area`+`Na-Ka Area`+`Rh-La Area`+`V -Ka Area`+depth_bin, 
          data=wide_joined)#included with depth
summary(all_lm)

clean_df <- wide_joined %>%
  rename_with(~ gsub("[- ]", "_", .x))

#Sulfur:
ggplot(data=alex_fixed_dna_joined,aes(x=S__Ka_Area,y=Total_DNA))+
  geom_point()+
  geom_smooth(method="lm")
s_lm<-lm(Total_DNA~S__Ka_Area, data=alex_fixed_dna_joined)
summary(s_lm)

#Lead:
ggplot(data=alex_fixed_dna_joined,aes(x=Pb_La_Area,y=Total_DNA))+
  geom_point()+
  geom_smooth(method="lm")
pb_lm<-lm(Total_DNA~Pb_La_Area, data=alex_fixed_dna_joined)
summary(pb_lm)

ggplot(data=alex_fixed_dna_joined,aes(x=Hg_La_Area,y=Total_DNA))+
  geom_point()+
  geom_smooth(method="lm")
summary(lm(Total_DNA~Hg_La_Area, data=alex_fixed_dna_joined))
cor(alex_fixed_dna_joined$Total_DNA,alex_fixed_dna_joined$Hg_La_Area, use = "complete.obs")


#correlation plot:
correlem<-cor(clean_df[,c("depth_bin", "Ca_Ka_Area","Mg_Ka_Area","Fe_Ka_Area","P__Ka_Area","S__Ka_Area",
                          "Cl_Ka_Area","Al_Ka_Area","Mn_Ka_Area","Si_Ka_Area","Ba_La_Area","Rh_La_Area",
                          "V__Ka_Area")], use="complete.obs")
corrplot(correlem, method = "circle", type = "upper", tl.col = "black", addCoef.col = "black")

corr30kv<-cor(new_clean_df[,element_pred], use="complete.obs")
corrplot(corr30kv, method = "circle", type = "upper", tl.col = "black", addCoef.col = "black")


#compare sites
ExtractionRegister<-ExtractionRegister|>
  rename("Total_DNA"="total.DNA..ng.")
ExtractionRegister$Total_DNA<-as.numeric(ExtractionRegister$Total_DNA)
hasDNA<-ExtractionRegister|>
  filter(Total_DNA!=0&!is.na(Total_DNA))
hasDNA<-hasDNA|>
  filter(Site!="Do Not Sequence! Not outside EEZ!")

hasDNA|>
  group_by(Basin)|>
  summarize(MeanDNABySite=mean(Total_DNA),
            DNA_SD=sd(Total_DNA),
            counts=n())

hasDNA<-hasDNA|>
  mutate(oxygen_level=case_when(
    Basin=="Santa Barbara Basin"~"Low",
    Basin=="Santa Monica Basin"~"Moderate",
    Basin=="Santa Cruz Basin"~"High"
  ))

hasDNA<-hasDNA|>
  rename("depth_bin"="Core.Depth.cm.")
#inter-basin comparison
basinDNA_1way<-hasDNA %>%
  group_by(Basin) %>%
  filter(n() > 1&Basin!="Caymans") %>%
  ungroup() %>%
  oneway.test(Total_DNA~Basin, data=.)

summary(basinDNA_1way)

basinDNA_aov<-hasDNA%>%
  group_by(Basin) %>%
  filter(n() > 1&Basin!="Caymans") %>%
  ungroup() %>%
  aov(Total_DNA~Basin, data=.)

summary(basinDNA_aov)
TukeyHSD(basinDNA_aov)
plot(TukeyHSD(basinDNA_aov))

hasDNA$Basin_Site <- with(hasDNA, paste(Basin, Site, sep = "_"))
hasDNA$Basin_Site <- as.factor(hasDNA$Basin_Site)

hasDNA <- hasDNA %>%
  mutate(Basin_Site = ifelse(Basin == "Santa Barbara Basin",
                             paste(Basin, Site, sep = "_"),
                             Basin)) # All other basins keep their name only

noEmpty <- hasDNA %>% filter(Basin != "")

# Fit the model AFTER defining Basin_Site correctly
site_aov <- aov(Total_DNA ~ Basin_Site, data = noEmpty)

tukey_df <- as.data.frame(TukeyHSD(site_aov)$Basin_Site)
tukey_df$comparison <- rownames(tukey_df)

# Split into two site names
tukey_df <- separate(tukey_df, comparison, into = c("group1", "group2"), sep = "-")

# Reorder alphabetically so A-B and B-A become the same
tukey_df$pair <- apply(tukey_df[, c("group1", "group2")], 1, function(x) paste(sort(x), collapse = "-"))

# Remove duplicates
tukey_df_unique <- tukey_df[!duplicated(tukey_df$pair), ]

ggplot(tukey_df_unique, aes(x = diff, y = pair)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lwr, xmax = upr)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Difference in Mean Total DNA", y = "Site Comparisons")





hasDNA<-hasDNA|>
  rename("Water_Depth"="Water.Depth..meters.")

all_basins_var<-lm(Total_DNA~depth_bin+Water_Depth+Basin+oxygen_level, data=hasDNA)
summary(all_basins_var)
anova(all_basins_var)
aov_all_basin<-aov(all_basins_var)
TukeyHSD(aov_all_basin, "Water_Depth")


#Max's model:
MAX_model <- lmer(Total_DNA ~ Basin + depth_bin + (1|Site), data = hasDNA)
summary(MAX_model)

#AOV Tukey test
TukeyHSD(basinDNA_aov)
 
#comparing sites within SBB
onlySBB<-hasDNA|>
  filter(Basin=="Santa Barbara Basin")
onlySBB|>
  group_by(Site)|>
  summarise(MeanSiteDNA=mean(Total_DNA),
            count=n())#TODO: add Max's site 1 DNA
onlySBB%>%
  group_by(Site)%>%
  filter(n()>1)%>%
  ungroup()%>%
  oneway.test(Total_DNA~Site, data=.)

SBBSitesDiff<-onlySBB%>%
  group_by(Site)%>%
  filter(n()>1)%>%
  ungroup()%>%
  aov(Total_DNA~Site, data=.)
TukeyHSD(SBBSitesDiff)
plot(TukeyHSD(SBBSitesDiff))


#forest plot for elemental contribution to DNA
coef_elements<-tidy(all_lm,conf.int=TRUE)#95% level
coef_elements|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()


coef_30kv<-tidy(all_30kv_lm,conf.int=TRUE)
coef_30kv|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()

coef_all_energy<-tidy(all_ele_lm, conf.int=TRUE)
coef_all_energy|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()



# #Sulfur
# ggplot(data=all_energy_levels, aes(x=S__Ka_Area, y=Total_DNA))+
#   geom_point()+
#   geom_smooth(method="lm")
# s_lm<-lm(Total_DNA~S__Ka_Area, data=all_energy_levels)
# summary(s_lm)
# #Iron
# ggplot(data=all_energy_levels, aes(x=Fe_Ka_Area, y=Total_DNA))+
#   geom_point()
# fe_lm<-lm(Total_DNA~Fe_Ka_Area, data=all_energy_levels)
# summary(fe_lm)
# #lead
# ggplot(data=all_energy_levels, aes(x=Pb_La_Area.x, y=Total_DNA))+
#   geom_point()+
#   geom_smooth(method="lm")
# summary(lm(Total_DNA~Pb_La_Area.x, data=all_energy_levels))
# #potassium:
# ggplot(data=all_energy_levels, aes(x=K__Ka_Area, y=Total_DNA))+
#   geom_point()+
#   geom_smooth(method="lm")
# summary(lm(Total_DNA~K__Ka_Area, data=all_energy_levels))
# #barium:
# ggplot(data=all_energy_levels, aes(x=Ba_La_Area, y=Total_DNA))+
#   geom_point()



#all basin comparison:
hasDNA$depth_bin<-as.numeric(hasDNA$depth_bin)
hasDNA$Water.Depth..meters.<-as.numeric(hasDNA$Water.Depth..meters.)
hasDNA$QUBIT.DNA.conc..ng...uL.<-as.numeric(hasDNA$QUBIT.DNA.conc..ng...uL.)
hasDNA<-hasDNA|>
  rename("Qubit_Concentration"="QUBIT.DNA.conc..ng...uL.")
filtered_DNA<-hasDNA|>
  filter(Total_DNA>2000)


ggplot(data=filtered_DNA, aes(x=depth_bin, y=Total_DNA))+
  geom_point()

just_depth<-lm(Total_DNA~depth_bin, data=filtered_DNA)
summary(just_depth)

ggplot(data=onlySBB, aes(x=Site, y=Total_DNA))+
  geom_boxplot()

#TODO:

