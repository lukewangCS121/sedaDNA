library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(lme4)
library(stringr)
library(broom)
library(multcomp)
library(patchwork)
library(corrplot)

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
MaxSample<-read_csv("/Users/lukewang/Downloads/Luke Wang Copy of aces_cruise_data_SR2422 - extraction_data (3).csv")
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


# #LM of all elements and depth of 10kv
# element_pred<-alex_fixed_dna_joined|>
#   select(contains("_Area"))|>
#   colnames()
# all_pred<-c("depth_bin",element_pred)
# formula <- as.formula(paste("Total_DNA ~", paste(all_pred, collapse = " + ")))
#making a formula for all elements across all levels
all_ele_lm<-lm(formula,data=alex_fixed_dna_joined)
summary(all_ele_lm)


# all_30kv_lm<-lm(formula, data=new_clean_df)
# summary(all_30kv_lm)
# 
# all_lm<-lm(Total_DNA~`Ca-Ka Area`+`Mg-Ka Area`+`Fe-Ka Area`+`P -Ka Area`+`S -Ka Area`+`Cl-Ka Area`+`Al-Ka Area`+`K -Ka Area`+`Mn-Ka Area`+`Si-Ka Area`+`Ti-Ka Area`+`Ba-La Area`+`Na-Ka Area`+`Rh-La Area`+`V -Ka Area`+depth_bin, 
#           data=wide_joined)#included with depth
# summary(all_lm)

# clean_df <- wide_joined %>%
#   rename_with(~ gsub("[- ]", "_", .x))
# 
# 
zn<-ggplot(data=filter(alex_fixed_dna_joined, depth_bin<68), aes(x=depth_bin,y=Zn_Ka_Area))+
  geom_line()+
  labs(x="Core Depth(cm)", y="Zinc Relative Abundance %")
Y<-ggplot(data=filter(alex_fixed_dna_joined, depth_bin<68), aes(x=depth_bin,y=Y__Ka_Area))+
  geom_line()+
  labs(x="Core Depth(cm)", y="Yttrium Relative Abundance %")
Ti<-ggplot(data=filter(alex_fixed_dna_joined, depth_bin<68), aes(x=depth_bin,y=Ti_Ka_Area))+
  geom_line()+
  labs(x="Core Depth(cm)", y="Titanium Relative Abundance %")
rb<-ggplot(data=filter(alex_fixed_dna_joined, depth_bin<68), aes(x=depth_bin,y=Rb_Ka_Area))+
  geom_line()+
  labs(x="Core Depth(cm)", y="Rubidium Relative Abundance %")
mn<-ggplot(data=filter(alex_fixed_dna_joined, depth_bin<68), aes(x=depth_bin,y=Mn_Ka_Area))+
  geom_line()+
  labs(x="Core Depth(cm)", y="Manganese Relative Abundance %")
k_ka<-ggplot(data=filter(alex_fixed_dna_joined, depth_bin<68), aes(x=depth_bin,y=K__Ka_Area))+
  geom_line()+
  labs(x="Core Depth(cm)", y="Potassium Relative Abundance %")
fe<-ggplot(data=filter(alex_fixed_dna_joined, depth_bin<68), aes(x=depth_bin,y=Zn_Ka_Area))+
  geom_line()+
  labs(x="Core Depth(cm)", y="Iron Relative Abundance %")
zn+Y+Ti+rb+mn+k_ka+fe


# #correlation plot:
# correlem<-cor(clean_df[,c("depth_bin", "Ca_Ka_Area","Mg_Ka_Area","Fe_Ka_Area","P__Ka_Area","S__Ka_Area",
#                           "Cl_Ka_Area","Al_Ka_Area","Mn_Ka_Area","Si_Ka_Area","Ba_La_Area","Rh_La_Area",
#                           "V__Ka_Area")], use="complete.obs")
# corrplot(correlem, method = "circle", type = "upper", tl.col = "black", addCoef.col = "black")
# 
# corr30kv<-cor(new_clean_df[,element_pred], use="complete.obs")
# corrplot(corr30kv, method = "circle", type = "upper", tl.col = "black", addCoef.col = "black")

corr_final<-cor(alex_fixed_dna_joined[c("depth_bin","S__Ka_Area","Fe_Ka_Area", "Zn_Ka_Area")])
corrplot(corr_final, method = "circle", type = "upper", tl.col = "black", addCoef.col = "black")


#compare sites
ExtractionRegister<-ExtractionRegister|>
  rename("Total_DNA"="total.DNA..ng.")
ExtractionRegister$Total_DNA<-as.numeric(ExtractionRegister$Total_DNA)
hasDNA<-ExtractionRegister|>
  filter(Total_DNA!=0&!is.na(Total_DNA))
hasDNA<-hasDNA|>
  filter(Site!="Do Not Sequence! Not outside EEZ!")
hasDNA_combined|>
  group_by(Basin)|>
  summarize(MeanDNABySite=mean(Total_DNA),
            DNA_SD=sd(Total_DNA),
            counts=n())

hasDNA<-hasDNA|>
  rename("depth_bin"="Core.Depth.cm.")
hasDNA<-hasDNA|>
  rename("Extraction_Kit"="Extraction.kit")
#inter-basin comparison
hasDNA_combined %>%
  group_by(Basin) %>%
  filter(n() > 1&Basin!="Caymans"&!is.na(Basin)) %>%
  ungroup() %>%
  oneway.test(Total_DNA~Basin+depth_bin, data=.)
summary(basinDNA_1way)

interbasin<-lm(Total_DNA~Basin+depth_bin, 
               data=filter(hasDNA_combined, Basin!="Caymans"&Basin!=""))

anova(interbasin)
coef_inter<-tidy(interbasin,conf.int=TRUE)
coef_inter|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()#Relationship holds 

basinDNA_aov<-hasDNA_combined%>%
  group_by(Basin) %>%
  filter(n() > 1&Basin!="Caymans"&!is.na(Basin)) %>%
  ungroup() %>%
  aov(Total_DNA~Basin+depth_bin,data=.)

summary(basinDNA_aov)
TukeyHSD(basinDNA_aov)
plot(TukeyHSD(basinDNA_aov))

#Joining hasDNA with Max's Site 1 samples
MaxSample$depth_bin<-as.numeric(MaxSample$depth_bin)
hasDNA$depth_bin<-as.numeric(hasDNA$depth_bin)
MaxSample$Sample.mass <- as.numeric(MaxSample$Sample.mass)
MaxSample$Sample.mass.post.PBS.wash..added.by.Max. <- as.numeric(MaxSample$Sample.mass.post.PBS.wash..added.by.Max.)
MaxSample$Total.volume..µl.<-as.numeric(MaxSample$Total.volume..µl.)
hasDNA_combined<-bind_rows(hasDNA,MaxSample)
hasDNA_combined <- hasDNA_combined %>%
  mutate(
    Basin = ifelse(grepl("^SBB1", `Sample ID`), "Santa Barbara Basin", Basin),
    Site = ifelse(grepl("^SBB1", `Sample ID`), "1", Site)
  )
hasDNA_combined<-hasDNA_combined%>%
  mutate(Extraction_Kit = ifelse(grepl("POWERSOIL", Extraction.Box, ignore.case = TRUE),
                                 "PowerSoil", Extraction_Kit))
hasDNA_combined<-hasDNA_combined%>%
  mutate(Total_DNA = ifelse(is.na(Total_DNA), `total DNA (ng)`, Total_DNA))
hasDNA_combined<-hasDNA_combined%>%
  mutate(Water_Depth=ifelse(is.na(Water_Depth), `Water Depth (added by Alex) (meters)`,Water_Depth))
hasDNA_combined<-hasDNA_combined%>%
  mutate(CORE.ID..added.by.Max.<-
           ifelse(is.na(CORE.ID..added.by.Max.), `CORE ID (added by Max)`, CORE.ID..added.by.Max.))
hasDNA_combined$Total_DNA<-as.numeric(hasDNA_combined$Total_DNA)

# #Basin+Site comparison(ignore for now)
# hasDNA$Basin_Site <- with(hasDNA, paste(Basin, Site, sep = "_"))
# hasDNA$Basin_Site <- as.factor(hasDNA$Basin_Site)
# 
# hasDNA <- hasDNA %>%
#   mutate(Basin_Site = ifelse(Basin == "Santa Barbara Basin",
#                              paste(Basin, Site, sep = "_"),
#                              Basin)) # All other basins keep their name only
# 
# noEmpty <- hasDNA %>% filter(Basin != "")

# Fit the model AFTER defining Basin_Site correctly
# site_aov <- aov(Total_DNA ~ Basin_Site, data = noEmpty)
# 
# tukey_df <- as.data.frame(TukeyHSD(site_aov)$Basin_Site)
# tukey_df$comparison <- rownames(tukey_df)
# 
# # Split into two site names
# tukey_df <- separate(tukey_df, comparison, into = c("group1", "group2"), sep = "-")
# 
# # Reorder alphabetically so A-B and B-A become the same
# tukey_df$pair <- apply(tukey_df[, c("group1", "group2")], 1, function(x) paste(sort(x), collapse = "-"))
# 
# # Remove duplicates
# tukey_df_unique <- tukey_df[!duplicated(tukey_df$pair), ]
# tukey_df_unique$label <- paste(tukey_df_unique$group1, "-", tukey_df_unique$group2)
# ggplot(tukey_df_unique, aes(x = diff, y = label)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin = lwr, xmax = upr)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   theme_minimal() +
#   labs(x = "Difference in Mean Total DNA (Group2 - Group1)",
#        y = "Site Comparisons")
#end of basin_site

hasDNA<-hasDNA|>
  rename("Water_Depth"="Water.Depth..meters.")


#Max's model:
MAX_model <- lmer(Total_DNA ~ Basin + depth_bin + (1|Site), data = hasDNA)
summary(MAX_model)

 
#comparing sites within SBB
onlySBB<-hasDNA_combined|>
  filter(Basin=="Santa Barbara Basin")
onlySBB|>
  group_by(Extraction_Kit)|>
  summarise(MeanSiteDNA=mean(Total_DNA),
            count=n())
onlySBB%>%
  group_by(Site)%>%
  filter(n()>1)%>%
  ungroup()%>%
  oneway.test(Total_DNA~Site, data=.)

intrabasin<-lm(Total_DNA~Site+depth_bin, data=filter(onlySBB, Extraction_Kit=="PowerSoil"))
anova(intrabasin)
coef_intra<-tidy(intrabasin,conf.int=TRUE)
coef_intra|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()
table(filter(onlySBB, Extraction_Kit=="PowerSoil")$Site)


#sites 2 and 3 are deeper and more anoxic
#site 1 kind of on an ledge 
#TODO: add this to presentation + extraction_kit
SBBSitesDiff<-onlySBB%>%
  group_by(Site)%>%
  filter(n()>1)%>%
  ungroup()%>%
  aov(Total_DNA~Site+depth_bin, data=.)
TukeyHSD(SBBSitesDiff)
plot(TukeyHSD(SBBSitesDiff))

#QQ:
par(mfrow = (c(2, 2)))
plot(SBBSitesDiff)


#forest plot for elemental contribution to DNA
# coef_elements<-tidy(all_lm,conf.int=TRUE)#95% level
# coef_elements|>
#   filter(term!="(Intercept)")|>
#   ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
#   geom_pointrange()+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   coord_flip()
# 
# 
# coef_30kv<-tidy(all_30kv_lm,conf.int=TRUE)
# coef_30kv|>
#   filter(term!="(Intercept)")|>
#   ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
#   geom_pointrange()+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   coord_flip()
# 
coef_all_energy<-tidy(all_ele_lm, conf.int=TRUE)
coef_all_energy|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()

selected_elem_lm<-lm(Total_DNA~Fe_Ka_Area+S__Ka_Area+Zn_Ka_Area+depth_bin, 
                     data=alex_fixed_dna_joined)
summary(selected_elem_lm)
coef_selected_elem<-tidy(selected_elem_lm,conf.int=TRUE)
coef_selected_elem|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()

m<-cor(alex_fixed_dna_joined[,c("Fe_Ka_Area","S__Ka_Area","Zn_Ka_Area","depth_bin","Total_DNA", "Cu_Ka_Area",
                                "Mn_Ka_Area", "Ca_Ka_Area")], 
       use="complete.obs")

corrplot(m, method = "circle", type = "upper", tl.col = "black", addCoef.col = "black")


# library(car)
# vif(selected_elem_lm)
# par(mfrow = c(2, 2))e
# plot(selected_elem_lm)
# anova(selected_elem_lm)



#PCA:

selected_elements<-alex_fixed_dna_joined[,c("Al_Ka_Area", "As_Ka_Area", "Ba_Ka_Area", "Ba_La_Area", "Br_Ka_Area", "Ca_Ka_Area", "Cd_Ka_Area", "Cl_Ka_Area", "Co_Ka_Area", "Cu_Ka_Area", "Fe_Ka_Area", "Ga_Ka_Area", "Hg_La_Area", "I__Ka_Area", "K__Ka_Area", "Mg_Ka_Area", "Mn_Ka_Area", "Mo_Ka_Area", 
                                            "Na_Ka_Area", "Nb_Ka_Area", "Ni_Ka_Area", "P__Ka_Area", "Pb_La_Area", "Rb_Ka_Area", "S__Ka_Area", "Sb_Ka_Area", "Si_Ka_Area", "Sn_Ka_Area", "Sr_Ka_Area", "Te_Ka_Area", "Ti_Ka_Area", "U__La_Area", "V__Ka_Area", "Y__Ka_Area", "Zn_Ka_Area", "Zr_Ka_Area")]
pca_result <- prcomp(selected_elements, center = TRUE, scale. = TRUE)#metal ions
summary(pca_result)

p<-autoplot(pca_result, data = alex_fixed_dna_joined, 
         loadings = TRUE, loadings.label = TRUE)

loadings_df <- as.data.frame(pca_result$rotation[, 1:2])
loadings_df$Element <- rownames(loadings_df)

lm_smiple<-glm(Total_DNA~Al_Ka_Area+As_Ka_Area, data=alex_fixed_dna_joined)
summary(lm_smiple)

#Filter out PC1<-0.1, we select the elements vary the most, using covariates to see if they are correlated. 

loadings_filtered<-loadings_df|>
  filter(PC1<(-0.2))#Fe,K, Mn, Rb, Ti,Y, Zn
new_PC_elements<-lm(Total_DNA~Fe_Ka_Area+K__Ka_Area+Mn_Ka_Area+Rb_Ka_Area+Ti_Ka_Area+Y__Ka_Area+Zn_Ka_Area+depth_bin, 
                                      data=alex_fixed_dna_joined)
coef_PCA<-tidy(new_PC_elements,conf.int=TRUE)
coef_PCA|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()

# Add a highlight layer for Fe, Cu, Ca, Mn
highlight_vars <- c("Fe_Ka_Area", "Zn_Ka_Area", "Y__Ka_Area", "Mn_Ka_Area", "Ti_Ka_Area", "Rb_Ka_Area",
                    "K__Ka_Area","Rb_Ka_Area")

p + 
  geom_point(data = subset(loadings_df, Element %in% highlight_vars), 
             aes(x = PC1, y = PC2), 
             color = "black", size = 3) +
  geom_text(data = subset(loadings_df, Element %in% highlight_vars), 
            aes(x = PC1, y = PC2, label = Element), 
            color = "black", vjust = -0.5, fontface = "bold")



#all basin comparison:
hasDNA_combined$depth_bin<-as.numeric(hasDNA_combined$depth_bin)
hasDNA_combined$Water.Depth..meters.<-as.numeric(hasDNA_combined$Water.Depth..meters.)
hasDNA_combined$QUBIT.DNA.conc..ng...uL.<-as.numeric(hasDNA_combined$QUBIT.DNA.conc..ng...uL.)
hasDNA<-hasDNA|>
  rename("Qubit_Concentration"="QUBIT.DNA.conc..ng...uL.")
filtered_DNA<-hasDNA|>
  filter(Total_DNA>2000)
hasDNA<-hasDNA|>
  rename("Water_Depth"="Water.Depth..meters.")
hasDNA_combined<-hasDNA_combined|>
  rename("Water_Depth"="Water.Depth..meters.")

model_overall<-lm(Total_DNA~Basin, data=filter(hasDNA_combined, Basin!="Caymans"&Basin!=""&Extraction_Kit=="PowerMax"))
summary(model_overall)
coef_overall_model<-tidy(model_overall, conf.int=TRUE)
coef_overall_model|>
  filter(term!="(Intercept)")|>
  ggplot(aes(x=term,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 0,linetype="dashed")+
  coord_flip()

anova(model_overall)

#TODO: most important. Pearson's plot for water depth and basin
#random graphs
hasDNA_filtered <- hasDNA_combined %>%
  filter(Total_DNA>1000)

ggplot(data=filter(hasDNA_combined, Basin!="Caymans"&Basin!=""&Extraction_Kit=="PowerMax"&Basin!="Santa Monica Basin"),aes(x=depth_bin,y=log(Total_DNA), color=Basin))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Core Depth (cm)", y="Log(Total DNA)")

ggplot(data=hasDNA_filtered, aes(x=depth_bin, y=Total_DNA,color=Extraction_Kit))+
  geom_point()

just_depth<-lm(Total_DNA~depth_bin, data=hasDNA)
summary(just_depth)


ggplot(data=filter(onlySBB, Total_DNA>500), aes(x=Site, y=log(Total_DNA)))+
  geom_violin()

onlySBB<-filter(hasDNA_combined, Basin=="Santa Barbara Basin")
onlySBB$depth_bin<-as.numeric(onlySBB$depth_bin)

ggplot(data=filter(hasDNA_combined, Basin=="Santa Barbara Basin"|Basin=="Santa Cruz Basin"|Basin=="Santa Monica Basin"), 
       aes(x=depth_bin, y=log(Total_DNA),color=Basin))+
  geom_point()+
  geom_smooth(method="lm")


MaxSample$`total DNA (ng)`<-as.numeric(MaxSample$`total DNA (ng)`)
ggplot(data=filter(MaxSample, `total DNA (ng)`>100), aes(x=depth_bin, y=`total DNA (ng)`))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Depth(cm)", title="DNA across Depth within SBB1")

ggplot(data=filter(hasDNA_combined, Basin!=""),aes(x=Basin, y=log(Total_DNA)))+
  geom_violin()


filtered_SBB <- onlySBB %>%
  filter(Total_DNA < quantile(Total_DNA, 0.9, na.rm = TRUE)) 

ggplot(data=filter(onlySBB, depth_bin<=25&Extraction_Kit=="PowerSoil"), aes(x=depth_bin, y=log(Total_DNA),color=Site))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)+
  facet_wrap(~Site)+
  labs(x="Core Depth(cm)", y="log(Total DNA(ng))")

#Keep

