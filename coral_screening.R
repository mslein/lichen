install.packages("revtools")
library(revtools)
library(tidyverse)
pacman::p_load(viridis, revtools, nlme, lme4, MuMIn, patchwork)

# 18 jan 2024
# ALL FIELDS = ("temperature" AND ("symbio*" OR "holobiont") AND ("respiration" OR "photosynthesis")) AND "coral"
# 384 papers ----> coral.bib

#coral1 <- read_bibliography("lit search files/coral.bib")
#screen_abstracts(coral1)

#72 candidate studies, 24 only have 2 temps
#17 isolated algae studies, 4 have 2 temps		



#29 jan 2024
# ABSTRACT = "zooxan*" AND "temperature" AND ("photosynthesis" OR "respiration")
#65 papers -----> zoox.bib
#coral2 <- read_bibliography("lit search files/zoox.bib")
#lichen2 <- read_bibliography("wos155.bib")
#screen_abstracts(coral2)

#### mixotroph search 
# 18 jan 
#ALL FIELDS = "mixotroph*" AND "photosynthesis" AND "temperature"
# 71 papers ----> mixotroph.bib
#mixo1 <- read_bibliography("lit search files/mixotroph.bib")
#3screen_abstracts(mixo1)



#### symbiodinium search
# 13 march 
# ABSTRACT = "symbiodinium" AND "temperature" and "photosynthesis"
# 41 papers ----> symbiodinium.bib
#symbio1 <- read_bibliography("lit search files/symbiodinium.bib")
#screen_abstracts(symbio1)



########## 16 April 2024


#so this study already log transformed these, i tried transforming them back but it doesn't
#seem like it worked

#corals
coral <- read_csv("extraction/coral dataset/coral/coral_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("log umol_o2") ~ (exp(response_value))*(1/10^6),
                         gas_units %in% c("mg_o2") ~ response_value*(1/32000),
                         gas_units %in% c("ug_o2") ~ response_value*(1/32000000),
                         gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                         gas_units %in% c("umol_o2") ~ response_value*(1/10^6),
                         gas_units %in% c("log umol_o2") ~ (response_value^10)/10^-6, 
                         gas_units %in% c("mmol_o2") ~ response_value*(1/10^3),
                         gas_units %in% c("pg_DIC") ~ response_value*(1/(1.2011*10^13)),
                         gas_units %in% c("uL_o2") ~ response_value*(4*10^-14), 
                         gas_units %in% c("umol_c") ~ response_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("468 mm^2") ~ 10^3/(468), 
                          `mass/area` %in% c("l") ~ 1/1000, 
                          `mass/area` %in% c("cm^2") ~ 1, 
                          `mass/area` %in% c("chla") ~ 1, 
                          `mass/area` %in% c("larvae^-1") ~ 1, 
                          `mass/area` %in% c("cell") ~ 1, 
                          `mass/area` %in% c("chla") ~ 1, 
                          `mass/area` %in% c("g ADFW") ~ 1, 
                          `mass/area` %in% c("g CaCO3") ~ 1,
                          `mass/area` %in% c("mg") ~ 10^3,
                          `mass/area` %in% c("zoox x 10^6") ~ 10^6),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60, 
                         time %in% c("day") ~ 1440),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral))




npp_coral <- coral %>%
  filter(metabolic_category == "npp") %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  mutate(centre_temp = I(inv_T-mean(inv_T))) %>%
  filter(case_when(study_id %in% c("castillohelmuth2005") ~ inv_T < 38, 
                   study_id %in% c("jiang2021") ~ inv_T < 38.25, 
                   TRUE ~ inv_T < 50))


r_coral <- coral %>%
  filter(metabolic_category == "r") %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  mutate(centre_temp = I(inv_T-mean(inv_T)))%>%
  filter(case_when(study_id %in% c("castillohelmuth2005") ~ inv_T < 38, 
                   study_id %in% c("dorey2020") ~ inv_T < 41, 
                   study_id %in% c("kemp2011") ~ inv_T < 40.5,
                   TRUE ~ inv_T < 50))

gpp_coral <- coral %>%
  filter(metabolic_category == "gpp")


npp_coral$depth_broad<- relevel(npp_coral$depth_broad, "shallow")
npp_coral$broad_coral<- relevel(npp_coral$broad_coral, "hard")
r_coral$depth_broad<- relevel(r_coral$depth_broad, "shallow")
r_coral$broad_coral<- relevel(r_coral$broad_coral, "hard")


npp_coral 


#random effects for npp
model00npp_coral<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|species), data=npp_coral, REML=FALSE)
model0npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id), data=npp_coral, REML=FALSE)
model1npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp + (1|study_id) , data=npp_coral, REML=FALSE)
model2npp_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp, data=npp_coral, na.action=na.exclude)
model.sel(model00npp_coral, model0npp_coral, model1npp_coral, model2npp_coral)
#model0npp_coral is best by > 2 AIC
#random effects for respiration
#this one has a singular fit
model00r_coral<- lmer(log(abs(mol_grammin)) ~abs(latitude) +centre_temp  + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|species) , data=r_coral, REML=FALSE)
model0r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=r_coral, REML=FALSE)
model1r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp + (1|study_id) , data=r_coral, REML=FALSE)
model2r_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp, data=r_coral, na.action=na.exclude)
model.sel(model00r_coral, model0r_coral, model1r_coral, model2r_coral)
#model0r_coral is best by > 2 AIC



#model0npp == best by 2 AIC
#model0r == best by 2 AIC



#model selection for photosynthesis
#model selection for photosynthesis
#missing depth
model3npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +centre_temp +  (1+centre_temp|study_id), data=npp_coral)
#missing latitude
model4npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + depth_broad +centre_temp + (1+centre_temp|study_id) , data=npp_coral)
#missing both latitude and depth
model5npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + (1+centre_temp|study_id), data=npp_coral)
#missing interaction between lichen type + centre temp
#missing temp
model7npp_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + depth_broad  + (1+centre_temp|study_id), data=npp_coral)
#full model
model9npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad + (1+centre_temp|study_id) , data=npp_coral)
model.sel(model3npp_coral, model4npp_coral, model5npp_coral, model7npp_coral, model9npp_coral)
#model 7, 9 are tied, going with model 9 because it has everything
summary(model9npp_coral)
confint(model9npp_coral)



#model selection for respiration
#missing depth
model3r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +centre_temp +  (1+centre_temp|study_id), data=r_coral)
#missing latitude
model4r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + depth_broad +centre_temp + (1+centre_temp|study_id) , data=r_coral)
#missing both latitude and depth
model5r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp +centre_temp + (1+centre_temp|study_id), data=r_coral)
#missing interaction between lichen type + centre temp
model6r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +broad_coral + (1+centre_temp|study_id), data=r_coral)
#missing temp
model7r_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + depth_broad +broad_coral + (1+centre_temp|study_id), data=r_coral)
#missing cnidarian type
model8r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad + (1+centre_temp|study_id), data=r_coral)
#full model
model9r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +centre_temp + (1+centre_temp|study_id) , data=r_coral)
model.sel(model3r_coral, model4r_coral, model5r_coral, model6r_coral, model7r_coral, model8r_coral, model9r_coral)


#model 6r is best by more than 2 AIC
summary(model6r_coral)
confint(model6r_coral)






#plots
#install.packages("remotes")
remotes::install_github("palaeoverse/rphylopic")
library(rphylopic)
#code to get the uuid's for the phylopic pngs
uuid <- rphylopic::get_uuid(name = "Acropora florida")



npp_coral_df<- count(npp_coral, study_id, latitude, depth_broad)

randslope_npp_coral<- random.effects(model9npp_coral) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope = centre_temp)

npp_joined_coral<- left_join(randslope_npp_coral,npp_coral_df, by="study_id") %>%
  mutate(c_intercept=case_when(depth_broad == "shallow" ~ -19.85209,
                               depth_broad == "deep" ~ -19.85209+-8.18040),
         latitudeii=abs(latitude)*-0.03817, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope+0.30833) %>%
  as.data.frame()



npp_coral_plot<- ggplot(data=npp_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="#abbb80")+
  geom_point(colour="#abbb80", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=npp_joined_coral, aes(slope=r_slope, intercept=f_intercept), colour="#abbb80", size=1, alpha=0.5)+
  geom_abline(aes(slope=-0.201, intercept=-20.21535), colour="#99a873", size=2)+
  #geom_hline(yintercept=-9.1055, size=1.5)+  
  xlab("Temperature (1/kT)")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("NPP")+
  ylab("")+
  xlab("Temperature (1/kT)")+
  ylim(-55,-5)+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=-3.7, y=-50, ysize=8, alpha=1,fill = "#abbb80")

  

r_coral_df<- count(r_coral, study_id, broad_coral, depth_broad, latitude)


randslope_r_coral<- random.effects(model6r_coral) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

r_joined_coral<- left_join(randslope_r_coral,r_coral_df, by="study_id") %>%
  as.data.frame() %>%
  mutate(c_intercept=case_when(broad_coral == "hard" & depth_broad == "shallow" ~ -20.721147,
                               broad_coral == "soft" & depth_broad == "shallow"~ -20.721147 +0.571205,
                               broad_coral == "hard" & depth_broad == "deep"~ -20.721147 +-2.57398 ,
                               broad_coral == "soft" & depth_broad == "deep"~ -20.721147 +0.571205+-2.57398), 
         latitudeii=abs(latitude)*-0.02755, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope+-0.552273) %>%
  as.data.frame() 

#double check the math on these

r_coral <- r_coral %>%
  mutate(temp_plus = temp+40, 
         inv_T_plus=(1/((8.617333262145*10^-5)*(temp_plus+273.15))),
    centre_temp_plus = I(inv_T_plus-mean(inv_T_plus)))
  
r_coral_plot <- ggplot(data=r_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="grey")+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_point(colour="grey", alpha=0.3, size=4)+
  geom_abline(data=r_joined_coral, aes(slope=f_slope, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=-0.552273, intercept = -20.721147, size=2.5, colour="grey65")+  
  xlab("Temperature (1/kT)")+
  ylab("")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("R")+
  ylim(-55,-5)+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=-3.7, y=-50, ysize=8, alpha=1,fill = "grey")
  #add_phylopic(uuid = "32eeaa29-2e90-40b5-92fb-eb8d7fcab9ab", x=2.5, y=-40, ysize=1.2, alpha=1,fill = "grey")
#scale_x_continuous(sec.axis = sec_axis(~.*(1/(8.61*10^-5))-273.15), name = "Temperature (°C)")
#scale_x_continuous(sec.axis = sec_axis(~ . * 20, name = "Temperature (°C)"))


gpp_coral_df<- count(gpp_comb_coral, study_id, depth_broad, latitude) %>% mutate(latitude=abs(latitude))


randslope_gpp_coral<- random.effects(model6gpp_coral) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

gpp_joined_coral<- left_join(randslope_gpp_coral,gpp_coral_df, by="study_id") %>%
  as.data.frame() %>%
  mutate(c_intercept=case_when(depth_broad == "shallow" ~ -26.194494,
                               depth_broad == "deep"~ -26.194494 +5.190041),
         latitudeii=abs(latitude)*-0.044934,
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope) %>%
  as.data.frame() 

gpp_comb_coral$gmpred <- predict(model4gpp_coral)

gpp_coral_plot <- ggplot(data=gpp_comb_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="olivedrab2")+
  geom_point(colour="olivedrab1", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=gpp_joined_coral, aes(slope=f_slope*-1, intercept=f_intercept), colour="olivedrab2", size=1, alpha=0.5)+
  geom_abline(slope=0, intercept = -26.194494, size=2.5, colour="olivedrab2")+  
  xlab("Temperature (1/kT)")+
  ylab("")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("GPP")+
  ylab("Metabolic rate (log(mol O2 /mg / min))")+
  ylim(-55,-5)+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=-3.7, y=-50, ysize=8, alpha=1,fill = "olivedrab2")



coral_plots <- (gpp_coral_plot + npp_coral_plot + r_coral_plot)

fig_c<- coral_plots + plot_layout(ncol=3) 



ggsave(fig_c, filename = "./figures/fig_c.png", dpi=700, width=16, height=5.5)







######## 
View(count(coral, study_id, metabolic_category))


#npp_and_r <- coral %>% 
  #filter(study_id %in% c("bahar2018", "godefroid2023", 
      #                  "hadjioannou2019", "hill2014", "howe2001", 
      #                  "jiang2021", "juillet-leclerc2014", "rodolfo-metalpa2006", 
      #                  "samiel2015"))

#write_csv(npp_and_r, "extraction/coral dataset/coral/gpp_arithmetic_coral1jun2024.csv")


gpp_arithmetic_coral <- read_csv("extraction/coral dataset/coral/gpp_arithmetic_coral1jun2024.csv") %>%
  mutate(r_value = abs(r_value)*-1,
         gpp = npp_value + r_value) %>%
  mutate(mol = case_when(gas_units %in% c("mg_o2") ~ gpp*(1/32000), 
                         gas_units %in% c("mmol_o2") ~ gpp*(1/10^9), 
                         gas_units %in% c("ug_o2") ~ gpp*(1/32000000),
                         gas_units %in% c("umol_o2") ~ gpp*(1/10^6)),
         gram = case_when(`mass/area` %in% c("468 mm^2") ~ 10^3/(468), 
                          `mass/area` %in% c("l") ~ 1/1000, 
                          `mass/area` %in% c("cm^2") ~ 1000, 
                          `mass/area` %in% c("chla") ~ 1, 
                          `mass/area` %in% c("g CaCO3") ~ 1),
         min = case_when(time %in% c("hr") ~ 60,
                         time %in% c("sec") ~ 1/60),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min), 
         centre_celsius = I(temp-mean(temp)), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral))%>% 
  select(study_id, temp, broad_coral, depth_broad, latitude, species, mol_grammin)


gpp_coral_collect <- gpp_coral %>% select(study_id, temp, broad_coral, depth_broad, latitude, species, mol_grammin)


gpp_comb_coral <- rbind(gpp_coral_collect, gpp_arithmetic_coral) %>%
  filter(study_id != "becker2021") %>% 
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))), 
                centre_temp = I(inv_T-mean(inv_T))) %>%
  filter(case_when(study_id %in% c("banc-prandi2022", "hill2014") ~ inv_T < 38.5, 
                   study_id %in% c("castillohelmuth2005") ~ inv_T < 38, 
                   study_id %in% c("higuchi2015") ~ inv_T < 40, 
                   study_id %in% c("howe2001", "kemp2011") ~ inv_T < 40.5, 
                   study_id %in% c("jiang2021", "jones1998") ~ inv_T < 38.25, 
                   study_id %in% c("juillet-leclerc2014") ~ inv_T < 39, 
                   TRUE ~ inv_T < 50))



gpp_comb_coral$depth_broad<- relevel(gpp_comb_coral$depth_broad, "shallow")

model00gpp<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +depth_broad + (1+centre_temp|species) , data=gpp_comb_coral, REML=FALSE)
model0gpp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +depth_broad + (1+centre_temp|study_id) , data=gpp_comb_coral, REML=FALSE)
model1gpp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +depth_broad + (1|study_id) , data=gpp_comb_coral, REML=FALSE)
model2gpp <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +depth_broad, data=gpp_comb_coral, na.action=na.exclude)

model.sel(model00gpp, model0gpp, model1gpp, model2gpp)

#model0gpp is best by more than 2 AIC


model3gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +depth_broad + (1+centre_temp|study_id), data=gpp_comb_coral)
#missing latitude
model4gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + depth_broad + (1+centre_temp|study_id) , data=gpp_comb_coral)
#missing both latitude and elevation
model5gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp  + (1+centre_temp|study_id), data=gpp_comb_coral)
#missing temp
model6gpp_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + depth_broad + (1+centre_temp|study_id), data=gpp_comb_coral)
model.sel(model3gpp_coral, model4gpp_coral, model5gpp_coral, model6gpp_coral)

#model6gpp_coral is tied w 3, going w/ 6
summary(model6gpp_coral)
confint(model6gpp_coral)





#creating the friedmann and sun figure 



#has npp, gpp, and r
#all_coral <- coral %>%
  #filter(study_id %in% c("nakamura2003","jiang2021", "castillohelmuth2005", 
                        # "silbiger2019", "hadjioannou2019")) 
#write_csv(all_coral, "f+s_npp_gpp_r.csv")

#has npp and r
#npp_and_r_coral <- coral %>%
  #filter(study_id %in% c("bahr2018", "godefroid2023", "hill2014", "howe2001", 
                         #"juillet-leclerc2014", "rodolfo-metalpa2006", "samiel2015"))
#write_csv(npp_and_r_coral, "f+s_npp_r.csv")
  

#has gpp and r 
#omitting banc-prandi for now because their gpp data literally makes no fucking sense
#gpp_and_r_coral <- coral %>%
  #filter(study_id %in% c("higuchi2015", "kemp2011"))
#write_csv(gpp_and_r_coral, "f+s_gpp_r.csv")



f_s_1 <- read_csv("f+s_coral_gpp_r.csv") %>% 
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         npp_mol = case_when(gas_units %in% c("mg_o2") ~ npp*(1/32000),
                         gas_units %in% c("nmol_o2") ~ npp*(1/10^9),
                         gas_units %in% c("umol_o2") ~ npp*(1/10^6)),
         r_mol = case_when(gas_units %in% c("mg_o2") ~ r*(1/32000),
                             gas_units %in% c("nmol_o2") ~ r*(1/10^9),
                             gas_units %in% c("umol_o2") ~ r*(1/10^6)),
         gpp_mol = case_when(gas_units %in% c("mg_o2") ~ gpp*(1/32000),
                             gas_units %in% c("nmol_o2") ~ gpp*(1/10^9),
                             gas_units %in% c("umol_o2") ~ gpp*(1/10^6)),
         gram = case_when(`mass/area` %in% c("cm^2") ~ 1),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1),
         npp_mol_grammin= npp_mol/(gram*min),
         r_mol_grammin= r_mol/(gram*min),
         gpp_mol_grammin= gpp_mol/(gram*min),
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  select(study_id, inv_T, broad_coral, depth_broad, latitude, species, npp_mol_grammin, gpp_mol_grammin, r_mol_grammin)
  
  
f_s_2 <- read_csv("f+s_coral_npp_gpp_r.csv")%>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         npp_mol = case_when(gas_units %in% c("ug_o2") ~ npp*(1/32000000),
                         gas_units %in% c("nmol_o2") ~ npp*(1/10^9),
                         gas_units %in% c("umol_o2") ~ npp*(1/10^6)),
         r_mol = case_when(gas_units %in% c("ug_o2") ~ r*(1/32000000),
                           gas_units %in% c("nmol_o2") ~ r*(1/10^9),
                           gas_units %in% c("umol_o2") ~ r*(1/10^6)), 
         gpp_mol = case_when(gas_units %in% c("ug_o2") ~ gpp*(1/32000000),
                             gas_units %in% c("nmol_o2") ~ gpp*(1/10^9),
                             gas_units %in% c("umol_o2") ~ gpp*(1/10^6)),
         gram = case_when(`mass/area` %in% c("cm^2") ~ 1, 
                          `mass/area` %in% c("larvae^-1") ~ 1),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1),
         npp_mol_grammin= npp_mol/(gram*min),
         r_mol_grammin= r_mol/(gram*min), 
         gpp_mol_grammin=gpp_mol/(gram*min),
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  select(study_id, inv_T, broad_coral, depth_broad, latitude, species, npp_mol_grammin, gpp_mol_grammin, r_mol_grammin)



f_s_3 <- read_csv("f+s_coral_npp_r.csv") %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         npp_mol = case_when(gas_units %in% c("mg_o2") ~ npp*(1/32000),
                         gas_units %in% c("ug_o2") ~ npp*(1/32000000),
                         gas_units %in% c("umol_o2") ~ npp*(1/10^6), 
                         gas_units %in% c("mmol_o2") ~ (npp*(10^3))*(1/10^3)),
         r_mol = case_when(gas_units %in% c("mg_o2") ~ r*(1/32000),
                             gas_units %in% c("ug_o2") ~ r*(1/32000000),
                             gas_units %in% c("umol_o2") ~ r*(1/10^6), 
                             gas_units %in% c("mmol_o2") ~ r*(1/10^3)),
         gpp_mol = case_when(gas_units %in% c("mg_o2") ~ gpp*(1/32000),
                             gas_units %in% c("ug_o2") ~ gpp*(1/32000000),
                             gas_units %in% c("umol_o2") ~ gpp*(1/10^6), 
                             gas_units %in% c("mmol_o2") ~ gpp*(1/10^3)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("468 mm^2") ~ 0.002136752, 
                          `mass/area` %in% c("l") ~ 1/1000, 
                          `mass/area` %in% c("cm^2") ~ 1, 
                          `mass/area` %in% c("chla") ~ 1,   
                          `mass/area` %in% c("g CaCO3") ~ 1),
         min = case_when(time %in% c("hr") ~ 60,
                         time %in% c("sec") ~ 1/60),
         npp_mol_grammin= npp_mol/(gram*min),
         r_mol_grammin = r_mol/(gram*min),
         gpp_mol_grammin = gpp_mol/(gram*min),
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  select(study_id, inv_T, broad_coral, depth_broad, latitude, species, npp_mol_grammin, gpp_mol_grammin, r_mol_grammin)

f_s_comb <- rbind(f_s_1, f_s_2, f_s_3) 

f_s_comb_long <- f_s_comb %>%
  pivot_longer(cols=c(npp_mol_grammin, gpp_mol_grammin, r_mol_grammin), names_to = "metabolic_category", 
               values_to = "response") %>%
  mutate(centre_temp = I(inv_T -mean(inv_T))) %>%
  mutate(latitude=as.numeric(latitude)) %>%
  mutate(response_add = response+1)
  

f_s_comb_long %>%
  ggplot(aes(x=response))+
  geom_histogram()+
  facet_wrap(~study_id)


hist(f_s_comb_long$response)
max(f_s_comb_long$response)

min(f_s_comb_long$response)

f_s_comb_long %>%
  ggplot()+
  geom_point(aes(y=log(abs(response)), x=centre_temp,  colour=metabolic_category), alpha=0.5)+
  facet_wrap(~metabolic_category)+
  geom_abline(aes(intercept=i, slope=s, colour=metabolic_category), 
            data=data.frame(metabolic_category=c("gpp_mol_grammin","npp_mol_grammin", 
                                                 "r_mol_grammin"), i=c(-15.28187, -15.28187+-0.61731, -15.28187+-0.50402), 
                            s=c(-0.26077, 0.25232+-0.26077, -0.50402+-0.14885)))+
  theme_bw() + 
  theme(aspect.ratio = 1)






model00fs<- lmer(log(response_add) ~centre_temp*metabolic_category + abs(latitude)  + (1+centre_temp|species) , data=f_s_comb_long, REML=FALSE)
model0fs <- lmer(log(response_add) ~centre_temp*metabolic_category + abs(latitude)   + (1+centre_temp|study_id) , data=f_s_comb_long, REML=FALSE)
model1fs <- lmer(log(response_add) ~centre_temp*metabolic_category + abs(latitude)  + (1|study_id) , data=f_s_comb_long, REML=FALSE)
model2fs <- gls(log(response_add) ~centre_temp*metabolic_category + abs(latitude) , data=f_s_comb_long, na.action=na.exclude)

model.sel(model00fs, model0fs, model1fs, model2fs)

#model0fs is best by more than 2 AIC

#missing metabolic category
model3fs_coral <- lmer(log(abs(response)) ~centre_temp + abs(latitude) +depth_broad + (1+centre_temp|study_id), data=f_s_comb_long)
#missing latitude
model4fs_coral <- lmer(log(abs(response)) ~centre_temp + depth_broad + (1+centre_temp|study_id) , data=f_s_comb_long)
#missing both latitude and depth
model5fs_coral <- lmer(log(abs(response)) ~centre_temp  + (1+centre_temp|study_id), data=f_s_comb_long)
#missing temp
model6fs_coral <- lmer(log(abs(response)) ~metabolic_category + depth_broad + abs(latitude)  + (1+centre_temp|study_id) , data=f_s_comb_long)
#missing interaction between temp and metabolic category
model7fs_coral <- lmer(log(abs(response)) ~centre_temp + metabolic_category + abs(latitude) + depth_broad + (1+centre_temp|study_id), data=f_s_comb_long)
#full model
model8fs_coral <- lmer(log(abs(response)) ~centre_temp*metabolic_category + abs(latitude) + depth_broad + (1+centre_temp|study_id), data=f_s_comb_long)

model.sel(model3fs_coral, model4fs_coral, model5fs_coral, model6fs_coral, model7fs_coral, model8fs_coral)

#model8fs is best > 2 AIC 
summary(model8fs_coral)

###
#gpp = -0.260772
#npp = -0.008454
#r = -0.40

expectations<-ggplot() +
  geom_function(fun = ~ 0.3^.x, colour="olivedrab1", size=2) +
  geom_function(fun = ~ -0.6^.x, color = "grey", size=2) +
  geom_function(fun = ~ -0.3^.x, colour="#abbb80", size=2) +
  xlim(5,0)+
  ggtitle("Expectations")+
  xlab("Temperature (1/kT)")+
  ylab("Metabolic rate")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))
  
  coral<- ggplot() +
    geom_function(fun = ~ 0.260772^.x, colour="olivedrab1", size=2) +
    geom_function(fun = ~ -0.40^.x, color = "grey", size=2) +
    geom_function(fun = ~ -0.008454^.x, colour="#abbb80", size=2) +
    xlim(5,0)+
    ggtitle("Coral")+
  xlab("Temperature (1/kT)")+
    ylab("Metabolic rate")+
    theme_bw()+
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          axis.text.x = element_text( size = 20),
          legend.position="none",
          # The new stuff
          strip.text = element_text(size = 20), 
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))
  
  
expectations+ coral




