install.packages("revtools")
library(revtools)
library(tidyverse)
pacman::p_load(viridis, revtools, nlme, lme4, MuMIn, patchwork)

# 18 jan 2024
# ALL FIELDS = ("temperature" AND ("symbio*" OR "holobiont") AND ("respiration" OR "photosynthesis")) AND "coral"
# 384 papers ----> coral.bib

coral1 <- read_bibliography("lit search files/coral.bib")
#lichen2 <- read_bibliography("wos155.bib")
screen_abstracts(coral1)

#72 candidate studies, 24 only have 2 temps
#17 isolated algae studies, 4 have 2 temps		



#29 jan 2024
# ABSTRACT = "zooxan*" AND "temperature" AND ("photosynthesis" OR "respiration")
#65 papers -----> zoox.bib
coral2 <- read_bibliography("lit search files/zoox.bib")
#lichen2 <- read_bibliography("wos155.bib")
screen_abstracts(coral2)

#### mixotroph search 
# 18 jan 
#ALL FIELDS = "mixotroph*" AND "photosynthesis" AND "temperature"
# 71 papers ----> mixotroph.bib
mixo1 <- read_bibliography("lit search files/mixotroph.bib")
#lichen2 <- read_bibliography("wos155.bib")
screen_abstracts(mixo1)



#### symbiodinium search
# 13 march 
# ABSTRACT = "symbiodinium" AND "temperature" and "photosynthesis"
# 41 papers ----> symbiodinium.bib
symbio1 <- read_bibliography("lit search files/symbiodinium.bib")
#lichen2 <- read_bibliography("wos155.bib")
screen_abstracts(symbio1)



########## 16 April 2024

#corals
coral <- read_csv("extraction/coral dataset/coral_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol = case_when(gas_units %in% c("mg_o2") ~ response_value*(1/32000),
                         gas_units %in% c("ug_o2") ~ response_value*(1/32000000),
                         gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                         gas_units %in% c("umol_o2") ~ response_value*(1/10^6),
                         gas_units %in% c("log umol_o2") ~ (response_value^10)*(1/10^6),
                         gas_units %in% c("mmol_o2") ~ response_value*(1/10^3)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("468 mm^2") ~ 10^6/(468), 
                          `mass/area` %in% c("l") ~ 1/1000, 
                          `mass/area` %in% c("cm^2") ~ 1000, 
                          `mass/area` %in% c("chla") ~ 1, 
                          `mass/area` %in% c("larvae^-1") ~ 1),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min),
         centre_temp = I(inv_T-mean(inv_T)), 
         centre_celsius = I(temp-mean(temp)))

hist(log(coral$mol_grammin))



npp_coral <- coral %>%
  filter(metabolic_category == "npp") %>%
  filter_all(all_vars(!is.infinite(.)))

r_coral <- coral %>%
  filter(metabolic_category == "r") %>%
  filter_all(all_vars(!is.infinite(.)))

gpp_coral <- coral %>%
  filter(metabolic_category == "gpp")



#random effects for npp
#singular fit?
model00npp_coral<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|species), data=npp_coral, REML=FALSE)
model0npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|study_id) , data=npp_coral, REML=FALSE)
model1npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + (1|study_id) , data=npp_coral, REML=FALSE)
model2npp_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude), data=npp_coral)

model.sel(model00npp_coral, model0npp_coral, model1npp_coral, model2npp_coral)
#random effects for respiration
#this has a singular fit
model00r_coral<- lmer(log(abs(mol_grammin)) ~abs(latitude) +centre_temp + (1+centre_temp|species) , data=r_coral, REML=FALSE)
#this also has a singular fit 
model0r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|study_id) , data=r_coral, REML=FALSE)
model1r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=r_coral, REML=FALSE)
model2r_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude), data=r_coral)
model.sel(model00r_coral, model0r_coral, model1r_coral, model2r_coral)
#random effects for gpp
#singular fit?
#this has a singular fit
model00gpp_coral<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|species) , data=gpp_coral, REML=FALSE)
#this also has a singular fit
model0gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + (1+centre_temp|study_id) , data=gpp_coral, REML=FALSE)
model1gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=gpp_coral, REML=FALSE)
model2gpp_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude) , data=gpp_coral)
model.sel(model00gpp_coral, model0gpp_coral, model1gpp_coral, model2gpp_coral)

#model0npp and model00npp are tied, going with model 0npp to standardize
#model0r == best by 2 AIC
#model0gpp == best by 2 AIC

#model selection for photosynthesis
#missing latitude
model4npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + (1|study_id) , data=npp_coral)
#missing temp
model7npp_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + (1|study_id), data=npp_coral)
#full model
model9npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=npp_coral)
model.sel(model4npp_coral, model7npp_coral, model9npp_coral)
#model 7 is best

#model selection for respiration
#missing latitude
model4r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + (1|study_id) , data=r_coral)
#missing temp
model7r_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + (1|study_id), data=r_coral)
#full model
model9r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=r_coral)
model.sel(model4r_coral, model7r_coral, model9r_coral)

#model 4r is best
summary(model4r_coral)
confint(model4r)

#model selection for gpp --> these all get singular fits
#missing latitude
model4gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + (1|study_id) , data=gpp_coral)
#missing temp
model7gpp_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + (1|study_id), data=gpp_coral)
#full model --> gets a singular fit
model9gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=gpp_coral)
model.sel(model4gpp_coral, model7gpp_coral, model9gpp_coral)

#model 4gpp is best
summary(model4gpp)

summary(model7npp_coral)
confint(model7npp_coral)

summary(model4r_coral)
confint(model4r_coral)

summary(model4gpp_coral)
confint(model4gpp_coral)

npp_coral_df<- count(npp, study_id, latitude, centre_temp)

randslope_p_coral<- random.effects(model7npp_coral) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp)

photosynthesis_joined_coral<- left_join(randslope_p_coral,npp_coral_df, by="study_id") %>%
  mutate(latitudeii=abs(latitude)*-0.03817, 
         f_intercept=r_intercept+latitudeii) %>%
  as.data.frame()

npp_coral_plot<- ggplot(data=npp_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="#abbb80")+
  geom_point(colour="#abbb80", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(2, -2.5))+
  #geom_abline(data=photosynthesis_joined_coral, aes(slope=0, intercept=f_intercept), colour="#abbb80", size=1, alpha=0.5)+
  geom_abline(aes(slope=0, intercept=-19.89176), colour="#99a873", size=2)+
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
  xlab("Temperature (1/kT)")
  #add_phylopic(uuid = "32eeaa29-2e90-40b5-92fb-eb8d7fcab9ab", x=2.5, y=-25, ysize=1.2, alpha=1,fill = "#abbb80")
  

respiration_df<- count(r, study_id, lichen_type, elevation_broad, latitude, centre_temp)


randslope_r<- random.effects(model8r) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

#respiration_joined<- left_join(randslope_r,respiration_df, by="study_id") %>%
#  mutate(c_intercept=case_when(lichen_type == "cyanobacteria" & elevation_broad == "high" ~0.7043+-0.1195+-15.34,
#                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~0.7043+-15.34, 
#                               lichen_type == "green algae" & elevation_broad == "high" ~-15.34+-0.1195,
#                               lichen_type == "green algae" & elevation_broad == "neutral" ~-15.34,
#                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~0.3837+-15.34+-0.1195,
#                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~0.3837+-15.34),
#         latitudeii=abs(latitude)*0.008409, 
#         f_intercept=r_intercept+c_intercept+latitudeii, 
#         f_slope=r_slope+-0.60) %>%
#  as.data.frame()

#double check the math on these
r_coral_plot <- ggplot(data=r_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="grey")+
  geom_point(colour="grey", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(2, -2.5))+
  #geom_abline(data=respiration_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=1.3082, intercept = -20.9268, size=2.5, colour="grey65")+  
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
  ggtitle("R")
  #add_phylopic(uuid = "32eeaa29-2e90-40b5-92fb-eb8d7fcab9ab", x=2.5, y=-40, ysize=1.2, alpha=1,fill = "grey")
#scale_x_continuous(sec.axis = sec_axis(~.*(1/(8.61*10^-5))-273.15), name = "Temperature (°C)")
#scale_x_continuous(sec.axis = sec_axis(~ . * 20, name = "Temperature (°C)"))

#how do i plot this if this is the centred term?


gpp_coral_plot <- ggplot(data=gpp_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="grey")+
  geom_point(colour="olivedrab1", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(2, -2.5))+
  #geom_abline(data=respiration_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=-1.1397, intercept = -24.2411, size=2.5, colour="olivedrab2")+  
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
  ylab("Metabolic rate (log(mmol O2 /mg / min))")
  #add_phylopic(uuid = "32eeaa29-2e90-40b5-92fb-eb8d7fcab9ab", x=2.5, y=-50, ysize=1.2, alpha=1,fill = "olivedrab2")



gpp_coral_plot + npp_coral_plot + r_coral_plot


install.packages("remotes")
remotes::install_github("palaeoverse/rphylopic")
library(rphylopic)
#code to get the uuid's for the phylopic pngs
uuid <- rphylopic::get_uuid(name = "Acropora palmata")





























