install.packages("revtools")
library(revtools)
library(tidyverse)
pacman::p_load(viridis, revtools, nlme, lme4, MuMIn, patchwork, tidyverse)


#######screening abstracts
#lichen1 <- read_bibliography("lichen24sept21.bib")
#lichen2 <- read_bibliography("wos155.bib")
#screen_abstracts(lichen2)

##ask mary about the difference in units between the spreadsheet 
## and the plos one figure 1....


#i actually think i forgot to differentiate the molar mass shit fuck between co2 and o2

library(tidyverse)
lichen <- read_csv("extraction/lichen dataset/lichen_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol = case_when(gas_units %in% c("mg_co2") ~ response_value*(1/44010), 
                         gas_units %in% c("mg_o2") ~ response_value*(1/32000), 
                         gas_units %in% c("nmol_co2") ~ response_value*(1/10^9), 
                         gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                         gas_units %in% c("umol_co2") ~ response_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("kg") ~ 1000, 
                          `mass/area` %in% c("m") ~ 1000, 
                          `mass/area` %in% c("mg") ~ 1/1000, 
                          `mass/area` %in% c("mg chla") ~ 1/1000, 
                          `mass/area` %in% c("mg chl") ~ 1/1000),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min),
         latitude_broad=as.factor(latitude_broad), 
         elevation_broad=as.factor(elevation_broad), 
         lichen_type=as.factor(lichen_type), 
         centre_temp = I(inv_T-mean(inv_T)), 
         centre_celsius = I(temp-mean(temp)))



lichen_max_ppfd <- lichen %>%
  group_by(study_id) %>%
  filter(avg_ppfd == max(avg_ppfd)) %>%
  ungroup()

npp <- lichen_max_ppfd %>%
  filter(metabolic_category == "npp") %>%
  filter_all(all_vars(!is.infinite(.)))

r <- lichen %>%
  filter(metabolic_category == "r") %>%
  filter_all(all_vars(!is.infinite(.)))

gpp <- lichen_max_ppfd %>%
  filter(metabolic_category == "gpp")




#subset_lichen$broad_responses<- relevel(subset_lichen$broad_responses, "respiration")
npp$latitude_broad<- relevel(npp$latitude_broad, "temperate")
npp$elevation_broad<- relevel(npp$elevation_broad, "neutral")
npp$lichen_type<- relevel(npp$lichen_type, "green algae")
r$latitude_broad<- relevel(r$latitude_broad, "tropical")
r$elevation_broad<- relevel(r$elevation_broad, "neutral")
r$lichen_type<- relevel(r$lichen_type, "green algae")
gpp$latitude_broad<- relevel(gpp$latitude_broad, "temperate")

#random effects for npp
#singular fit?
model00npp<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad  + (1+centre_temp|species), data=npp, REML=FALSE)
model0npp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad  + (1+centre_temp|study_id) , data=npp, REML=FALSE)
model1npp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad  + (1|study_id) , data=npp, REML=FALSE)
model2npp <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad, data=npp, na.action=na.exclude)
#random effects for respiration
model00r<- lmer(log(abs(mol_grammin)) ~abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|species) , data=r, REML=FALSE)
model0r <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=r, REML=FALSE)
model1r <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|study_id) , data=r, REML=FALSE)
model2r <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad, data=r, na.action=na.exclude)
#random effects for gpp
#singular fit?
model00gpp<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|species) , data=gpp, REML=FALSE)
model0gpp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + (1+centre_temp|study_id) , data=gpp, REML=FALSE)
model1gpp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=gpp, REML=FALSE)
model2gpp <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude) , data=gpp, na.action=na.exclude)

#model0npp == best by 2 AIC
#model0r == best by 2 AIC
#model1gpp <- different random effects structure


#model selection -- random effects --> all the data togetehr (not relevant anymore)
#model00 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|centre_temp) , data=subset_lichen, REML=FALSE)
#model0 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=subset_lichen, REML=FALSE)
#model1 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|centre_temp) + (1 |study_id) , data=subset_lichen, REML=FALSE)
#model1 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|study_id) , data=subset_lichen, REML=FALSE)
#model2 <- gls(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad, data=subset_lichen)

#model selection for npp random effects
model.sel(model00npp, model0npp, model1npp, model2npp)

#model selection for respiration random effects
model.sel(model00r, model0r, model1r, model2r)

#model selection for gpp random effects
#model.sel(model1gpp, model2gpp)

#model0 prevails for npp and r 
#model 1 prevails for gpp


#model selection for photosynthesis
#missing elevation
model3npp <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp +abs(latitude)+  (1+centre_temp|study_id), data=npp)
#missing latitude
model4npp <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=npp)
#missing both latitude and elevation
model5npp <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=npp)
#missing interaction between lichen type + centre temp
model6npp <- lmer(log(abs(mol_grammin)) ~lichen_type +centre_temp + abs(latitude)+ elevation_broad + (1+centre_temp|study_id), data=npp)
#missing temp
model7npp <- lmer(log(abs(mol_grammin)) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=npp)
#missing lichen type
model8npp <- lmer(log(abs(mol_grammin)) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=npp)
#full model
model9npp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=npp)
model.sel(model3npp, model4npp, model5npp, model6npp, model7npp, model8npp, model9npp)

#model3npp, model7npp, model9npp are tied


#model selection for respiration
#missing elevation
model3r <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=r)
#missing latitude
model4r <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=r)
#missing both latitude and elevation
model5r <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=r)
#missing temp
model6r <- lmer(log(abs(mol_grammin)) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=r)
#missing lichen type
model7r <- lmer(log(abs(mol_grammin)) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=r)
#missing interaction between lichen type + centre temp
model8r <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=r)
model9r <- lmer(log(abs(mol_grammin))~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=r)
model.sel(model3r, model4r, model5r, model6r, model7r, model8r, model9r)

#model8r is best by < 2 AIC


summary(model8r)
confint(model8r)


summary(model7npp)
confint(model7npp)




npp_df<- count(npp, study_id, lichen_type, elevation_broad, latitude)


randslope_p<- random.effects(model7npp) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

photosynthesis_joined<- left_join(randslope_p,npp_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "green algae" & elevation_broad == "neutral" ~-11.933541,
                               lichen_type == "green algae" & elevation_broad == "high" ~-11.933541+0.028084, 
                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~-11.933541+0.683963,
                               lichen_type == "cyanobacteria" & elevation_broad == "high" ~-11.933541+0.683963+0.028084, 
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~ -11.933541+1.222073, 
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~ -11.933541+1.222073+0.028084),
         latitudeii=abs(latitude)*-0.052765, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope) %>%
  as.data.frame()

p_plot<- ggplot(data=npp, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="#abbb80")+
  geom_point(colour="#abbb80", alpha=0.3, size=4)+
  #scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=photosynthesis_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="#abbb80", size=1, alpha=0.5)+
  geom_abline(aes(slope=0, intercept=-11.933541), colour="#99a873", size=2)+
  #geom_hline(yintercept=-9.1055, size=1.5)+  
  xlab("Temperature (1/kT)")+
  #ylab("Metabolic rate (log(mmol CO2 per mg per min))")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("NPP")+
  xlab("Temperature (1/kT)")
  #ylim(-31, -4)

respiration_df<- count(r, study_id, lichen_type, elevation_broad, latitude, centre_temp)


randslope_r<- random.effects(model8r) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

respiration_joined<- left_join(randslope_r,respiration_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "cyanobacteria" & elevation_broad == "high" ~0.708692+-0.108135+-15.304981,
                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~0.708692+-15.304981, 
                               lichen_type == "green algae" & elevation_broad == "high" ~-15.304981+-0.108135,
                               lichen_type == "green algae" & elevation_broad == "neutral" ~-15.304981,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~0.379030+-15.304981+-0.108135,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~0.379030+-15.304981),
         latitudeii=abs(latitude)*0.008146, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope+-0.584779) %>%
  as.data.frame()

#double check the math on these
r_plot <- ggplot(data=r, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="grey")+
  geom_point(colour="grey", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=respiration_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=0.6018, intercept = -15.34, size=2.5, colour="grey65")+  
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
  #scale_x_continuous(sec.axis = sec_axis(~.*(1/(8.61*10^-5))-273.15), name = "Temperature (°C)")
  #scale_x_continuous(sec.axis = sec_axis(~ . * 20, name = "Temperature (°C)"))

#how do i plot this if this is the centred term?


ggplot(data=r, aes(x=temp, y=log(abs(mol_grammin))), colour="grey")+
  geom_point(colour="grey", alpha=0.3, size=4)+
  geom_abline(data=respiration_joined, aes(slope=f_slope, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=-0.6018, intercept = -15.34, size=2.5, colour="grey65")+  
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
  ggtitle("Respiration")




p_plot + r_plot


### making the friedman and sun figure

  
ggplot() +
  geom_function(fun = ~ -0.3^.x, colour="olivedrab1", size=2) +
  geom_function(fun = ~ 0.6^.x, color = "grey", size=2) +
  geom_function(fun = ~ 0.3^.x, colour="#abbb80", size=2) +
  xlim(2.5,-2.5)+
  xlab("Temperature (1/kT)")+
  ylab("Metabolic rate (mmol CO2 per mg per min)")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))


ggplot() +
  geom_function(fun = ~ .x^3/4, colour="black", size=2)+
  xlim(0,20)
  
ggplot() +
  geom_function(fun = ~ .x^-1/4, color = "grey", size=2) +
  #geom_function(fun = ~ 0.3^.x, colour="#abbb80", size=2) +
  xlim(0,20)+
  xlab("Temperature (1/kT)")+
  ylab("Metabolic rate (mmol CO2 per mg per min)")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))



# trying to back calculate GPP 

gpp_arithmetic <- read_csv("extraction/lichengpp_arithmetic3oct2023.csv") %>%
  mutate(r_value = r_value*-1,
    gpp = npp_value + r_value) %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol_gpp = case_when(gas_units %in% c("mg_co2") ~ gpp*(1/44010), 
                         gas_units %in% c("mg_o2") ~ gpp*(1/32000), 
                         gas_units %in% c("nmol_co2") ~ gpp*(1/10^9), 
                         gas_units %in% c("nmol_o2") ~ gpp*(1/10^9),
                         gas_units %in% c("umol_co2") ~ gpp*(1/10^6)),
         mol_r = case_when(gas_units %in% c("mg_co2") ~ r_value*(1/44010), 
                           gas_units %in% c("mg_o2") ~ r_value*(1/32000), 
                           gas_units %in% c("nmol_co2") ~ r_value*(1/10^9), 
                           gas_units %in% c("nmol_o2") ~ r_value*(1/10^9),
                           gas_units %in% c("umol_co2") ~ r_value*(1/10^6)),
         mol_npp = case_when(gas_units %in% c("mg_co2") ~ npp_value*(1/44010), 
                             gas_units %in% c("mg_o2") ~ npp_value*(1/32000), 
                             gas_units %in% c("nmol_co2") ~ npp_value*(1/10^9), 
                             gas_units %in% c("nmol_o2") ~ npp_value*(1/10^9),
                             gas_units %in% c("umol_co2") ~ npp_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("kg") ~ 1000, 
                          `mass/area` %in% c("m") ~ 1000, 
                          `mass/area` %in% c("mg") ~ 1/1000, 
                          `mass/area` %in% c("mg chla") ~ 1/1000),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60),
         gpp_mol_grammin= mol_gpp/(gram*min),
         r_mol_grammin = mol_r/(gram*min),
         npp_mol_grammin= mol_npp/(gram*min),
         latitude_broad=as.factor(latitude_broad), 
         elevation_broad=as.factor(elevation_broad), 
         lichen_type=as.factor(lichen_type), 
         centre_temp = I(inv_T-mean(inv_T)))

gpp_arithmetic$latitude_broad<- relevel(gpp_arithmetic$latitude_broad, "temperate")
gpp_arithmetic$elevation_broad<- relevel(gpp_arithmetic$elevation_broad, "neutral")
gpp_arithmetic$lichen_type<- relevel(gpp_arithmetic$lichen_type, "green algae")


model3gpp <- lmer(log(abs(gpp_mol_grammin)) ~centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=gpp_arithmetic)
#missing latitude
model4gpp <- lmer(log(abs(gpp_mol_grammin)) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=gpp_arithmetic)
#missing both latitude and elevation
model5gpp <- lmer(log(abs(gpp_mol_grammin)) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing temp
model6gpp <- lmer(log(abs(gpp_mol_grammin)) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing lichen type
model7gpp <- lmer(log(abs(gpp_mol_grammin)) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing interaction between lichen type + centre temp
model8gpp <- lmer(log(abs(gpp_mol_grammin)) ~centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
model9gpp <- lmer(log(abs(gpp_mol_grammin))~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
model.sel(model3gpp, model4gpp, model5gpp, model6gpp, model7gpp, model8gpp, model9gpp)


#model6gpp is best 
summary(model6gpp)


model3npp2 <- lmer(log(abs(npp_mol_grammin)) ~centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=gpp_arithmetic)
#missing latitude
model4npp2 <- lmer(log(abs(npp_mol_grammin)) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=gpp_arithmetic)
#missing both latitude and elevation
model5npp2 <- lmer(log(abs(npp_mol_grammin)) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing temp
model6npp2 <- lmer(log(abs(npp_mol_grammin)) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing lichen type
model7npp2 <- lmer(log(abs(npp_mol_grammin)) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing interaction between lichen type + centre temp
model8npp2 <- lmer(log(abs(npp_mol_grammin)) ~centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
model9npp2 <- lmer(log(abs(npp_mol_grammin))~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
model.sel(model3npp2, model4npp2, model5npp2, model6npp2, model7npp2, model8npp2, model9npp2)



#model8npp2 is best model w temperature
summary(model6npp2)


model3r2 <- lmer(log(abs(r_mol_grammin)) ~centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=gpp_arithmetic)
#missing latitude
model4r2 <- lmer(log(abs(r_mol_grammin)) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=gpp_arithmetic)
#missing both latitude and elevation
model5r2 <- lmer(log(abs(r_mol_grammin)) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing temp
model6r2 <- lmer(log(abs(r_mol_grammin)) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing lichen type
model7r2 <- lmer(log(abs(r_mol_grammin)) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
#missing interaction between lichen type + centre temp
model8r2 <- lmer(log(abs(r_mol_grammin)) ~centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
model9r2 <- lmer(log(abs(r_mol_grammin))~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_arithmetic)
model.sel(model3r2, model4r2, model5r2, model6r2, model7r2, model8r2, model9r2)


#model4r2 is best
summary(model4r2)

g_df<- count(gpp_arithmetic, study_id, lichen_type, elevation_broad, latitude)


randslope_g<- random.effects(model6gpp) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

g_joined<- left_join(randslope_g,g_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "cyanobacteria" & elevation_broad == "high" ~3.007445+0.560869+-14.278130,
                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~3.007445+-14.278130, 
                               lichen_type == "green algae" & elevation_broad == "high" ~-14.278130+0.560869,
                               lichen_type == "green algae" & elevation_broad == "neutral" ~-14.278130,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~1.272075+-14.278130+0.560869,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~1.272075+-14.278130),
         latitudeii=abs(latitude)*-0.032704, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope) %>%
  as.data.frame()





gpp_plot<-ggplot(data=gpp_arithmetic, aes(x=centre_temp, y=log(abs(gpp_mol_grammin))), colour="olivedrab1")+
  geom_point(colour="olivedrab1", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=g_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="olivedrab1", size=1, alpha=0.5)+
  geom_abline(slope=0.174060, intercept = -13.772024, size=1.5, colour="olivedrab2")+ 
  xlab("Temperature (1/kT)")+
  ylab("Metabolic rate (log(mmol CO2 / mg / min))")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("GPP")





randslope_n<- random.effects(model6npp2) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

n_joined<- left_join(randslope_n,g_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "cyanobacteria" & elevation_broad == "high" ~0.599995+0.754487+-11.772943,
                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~0.599995+-11.772943, 
                               lichen_type == "green algae" & elevation_broad == "high" ~-11.772943+0.754487,
                               lichen_type == "green algae" & elevation_broad == "neutral" ~-11.772943,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~1.206334+-11.772943+0.754487,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~1.206334+-11.772943),
         latitudeii=abs(latitude)*-0.058892, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope+-0.049502) %>%
  as.data.frame()


npp_plot <- ggplot(data=gpp_arithmetic, aes(x=centre_temp, y=log(abs(npp_mol_grammin))), colour="#abbb80")+
  geom_point(colour="#abbb80", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=n_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="#abbb80", size=1, alpha=0.5)+
  geom_abline(slope=0.049502, intercept = -11.772943, size=1.5, colour="#99a873")+ 
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
  ggtitle("NPP")



randslope_r<- random.effects(model4r2) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

r_joined<- left_join(randslope_r,g_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "cyanobacteria" & elevation_broad == "high" ~4.06343+0.50715+-16.32424,
                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~4.06343+-16.32424, 
                               lichen_type == "green algae" & elevation_broad == "high" ~-16.32424+0.50715,
                               lichen_type == "green algae" & elevation_broad == "neutral" ~-16.32424,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~1.32314+-16.32424+0.50715,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~1.32314+-16.32424),
         f_intercept=r_intercept+c_intercept, 
         p_slope= case_when(lichen_type== "cyanobacteria" ~ -0.10089+-0.52651, 
                            lichen_type== "cyanobacteria + green algae" ~ -0.36705+-0.52651, 
                            lichen_type== "green algae" ~ -0.52651),
         f_slope=r_slope+p_slope) %>%
  as.data.frame()

resp_plot<- ggplot(data=gpp_arithmetic, aes(x=centre_temp, y=log(abs(r_mol_grammin))), colour="grey")+
  geom_point(colour="grey", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=r_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=0.52651, intercept = -16.32424, size=1.5, colour="grey65")+ 
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




gpp_plot + npp_plot + resp_plot













