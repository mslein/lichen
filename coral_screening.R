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
coral <- read_csv("extraction/coral dataset/coral/coral_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol = case_when(gas_units %in% c("log umol_o2") ~ (response_value^10)*(1/10^6),
                         gas_units %in% c("mg_o2") ~ response_value*(1/32000),
                         gas_units %in% c("ug_o2") ~ response_value*(1/32000000),
                         gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                         gas_units %in% c("umol_o2") ~ response_value*(1/10^6),
                         gas_units %in% c("log umol_o2") ~ (response_value^10)*(1/10^6),
                         gas_units %in% c("mmol_o2") ~ response_value*(1/10^3),
                         gas_units %in% c("pg_DIC") ~ response_value*(1/(1.2011*10^13)),
                         gas_units %in% c("uL_o2") ~ response_value*(4*10^-14), 
                         gas_units %in% c("umol_c") ~ response_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("468 mm^2") ~ 10^3/(468), 
                          `mass/area` %in% c("l") ~ 1/1000, 
                          `mass/area` %in% c("cm^2") ~ 1000, 
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
         centre_temp = I(inv_T-mean(inv_T)), 
         centre_celsius = I(temp-mean(temp)), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type))




npp_coral <- coral %>%
  filter(metabolic_category == "npp") %>%
  filter_all(all_vars(!is.infinite(.)))

r_coral <- coral %>%
  filter(metabolic_category == "r") %>%
  filter_all(all_vars(!is.infinite(.)))

gpp_coral <- coral %>%
  filter(metabolic_category == "gpp")


npp_coral$depth_broad<- relevel(npp_coral$depth_broad, "shallow")
npp_coral$cnidarian_type<- relevel(npp_coral$cnidarian_type, "hermatypic coral")
r_coral$depth_broad<- relevel(r_coral$depth_broad, "shallow")
r_coral$cnidarian_type<- relevel(r_coral$cnidarian_type, "hermatypic coral")
gpp_coral$cnidarian_type<- relevel(gpp_coral$cnidarian_type, "hermatypic coral")


#random effects for npp
model00npp_coral<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|species), data=npp_coral, REML=FALSE)
model0npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=npp_coral, REML=FALSE)
model1npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp + (1|study_id) , data=npp_coral, REML=FALSE)
model2npp_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp, data=npp_coral, na.action=na.exclude)
model.sel(model00npp_coral, model0npp_coral, model1npp_coral, model2npp_coral)
#random effects for respiration
model00r_coral<- lmer(log(abs(mol_grammin)) ~abs(latitude) +centre_temp  + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|species) , data=r_coral, REML=FALSE)
model0r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=r_coral, REML=FALSE)
model1r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp + (1|study_id) , data=r_coral, REML=FALSE)
model2r_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  + depth_broad +cnidarian_type*centre_temp, data=r_coral, na.action=na.exclude)
model.sel(model00r_coral, model0r_coral, model1r_coral, model2r_coral)
#random effects for gpp
#singular fit?
#this has a singular fit
model00gpp_coral<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|species) , data=gpp_coral, REML=FALSE)
#this also has a singular fit
model0gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|study_id) , data=gpp_coral, REML=FALSE)
model1gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=gpp_coral, REML=FALSE)
model2gpp_coral <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude), data=gpp_coral, , na.action=na.exclude)
model.sel(model00gpp_coral, model0gpp_coral, model1gpp_coral, model2gpp_coral)

#model0npp == best by 2 AIC
#model0r == best by 2 AIC
#model0gpp == best by 2 AIC ---> has a singular fit though


#model selection for photosynthesis
#model selection for photosynthesis
#missing depth
model3npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp +  (1+centre_temp|study_id), data=npp_coral)
#missing latitude
model4npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=npp_coral)
#missing both latitude and depth
model5npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp +cnidarian_type*centre_temp + (1+centre_temp|study_id), data=npp_coral)
#missing interaction between lichen type + centre temp
model6npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type + (1+centre_temp|study_id), data=npp_coral)
#missing temp
model7npp_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + depth_broad +cnidarian_type + (1+centre_temp|study_id), data=npp_coral)
#missing cnidarian type
model8npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad + (1+centre_temp|study_id), data=npp_coral)
#full model
model9npp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=npp_coral, REML=FALSE)
model.sel(model3npp_coral, model4npp_coral, model5npp_coral, model6npp_coral, model7npp_coral, model8npp_coral, model9npp_coral)
#model 7 is best
summary(model7npp_coral)

final_npp_coral<- model.avg(model6npp_coral, model7npp_coral)
random.effects(model7npp_coral)

#model selection for respiration
#missing depth
model3r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp +  (1+centre_temp|study_id), data=r_coral)
#missing latitude
model4r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=r_coral)
#missing both latitude and depth
model5r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp +cnidarian_type*centre_temp + (1+centre_temp|study_id), data=r_coral)
#missing interaction between lichen type + centre temp
model6r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type + (1+centre_temp|study_id), data=r_coral)
#missing temp
model7r_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + depth_broad +cnidarian_type + (1+centre_temp|study_id), data=r_coral)
#missing cnidarian type
model8r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad + (1+centre_temp|study_id), data=r_coral)
#full model
model9r_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + depth_broad +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=r_coral, REML=FALSE)
model.sel(model3r_coral, model4r_coral, model5r_coral, model6r_coral, model7r_coral, model8r_coral, model9r_coral)

#model 4r is best, but failed to converge, going with model 5r
summary(model5r_coral)
confint(model5r_coral)

#model selection for gpp -->  these are all singular fits
#missing temp
model7gpp_coral <- lmer(log(abs(mol_grammin)) ~abs(latitude) + (1|study_id), data=gpp_coral)
#missing latitude
model8gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp  + (1|study_id), data=gpp_coral)
#full model
model9gpp_coral <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id) , data=gpp_coral, REML=FALSE)
model.sel(model7gpp_coral, model8gpp_coral, model9gpp_coral)

#model 8gpp and 9 are tied, going w/ model 9 
summary(model9gpp_coral)



npp_coral_df<- count(npp_coral, study_id, latitude, cnidarian_type)

randslope_npp_coral<- random.effects(model7npp_coral) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope = centre_temp)

npp_joined_coral<- left_join(randslope_npp_coral,npp_coral_df, by="study_id") %>%
  mutate(c_intercept=case_when(cnidarian_type == "hermatypic coral" & depth_broad == "shallow" ~ -17.74629,
                               cnidarian_type == "hermatypic coral" & depth_broad == "deep" ~ -17.74629+-8.84861,
                               cnidarian_type == "brooding coral" & depth_broad == "deep" ~ -17.74629+-4.02401+-8.84861,
                               cnidarian_type == "brooding coral" & depth_broad == "shallow" ~ -17.74629+-4.02401),
         latitudeii=abs(latitude)*-0.03814, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope) %>%
  as.data.frame()



npp_coral_plot<- ggplot(data=npp_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="#abbb80")+
  geom_point(colour="#abbb80", alpha=0.3, size=4)+
  #scale_x_reverse(limits=c(2, -2.5))+
  geom_abline(data=npp_joined_coral, aes(slope=r_slope, intercept=f_intercept), colour="#abbb80", size=1, alpha=0.5)+
  geom_abline(aes(slope=0, intercept=-18.95686), colour="#99a873", size=2)+
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
  

r_coral_df<- count(r_coral, study_id, cnidarian_type)


randslope_r_coral<- random.effects(model5r_coral) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

r_joined_coral<- left_join(randslope_r_coral,r_coral_df, by="study_id") %>%
  as.data.frame() %>%
  mutate(c_intercept=case_when(cnidarian_type == "hermatypic coral"~ -22.2736,
                               cnidarian_type == "soft coral" ~ -22.2736 +0.5680,
                               cnidarian_type == "brooding coral" ~ -22.2736+-1.0658), 
         f_intercept=r_intercept+c_intercept, 
         c_slope = case_when(cnidarian_type == "hermatypic coral"~ -0.7772,
                             cnidarian_type == "soft coral" ~ -0.7772 + 0.2312,
                             cnidarian_type == "brooding coral" ~ -0.7772+0.5261), 
         f_slope=r_slope+c_slope) %>%
  as.data.frame() 

#double check the math on these
r_coral_plot <- ggplot(data=r_coral, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="grey")+
  geom_point(colour="grey", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(3.5, -2.5))+
  geom_abline(data=r_joined_coral, aes(slope=f_slope*-1, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=0.7772, intercept = -22.2736, size=2.5, colour="grey65")+  
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













###### ISOLATED ALGAE ANALYSIS
isolated_algae <- read_csv("extraction/coral dataset/isolated algae/isolated_algae.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol = case_when(gas_units %in% c("ug_o2") ~ response_value*(1/32000000),
                         gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                         gas_units %in% c("umol_o2") ~ response_value*(1/10^6),
                         gas_units %in% c("mmol_o2") ~ response_value*(1/10^3),
                         gas_units %in% c("fmol_o2") ~ response_value*(1/10^15),
                         gas_units %in% c("pmol_o2") ~ response_value*(1/10^12),
                         gas_units %in% c("pg_c") ~ response_value*(1/(1.2011*10^13)),
                         gas_units %in% c("ug_C") ~ response_value*(1/(1.2011*10^6)),
                         gas_units %in% c("ug_o2") ~ response_value*(3.1999*10^-7)),
         gram = case_when(`mass/area` %in% c("cell") ~ 1, 
                          `mass/area` %in% c("mg chla") ~ 10^3,
                          `mass/area` %in% c("ug chl a") ~ 10^6,
                          `mass/area` %in% c("ug chla") ~ 10^6,
                          `mass/area` %in% c("cell x10^6") ~ 10^6,
                          `mass/area` %in% c("cell x10^9") ~ 10^9,
                          `mass/area` %in% c("cell x 10^9") ~ 10^9),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("day") ~ 1440),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min),
         centre_temp = I(inv_T-mean(inv_T)), 
         centre_celsius = I(temp-mean(temp)), 
         cnidarian_type = as.factor(cnidarian_type))


npp_algae <- isolated_algae %>%
  filter(metabolic_category == "npp") %>%
  filter_all(all_vars(!is.infinite(.)))

r_algae <- isolated_algae %>%
  filter(metabolic_category == "r") %>%
  filter_all(all_vars(!is.infinite(.)))

gpp_algae <- isolated_algae %>%
  filter(metabolic_category == "gpp")



npp_algae$cnidarian_type<- relevel(npp_algae$cnidarian_type, "anemone")
r_algae$cnidarian_type<- relevel(r_algae$cnidarian_type, "anemone")
gpp_algae$cnidarian_type<- relevel(gpp_algae$cnidarian_type, "hermatypic coral")


#random effects for npp -- these are all singular fits except for model 1, 2
model00npp_algae<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  +cnidarian_type*centre_temp + (1+centre_temp|species), data=npp_algae, REML=FALSE)
model0npp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=npp_algae, REML=FALSE)
model1npp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)   +cnidarian_type*centre_temp + (1|study_id) , data=npp_algae, REML=FALSE)
model2npp_algae <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  +cnidarian_type*centre_temp, data=npp_algae, na.action=na.exclude)
model.sel(model00npp_algae, model0npp_algae, model1npp_algae, model2npp_algae)
#random effects for respiration
#this one is a singular fit
model00r_algae<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  +cnidarian_type*centre_temp + (1+centre_temp|species), data=r_algae, REML=FALSE)
model0r_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=r_algae, REML=FALSE)
model1r_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)   +cnidarian_type*centre_temp + (1|study_id) , data=r_algae, REML=FALSE)
model2r_algae <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  +cnidarian_type*centre_temp, data=r_algae, na.action=na.exclude)
model.sel(model00r_algae, model0r_algae, model1r_algae, model2r_algae)
#model 0 is best


#random effects for gpp
#these are all singular fits except for model 1,2
model00gpp_algae<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  +cnidarian_type*centre_temp + (1+centre_temp|species), data=gpp_algae, REML=FALSE)
model0gpp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=gpp_algae, REML=FALSE)
model1gpp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude)   +cnidarian_type*centre_temp + (1|study_id) , data=gpp_algae, REML=FALSE)
model2gpp_algae <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude)  +cnidarian_type*centre_temp, data=gpp_algae, na.action=na.exclude)
model.sel(model00gpp_algae, model0gpp_algae, model1gpp_algae, model2gpp_algae)

#model selection for npp
#missing latitude
model4npp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp +cnidarian_type*centre_temp + (1|study_id) , data=npp_algae)
#missing interaction between lichen type + centre temp
model6npp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type + (1|study_id), data=npp_algae)
#missing temp
model7npp_algae <- lmer(log(abs(mol_grammin)) ~abs(latitude) +cnidarian_type + (1|study_id), data=npp_algae)
#missing cnidarian type
model8npp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id), data=npp_algae)
#full model
model9npp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp + (1|study_id) , data=npp_algae, REML=FALSE)
model.sel(model4npp_algae, model6npp_algae, model7npp_algae, model8npp_algae, model9npp_algae)

#model9npp_algae is best

#model selection for respiration
#missing latitude
model4r_algae <- lmer(log(abs(mol_grammin)) ~centre_temp +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=r_algae)
#missing interaction between lichen type + centre temp
model6r_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type + (1+centre_temp|study_id), data=r_algae)
#missing temp
model7r_algae <- lmer(log(abs(mol_grammin)) ~abs(latitude) +cnidarian_type + (1+centre_temp|study_id), data=r_algae)
#missing cnidarian type
model8r_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1+centre_temp|study_id), data=r_algae)
#full model
model9r_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp + (1+centre_temp|study_id) , data=r_algae, REML=FALSE)
model.sel(model4r_algae, model6r_algae, model7r_algae, model8r_algae, model9r_algae)

#model9r_algae is best


#model selection for gpp
#these are all singular fits
#missing latitude
model4gpp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp +cnidarian_type*centre_temp + (1|study_id) , data=gpp_algae)
#missing interaction between lichen type + centre temp
model6gpp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type + (1|study_id), data=gpp_algae)
#missing temp
model7gpp_algae <- lmer(log(abs(mol_grammin)) ~abs(latitude) +cnidarian_type + (1|study_id), data=gpp_algae)
#missing cnidarian type
model8gpp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) + (1|study_id), data=gpp_algae)
#full model
model9gpp_algae <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +cnidarian_type*centre_temp + (1|study_id) , data=gpp_algae, REML=FALSE)
model.sel(model4gpp_algae, model6gpp_algae, model7gpp_algae, model8gpp_algae, model9gpp_algae)


#model 6 and 4 are tied
model.avg(model4gpp_algae, model6gpp_algae)





npp_algae_df<- count(npp_algae, study_id, latitude, cnidarian_type)

randslope_npp_algae<- random.effects(model9npp_algae) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp)

npp_joined_algae<- left_join(randslope_npp_algae, npp_algae_df, by="study_id") %>%
  mutate(c_intercept=case_when(cnidarian_type == "anemone" ~ -28.379141,
                               cnidarian_type == "bivalve" ~ -1.141867+-28.379141,
                               cnidarian_type == "free-living" ~ 0.081163+-28.379141,
                               cnidarian_type == "hermatypic coral"~ -0.273064+-28.379141,
                               cnidarian_type == "jellyfish"  ~ -1.180844+-28.379141, 
                               cnidarian_type == "octocoral"  ~ -0.160391 +-28.379141),
         latitudeii=abs(latitude)*-0.007744, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope= case_when(cnidarian_type == "anemone" ~ 0.038349,
                            cnidarian_type == "bivalve" ~ 0.038349+1.642649,
                            cnidarian_type == "free-living" ~ 0.038349+-0.101801,
                            cnidarian_type == "hermatypic coral"~ 0.038349+-0.086094,
                            cnidarian_type == "jellyfish"  ~ 0.038349+0.630880, 
                            cnidarian_type == "octocoral"  ~ 0.038349 +-0.460304)) %>%
  as.data.frame()





npp_algae_plot<- ggplot(data=npp_algae, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="#abbb80")+
  geom_point(colour="#abbb80", alpha=0.3, size=4)+
  #scale_x_reverse(limits=c(2, -2.5))+
  geom_abline(data=npp_joined_algae, aes(slope=f_slope*-1, intercept=f_intercept), colour="#abbb80", size=1, alpha=0.5)+
  geom_abline(aes(slope=-0.038349, intercept=-28.379141), colour="#99a873", size=2)+
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


r_algae_df<- count(r_algae, study_id, latitude, cnidarian_type)


randslope_r_algae<- random.effects(model9r_algae) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)


r_joined_algae<- left_join(randslope_r_algae, r_algae_df, by="study_id") %>%
  mutate(c_intercept=case_when(cnidarian_type == "anemone" ~ -30.626787,
                               cnidarian_type == "bivalve" ~ 1.763309+-30.626787,
                               cnidarian_type == "free-living" ~ 0.011017+-30.626787,
                               cnidarian_type == "hermatypic coral"~ -0.370946+-30.626787,
                               cnidarian_type == "jellyfish"  ~ -0.846593+-30.626787, 
                               cnidarian_type == "octocoral"  ~ -0.483789  +-30.626787),
         latitudeii=abs(latitude)*-0.002289, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         c_slope= case_when(cnidarian_type == "anemone" ~ -0.041305,
                            cnidarian_type == "bivalve" ~ -0.041305+0.420322,
                            cnidarian_type == "free-living" ~ -0.041305+-0.080356,
                            cnidarian_type == "hermatypic coral"~ -0.041305+0.021492,
                            cnidarian_type == "jellyfish"  ~ -0.041305+0.309925, 
                            cnidarian_type == "octocoral"  ~ -0.041305 +0.019730), 
         f_slope=r_slope+c_slope) %>%
  as.data.frame()





#double check the math on these
r_algae_plot <- ggplot(data=r_algae, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="grey")+
  geom_point(colour="grey", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(3.5, -2.5))+
  geom_abline(data=r_joined_algae, aes(slope=f_slope*-1, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=0.041305, intercept = -30.626787, size=2.5, colour="grey65")+  
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




full_plot <- (p_plot + r_plot) / (npp_coral_plot + r_coral_plot) / (npp_algae_plot + r_algae_plot)

ggsave(filename = "fullplot.png", plot=full_plot, width = 10, height = 12, dpi = 300)



