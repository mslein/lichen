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
  filter_all(all_vars(!is.infinite(.))) %>%
  filter(case_when(study_id %in% c("adams_1971", "delprado_1995", "domaschke_2013", 
                                   "eickmeieradams1973", "kershaw1983", 
                                   "macfarlane1983", "aubert_2007", 
                                   "sanchokappen1989") ~ inv_T < 41, 
                   study_id %in% c("coxson_2003", "reiter2000", 
                                   "sancho1997",
                                   "kappen_2000") ~ inv_T < 42, 
                   study_id %in% c("brownkershaw1984", "lange2000", 
                                   "zotz2003", "lechowicz1973", "lechowiczadams1974", 
                                   "sonesson1989", "sundberg1997", "tretiach1997") ~ inv_T < 40,
                   study_id %in% c() ~ inv_T < 43,
                   study_id %in% c("kappen1983") ~ inv_T < 41.5,
                   study_id %in% c("schroeter1995") ~ inv_T < 42.5,
                   study_id %in% c("zotz1998", "lange2004") ~ inv_T < 39.5,
                   TRUE ~ inv_T < 50))
r <- lichen %>%
  filter(metabolic_category == "r") %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  filter(case_when(study_id %in% c("castillohelmuth2005") ~ inv_T < 38, 
                   study_id %in% c("dorey2020") ~ inv_T < 41, 
                   study_id %in% c("kemp2011") ~ inv_T < 40.5,
                   TRUE ~ inv_T < 50))

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
model00npp<- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad  + (1+centre_temp|species), data=npp, REML=FALSE)
model0npp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad  + (1+centre_temp|study_id) , data=npp, REML=FALSE)
model1npp <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad  + (1|study_id) , data=npp, REML=FALSE)
model2npp <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad, data=npp, na.action=na.exclude)
#random effects for respiration
model00r<- lmer(log(abs(mol_grammin)) ~abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|species) , data=r, REML=FALSE)
model0r <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=r, REML=FALSE)
model1r <- lmer(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|study_id) , data=r, REML=FALSE)
model2r <- gls(log(abs(mol_grammin)) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad, data=r, na.action=na.exclude)

#model selection for npp random effects
model.sel(model00npp, model0npp, model1npp, model2npp)

#model selection for respiration random effects
model.sel(model00r, model0r, model1r, model2r)

#model0 prevails for npp and r 



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

#model7npp
#going with model9npp
summary(model6npp)
confint(model6npp)

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




#plots
#install.packages("remotes")
remotes::install_github("palaeoverse/rphylopic")
library(rphylopic)
#code to get the uuid's for the phylopic pngs
uuid <- rphylopic::get_uuid(name = "Hypogymnia physodes")





npp_df<- count(npp, study_id, lichen_type, elevation_broad, latitude)


randslope_p<- random.effects(model7npp) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

photosynthesis_joined<- left_join(randslope_p,npp_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "green algae" & elevation_broad == "neutral" ~-11.626354,
                               lichen_type == "green algae" & elevation_broad == "high" ~-11.626354+0.110541, 
                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~-11.6263541+0.625760,
                               lichen_type == "cyanobacteria" & elevation_broad == "high" ~-11.626354+0.625760+0.110541, 
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~ -11.626354+1.216862, 
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~ -11.626354+1.216862+0.110541),
         latitudeii=abs(latitude)*-0.056548,
         f_intercept=r_intercept+c_intercept+latitudeii,
         f_slope=r_slope) %>%
  as.data.frame()

npp_plot<- ggplot(data=npp, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="#abbb80")+
  geom_point(colour="#abbb80", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=photosynthesis_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="#abbb80", size=1, alpha=0.5)+
  geom_abline(aes(slope=0, intercept=-11.626354), colour="#99a873", size=2)+
  #geom_hline(yintercept=-9.1055, size=1.5)+  
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
  ggtitle("NPP")+
  xlab("Temperature (1/kT)")+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=-4, y=-20, ysize=3, alpha=1,fill = "#99a873")+
  ylim(-22,-5)



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
  geom_abline(slope=0.584779, intercept = -15.304981, size=2.5, colour="grey65")+  
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
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=-4, y=-20, ysize=3, alpha=1, fill = "grey65")+
  ylim(-22,-5)
  #scale_x_continuous(sec.axis = sec_axis(~.*(1/(8.61*10^-5))-273.15), name = "Temperature (°C)")
  #scale_x_continuous(sec.axis = sec_axis(~ . * 20, name = "Temperature (°C)"))

#how do i plot this if this is the centred term?



gpp_df<- count(gpp_comb, study_id, lichen_type, elevation_broad, centre_temp, latitude)

randslope_gpp<- random.effects(model6gpp) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

gpp_joined<- left_join(randslope_gpp,gpp_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "cyanobacteria"  ~-13.19796+1.34267,
                               lichen_type == "green algae"  ~-13.19796,
                               lichen_type == "cyanobacteria + green algae" ~-13.19796+0.60857),
         latitudeii=abs(latitude)*-0.03066,
         f_intercept=r_intercept+c_intercept+latitudeii,
         f_slope=r_slope) %>%
  as.data.frame()


gpp_plot <- ggplot(data=gpp_comb, aes(x=centre_temp, y=log(abs(mol_grammin))), colour="olivedrab1")+
  geom_point(colour="olivedrab1", alpha=0.3, size=4)+
  scale_x_reverse(limits=c(4, -4.5))+
  geom_abline(data=gpp_joined, aes(slope=f_slope*-1, intercept=f_intercept), colour="olivedrab1", size=1, alpha=0.5)+
  geom_abline(slope=0, intercept = -11.71171, size=2.5, colour="olivedrab2")+  
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
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=-4, y=-20, ysize=3, alpha=1, fill = "olivedrab2")+
  ylim(-22,-5)+
  ylab("Metabolic rate (log(mmol O2 /mg / min))")

lichen_plots <- (gpp_plot + npp_plot + r_plot)
full_plots <- lichen_plots / coral_plots

full_fig<- full_plots  



ggsave(full_fig, filename = "./figures/fig_c.png", dpi=700, width=20, height=15)


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


#figuring out which studies have both npp and r
count(lichen, study_id, metabolic_category) 

gpp_arithmetic <- read_csv("extraction/lichen dataset/gpp_arithmetic2jun2024.csv") %>%
  mutate(r_value = abs(r_value)*-1,
    gpp = npp_value + r_value) %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol = case_when(gas_units %in% c("mg_co2") ~ gpp*(1/44010), 
                         gas_units %in% c("nmol_co2") ~ gpp*(1/10^9), 
                         gas_units %in% c("umol_co2") ~ gpp*(1/10^6)),
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
         latitude_broad=as.factor(latitude_broad), 
         elevation_broad=as.factor(elevation_broad), 
         lichen_type=as.factor(lichen_type), 
         centre_temp = I(inv_T-mean(inv_T))) %>%
  select(study_id, centre_temp, lichen_type, elevation_broad, latitude, species, mol_grammin, inv_T)

gpp_collect <- gpp %>% select(study_id, centre_temp, lichen_type, elevation_broad, latitude, species, mol_grammin, inv_T) %>%
  mutate(elevation_broad=as.factor(elevation_broad), 
         lichen_type = as.factor(lichen_type))


gpp_comb <- rbind(gpp_collect, gpp_arithmetic) %>%
  filter(case_when(study_id %in% c("adams_1971", "delprado_1995", "domanschke_2013", 
                                   "eickermeieradams1973", "kershaw1983", 
                                   "lechowiczadams1974", "macfarlane1983", "kershaw1977", 
                                   "lange_1980", "lange_1991", "lechowicz1973", "sancho1997",
                                   "smith2001", "sundberg1997", "treitach1997") ~ inv_T < 41, 
                   study_id %in% c("aubert_2007", "coxson_2003", "reiter2000",
                                    "sanchokappen1989", "sonesson1989", 
                                   "green1998") ~ inv_T < 42, 
                   study_id %in% c("brownkershaw1984", "lange2004", "lange2000", 
                                   "zotz2003") ~ inv_T < 40,
                   study_id %in% c("kappen_2000", "schroeter1995") ~ inv_T < 43,
                   study_id %in% c("kappen1983") ~ inv_T < 41.5,
                   study_id %in% c("schroeter1995") ~ inv_T < 42.5,
                   study_id %in% c("zotz1998") ~ inv_T < 39.5,
                   TRUE ~ inv_T < 50))

gpp_comb$elevation_broad<- relevel(gpp_comb$elevation_broad, "neutral")
gpp_comb$lichen_type<- relevel(gpp_comb$lichen_type, "green algae")

model00gpp<- lmer(log(abs(mol_grammin)) ~lichen_type*centre_temp +  abs(latitude) + elevation_broad  + (1+centre_temp|species) , data=gpp_comb, REML=FALSE)
model0gpp <- lmer(log(abs(mol_grammin)) ~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id) , data=gpp_comb, REML=FALSE)
model1gpp <- lmer(log(abs(mol_grammin)) ~lichen_type*centre_temp +  abs(latitude) + elevation_broad+ (1|study_id) , data=gpp_comb, REML=FALSE)
model2gpp <- gls(log(abs(mol_grammin)) ~lichen_type*centre_temp +  abs(latitude) + elevation_broad , data=gpp_comb, na.action=na.exclude)

model.sel(model00gpp, model0gpp, model1gpp, model2gpp)
#model0gpp is best


model3gpp <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=gpp_comb)
#missing latitude
model4gpp <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=gpp_comb)
#missing both latitude and elevation
model5gpp <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=gpp_comb)
#missing temp
model6gpp <- lmer(log(abs(mol_grammin)) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=gpp_comb)
#missing lichen type
model7gpp <- lmer(log(abs(mol_grammin)) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_comb)
#missing interaction between lichen type + centre temp
model8gpp <- lmer(log(abs(mol_grammin)) ~centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_comb)
model9gpp <- lmer(log(abs(mol_grammin))~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=gpp_comb)
model.sel(model3gpp, model4gpp, model5gpp, model6gpp, model7gpp, model8gpp, model9gpp)

#model 6 is best
summary(model6gpp)
confint(model6gpp)



#creating the friedmann and sun figure 

count(lichen_max_ppfd, study_id, metabolic_category)

#has npp, gpp, and r
all_lichen <- lichen_max_ppfd  %>%
  filter(study_id %in% c("nakamura2003","jiang2021", "castillohelmuth2005", 
                         "silbiger2019", "hadjioannou2019")) 
write_csv(all_coral, "f+s_npp_gpp_r.csv")

#has npp and r
npp_and_r_coral <- coral %>%
  filter(study_id %in% c("bahr2018", "godefroid2023", "hill2014", "howe2001", 
                         "juillet-leclerc2014", "rodolfo-metalpa2006", "samiel2015"))
write_csv(npp_and_r_coral, "f+s_npp_r.csv")


#has gpp and r 
#omitting banc-prandi for now because their gpp data literally makes no fucking sense
gpp_and_r_coral <- coral %>%
  filter(study_id %in% c("higuchi2015", "kemp2011"))
write_csv(gpp_and_r_coral, "f+s_gpp_r.csv")

















