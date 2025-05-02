library(revtools)
library(tidyverse)
pacman::p_load(viridis, revtools, nlme, lme4, MuMIn, patchwork, tidyverse)


lichen_tpc_Eas <- read_csv("analysis/tidy data/lichen_tpc_Eas.csv")
lichen_raw_npp <- read_csv("analysis/tidy data/lichen_npp_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, inv_T, mol_gCmin)
lichen_raw_gpp <- read_csv("analysis/tidy data/lichen_gpp_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, inv_T, mol_gCmin)
lichen_raw_r <- read_csv("analysis/tidy data/lichen_r_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, inv_T, mol_gCmin)

gpp_lichen_all <- lichen_tpc_Eas %>% filter(metabolic_category == "gpp") %>%
  left_join(lichen_raw_gpp, by="tpc_grp") %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE)
npp_lichen_all <- lichen_tpc_Eas %>% filter(metabolic_category == "npp") %>%
  left_join(lichen_raw_npp, by="tpc_grp")  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE)
r_lichen_all <- lichen_tpc_Eas %>% filter(metabolic_category == "r") %>%
  left_join(lichen_raw_r, by="tpc_grp")  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE)


#full model --> assessing collinearity 
all_lichen_mod <- lmer(e ~ abs(latitude) +lichen_type + elevation_broad  +(1|study_id), data=npp_lichen_all,  REML=FALSE) 
car::vif(all_lichen_mod) #removing lichen type as this is causing problems

#model selection for photosynthesis
#missing elevation
model3npp_lichen_final <- lmer(e~ abs(latitude) +(1|study_id), data=npp_lichen_all, REML=FALSE)
#missing latitude
model4npp_lichen_final <- lmer(e ~ elevation_broad  +(1|study_id), data=npp_lichen_all, REML=FALSE)
#missing both latitude and elevation
model5npp_lichen_final <- lmer(e ~ 1 +(1|study_id), data=npp_lichen_all, REML=FALSE) 
#missing interaction between lichen type + centre temp
model6npp_lichen_final <- lmer(e ~ abs(latitude)+ elevation_broad +(1|study_id), data=npp_lichen_all,  REML=FALSE)
model.sel(model3npp_lichen_final, model4npp_lichen_final, model5npp_lichen_final, model6npp_lichen_final)
#model6npp, model3npp are tied, going with model5npp to describe overall temperature dependence
summary(model5npp_lichen_final) #Ea= 0.42304
confint(model5npp_lichen_final) #0.34816037 0.49782240


#model selection for respiration
#missing elevation
model3r_lichen_final <- lmer(e~ abs(latitude) +(1|study_id), data=r_lichen_all, REML=FALSE)
#missing latitude
model4r_lichen_final <- lmer(e ~ elevation_broad  +(1|study_id), data=r_lichen_all, REML=FALSE)
#missing both latitude and elevation
model5r_lichen_final <- lmer(e ~ 1 +(1|study_id), data=r_lichen_all, REML=FALSE) 
#missing interaction between lichen type + centre temp
model6r_lichen_final <- lmer(e ~ abs(latitude)+ elevation_broad +(1|study_id), data=r_lichen_all,  REML=FALSE)
model.sel(model3r_lichen_final, model4r_lichen_final, model5r_lichen_final, model6r_lichen_final)
#model6r is best going with model5r to describe overall temperature dependence
summary(model5r_lichen_final) #Ea= 0.5720
confint(model5r_lichen_final) #0.45030538 0.6936945



#arrhenius modelling for respiration
lichen_r_linear <- read_csv("analysis/tidy data/lichen_r_final.csv") %>%
  separate(tpc_grp, into = c("response_id", "study_id"), remove=FALSE) %>%
  mutate(mol_gCmin=abs(mol_gCmin), 
         elevation_broad=as.factor(elevation_broad), 
         lichen_type=as.factor(lichen_type)) %>%
  group_by(tpc_grp) %>%
  # Step 1: Find the temperature at which rate is maximized
  mutate(max_temp = temp[which.max(mol_gCmin)]) %>%
  # Step 2: Keep only temperatures below or equal to that max point
  filter(temp <= max_temp) %>%
  filter(mol_gCmin > 0)

lichen_r_linear$elevation_broad<- relevel(lichen_r_linear$elevation_broad, "neutral")
lichen_r_linear$lichen_type<- relevel(lichen_r_linear$lichen_type, "green algae")


all_lichen_mod <- lmer(log(mol_gCmin) ~inv_T + abs(latitude) +inv_T*lichen_type + elevation_broad  +(1+inv_T|study_id:response_id), data=lichen_r_linear,  REML=FALSE) 
car::vif(all_lichen_mod) #removing lichen type

optCtrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))


#model selection for respiration
#missing elevation
model3r_lichen_final <- lmer(log(mol_gCmin) ~inv_T +abs(latitude)+ (1+inv_T|study_id:response_id), data=lichen_r_linear, REML=FALSE, control=lmerControl(optimizer = "bobyqa",
                                                                                                                                                         optCtrl = list(maxfun = 2e5)))
#missing latitude
model4r_lichen_final <- lmer(log(mol_gCmin)  ~inv_T  + elevation_broad + (1+inv_T|study_id:response_id), data=lichen_r_linear, REML=FALSE, control=lmerControl(optimizer = "bobyqa",
                                                                                                                                                               optCtrl = list(maxfun = 2e5))) 
#missing both latitude and elevation
model5r_lichen_final <- lmer(log(mol_gCmin)  ~inv_T + (1+inv_T|study_id:response_id), data=lichen_r_linear, REML=FALSE, control=lmerControl(optimizer = "bobyqa",
                                                                                                                                           optCtrl = list(maxfun = 2e5)))
#missing temp
model6r_lichen_final <- lmer(log(mol_gCmin)  ~abs(latitude)  + elevation_broad +(1+inv_T|study_id:response_id), data=lichen_r_linear, REML=FALSE, control=lmerControl(optimizer = "bobyqa",
                                                                                                                                                                     optCtrl = list(maxfun = 2e5)))
#missing lichen type
model7r_lichen_final <- lmer(log(mol_gCmin) ~inv_T +  abs(latitude) + elevation_broad +(1+inv_T|study_id:response_id), data=lichen_r_linear, REML=FALSE,  control=lmerControl(optimizer = "bobyqa",
                                                                                                                                                                              optCtrl = list(maxfun = 2e5)))
#missing interaction between lichen type + centre temp
model.sel(model3r_lichen_final, model4r_lichen_final, model5r_lichen_final, model6r_lichen_final, model7r_lichen_final)
#model7r and model3r are tied, model3r is simpler, going w/ that one

summary(model3r_lichen_final) #Ea = -0.51176
confint(model3r_lichen_final) #-0.551134373 -0.47249748





#model selection for GPP
#missing elevation
model3gpp_lichen_final <- lmer(e~ abs(latitude) +(1|study_id), data=gpp_lichen_all, REML=FALSE)
#missing latitude
model4gpp_lichen_final <- lmer(e ~ elevation_broad  +(1|study_id), data=gpp_lichen_all, REML=FALSE)
#missing both latitude and elevation
model5gpp_lichen_final <- lmer(e ~ 1 +(1|study_id), data=gpp_lichen_all, REML=FALSE) 
#missing interaction between lichen type + centre temp
model6gpp_lichen_final <- lmer(e ~ abs(latitude)+ elevation_broad +(1|study_id), data=gpp_lichen_all,  REML=FALSE)
model.sel(model3gpp_lichen_final, model4gpp_lichen_final, model5gpp_lichen_final, model6gpp_lichen_final)
#model6r is best going with model5r to describe overall temperature dependence
summary(model5gpp_lichen_final) #Ea= 0.34629
confint(model5gpp_lichen_final) #0.2772515 0.41520588





################# CORAL ################# 

coral_tpc_Eas <- read_csv("analysis/tidy data/coral_tpc_Eas.csv") 
coral_raw_npp <- read_csv("analysis/tidy data/coral_npp_final.csv") %>% select(tpc_grp, depth_broad, broad_coral, latitude, species)
coral_raw_gpp <- read_csv("analysis/tidy data/coral_gpp_final.csv") %>% select(tpc_grp, depth_broad, broad_coral, latitude, species) 
coral_raw_r <- read_csv("analysis/tidy data/coral_r_final.csv") %>% select(tpc_grp, depth_broad, broad_coral, latitude, species)

gpp_coral_all <- coral_tpc_Eas %>% filter(metabolic_category == "gpp") %>%
  left_join(coral_raw_gpp, by="tpc_grp") %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE)
npp_coral_all <- coral_tpc_Eas %>% filter(metabolic_category == "npp") %>%
  left_join(coral_raw_npp, by="tpc_grp")  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE)
r_coral_all <- coral_tpc_Eas %>% filter(metabolic_category == "r") %>%
  left_join(coral_raw_r, by="tpc_grp")  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE)

#model selection for npp
#missing depth
model3npp_coral_final <- lmer(e ~ abs(latitude) + (1|study_id), data=npp_coral_all, REML=FALSE)
#missing latitude
model4npp_coral_final <- lmer(e ~ 1 + (1|study_id), data=npp_coral_all, REML=FALSE)
model.sel(model3npp_coral_final, model4npp_coral_final)
#model3npp_coral_final is best
#using model4npp to describe the overall Ea for NPP
summary(model4npp_coral_final) # Ea= 0.8669 
confint(model4npp_coral_final) #0.4336468 1.3065570

#model selection for respiration
#missing depth
model3r_coral_final <- lmer(e ~ abs(latitude) +  (1|study_id), data=r_coral_all, REML=FALSE)
#missing latitude
model4r_coral_final <- lmer(e ~depth_broad + (1|study_id), data=r_coral_all, REML=FALSE)
#missing both latitude and depth
model5r_coral_final <- lmer(e ~1 + (1|study_id), data=r_coral_all, REML=FALSE)
#missing interaction between lichen type + centre temp
model6r_coral_final <- lmer(e ~ abs(latitude) + depth_broad + (1|study_id), data=r_coral_all,  REML=FALSE)
model.sel(model3r_coral_final, model4r_coral_final, model5r_coral_final, model6r_coral_final)
#model6r_coral_final and model3r are tied
#using model5r_coral_final to describe the overall Ea for NPP

summary(model5r_coral_final) #e= 1.0218 
confint(model5r_coral_final) #0.6470683 1.3974794



#model selection for gpp
#missing depth
model3gpp_coral_final <- lmer(e ~abs(latitude) +  (1|study_id), data=gpp_coral_all, REML=FALSE)
#missing latitude
model4gpp_coral_final <- lmer(e ~ depth_broad +(1|study_id) , data=gpp_coral_all, REML=FALSE)
#missing both latitude and depth
model5gpp_coral_final <- lmer(e ~ 1+  (1|study_id), data=gpp_coral_all, REML=FALSE)
#missing temp
model6gpp_coral_final <- lmer(e ~abs(latitude) + depth_broad  + (1|study_id), data=gpp_coral_all, REML=FALSE)
#missing cnidarian type
#full model
model.sel(model3gpp_coral_final, model4gpp_coral_final, model5gpp_coral_final,  model6gpp_coral_final)
#model3gpp_coral_final, model6gpp_coral_final are tied
#going with model5gpp_coral_final to understand overal Ea
summary(model5gpp_coral_final) #e= 0.9210
confint(model5gpp_coral_final) #e=0.5467659 1.2965725




##################### Figures ######################
#plots
#install.packages("remotes")
pacman::p_load(rphylopic)
#code to get the uuid's for the phylopic pngs
uuid <- rphylopic::get_uuid(name = "Hypogymnia physodes")

lichen_randslope_npp<- random.effects(model5npp_lichen_final) %>%
  as.data.frame() %>%
  select(-condsd, -grpvar) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp) %>%
  mutate(f_intercept=r_intercept+0.423)

lichen_npp_plot<- ggplot()+
  scale_x_reverse(limits=c(0, -5))+
  scale_y_continuous(limits=c(0,7))+
  geom_abline(data=lichen_randslope_npp, aes(slope=f_intercept, intercept=0), colour="#abbb80", size=1, alpha=0.25)+
  geom_abline(aes(slope=0.423, intercept=0), colour="#abbb80", size=3)+
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
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=-4.7, y=6, height=2, alpha=1,fill = "#99a873")+
  geom_text(aes(x=-1, y=6.7,label = "Ea = -0.423±0.07"), size=6)


lichen_randslope_r<- random.effects(model3r_lichen_final) %>%
  as.data.frame() %>%
  select(-condsd, -grpvar) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=inv_T) %>%
  mutate(f_slope=r_slope+-0.51176, 
         f_intercept=r_intercept+4.79692)



#not sure how to plot this be
lichen_r_plot <- ggplot()+
  scale_x_reverse(limits=c(0, -5))+
  scale_y_continuous(limits=c(0,7))+
  geom_abline(data=lichen_randslope_r, aes(slope=f_slope*-1, intercept=f_intercept), colour="grey", size=1, alpha=0.5)+
  geom_abline(slope=0.51176, intercept = 4.79692, size=2.5, colour="grey65")+  
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
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=-4.7, y=6, height=2, alpha=1,fill = "grey65")+
  geom_text(aes(x=-0.8, y=6.4,label = "Ea = -0.51176"), size=8)

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






