library(revtools)
library(tidyverse)
pacman::p_load(viridis, revtools, nlme, lme4, MuMIn, patchwork, tidyverse, ggh4x, dplyr, purrr, broom)

lichen_tpc_Eas <- read_csv("analysis/tidy data/lichen_tpc_Eas.csv")
lichen_raw_npp <- read_csv("analysis/tidy data/lichen_npp_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, inv_T, mol_gCmin)
lichen_raw_gpp <- read_csv("analysis/tidy data/lichen_gpp_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, inv_T, mol_gCmin)
lichen_raw_r <- read_csv("analysis/tidy data/lichen_r_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, inv_T, mol_gCmin)

count(lichen_raw_r, species)


gpp_lichen_Eas <- lichen_tpc_Eas %>% filter(metabolic_category == "gpp") %>%
  left_join(lichen_raw_gpp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2))
npp_lichen_Eas <- lichen_tpc_Eas %>% filter(metabolic_category == "npp")  %>%
  left_join(lichen_raw_npp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2)) %>% drop_na(e)
r_lichen_Eas <- lichen_tpc_Eas %>% filter(metabolic_category == "r")  %>%
  left_join(lichen_raw_r, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2))






#### lichen 

#model selection for NPP
lichen_npp <- lme(e ~ abs(latitude) +lichen_type + elevation_broad +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=npp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_npp)

lichen_npp_model_set <- dredge(lichen_npp, trace = TRUE)
subset(lichen_npp_model_set, delta < 2) #tied w/ 2 best, simplest model includes only breadth

best_lichen_npp_model <- get.models(lichen_npp_model_set, 2)[[1]]
best_lichen_npp_model <- update(best_lichen_npp_model, method = "REML")
summary(best_lichen_npp_model) #0.3926778
intervals(best_lichen_npp_model) #0.34549059 0.43986506

#model selection for GPP
lichen_gpp <- lme(e ~ abs(latitude) +lichen_type + elevation_broad +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=gpp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_gpp)

lichen_gpp_model_set <- dredge(lichen_gpp, trace = TRUE)
subset(lichen_gpp_model_set, delta < 2) #best simplest model has breadth and topt

best_lichen_gpp_model <- get.models(lichen_gpp_model_set, 1)[[1]]
best_lichen_gpp_model <- update(best_lichen_gpp_model, method = "REML")
summary(best_lichen_gpp_model) #0.538573046
intervals(best_lichen_gpp_model) #0.379509242 0.6976368495


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


# Nest data by group
nested <- lichen_r_linear %>%
  group_by(tpc_grp) %>%
  nest()

# Fit linear model and extract slope + SE
group_ea <- nested %>%
  mutate(model = map(data, ~ lm(log(mol_gCmin) ~ inv_T, data = .x)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "inv_T") %>%
  separate(tpc_grp, into = c("response_id","study_id"), sep = "_", remove=FALSE) %>%
  rename(Ea = estimate, Ea_SE = std.error) %>%
  select(study_id, response_id, Ea, Ea_SE)


group_ea <- left_join(group_ea, lichen_r_linear, by = c("study_id", "response_id"), relationship = "many-to-many") %>%
  select(-tpc_grp.y) %>% rename(tpc_grp = tpc_grp.x) %>% distinct(Ea, .keep_all = TRUE) %>% drop_na(Ea, Ea_SE) 


group_ea$elevation_broad<- relevel(group_ea$elevation_broad, "neutral")
group_ea$lichen_type<- relevel(group_ea$lichen_type, "green algae")
# 4: Mixed model with weights
ea_model <- lme(Ea ~ abs(latitude) + lichen_type + elevation_broad,
                random = ~1 | study_id,
                data = group_ea,
                weights = varFixed(~ Ea_SE^2),
                method = "ML")


lichen_r_model_set <- dredge(ea_model, trace = TRUE)
subset(lichen_r_model_set, delta < 2) #tied, going with model 5

best_lichen_r_model <- get.models(lichen_r_model_set, 1)[[1]]
best_lichen_r_model <- update(best_lichen_r_model, method = "REML")
summary(best_lichen_r_model) #-0.5894234
intervals(best_lichen_r_model) #-0.6041062 -0.5747406








################# CORAL ################# 

coral_tpc_Eas <- read_csv("analysis/tidy data/coral_tpc_Eas.csv") 
coral_raw_npp <- read_csv("analysis/tidy data/coral_npp_final.csv") %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin)
coral_raw_gpp <- read_csv("analysis/tidy data/coral_gpp_final.csv") %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin) 
coral_raw_r <- read_csv("analysis/tidy data/coral_r_final.csv") %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin)

lichen_raw_r %>% left_join(lichen_raw_gpp) %>% left_join(lichen_raw_npp) %>% count(species)

gpp_coral_Eas <- coral_tpc_Eas %>% filter(metabolic_category == "gpp") %>%
  left_join(coral_raw_gpp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2)) %>% drop_na(e, e_se)
npp_coral_Eas <- coral_tpc_Eas %>% filter(metabolic_category == "npp") %>%
  left_join(coral_raw_npp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth,.keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2)) %>% drop_na(e, e_se)
r_coral_Eas <- coral_tpc_Eas %>% filter(metabolic_category == "r") %>%
  left_join(coral_raw_r, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2)) %>% drop_na(e, e_se)




#model selection for NPP
coral_npp <- lme(e ~ abs(latitude) +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=npp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_npp)

coral_npp_model_set <- dredge(coral_npp, trace = TRUE)
subset(coral_npp_model_set, delta < 2) #tied, choosing model includes only breadth

best_coral_npp_model <- get.models(coral_npp_model_set, 2)[[1]]
best_coral_npp_model <- update(best_coral_npp_model, method = "REML")
summary(best_coral_npp_model) #0.5937617
intervals(best_coral_npp_model) #0.2367783 0.77740006

#model selection for GPP
coral_gpp <- lme(e ~ abs(latitude) +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=gpp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_gpp)

coral_gpp_model_set <- dredge(coral_gpp, trace = TRUE)
subset(coral_gpp_model_set, delta < 2) 

best_coral_gpp_model <- get.models(coral_gpp_model_set, 1)[[1]]
best_coral_gpp_model <- update(best_coral_gpp_model, method = "REML")
summary(best_coral_gpp_model) #0.6564679
intervals(best_coral_gpp_model) #0.311348240.1.00158748

#model selection for R
coral_r <- lme(e ~ abs(latitude) +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=r_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_r) #removing latitude because of correlation
coral_r <- lme(e ~ I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=r_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 



coral_r_model_set <- dredge(coral_r, trace = TRUE)
subset(coral_r_model_set, delta < 2) 

best_coral_r_model <- get.models(coral_r_model_set, 1)[[1]]
best_coral_r_model <- update(best_coral_r_model, method = "REML")
summary(best_coral_r_model) #0.5758236
intervals(best_coral_r_model, which = "fixed") #0.5192733   0.632373938


#### mary's suggestion 

lme(log(rmax) ~ topt, random= ~1|study_id, data=gpp_coral_Eas, method = "REML") #0.04148732
lme(log(rmax) ~ topt, random= ~1|study_id, data=npp_coral_Eas, method = "REML") #0.0432332
lme(log(rmax) ~ topt, random= ~1|study_id, data=r_coral_Eas, method = "REML") #0.0105437
lme(log(rmax) ~ topt, random= ~1|study_id, data=gpp_lichen_Eas, method = "REML") #0.01736202
lme(log(rmax) ~ topt, random= ~1|study_id, data=npp_lichen_Eas, method = "REML") #0.06472847 






##### david's suggestion

npp_lichen_topt  <- npp_lichen_Eas %>% select(topt, study_id) %>% mutate(metabolic_category = "npp")
gpp_lichen_topt  <- gpp_lichen_Eas %>% select(topt, study_id) %>% mutate(metabolic_category = "gpp")
r_lichen_topt <- r_lichen_Eas %>% select(topt, study_id) %>% mutate(metabolic_category = "r")

npp_coral_topt  <- npp_coral_Eas %>% select(topt, study_id) %>% mutate(metabolic_category = "npp") 
gpp_coral_topt  <- gpp_coral_Eas %>% select(topt, study_id) %>% mutate(metabolic_category = "gpp")
r_coral_topt  <- r_coral_Eas %>% select(topt, study_id) %>% mutate(metabolic_category = "r")

lichen_topt <- rbind(npp_lichen_topt, gpp_lichen_topt, r_lichen_topt) %>% mutate(dataset="Lichen")
coral_topt <- rbind(npp_coral_topt, gpp_coral_topt, r_coral_topt) %>% mutate(dataset="Coral")

full_topt <- rbind(coral_topt, lichen_topt) %>% mutate(metabolic_category=as.factor(metabolic_category))


full_topt$metabolic_category <- relevel(full_topt$metabolic_category, ref = "npp")

topt_mod <-lme(topt ~ metabolic_category*dataset, random=~1|study_id, data=full_topt, method="REML")

summary(topt_mod)
intervals(topt_mod)


##thermal breadth

npp_lichen_tbr  <- npp_lichen_Eas %>% select(breadth, study_id) %>% mutate(metabolic_category = "npp")
gpp_lichen_tbr  <- gpp_lichen_Eas %>% select(breadth, study_id) %>% mutate(metabolic_category = "gpp")
r_lichen_tbr <- r_lichen_Eas %>% select(breadth, study_id) %>% mutate(metabolic_category = "r")

npp_coral_tbr  <- npp_coral_Eas %>% select(breadth, study_id) %>% mutate(metabolic_category = "npp") 
gpp_coral_tbr  <- gpp_coral_Eas %>% select(breadth, study_id) %>% mutate(metabolic_category = "gpp")
r_coral_tbr  <- r_coral_Eas %>% select(breadth, study_id) %>% mutate(metabolic_category = "r")

lichen_tbr <- rbind(npp_lichen_tbr, gpp_lichen_tbr, r_lichen_tbr) %>% mutate(dataset="Lichen")
coral_tbr <- rbind(npp_coral_tbr, gpp_coral_tbr, r_coral_tbr) %>% mutate(dataset="Coral")

full_tbr <- rbind(coral_tbr, lichen_tbr) %>% mutate(metabolic_category=as.factor(metabolic_category))


full_tbr$metabolic_category <- relevel(full_tbr$metabolic_category, ref = "npp")

tbr_mod <-lme(breadth ~ metabolic_category*dataset, random=~1|study_id, data=full_tbr, method="REML")

summary(tbr_mod)
intervals(tbr_mod)


#isolated algae
isl_algae_tpc_Eas <- read_csv("analysis/tidy data/isl_alg_tpc_Eas.csv") 
isl_algae_raw_npp <- read_csv("analysis/tidy data/isl_algae_npp_final.csv") %>% select(tpc_grp, latitude, species, mol_gCmin)
isl_algae_raw_gpp <- read_csv("analysis/tidy data/isl_algae_gpp_final.csv") %>% select(tpc_grp, latitude, species, mol_gCmin) 
isl_algae_raw_r <- read_csv("analysis/tidy data/isl_algae_r_final.csv") %>% select(tpc_grp, latitude, species, mol_gCmin)

gpp_isl_algae_Eas <- isl_algae_tpc_Eas %>% filter(metabolic_category == "gpp") %>%
  left_join(isl_algae_raw_gpp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2)) %>% drop_na(e, e_se)
npp_isl_algae_Eas <- isl_algae_tpc_Eas %>% filter(metabolic_category == "npp") %>%
  left_join(isl_algae_raw_npp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth,.keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2)) %>% drop_na(e, e_se)
r_isl_algae_Eas <- isl_algae_tpc_Eas %>% filter(metabolic_category == "r") %>%
  left_join(isl_algae_raw_r, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>%
  mutate(weight = 1 / (e_se^2)) %>% drop_na(e, e_se)




#model selection for NPP
isl_alg_npp <- lme(e ~ abs(latitude) +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=npp_isl_algae_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_npp)

isl_alg_npp_model_set <- dredge(isl_alg_npp, trace = TRUE)
subset(isl_alg_npp_model_set, delta < 2) #latitude is best by more than 2 AIC

best_isl_alg_npp_model <- get.models(isl_alg_npp_model_set, 1)[[1]]
best_isl_alg_npp_model <- update(best_isl_alg_npp_model, method = "REML")
summary(best_isl_alg_npp_model) #1.9009684
intervals(best_isl_alg_npp_model) #0.27472413 3.52721263

#model selection for GPP
isl_alg_gpp <- lme(e ~ abs(latitude) +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=gpp_isl_algae_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(isl_alg_gpp)

isl_alg_gpp_model_set <- dredge(isl_alg_gpp, trace = TRUE)
subset(isl_alg_gpp_model_set, delta < 2) # tied, keeping model with latitude, breadth, topt

best_isl_alg_gpp_model <- get.models(isl_alg_gpp_model_set, 1)[[1]]
best_isl_alg_gpp_model <- update(best_isl_alg_gpp_model, method = "REML")
summary(best_isl_alg_gpp_model) #0.3808134
intervals(best_isl_alg_gpp_model, which = "fixed") #0.3221303944 0.439496314

#model selection for R
isl_alg_r <- lme(e ~ abs(latitude) +I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=r_isl_algae_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(isl_alg_r) #removing latitude because of correlation
isl_alg_r <- lme(e ~ I(breadth - mean(breadth, na.rm = TRUE)) +I(topt- mean(topt, na.rm = TRUE)), random= ~1|study_id, data=r_isl_algae_Eas, weights = varFixed(~ e_se^2), method = "ML") 



isl_alg_r_model_set <- dredge(isl_alg_r, trace = TRUE)
subset(isl_alg_r_model_set, delta < 2) #best model includes breath and topt 

best_isl_alg_r_model <- get.models(isl_alg_r_model_set, 1)[[1]]
best_isl_alg_r_model <- update(best_isl_alg_r_model, method = "REML")
summary(best_isl_alg_r_model) #0.3231435
intervals(best_isl_alg_r_model, which = "fixed") #0.04542771 0.6008593














##################### Figures ######################
#plots
#install.packages("remotes")
pacman::p_load(rphylopic)
#code to get the uuid's for the phylopic pngs
uuid <- rphylopic::get_uuid(name = "Hypogymnia physodes")


npp_coral_breadth <-npp_coral_Eas %>%
  mutate(centred_breadth=I(breadth-mean(breadth))) %>%
  ggplot(aes(x=as.numeric(centred_breadth), y=e))+
  geom_point(colour="#abbb80", alpha=0.5, size=2)+
  geom_abline(slope =-0.0938028, intercept = 0.59, colour="#abbb80")+
  theme_bw()+
  xlab("Centred Breadth")+
  ylab("Ea")+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=13.25, y=3.5, height=0.75, alpha=1,fill = "black")+
  ylim(-0.2,4)+
  xlim(-15,15)

gpp_coral_breadth <-gpp_coral_Eas %>%
  mutate(centred_breadth=I(breadth-mean(breadth))) %>%
  ggplot(aes(x=as.numeric(centred_breadth), y=e))+
  geom_point(colour="olivedrab2", alpha=0.5, size=2)+
  geom_abline(slope =-0.0493594, intercept = 0.65, colour="olivedrab2")+
  theme_bw()+
  xlab("Centred Breadth")+
  ylab("Ea")+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=13.25, y=3.5, height=0.75, alpha=1,fill = "black")+
  ylim(-0.2,4)+
  xlim(-15,15)

r_coral_breadth <-r_coral_Eas %>%
  mutate(centred_breadth=I(breadth-mean(breadth))) %>%
  ggplot(aes(x=as.numeric(centred_breadth), y=e))+
  geom_point(colour="grey", alpha=0.5, size=2)+
  geom_abline(slope =-0.0636694, intercept = 0.5210243, colour="grey")+
  theme_bw()+
  xlab("Centred Breadth")+
  ylab("Ea")+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=13.25, y=3.5, height=0.75, alpha=1,fill = "black")+
  ylim(-0.2,4)+
  xlim(-15,15)


npp_lichen_breadth <-npp_lichen_Eas %>%
  mutate(centred_breadth=I(breadth-mean(breadth))) %>%
  ggplot(aes(x=as.numeric(centred_breadth), y=e))+
  geom_point(colour="#abbb80", alpha=0.5, size=2)+
  geom_abline(slope =-0.0198103, intercept = 0.3926762, colour="#abbb80")+
  theme_bw()+
  xlab("Centred Breadth")+
  ylab("Ea")+
  ylim(-0.2,4)+
  xlim(-15,15)+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=13.25, y=3.5, height=1, alpha=1,fill = "black")


gpp_lichen_breadth <-gpp_lichen_Eas %>%
  mutate(centred_breadth=I(breadth-mean(breadth))) %>%
  ggplot(aes(x=as.numeric(centred_breadth), y=e))+
  geom_point(colour="olivedrab2", alpha=0.5, size=2)+
  geom_abline(slope =-0.0205494, intercept = 0.5389686, colour="olivedrab2")+
  theme_bw()+
  xlab("Centred Breadth")+
  ylab("Ea")+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=13.25, y=3.5, height=1, alpha=1,fill = "black")+
  ylim(-0.2,4)+
  xlim(-15,15)


(gpp_coral_breadth + gpp_lichen_breadth) / (npp_coral_breadth + npp_lichen_breadth) / (r_coral_breadth + plot_spacer())
  
  




#trying to just plot the Eas

npp_lichen_ea  <- npp_lichen_Eas %>% select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "npp")
gpp_lichen_ea  <- gpp_lichen_Eas %>% select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "gpp")
r_lichen_linear_ea  <- random.effects(best_lichen_r_model) %>% as.data.frame() %>% 
  select(`(Intercept)`) %>% rename(e=`(Intercept)`) %>% mutate(metabolic_category = "r", eh= c(""), topt= c("")) %>%    
  mutate(E=e+-0.5894235, e=abs(E)) %>% select(-E)
r_lichen_ea <- r_lichen_Eas %>% select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "r")

npp_coral_ea  <- npp_coral_Eas %>% select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "npp") 
gpp_coral_ea  <- gpp_coral_Eas %>% select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "gpp")
r_coral_ea  <- r_coral_Eas %>% select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "r")


lichen_all <- rbind(npp_lichen_ea, gpp_lichen_ea, r_lichen_ea) %>% mutate(dataset="Lichen")
coral_all <- rbind(npp_coral_ea, gpp_coral_ea, r_coral_ea) %>% mutate(dataset="Coral")
lit <- data.frame(e=c("", "", ""), eh=c("", "", ""), topt=c("", "", ""), breadth =c("", "", ""), rmax =c("", "", ""), study_id=c("", "", ""), metabolic_category=c("gpp", "r", "npp"), dataset=c("López-Urrutia et al. 2006", "López-Urrutia et al. 2006", "López-Urrutia et al. 2006"))

full_Ea <- rbind(lichen_all, coral_all, lit) %>% mutate(e=as.numeric(e), topt=as.numeric(topt))



full_topt_breadth <- rbind(lichen_all, coral_all) %>% mutate(e=as.numeric(e), topt=as.numeric(topt), breadth=as.numeric(breadth), 
                                                             metabolic_category=as.factor(metabolic_category), dataset=as.factor(dataset))


mod_npp_l<-coef(summary(best_lichen_npp_model)) %>% as.data.frame() %>% slice(1) %>% select(Value, `Std.Error`) %>%
  mutate(dataset="Lichen", metabolic_category="npp")  %>% rename(se=`Std.Error`)

mod_r_l<-coef(summary(best_lichen_r_model)) %>% as.data.frame() %>% slice(1) %>% select(Value, `Std.Error`) %>%
  mutate(dataset="Lichen", metabolic_category="r", Value=abs(Value)) %>% rename(se=`Std.Error`)

mod_gpp_l<-coef(summary(best_lichen_gpp_model)) %>% as.data.frame() %>% slice(1) %>% select(Value, `Std.Error`) %>%
    mutate(dataset="Lichen", metabolic_category="gpp")   %>% rename(se=`Std.Error`)
  
mod_npp_c<-coef(summary(best_coral_npp_model)) %>% as.data.frame() %>% slice(1) %>% select(Value, `Std.Error`) %>%
  mutate(dataset="Coral", metabolic_category="npp")  %>% rename(se=`Std.Error`)

mod_r_c<-coef(summary(best_coral_r_model)) %>% as.data.frame() %>% slice(1) %>% select(Value, `Std.Error`) %>%
  mutate(dataset="Coral", metabolic_category="r")  %>% rename(se=`Std.Error`)

mod_gpp_c<-coef(summary(best_coral_gpp_model)) %>% as.data.frame() %>% slice(1) %>% select(Value, `Std.Error`) %>%
  mutate(dataset="Coral", metabolic_category="gpp") %>% rename(se=`Std.Error`)



mod_lit <- data.frame(Value=c(0.33, 0.56, 0.29), 
                      metabolic_category=c("gpp", "r", "npp"), 
                      se=c(0.089,0.024,0.036), 
                      dataset=c("López-Urrutia et al. 2006", "López-Urrutia et al. 2006", "López-Urrutia et al. 2006"))

full_mod_Ea <- rbind(mod_npp_l, mod_r_l, mod_gpp_l, mod_npp_c, mod_r_c, mod_gpp_c,
                     mod_lit)





######topt/tbr plot

full_tbr <- rbind(coral_tbr, lichen_tbr) %>% mutate(metabolic_category=as.factor(metabolic_category))


full_tbr$metabolic_category <- relevel(full_tbr$metabolic_category, ref = "npp")
full_tbr$dataset <- relevel(full_tbr$dataset, ref = "coral")

tbr_mod <-lme(breadth ~ metabolic_category*dataset, random=~1|study_id, data=full_tbr, method="REML")
summary(tbr_mod)

full_topt_breadth$metabolic_category <- relevel(full_topt_breadth$metabolic_category, ref = "npp")
topt_mod <-lme(topt ~ metabolic_category*dataset, random=~1|study_id, data=full_topt_breadth, method="REML")
summary(topt_mod)



topt_tbr_modeloutput <- read_csv("topt_tbr_modeloutput.csv")
topt_mod <- topt_tbr_modeloutput %>% filter(response=="topt")
tbr_mod <- topt_tbr_modeloutput %>% filter(response=="tbr")


main_topt_plot<-ggplot()+
  geom_jitter(data=full_topt, aes(x=dataset, y=topt, colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=topt_mod , aes(x=dataset, y=estimate, ymin=estimate-se, ymax=estimate+se, colour=metabolic_category), size=0.8, linewidth=0.8)+
  facet_wrap2(~metabolic_category, ncol=1, strip = strip_themed(
    background_x = facet_colors))+
  coord_flip()+
 # geom_text(data = labels_df2,
          #  aes(x = x, y = y, label = label),
          #  color = "black", hjust = 0, vjust=1, size = 2.5) +
  xlab("")+
  ylab("Thermal optimum (Topt)")+
  scale_colour_manual(values=c("olivedrab2","#abbb80", "grey"))+
  theme_classic()+
  theme(legend.position="none")+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=2, y=40, height=0.6, alpha=1,fill = "black")+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=1, y=40, height=0.5, alpha=1,fill = "black")


#ggsave(main_topt_plot, filename = "./figures/main_topt_plot.png", dpi=700, width=5, height=6)

main_tbr_plot<-ggplot()+
  geom_jitter(data=full_topt, aes(x=dataset, y=breadth, colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=tbr_mod , aes(x=dataset, y=estimate, ymin=estimate-se, ymax=estimate+se, colour=metabolic_category), size=0.8, linewidth=0.8)+
  facet_wrap2(~metabolic_category, ncol=1, strip = strip_themed(
    background_x = facet_colors))+
  coord_flip()+
  # geom_text(data = labels_df2,
  #  aes(x = x, y = y, label = label),
  #  color = "black", hjust = 0, vjust=1, size = 2.5) +
  xlab("")+
  ylab("Thermal breadth (Tbr)")+
  scale_colour_manual(values=c("olivedrab2","#abbb80", "grey"))+
  theme_classic()+
  theme(legend.position="none")+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=2, y=30, height=0.6, alpha=1,fill = "black")+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=1, y=30, height=0.5, alpha=1,fill = "black")




##### figure 3




lu <- full_mod_Ea %>% filter(dataset=="López-Urrutia et al. 2006") 

fig3a <- ggplot()+
  #geom_jitter(aes(x=metabolic_category, y=e, colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=lu, aes(x=metabolic_category, y=Value, ymin=Value-se, ymax=Value+se, colour=metabolic_category), size=1, linewidth=2.4,lineend='round')+
  coord_flip()+
  scale_x_discrete(limits=c("gpp", "npp", "r"), labels=c("GPP", "NPP", "R"))+
  scale_y_continuous(limits = c(0, 2))+
  scale_colour_manual(values=c("olivedrab2","#abbb80", "grey"), )+
  labs( x="",
    y = "Temperature dependence (Ea)",
    title = "Non-symbiotic") +
  theme_bw() +
  theme(legend.position="none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

lichen_data <- full_Ea %>% filter(dataset=="Lichen") 
lichen_mod <- full_mod_Ea %>% filter(dataset=="Lichen") 

labels_df_l <- data.frame(
  dataset = c("Lichen", "Lichen", "Lichen"),           
  metabolic_category = c("gpp", "npp", "r"),
  y = c(1.75, 1.75, 1.75),  
  x= c(1, 2, 3),
  label = c("n = 43", "n = 53", "n = 8"))



fig3b <- ggplot()+
  geom_jitter(data=lichen_data, aes(x=metabolic_category, y=e, colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=lichen_mod, aes(x=metabolic_category, y=Value, ymin=Value-se, ymax=Value+se, colour=metabolic_category), size=1, linewidth=2.4,lineend='round')+
  coord_flip()+
  scale_x_discrete(limits=c("gpp", "npp", "r"), labels=c("GPP", "NPP", "R"))+
  scale_y_continuous(limits = c(0, 2))+
  scale_colour_manual(values=c("olivedrab2","#abbb80", "grey"), )+
  labs( x="",
        y = "Temperature dependence (Ea)",
        title = "Lichen") +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=3.3, y=1.75, height=0.5, alpha=1)+
   geom_text(data = labels_df_l,
   aes(x = x, y = y, label = label),
   color = "black", hjust = 0, vjust=1, size = 4)
  

coral_data <- full_Ea %>% filter(dataset=="Coral") 
coral_mod <- full_mod_Ea %>% filter(dataset=="Coral") 

labels_df_c <- data.frame(
  dataset = c("Coral", "Coral", "Coral"),           
  metabolic_category = c("gpp", "npp", "r"),
  y = c(1.75, 1.75, 1.75),  
  x= c(1, 2, 3),
  label = c("n = 46", "n = 40", "n = 57"))


fig3c <- ggplot()+
  geom_jitter(data=coral_data, aes(x=metabolic_category, y=e, colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=coral_mod, aes(x=metabolic_category, y=Value, ymin=Value-se, ymax=Value+se, colour=metabolic_category), size=1, linewidth=2.4,lineend='round')+
  coord_flip()+
  geom_text(data = labels_df_c,
            aes(x = x, y = y, label = label),
            color = "black", hjust = 0, vjust=1, size = 4)+
  scale_x_discrete(limits=c("gpp", "npp", "r"), labels=c("GPP", "NPP", "R"))+
  scale_y_continuous(limits = c(0, 2))+
  scale_colour_manual(values=c("olivedrab2","#abbb80", "grey"), )+
  labs( x="",
        y = "Temperature dependence (Ea)",
        title = "Coral") +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=3.3, y=1.75, height=0.5, alpha=1,fill = "black")



figure3 <- fig3a + fig3b + fig3c 

ggsave(figure3, filename = "./figures/figure3.png", dpi=700, width=12, height=4)






### making the friedman and sun figure



# Constants
k <- 8.617e-5  # eV/K

# Reference temperature for scaling (e.g. 20°C)
T_ref <- 20 + 273.15

my_colours <- c("#abbb80", "olivedrab1", "grey")

labels_df_4a <- data.frame(          
  metabolic_category = c("gpp", "npp", "r"),
  y = c(-7.9, -6.9, -5.9),  
  x= c(0, 0, 0),
  label = c("Ea = 0.32", "Ea = 0.36", "Ea = 0.65"),
  my_colours=  c("olivedrab1","#abbb80", "grey"))

# Build the plot
figure4a <- ggplot(data.frame(T_C = c(0, 45)), aes(x = T_C)) +
  # First function (Ea = 0.3)
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.3 / (k * T_ref))
    - exp(-0.3 / (k * T_K)) / ref_rate
  }, color = my_colours[1], size = 1.5, linetype="longdash") +
  
  # GPP (Ea = 0.32, positive)
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.32 / (k * T_ref))
    exp(-0.32 / (k * T_K)) / ref_rate
  }, color = my_colours[2], size = 1.5) +
  
  # Respiration (Ea = 0.65, negative flux)
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.65 / (k * T_ref))
    - exp(-0.65 / (k * T_K)) / ref_rate  # negative sign here
  }, color = my_colours[3], size = 1.5) +
  
  labs(
    x = "Temperature (°C)",
    y = "Metabolic rate (normalized at 20°C)",
    title = "Non-symbiotic"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )+
  ylim(-8,8)+
  geom_text(data = labels_df_4a, aes(x = x, y = y, label = label, color = my_colours),
            hjust = 0, vjust = 1, size = 4, show.legend = FALSE, fontface = "bold") +
  scale_color_identity()







# Build the plot

labels_df_4b <- data.frame(          
  metabolic_category = c("gpp", "npp", "r"),
  y = c(-7.9, -6.9, -5.9),  
  x= c(0, 0, 0),
  label = c("Ea = 0.54", "Ea = 0.39", "Ea = 0.59"),
  my_colours=  c("olivedrab1","#abbb80", "grey"))



figure4b<- ggplot(data.frame(T_C = c(0, 45)), aes(x = T_C)) +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.0504/ (k * T_ref))  # normalize to 1 at 20°C
    - exp(-0.0504 / (k * T_K)) / ref_rate
  }, color = my_colours[1], size = 1.5, linetype="longdash") +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.3936/ (k * T_ref))  # normalize to 1 at 20°C
    exp(-0.3936 / (k * T_K)) / ref_rate
  }, color = my_colours[1], size = 1.5) +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.5390 / (k * T_ref))
    exp(-0.5390 / (k * T_K)) / ref_rate
  }, color = my_colours[2], size = 1.5) +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.5894 / (k * T_ref))
    - exp(-0.5894 / (k * T_K)) / ref_rate
  }, color = my_colours[3], size = 1.5) +
  labs(
    x = "Temperature (°C)",
    y = "",
    title = "Lichen"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )+
  #ylim(0,8)+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=3, y=6.5, height=3.3, alpha=1)+
  ylim(-8,8)+
  geom_text(data = labels_df_4b, aes(x = x, y = y, label = label, color = my_colours),
            hjust = 0, vjust = 1, size = 4, show.legend = FALSE, fontface = "bold") +
  scale_color_identity()



labels_df_4c <- data.frame(          
  metabolic_category = c("gpp", "npp", "r"),
  y = c(-7.9, -6.9, -5.9),  
  x= c(0, 0, 0),
  label = c("Ea = 0.66", "Ea = 0.59", "Ea = 0.58"),
  my_colours=  c("olivedrab1","#abbb80", "grey"))



figure4c<- ggplot(data.frame(T_C = c(0, 45)), aes(x = T_C)) +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.13/ (k * T_ref))  # normalize to 1 at 20°C
    exp(-0.13 / (k * T_K)) / ref_rate
  }, color = my_colours[1], size = 1.5, linetype="longdash") +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.5947 / (k * T_ref))  # normalize to 1 at 20°C
    exp(-0.5947 / (k * T_K)) / ref_rate
  }, color = my_colours[1], size = 1.5) +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.6567 / (k * T_ref))
    exp(-0.6567 / (k * T_K)) / ref_rate
  }, color = my_colours[2], size = 1.5) +
  stat_function(fun = function(T_C) {
    T_K <- T_C + 273.15
    ref_rate <- exp(-0.5210/ (k * T_ref))
    - exp(-0.5210 / (k * T_K)) / ref_rate
  }, color = my_colours[3], size = 1.5) +
  labs(
    x = "Temperature (°C)",
    y = "",
    title = "Coral"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",x=5.5, y=6.8, height=3, alpha=1,fill = "black")+
  ylim(-8,8)+
  geom_text(data = labels_df_4c, aes(x = x, y = y, label = label, color = my_colours),
            hjust = 0, vjust = 1, size = 4, show.legend = FALSE, fontface = "bold") +
  scale_color_identity()



figure4 <- figure4a + figure4b + figure4c 

ggsave(figure4, filename = "./figures/figure4.png", dpi=700, width=11, height=4)



#topt vs Ea

lichen_topt_ea<- full_Ea %>%
  filter(dataset == "Lichen") %>%
  ggplot(aes(x = topt, y = e, colour = metabolic_category)) +
  geom_point() +
  theme_bw()+
  geom_smooth(method="lm")+
  scale_colour_manual(values=c("olivedrab1","#abbb80", "grey"))+
  facet_wrap(~metabolic_category)+
  theme(legend.position="none")+
  add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76",x=40, y=0.8, height=0.2, alpha=1,fill = "black")+
  xlab("Thermal optimum")+
  ylab("Ea")






coral_topt_ea<- full_Ea %>%
  filter(dataset == "Coral") %>%
  ggplot(aes(x = topt, y = e, colour = metabolic_category)) +
  geom_point() +
  scale_colour_manual(values=c("olivedrab1","#abbb80", "grey"))+
  theme_bw()+
  geom_smooth(method="lm")+
  facet_wrap(~metabolic_category) +
  theme(legend.position="none")+
  add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",x=35, y=3, height=0.5, alpha=1,fill = "black")+
  xlab("Thermal optimum")+
  ylab("Ea")
  



lichen_topt_ea / coral_topt_ea






