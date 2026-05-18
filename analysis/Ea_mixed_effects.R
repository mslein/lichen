
#loading r packages needed
#install.packages("pacman")
pacman::p_load(viridis, revtools, nlme, lme4, MuMIn, patchwork, tidyverse, ggh4x, dplyr, purrr, broom, tibble, rphylopic, officer, flextable, systemfonts, textshaping)


#reading in the data 
##lichen
lichen_tpc_Eas <- read_csv("analysis/tidy data/lichen_tpc_Eas.csv")
lichen_raw_npp <- read_csv("analysis/tidy data/lichen_npp_final.csv") %>% dplyr::select(tpc_grp, elevation_broad, lichen_type, latitude, species, temp, mol_gCmin)
lichen_raw_gpp <- read_csv("analysis/tidy data/lichen_gpp_final.csv") %>% dplyr::select(tpc_grp, elevation_broad, lichen_type, latitude, species, temp, mol_gCmin)
lichen_raw_r <- read_csv("analysis/tidy data/lichen_r_final.csv") %>% dplyr::select(tpc_grp, elevation_broad, lichen_type, latitude, species, temp, mol_gCmin)
#join the tpc parameter estimates w/ the raw data to get the additional information
gpp_lichen_Eas <- lichen_tpc_Eas %>% filter(metabolic_category == "gpp") %>%
  left_join(lichen_raw_gpp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE)  %>% drop_na(e) %>% 
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)),
         lichen_type=as.factor(lichen_type),
         elevation_broad=as.factor(elevation_broad))
npp_lichen_Eas <- lichen_tpc_Eas %>% filter(metabolic_category == "npp") %>%
  left_join(lichen_raw_npp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)),
         lichen_type=as.factor(lichen_type),
         elevation_broad=as.factor(elevation_broad))
r_lichen_Eas <- lichen_tpc_Eas %>% filter(metabolic_category == "r")  %>%
  left_join(lichen_raw_r, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)),
         lichen_type=as.factor(lichen_type),
         elevation_broad=as.factor(elevation_broad))

##coral
coral_tpc_Eas <- read_csv("analysis/tidy data/coral_tpc_Eas.csv") 
coral_raw_npp <- read_csv("analysis/tidy data/coral_npp_final.csv") %>% dplyr::select(tpc_grp, depth_broad, broad_coral, latitude, species, temp, mol_gCmin)
coral_raw_gpp <- read_csv("analysis/tidy data/coral_gpp_final.csv") %>% dplyr::select(tpc_grp, depth_broad, broad_coral, latitude, species, temp, mol_gCmin) 
coral_raw_r <- read_csv("analysis/tidy data/coral_r_final.csv") %>% dplyr::select(tpc_grp, depth_broad, broad_coral, latitude, species, temp, mol_gCmin)
#join the tpc parameter estimates w/ the raw data to get the additional information
gpp_coral_Eas <- coral_tpc_Eas %>% filter(metabolic_category == "gpp") %>% 
  left_join(coral_raw_gpp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)))
npp_coral_Eas <- coral_tpc_Eas %>% filter(metabolic_category == "npp") %>%
  left_join(coral_raw_npp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth,.keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)))
r_coral_Eas <- coral_tpc_Eas %>% filter(metabolic_category == "r") %>%
  left_join(coral_raw_r, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)))


# foraminifiera -- also known as cell
cell_tpc_Eas <- read_csv("analysis/tidy data/cell_tpc_Eas.csv") 
cell_raw_npp <- read_csv("analysis/tidy data/cell_npp_final.csv") %>% dplyr::select(tpc_grp, latitude, species, temp, mol_gCmin)
cell_raw_gpp <- read_csv("analysis/tidy data/cell_gpp_final.csv") %>% dplyr::select(tpc_grp, latitude, species, temp, mol_gCmin) 
cell_raw_r <- read_csv("analysis/tidy data/cell_r_final.csv") %>% dplyr::select(tpc_grp, latitude, species, temp, mol_gCmin)
#join the tpc parameter estimates w/ the raw data to get the additional information
gpp_cell_Eas <- cell_tpc_Eas %>% filter(metabolic_category == "gpp") %>% 
  left_join(cell_raw_gpp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)))
npp_cell_Eas <- cell_tpc_Eas %>% filter(metabolic_category == "npp") %>%
  left_join(cell_raw_npp, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth,.keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)))
r_cell_Eas <- cell_tpc_Eas %>% filter(metabolic_category == "r") %>%
  left_join(cell_raw_r, by="tpc_grp") %>% distinct(e, tpc_grp, topt, breadth, .keep_all = TRUE)  %>% separate(tpc_grp, into=c("response_id", "study_id"), remove = FALSE) %>% drop_na(e) %>%
  mutate(scaled_latitude = I(abs(latitude)- mean(abs(latitude))), 
         scaled_breadth=I(breadth - mean(breadth, na.rm = TRUE)), 
         scaled_topt=I(topt- mean(topt, na.rm = TRUE)))



#################################
# mixed effects models for lichen
#################################

#######
# temperature dependence 
#######
#setting correct reference levels
npp_lichen_Eas$lichen_type <- relevel(npp_lichen_Eas$lichen_type, ref = "green algae")
npp_lichen_Eas$elevation_broad <- relevel(npp_lichen_Eas$elevation_broad, ref = "neutral")
#model selection for NPP
#full model
lichen_npp <- lme(e ~ scaled_latitude + lichen_type + elevation_broad + scaled_breadth + scaled_topt, random= ~1|study_id, data=npp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_npp) #low correlation
#getting all possible combinations of models + comparing
lichen_npp_model_set <- dredge(lichen_npp, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
lichen_npp_model_top_reml <- lapply(get.models(lichen_npp_model_set, subset = delta < 2),
                          function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
lichen_npp_model_top_avg <- model.avg(lichen_npp_model_top_reml, full=TRUE)
summary(lichen_npp_model_top_avg)


#model selection for GPP
#setting correct reference levels
gpp_lichen_Eas$lichen_type <- relevel(gpp_lichen_Eas$lichen_type, ref = "green algae")
gpp_lichen_Eas$elevation_broad <- relevel(gpp_lichen_Eas$elevation_broad, ref = "neutral")
lichen_gpp <- lme(e ~ scaled_latitude +lichen_type + elevation_broad + scaled_breadth +scaled_topt, random= ~1|study_id, data=gpp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_gpp) # low correlation

#getting all possible combinations of models + comparing
lichen_gpp_model_set <- dredge(lichen_gpp, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
lichen_gpp_model_top_reml <- lapply(get.models(lichen_gpp_model_set, subset = delta < 2),
                                    function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
lichen_gpp_model_top_avg <- model.avg(lichen_gpp_model_top_reml, full=TRUE)
summary(lichen_gpp_model_top_avg)

#model selection for R (Sharpe Schoolfield)
#setting correct reference levels
r_lichen_Eas$lichen_type <- relevel(r_lichen_Eas$lichen_type, ref = "green algae")
r_lichen_Eas$elevation_broad <- relevel(r_lichen_Eas$elevation_broad, ref = "neutral")
lichen_r <- lme(e ~ scaled_latitude, random= ~1|study_id, data=r_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_r) # removed topt and breadth because they were highly colinear

lichen_r_tpc_final <- lme(e ~ scaled_latitude, random= ~1|study_id, data=r_lichen_Eas, weights = varFixed(~ e_se^2), method = "REML") 
summary(lichen_r_tpc_final)
intervals(lichen_r_tpc_final)




#model selection for R -- linearized Arrhenius
r_lichen_Eas$lichen_type <- relevel(r_lichen_Eas$lichen_type, ref = "green algae")
r_lichen_Eas$elevation_broad <- relevel(r_lichen_Eas$elevation_broad, ref = "neutral")
#very simple model -> could not fit random effects
lichen_r <- lm(e ~ scaled_latitude + scaled_breadth +scaled_topt, data=r_lichen_Eas) 
summary(lichen_r)

#arrhenius modelling for respiration
lichen_r_linear <- read_csv("analysis/tidy data/lichen_r_final.csv") %>%
  separate(tpc_grp,
    into = c("response_id", "study_id"),
    sep = "_",
    extra = "merge",
    fill = "right",
    remove = FALSE) %>%
  mutate(respiration_abs = abs(mol_gCmin),
    elevation_broad = as.factor(elevation_broad),
    lichen_type = as.factor(lichen_type),
    inv_T = 1 / ((8.617333262145e-5) * (temp + 273.15))) %>%
  filter(is.finite(respiration_abs),
    respiration_abs > 0,
    is.finite(temp)) %>%
  group_by(tpc_grp) %>%
  arrange(temp, .by_group = TRUE) %>%
  mutate(temp_at_max = temp[which.max(respiration_abs)]) %>%
  filter(temp <= temp_at_max) %>%   # keeps only the rising limb
  ungroup()

# Nest data by group
nested <- lichen_r_linear %>%
  group_by(tpc_grp) %>%
  nest()

# Fit linear model and extract slope + SE
group_ea <- nested %>%
  mutate(model = purrr::map(data, ~ lm(log(respiration_abs) ~ inv_T, data = .x)),
    tidied = purrr::map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "inv_T") %>%
  separate(tpc_grp, into = c("response_id", "study_id"), sep = "_", remove = FALSE) %>%
  rename(slope = estimate, Ea_SE = std.error) %>%
  dplyr::select(study_id, response_id, slope, Ea_SE)

group_ea_r <- left_join(group_ea, lichen_r_linear, by = c("study_id", "response_id"), relationship = "many-to-many") %>%
  dplyr::select(-tpc_grp.y) %>%
  rename(tpc_grp = tpc_grp.x) %>%
  distinct(slope, .keep_all = TRUE) %>%
  drop_na(slope, Ea_SE) %>%
  mutate(Ea = slope * -1,
    scaled_latitude = abs(latitude) - mean(abs(latitude), na.rm = TRUE),
    elevation_broad=as.factor(elevation_broad),
    lichen_type=as.factor(lichen_type))

group_ea_r$elevation_broad <- relevel(group_ea_r$elevation_broad, "neutral")
group_ea_r$lichen_type <- relevel(group_ea_r$lichen_type, "green algae")

# Mixed model with weights
lichen_ea_model <- lme(Ea ~ scaled_latitude + lichen_type + elevation_broad, random = ~1 | study_id, data = group_ea_r, weights = varFixed(~ Ea_SE^2), method = "ML")

# model selection
lichen_r_model_set <- dredge(lichen_ea_model, trace = TRUE)

# refit all tied models (AIC within 2) with REML
lichen_r_model_top_reml <- lapply(get.models(lichen_r_model_set, subset = delta < 2), function(m) update(m, method = "REML"))

# average top REML models
lichen_r_model_top_avg <- model.avg(lichen_r_model_top_reml, full = TRUE)
summary(lichen_r_model_top_avg)
confint(lichen_r_model_top_avg)


#######
# thermal breadth
#######

#npp
lichen_npp_breadth <- lme(breadth ~ scaled_latitude +lichen_type + elevation_broad, random= ~1|study_id, data=npp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_npp_breadth) #low correlation
#getting all possible combinations of models + comparing
lichen_npp_breadth_set <- dredge(lichen_npp_breadth, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
lichen_npp_breadth_reml <- lapply(get.models(lichen_npp_breadth_set, subset = delta < 2),
                                    function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of
#only 1 model that is better by 2 AIC
lichen_npp_breadth_top <- lichen_npp_breadth_reml[[1]]
summary(lichen_npp_breadth_top)


#gpp
lichen_gpp_breadth <- lme(breadth ~ scaled_latitude +lichen_type + elevation_broad, random= ~1|study_id, data=gpp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_gpp_breadth) #low correlation
#getting all possible combinations of models + comparing
lichen_gpp_breadth_set <- dredge(lichen_gpp_breadth, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
lichen_gpp_breadth_reml <- lapply(get.models(lichen_gpp_breadth_set, subset = delta < 2),
                                  function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
lichen_gpp_breadth_top <- lichen_gpp_breadth_reml[[1]]
summary(lichen_gpp_breadth_top)

#r
lichen_r_breadth <- lm(breadth ~ scaled_latitude +lichen_type + elevation_broad, data=r_lichen_Eas, na.action = na.fail) 
performance::check_collinearity(lichen_r_breadth) #low correlation
#getting all possible combinations of models + comparing
lichen_r_breadth_set <- dredge(lichen_r_breadth, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
lichen_r_breadth_reml <- lapply(get.models(lichen_gpp_breadth_set, subset = delta < 2),
                                  function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
lichen_r_breadth_top <- lichen_r_breadth_reml[[1]]
summary(lichen_r_breadth_top)



#######
# thermal optimum
#######
#npp
lichen_npp_topt <- lme(topt ~ scaled_latitude +lichen_type + elevation_broad, random= ~1|study_id, data=npp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_npp_topt) #low correlation
#getting all possible combinations of models + comparing
lichen_npp_topt_set <- dredge(lichen_npp_topt, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
lichen_npp_topt_reml <- lapply(get.models(lichen_npp_topt_set, subset = delta < 2),
                                  function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
lichen_npp_topt_top <- lichen_npp_topt_reml[[1]]
summary(lichen_npp_topt)


#gpp
lichen_gpp_topt <- lme(topt ~ scaled_latitude +lichen_type + elevation_broad, random= ~1|study_id, data=gpp_lichen_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(lichen_gpp_topt) #low correlation
#getting all possible combinations of models + comparing
lichen_gpp_topt_set <- dredge(lichen_gpp_topt, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
lichen_gpp_topt_reml <- lapply(get.models(lichen_gpp_topt_set, subset = delta < 2),
                                  function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
lichen_gpp_topt_top_avg <- model.avg(lichen_gpp_topt_reml, full=TRUE)
summary(lichen_gpp_topt_top_avg)

#r
lichen_r_topt <- lm(topt ~ scaled_latitude +lichen_type + elevation_broad, data=r_lichen_Eas, na.action = na.fail) 
performance::check_collinearity(lichen_r_topt) #low correlation
#getting all possible combinations of models + comparing
lichen_r_topt_set <- dredge(lichen_r_topt, trace = TRUE)
#averaging the REML fitted top models to get the full model average of 
lichen_r_topt_top <- get.models(lichen_r_topt_set, subset = delta < 2)[[1]]
summary(lichen_r_topt_top)

#################################
# mixed effects models for coral
#################################

#model selection for NPP
coral_npp <- lme(e ~ scaled_latitude + scaled_breadth + scaled_topt, random= ~1|study_id, data=npp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_npp)
#getting all possible combinations of models + comparing
coral_npp_model_set <- dredge(coral_npp, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_npp_model_top_reml <- lapply(get.models(coral_npp_model_set, subset = delta < 2),
                                    function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
coral_npp_model_top_avg <- model.avg(coral_npp_model_top_reml, full=TRUE)
summary(coral_npp_model_top_avg)

#model selection for GPP
coral_gpp <- lme(e ~  scaled_latitude + scaled_breadth + scaled_topt, random= ~1|study_id, data=gpp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_gpp)

#getting all possible combinations of models + comparing
coral_gpp_model_set <- dredge(coral_gpp, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_gpp_model_top_reml <- lapply(get.models(coral_gpp_model_set, subset = delta < 2),
                                   function(m) update(m, method = "REML"))
#only 1 model that is better by 2 AIC
coral_gpp_model_top <- coral_gpp_model_top_reml[[1]]
summary(coral_gpp_model_top)

#model selection for R
coral_r <- lme(e ~  scaled_latitude + scaled_breadth + scaled_topt, random= ~1|study_id, data=r_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_r) 

#getting all possible combinations of models + comparing
coral_r_model_set <- dredge(coral_r, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_r_model_top_reml <- lapply(get.models(coral_r_model_set, subset = delta < 2),
                                   function(m) update(m, method = "REML"))
#only 1 model that is better by 2 AIC
coral_r_model_top <- coral_r_model_top_reml[[1]]
summary(coral_r_model_top)



#######
# thermal breadth
#######

#npp
coral_npp_breadth <- lme(breadth ~ scaled_latitude, random= ~1|study_id, data=npp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_npp_breadth) #low correlation
#getting all possible combinations of models + comparing
coral_npp_breadth_set <- dredge(coral_npp_breadth, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_npp_breadth_reml <- lapply(get.models(coral_npp_breadth_set, subset = delta < 2),
                                  function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of
#only 1 model that is better by 2 AIC
coral_npp_breadth_top <- coral_npp_breadth_reml[[1]]
summary(coral_npp_breadth_top)


#gpp
coral_gpp_breadth <- lme(breadth ~ scaled_latitude, random= ~1|study_id, data=gpp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_gpp_breadth) #low correlation
#getting all possible combinations of models + comparing
coral_gpp_breadth_set <- dredge(coral_gpp_breadth, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_gpp_breadth_reml <- lapply(get.models(coral_gpp_breadth_set, subset = delta < 2),
                               function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
coral_gpp_breadth_top <- coral_gpp_breadth_reml[[1]]
summary(coral_gpp_breadth_top)

#r
coral_r_breadth <- lme(breadth ~ scaled_latitude, random= ~1|study_id, data=r_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_r_breadth) #low correlation
#getting all possible combinations of models + comparing
coral_r_breadth_set <- dredge(coral_r_breadth, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_r_breadth_reml <- lapply(get.models(coral_r_breadth_set, subset = delta < 2),
                                       function(m) update(m, method = "REML"))
#only 1 model that is better by 2 AIC
coral_r_breadth_top <- coral_r_breadth_reml[[1]]
summary(coral_r_breadth_top)


#######
# thermal optimum
#######
#npp
coral_npp_topt <- lme(topt ~ scaled_latitude, random= ~1|study_id, data=npp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_npp_topt) #low correlation
#getting all possible combinations of models + comparing
coral_npp_topt_set <- dredge(coral_npp_topt, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_npp_topt_reml <- lapply(get.models(coral_npp_topt_set, subset = delta < 2),
                               function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
coral_npp_topt <- coral_npp_topt_reml[[1]]
summary(coral_npp_topt)


#gpp
coral_gpp_topt <- lme(topt ~ scaled_latitude, random= ~1|study_id, data=gpp_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_gpp_topt) #low correlation
#getting all possible combinations of models + comparing
coral_gpp_topt_set <- dredge(coral_gpp_topt, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_gpp_topt_reml <- lapply(get.models(coral_gpp_topt_set, subset = delta < 2),
                               function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
coral_gpp_topt <- coral_gpp_topt_reml[[1]]
summary(coral_gpp_topt)

#r
coral_r_topt <- lme(topt ~ scaled_latitude, random= ~1|study_id, data=r_coral_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(coral_r_topt) #low correlation
#getting all possible combinations of models + comparing
coral_r_topt_set <- dredge(coral_r_topt, trace = TRUE)
#refitting all tied models (AIC within 2) with REML
coral_r_topt_reml <- lapply(get.models(coral_r_topt_set, subset = delta < 2),
                              function(m) update(m, method = "REML"))
#averaging the REML fitted top models to get the full model average of 
coral_r_topt <- coral_r_topt_reml[[1]]
summary(coral_r_topt)

#########################
###### single cell ######
#########################

## gpp Sharpe-Schoolfield

cell_gpp <- lme(e ~ scaled_latitude + scaled_breadth, random= ~1|study_id, data=gpp_cell_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(cell_gpp) #remove Topt because it had moderate correlation

# compare candidate models using ML
cell_gpp_model_set <- dredge(cell_gpp, trace = TRUE)

# refit top model(s) with REML
cell_gpp_model_top_reml <- lapply( get.models(cell_gpp_model_set, subset = delta < 2), function(m) update(m, method = "REML")
)

# if only one top model, use it directly
if (length(cell_gpp_model_top_reml) == 1) {
  cell_gpp_ss_final <- cell_gpp_model_top_reml[[1]]
} else {
  cell_gpp_ss_final <- model.avg(cell_gpp_model_top_reml, full = TRUE)
}

summary(cell_gpp_ss_final)
intervals(cell_gpp_ss_final, which = "fixed")


#arrhenius modelling for gpp
cell_gpp_linear <- read_csv("analysis/tidy data/cell_gpp_final.csv") %>%
  separate(tpc_grp, into = c("response_id", "study_id"), sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
  mutate(mol_gCmin = abs(mol_gCmin), inv_T = 1 / ((8.617333262145e-5) * (temp + 273.15))) %>%
  filter(is.finite(mol_gCmin), is.finite(temp)) %>%
  group_by(tpc_grp) %>%
  arrange(temp, .by_group = TRUE) %>%
  mutate(temp_at_max = temp[which.max(mol_gCmin)]) %>%
  filter(temp <= temp_at_max) %>%
  ungroup()

nested <- cell_gpp_linear %>%
  group_by(tpc_grp) %>%
  nest()

group_ea <- nested %>%
  mutate(model = purrr::map(data, ~ lm(log(mol_gCmin) ~ inv_T, data = .x)),
    tidied = purrr::map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "inv_T") %>%
  separate(tpc_grp, into = c("response_id","study_id"), sep = "_", remove = FALSE) %>%
  rename(slope = estimate, Ea_SE = std.error) %>%
  dplyr::select(study_id, response_id, slope, Ea_SE)

group_ea_cell_gpp <- left_join(group_ea, cell_gpp_linear, by = c("study_id", "response_id"), relationship = "many-to-many") %>%
  dplyr::select(-tpc_grp.y, -temp_at_max) %>%
  rename(tpc_grp = tpc_grp.x) %>%
  distinct(slope, .keep_all = TRUE) %>%
  drop_na(slope, Ea_SE) %>%
  mutate(Ea = slope * -1, scaled_latitude = abs(latitude) - mean(abs(latitude), na.rm = TRUE))

cell_ea_model <- lme(Ea ~ scaled_latitude, random = ~1 | study_id, data = group_ea_cell_gpp, weights = varFixed(~ Ea_SE^2), method = "ML")

cell_gpp_model_set <- dredge(cell_ea_model, trace = TRUE)

cell_gpp_model_top_reml <- lapply(
  get.models(cell_gpp_model_set, subset = delta < 2),
  function(m) update(m, method = "REML")
)

if (length(cell_gpp_model_top_reml) == 1) {
  cell_gpp_final <- cell_gpp_model_top_reml[[1]]
} else {
  cell_gpp_final <- model.avg(cell_gpp_model_top_reml, full = TRUE)
}

summary(cell_gpp_final) #0.6395397
intervals(cell_gpp_final)

## npp Sharpe-Schoolfield

cell_npp <- lme(e ~ 1, random= ~1|study_id, data=npp_cell_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(cell_npp) # all terms were relatively correlated, move to intercept only

# compare candidate models using ML
cell_npp_model_set <- dredge(cell_npp, trace = TRUE)

# refit top model(s) with REML
cell_npp_model_top_reml <- lapply( get.models(cell_npp_model_set, subset = delta < 2), function(m) update(m, method = "REML"))

# if only one top model, use it directly
if (length(cell_npp_model_top_reml) == 1) {
  cell_npp_ss_final <- cell_npp_model_top_reml[[1]]
} else {
  cell_npp_ss_final <- model.avg(cell_npp_model_top_reml, full = TRUE)
}

summary(cell_npp_ss_final)
intervals(cell_npp_ss_final)

#arrhenius modelling for npp
cell_npp_linear <- read_csv("analysis/tidy data/cell_npp_final.csv") %>%
  separate(tpc_grp, into = c("response_id", "study_id"), sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
  mutate(mol_gCmin = abs(mol_gCmin), inv_T = 1 / ((8.617333262145e-5) * (temp + 273.15))) %>%
  filter(is.finite(mol_gCmin), is.finite(temp)) %>%
  group_by(tpc_grp) %>%
  arrange(temp, .by_group = TRUE) %>%
  mutate(temp_at_max = temp[which.max(mol_gCmin)]) %>%
  filter(temp <= temp_at_max) %>%
  ungroup()

nested <- cell_npp_linear %>%
  group_by(tpc_grp) %>%
  nest()

group_ea <- nested %>%
  mutate(model = purrr::map(data, ~ lm(log(mol_gCmin) ~ inv_T, data = .x)),
    tidied = purrr::map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "inv_T") %>%
  separate(tpc_grp, into = c("response_id","study_id"), sep = "_", remove = FALSE) %>%
  rename(slope = estimate, Ea_SE = std.error) %>%
  dplyr::select(study_id, response_id, slope, Ea_SE)

group_ea_cell_npp <- left_join(
  group_ea,
  cell_npp_linear,
  by = c("study_id", "response_id"),
  relationship = "many-to-many") %>%
  dplyr::select(-tpc_grp.y, -temp_at_max) %>%
  rename(tpc_grp = tpc_grp.x) %>%
  distinct(slope, .keep_all = TRUE) %>%
  drop_na(slope, Ea_SE) %>%
  mutate(Ea = slope * -1,
    scaled_latitude = abs(latitude) - mean(abs(latitude), na.rm = TRUE))

cell_ea_model <- lme(Ea ~ scaled_latitude, random = ~1 | study_id, data = group_ea_cell_npp, weights = varFixed(~ Ea_SE^2), method = "ML")

cell_npp_model_set <- dredge(cell_ea_model, trace = TRUE)

cell_npp_model_top_reml <- lapply(
  get.models(cell_npp_model_set, subset = delta < 2),
  function(m) update(m, method = "REML")
)

if (length(cell_npp_model_top_reml) == 1) {
  cell_npp_final <- cell_npp_model_top_reml[[1]]
} else {
  cell_npp_final <- model.avg(cell_npp_model_top_reml, full = TRUE)
}

summary(cell_npp_final) #1.086079
intervals(cell_npp_final, which = "fixed")



## r Sharpe-Schoolfield

cell_r <- lme(e ~ scaled_latitude + scaled_breadth, random= ~1|study_id, data=r_cell_Eas, weights = varFixed(~ e_se^2), method = "ML") 
performance::check_collinearity(cell_r) # removed topt because of moderate correlation

# compare candidate models using ML
cell_r_model_set <- dredge(cell_r, trace = TRUE)

# refit top model(s) with REML
cell_r_model_top_reml <- lapply( get.models(cell_r_model_set, subset = delta < 2), function(m) update(m, method = "REML"))

# if only one top model, use it directly
if (length(cell_r_model_top_reml) == 1) {
  cell_r_ss_final <- cell_r_model_top_reml[[1]]
} else {
  cell_r_ss_final <- model.avg(cell_r_model_top_reml, full = TRUE)
}

summary(cell_r_ss_final)
intervals(cell_r_ss_final)

#arrhenius modelling for respiration
cell_r_linear <- read_csv("analysis/tidy data/cell_r_final.csv") %>%
  separate(tpc_grp, into = c("response_id", "study_id"), sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
  mutate(mol_gCmin = abs(mol_gCmin),
    inv_T = 1 / ((8.617333262145e-5) * (temp + 273.15))) %>%
  filter(is.finite(mol_gCmin), is.finite(temp)) %>%
  group_by(tpc_grp) %>%
  arrange(temp, .by_group = TRUE) %>%
  mutate(temp_at_max = temp[which.max(mol_gCmin)]) %>%
  filter(temp <= temp_at_max) %>%
  ungroup()

nested <- cell_r_linear %>%
  group_by(tpc_grp) %>%
  nest()

group_ea <- nested %>%
  mutate(model = purrr::map(data, ~ lm(log(mol_gCmin) ~ inv_T, data = .x)),
    tidied = purrr::map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "inv_T") %>%
  separate(tpc_grp, into = c("response_id","study_id"), sep = "_", remove = FALSE) %>%
  rename(slope = estimate, Ea_SE = std.error) %>%
  dplyr::select(study_id, response_id, slope, Ea_SE)

group_ea_cell_r <- left_join(group_ea, cell_r_linear, by = c("study_id", "response_id"), relationship = "many-to-many") %>%
  dplyr::select(-tpc_grp.y, -temp_at_max) %>%
  rename(tpc_grp = tpc_grp.x) %>%
  distinct(slope, .keep_all = TRUE) %>%
  drop_na(slope, Ea_SE) %>%
  mutate(Ea = slope * -1,
    scaled_latitude = abs(latitude) - mean(abs(latitude), na.rm = TRUE))

cell_ea_model <- lme(Ea ~ scaled_latitude, random = ~1 | study_id, data = group_ea_cell_r, weights = varFixed(~ Ea_SE^2), method = "ML")

cell_r_model_set <- dredge(cell_ea_model, trace = TRUE)

refit_top_models_la <- function(model_set) {
  top_reml <- lapply(
    get.models(model_set, subset = delta < 2),
    function(m) update(m, method = "REML")
  )
  
  if (length(top_reml) == 1) {
    return(top_reml[[1]])
  } else {
    return(model.avg(top_reml, full = TRUE))
  }
}

cell_r_final <- refit_top_models_la(cell_r_model_set)
summary(cell_r_final) #0.6595894
intervals(cell_r_final)


#######
# thermal breadth
#######

# npp
cell_npp_breadth <- lme(
  breadth ~ scaled_latitude,
  random = ~1 | study_id,
  data = npp_cell_Eas,
  weights = varFixed(~ e_se^2),
  method = "ML"
)
performance::check_collinearity(cell_npp_breadth) # low correlation

# getting all possible combinations of models + comparing
cell_npp_breadth_set <- dredge(cell_npp_breadth, trace = TRUE)


# refit that exact fixed-effects structure from scratch with REML
cell_npp_breadth_top <- lm(breadth ~ scaled_latitude,
  data    = npp_cell_Eas,
  weights = 1 / (e_se^2),
  na.action = na.omit)

summary(cell_npp_breadth_top) #21.5580 


# gpp
cell_gpp_breadth <- lme(
  breadth ~ scaled_latitude,
  random = ~1 | study_id,
  data = gpp_cell_Eas,
  weights = varFixed(~ e_se^2),
  method = "ML"
)
performance::check_collinearity(cell_gpp_breadth) # low correlation

# getting all possible combinations of models + comparing
cell_gpp_breadth_set <- dredge(cell_gpp_breadth, trace = TRUE)

cell_gpp_breadth_top <- lm(breadth ~ scaled_latitude,
                           data    = gpp_cell_Eas,
                           weights = 1 / (e_se^2),
                           na.action = na.omit)

summary(cell_gpp_breadth_top) #21.4590




# r
cell_r_breadth <- lme(
  breadth ~ scaled_latitude,
  random = ~1 | study_id,
  data = r_cell_Eas,
  weights = varFixed(~ e_se^2),
  method = "ML"
)
performance::check_collinearity(cell_r_breadth) # low correlation

# getting all possible combinations of models + comparing
cell_r_breadth_set <- dredge(cell_r_breadth, trace = TRUE)

cell_r_breadth_top <- lm(breadth ~ scaled_latitude,
                           data    = r_cell_Eas,
                           weights = 1 / (e_se^2),
                           na.action = na.omit)

summary(cell_r_breadth_top) #22.8635


#######
# thermal optimum
#######

# npp
cell_npp_topt <- lme(
  topt ~ scaled_latitude,
  random = ~1 | study_id,
  data = npp_cell_Eas,
  weights = varFixed(~ e_se^2),
  method = "ML"
)
performance::check_collinearity(cell_npp_topt) # low correlation

# getting all possible combinations of models + comparing
cell_npp_topt_set <- dredge(cell_npp_topt, trace = TRUE)


cell_npp_topt_top <- lm(topt ~ scaled_latitude,
                         data    = npp_cell_Eas,
                         weights = 1 / (e_se^2),
                         na.action = na.omit)

summary(cell_npp_topt_top) #28.4026





# gpp
cell_gpp_topt <- lme(
  topt ~ scaled_latitude,
  random = ~1 | study_id,
  data = gpp_cell_Eas,
  weights = varFixed(~ e_se^2),
  method = "ML"
)
performance::check_collinearity(cell_gpp_topt) # low correlation

# getting all possible combinations of models + comparing
cell_gpp_topt_set <- dredge(cell_gpp_topt, trace = TRUE)

cell_gpp_topt_top <- lm(topt ~ scaled_latitude,
                        data    = gpp_cell_Eas,
                        weights = 1 / (e_se^2),
                        na.action = na.omit)

summary(cell_gpp_topt_top) #22.6596

# r
cell_r_topt <- lme(
  topt ~ scaled_latitude,
  random = ~1 | study_id,
  data = r_cell_Eas,
  weights = varFixed(~ e_se^2),
  method = "ML"
)
performance::check_collinearity(cell_r_topt) # low correlation

# getting all possible combinations of models + comparing
cell_r_topt_set <- dredge(cell_r_topt, trace = TRUE)

cell_r_topt_top <- lm(topt ~ scaled_latitude,
                        data    = r_cell_Eas,
                        weights = 1 / (e_se^2),
                        na.action = na.omit)

summary(cell_r_topt_top) #38.2453


#### figure 3 --> main figure (many parts)
# part 1: insets

#trying to just plot the Eas
npp_lichen_ea  <- npp_lichen_Eas %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "npp")
gpp_lichen_ea  <- gpp_lichen_Eas %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "gpp")
r_lichen_ea <- group_ea_r %>% mutate(e=slope*-1, eh="", topt="", breadth="", rmax="") %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "r")

npp_coral_ea  <- npp_coral_Eas %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "npp") 
gpp_coral_ea  <- gpp_coral_Eas %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "gpp")
r_coral_ea  <- r_coral_Eas %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "r")

npp_cell_ea <- group_ea_cell_npp %>% mutate(e=slope*-1, eh="", topt="", breadth="", rmax="") %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "npp")
gpp_cell_ea <- group_ea_cell_gpp %>% mutate(e=slope*-1, eh="", topt="", breadth="", rmax="") %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "gpp")
r_cell_ea <- group_ea_cell_r %>% mutate(e=slope*-1, eh="", topt="", breadth="", rmax="") %>% dplyr::select(e, eh, topt, breadth, rmax, study_id) %>% mutate(metabolic_category = "r")


lichen_all <- rbind(npp_lichen_ea, gpp_lichen_ea, r_lichen_ea) %>% mutate(dataset="lichen")
coral_all <- rbind(npp_coral_ea, gpp_coral_ea, r_coral_ea) %>% mutate(dataset="coral")
cell_all <- rbind(npp_cell_ea, gpp_cell_ea, r_cell_ea) %>% mutate(dataset="cell")
lit <- data.frame(e=c("", "", ""), eh=c("", "", ""), topt=c("", "", ""), breadth =c("", "", ""), rmax =c("", "", ""), study_id=c("", "", ""), metabolic_category=c("gpp", "r", "npp"), dataset=c("López-Urrutia et al. 2006", "López-Urrutia et al. 2006", "López-Urrutia et al. 2006"))

full_Ea <- rbind(lichen_all, coral_all, cell_all, lit) %>% mutate(e=as.numeric(e), topt=as.numeric(topt))



full_topt_breadth <- rbind(lichen_all, coral_all) %>% mutate(e=as.numeric(e), topt=as.numeric(topt), breadth=as.numeric(breadth), 
                                                             metabolic_category=as.factor(metabolic_category), dataset=as.factor(dataset))


#### figure 3 --> estimate + slope figures
# estimate figure: row 1 = GPP vs R, row 2 = NPP
# slope/curve figure: row 1 = GPP vs R, row 2 = NPP

# ============================================================
# 0) fitted cell model object names
# ============================================================

cell_gpp_model_obj <- cell_gpp_final
cell_npp_model_obj <- cell_npp_final
cell_r_model_obj   <- cell_r_final

# ============================================================
# 1) raw Ea values
# ============================================================

raw_ea_all <- bind_rows(
  npp_lichen_Eas %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "npp", dataset = "lichen"),
  
  gpp_lichen_Eas %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "gpp", dataset = "lichen"),
  
  group_ea_r %>%
    mutate(e = slope * -1, eh = NA, topt = NA, breadth = NA, rmax = NA) %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "r", dataset = "lichen"),
  
  npp_coral_Eas %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "npp", dataset = "coral"),
  
  gpp_coral_Eas %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "gpp", dataset = "coral"),
  
  r_coral_Eas %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "r", dataset = "coral"),
  
  group_ea_cell_npp %>%
    mutate(e = slope * -1, eh = NA, topt = NA, breadth = NA, rmax = NA) %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "npp", dataset = "cell"),
  
  group_ea_cell_gpp %>%
    mutate(e = slope * -1, eh = NA, topt = NA, breadth = NA, rmax = NA) %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "gpp", dataset = "cell"),
  
  group_ea_cell_r %>%
    mutate(e = slope * -1, eh = NA, topt = NA, breadth = NA, rmax = NA) %>%
    dplyr::select(e, eh, topt, breadth, rmax, study_id) %>%
    mutate(metabolic_category = "r", dataset = "cell")
) %>%
  mutate(
    e = as.numeric(e),
    topt = as.numeric(topt),
    breadth = as.numeric(breadth)
  )

# ============================================================
# 2) coefficient extraction
# ============================================================

get_coefs_CI <- function(mod, metabolic_category, dataset) {
  
  sm <- summary(mod)
  
  if (inherits(mod, "averaging")) {
    coefs <- as.data.frame(sm$coefmat.full)
    est_col <- "Estimate"
    se_col <- if ("Std. Error" %in% names(coefs)) "Std. Error" else "SE"
  } else if (inherits(mod, "lme")) {
    coefs <- as.data.frame(sm$tTable)
    est_col <- "Value"
    se_col <- "Std.Error"
  } else {
    stop("Unsupported model class: ", paste(class(mod), collapse = ", "))
  }
  
  coefs %>%
    tibble::rownames_to_column("term") %>%
    transmute(
      term,
      estimate = .data[[est_col]],
      se = .data[[se_col]],
      lower_95 = estimate - 1.96 * se,
      upper_95 = estimate + 1.96 * se,
      metabolic_category = metabolic_category,
      dataset = dataset
    )
}

all_coefs <- bind_rows(
  get_coefs_CI(lichen_gpp_model_top_avg, "gpp", "lichen"),
  get_coefs_CI(lichen_npp_model_top_avg, "npp", "lichen"),
  get_coefs_CI(lichen_r_model_top_avg,   "r",   "lichen"),
  
  get_coefs_CI(coral_gpp_model_top,      "gpp", "coral"),
  get_coefs_CI(coral_npp_model_top_avg,  "npp", "coral"),
  get_coefs_CI(coral_r_model_top,        "r",   "coral"),
  
  get_coefs_CI(cell_gpp_model_obj,       "gpp", "cell"),
  get_coefs_CI(cell_npp_model_obj,       "npp", "cell"),
  get_coefs_CI(cell_r_model_obj,         "r",   "cell"),
  
  tibble(
    estimate = c(0.33, 0.56, 0.29),
    se = c(0.089, 0.024, 0.036),
    metabolic_category = c("gpp", "r", "npp"),
    dataset = "literature"
  ) %>%
    mutate(
      term = "(Intercept)",
      lower_95 = estimate - 1.96 * se,
      upper_95 = estimate + 1.96 * se
    ) %>%
    dplyr::select(term, estimate, se, lower_95, upper_95, metabolic_category, dataset)
)

# ============================================================
# 3) sample size labels
# ============================================================

n_labels <- tibble::tribble(
  ~dataset,  ~metabolic_category, ~label,
  "lichen",  "gpp",               "n = 43",
  "lichen",  "npp",               "n = 53",
  "lichen",  "r",                 "n = 128",
  "coral",   "gpp",               "n = 46",
  "coral",   "npp",               "n = 38",
  "coral",   "r",                 "n = 56",
  "cell",    "gpp",               "n = 6",
  "cell",    "npp",               "n = 10",
  "cell",    "r",                 "n = 6"
)

# ============================================================
# 4) aesthetics
# ============================================================

proc_cols <- c(
  "gpp" = "yellowgreen",
  "r"   = "grey45",
  "npp" = "lightgoldenrod4"
)

proc_labels <- c(
  "gpp" = "GPP",
  "r"   = "R",
  "npp" = "NPP"
)

panel_titles <- c(
  "lichen" = "Lichen",
  "coral"  = "Coral",
  "cell"   = "Foraminifera"
)

fig_w <- 16
fig_h <- 10

phylo_uuids <- c(
  "lichen" = "a208bba4-f4bf-4810-bcc9-c5868836fc76",
  "coral"  = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",
  "cell"   = "c747a673-9e45-42bb-a657-7cbf20f0e848"
)

# Curve figure silhouettes
phylopic_specs <- tibble::tribble(
  ~system,   ~uuid,                                     ~x,   ~y,   ~height,
  "lichen",  "a208bba4-f4bf-4810-bcc9-c5868836fc76",   35, 0.30, 0.36,
  "coral",   "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",   35, 0.30, 0.36,
  "cell",    "c747a673-9e45-42bb-a657-7cbf20f0e848",   35, 0.30, 0.36
)

make_silhouette_inset <- function(system_name) {
  ggplot() +
    add_phylopic(
      uuid = phylo_uuids[[system_name]],
      x = 0.5,
      y = 0.5,
      height = 0.70,
      alpha = 1,
      fill = "black"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

# ============================================================
# 5) estimate panel builder
# ============================================================

make_ea_estimate_panel <- function(dataset_name,
                                   categories = c("gpp", "r"),
                                   y_limits = c(-0.5, 3.0),
                                   title = NULL,
                                   tag = NULL,
                                   show_ea_axis = TRUE,
                                   show_process_axis = TRUE,
                                   add_silhouette = TRUE) {
  
  cats <- categories
  cat_positions <- setNames(seq(1, by = 0.62, length.out = length(cats)), cats)
  
  raw_df <- raw_ea_all %>%
    filter(dataset == dataset_name, metabolic_category %in% cats, !is.na(e)) %>%
    mutate(
      metabolic_category = factor(metabolic_category, levels = cats),
      x_num = cat_positions[as.character(metabolic_category)] + 0.10
    )
  
  mod_df <- all_coefs %>%
    filter(dataset == dataset_name, term == "(Intercept)", metabolic_category %in% cats) %>%
    mutate(
      metabolic_category = factor(metabolic_category, levels = cats),
      x_num = cat_positions[as.character(metabolic_category)] + 0.10
    )
  
  lit_df <- all_coefs %>%
    filter(dataset == "literature", term == "(Intercept)", metabolic_category %in% cats) %>%
    mutate(
      metabolic_category = factor(metabolic_category, levels = cats),
      x_num = cat_positions[as.character(metabolic_category)] - 0.10
    )
  
  lab_df <- n_labels %>%
    filter(dataset == dataset_name, metabolic_category %in% cats) %>%
    mutate(
      metabolic_category = factor(metabolic_category, levels = cats),
      x_num = cat_positions[as.character(metabolic_category)] + 0.10,
      y = 2.52
    )
  
  p <- ggplot() +
    geom_jitter(
      data = raw_df,
      aes(x = x_num, y = e, colour = metabolic_category),
      alpha = 0.15,
      width = 0.018,
      height = 0.12,
      size = 3.0
    ) +
    geom_pointrange(
      data = lit_df,
      aes(x = x_num, y = estimate, ymin = lower_95, ymax = upper_95,
          colour = metabolic_category),
      size = 1.65,
      linewidth = 2.0,
      alpha = 0.55
    ) +
    geom_point(
      data = lit_df,
      aes(x = x_num, y = estimate),
      colour = "white",
      size = 3.3
    ) +
    geom_pointrange(
      data = mod_df,
      aes(x = x_num, y = estimate, ymin = lower_95, ymax = upper_95,
          colour = metabolic_category),
      size = 1.65,
      linewidth = 2.0
    ) +
    geom_label(
      data = lab_df,
      aes(x = x_num, y = y, label = label, fill = metabolic_category),
      colour = "white",
      size = 5.8,
      label.size = 0
    ) +
    coord_flip(clip = "off") +
    scale_x_continuous(
      breaks = unname(cat_positions[cats]),
      labels = proc_labels[cats],
      limits = c(min(cat_positions) - 0.35, max(cat_positions) + 0.35)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = c(0, 1, 2, 3)
    ) +
    scale_colour_manual(values = proc_cols[cats]) +
    scale_fill_manual(values = proc_cols[cats]) +
    labs(
      title = title,
      tag = tag,
      x = "",
      y = if (show_ea_axis) "Temperature dependence\n(Ea)" else ""
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 20),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold"),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 25, face = "bold"),
      plot.tag.position = c(0.02, 0.98),
      plot.margin = margin(8, 14, 8, 8)
    )
  
  if (!show_ea_axis) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  if (!show_process_axis) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      )
  }
  
  if (add_silhouette) {
    p <- p +
      patchwork::inset_element(
        make_silhouette_inset(dataset_name),
        left = 0.78,
        bottom = 0.07,
        right = 0.98,
        top = 0.27,
        align_to = "panel",
        clip = FALSE
      )
  }
  
  p
}

# ============================================================
# 5B) observed temperature bounds for trimming curves
# ============================================================

standardize_temp_bounds <- function(file_path, system_name, process_name) {
  read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::select(any_of(c("tpc_grp", "study_id", "temp"))) %>%
    mutate(
      system = system_name,
      process = process_name
    ) %>%
    mutate(
      tpc_grp = if ("tpc_grp" %in% names(.)) as.character(tpc_grp) else NA_character_,
      study_id = if ("study_id" %in% names(.)) as.character(study_id) else NA_character_,
      tpc_grp = dplyr::coalesce(tpc_grp, study_id)
    ) %>%
    filter(!is.na(tpc_grp), !is.na(temp)) %>%
    group_by(system, process, tpc_grp) %>%
    summarise(
      Tmin_obs = min(temp, na.rm = TRUE),
      Tmax_obs = max(temp, na.rm = TRUE),
      .groups = "drop"
    )
}

study_temp_bounds <- bind_rows(
  standardize_temp_bounds("analysis/tidy data/lichen_gpp_final.csv", "lichen", "gpp"),
  standardize_temp_bounds("analysis/tidy data/lichen_r_final.csv",   "lichen", "r"),
  standardize_temp_bounds("analysis/tidy data/lichen_npp_final.csv", "lichen", "npp"),
  standardize_temp_bounds("analysis/tidy data/coral_gpp_final.csv",  "coral", "gpp"),
  standardize_temp_bounds("analysis/tidy data/coral_r_final.csv",    "coral", "r"),
  standardize_temp_bounds("analysis/tidy data/coral_npp_final.csv",  "coral", "npp"),
  standardize_temp_bounds("analysis/tidy data/cell_gpp_final.csv",   "cell", "gpp"),
  standardize_temp_bounds("analysis/tidy data/cell_r_final.csv",     "cell", "r"),
  standardize_temp_bounds("analysis/tidy data/cell_npp_final.csv",   "cell", "npp")
)

system_temp_bounds <- study_temp_bounds %>%
  group_by(system, process) %>%
  summarise(
    Tmin_obs = min(Tmin_obs, na.rm = TRUE),
    Tmax_obs = max(Tmax_obs, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 6) normalized meta curves
# ============================================================

anchor_temps <- tibble::tribble(
  ~system,   ~Tref,
  "lichen",  15,
  "coral",   25,
  "cell",    25
)

topt_cutoffs <- tibble(
  system    = c("lichen", "lichen", "lichen",
                "coral",  "coral",  "coral",
                "cell",   "cell",   "cell"),
  process   = c("gpp",    "npp",    "r",
                "gpp",    "npp",    "r",
                "gpp",    "npp",    "r"),
  topt_mean = c(24.9,     22.4,     24.9,
                31.2,     29.8,     30.7,
                22.65,    28.40,    38.2453)
)

get_Ea <- function(model, term = "(Intercept)") {
  if (inherits(model, "averaging")) {
    cm <- summary(model)$coefmat.full
    est <- cm[term, "Estimate"]
    seCol <- grep("SE|Std", colnames(cm), value = TRUE)[1]
    se <- cm[term, seCol]
  } else if (inherits(model, "lme")) {
    tt <- summary(model)$tTable
    est <- tt[term, "Value"]
    se <- tt[term, "Std.Error"]
  } else {
    stop("Unsupported model class: ", paste(class(model), collapse = ", "))
  }
  
  list(Ea = as.numeric(est), Ea_se = as.numeric(se))
}

k_B <- 8.617e-5

arrhenius_anchor <- function(T, Ea, Tref) {
  T_K <- T + 273.15
  Tref_K <- Tref + 273.15
  exp(-Ea / (k_B * T_K)) / exp(-Ea / (k_B * Tref_K))
}

curve_spec <- tibble::tribble(
  ~system,  ~process, ~model,                    ~Ea_term,      ~Tmin, ~Tmax,
  "lichen", "gpp",    lichen_gpp_model_top_avg, "(Intercept)",  -8,    50,
  "lichen", "r",      lichen_r_model_top_avg,   "(Intercept)",  -8,    50,
  "lichen", "npp",    lichen_npp_model_top_avg, "(Intercept)",  -8,    50,
  "coral",  "gpp",    coral_gpp_model_top,      "(Intercept)",   6,    37,
  "coral",  "r",      coral_r_model_top,        "(Intercept)",   5,    39,
  "coral",  "npp",    coral_npp_model_top_avg,  "(Intercept)",  10,    37,
  "cell",   "gpp",    cell_gpp_model_obj,       "(Intercept)",   5,    45,
  "cell",   "r",      cell_r_model_obj,         "(Intercept)",   5,    45,
  "cell",   "npp",    cell_npp_model_obj,       "(Intercept)",   5,    45
)

make_meta_curve <- function(model, Ea_term, Tmin, Tmax, Tref, n = 200) {
  pars <- get_Ea(model, Ea_term)
  Ea <- pars$Ea
  Tseq <- seq(Tmin, Tmax, length.out = n)
  
  tibble(
    temp = Tseq,
    rate = arrhenius_anchor(Tseq, Ea, Tref = Tref),
    Ea = Ea
  )
}

curve_df_3x3 <- curve_spec %>%
  left_join(anchor_temps, by = "system") %>%
  left_join(topt_cutoffs, by = c("system", "process")) %>%
  left_join(system_temp_bounds, by = c("system", "process")) %>%
  mutate(
    topt_mean = ifelse(is.na(topt_mean), Inf, topt_mean),
    Tmin_obs = ifelse(is.na(Tmin_obs), -Inf, Tmin_obs),
    curve = pmap(
      list(model, Ea_term, Tmin, Tmax, Tref),
      ~ make_meta_curve(..1, ..2, ..3, ..4, ..5, n = 200)
    )
  ) %>%
  select(system, process, Tref, topt_mean, Tmin_obs, curve) %>%
  unnest(curve) %>%
  filter(temp >= Tmin_obs, temp <= topt_mean) %>%
  mutate(
    system = factor(system, levels = c("lichen", "coral", "cell")),
    process = factor(process, levels = c("gpp", "r", "npp"))
  )

# ============================================================
# 8) estimate figure
# ============================================================

estimate_lichen_gr <- make_ea_estimate_panel(
  "lichen", c("gpp", "r"),
  title = "Lichen",
  tag = "a",
  show_ea_axis = FALSE,
  show_process_axis = TRUE
)

estimate_coral_gr <- make_ea_estimate_panel(
  "coral", c("gpp", "r"),
  title = "Coral",
  tag = "b",
  show_ea_axis = FALSE,
  show_process_axis = TRUE
)

estimate_cell_gr <- make_ea_estimate_panel(
  "cell", c("gpp", "r"),
  title = "Foraminifera",
  tag = "c",
  show_ea_axis = FALSE,
  show_process_axis = TRUE
)

estimate_lichen_npp <- make_ea_estimate_panel(
  "lichen", c("npp"),
  title = NULL,
  tag = "d",
  show_ea_axis = TRUE,
  show_process_axis = TRUE
)

estimate_coral_npp <- make_ea_estimate_panel(
  "coral", c("npp"),
  title = NULL,
  tag = "e",
  show_ea_axis = TRUE,
  show_process_axis = TRUE
)

estimate_cell_npp <- make_ea_estimate_panel(
  "cell", c("npp"),
  title = NULL,
  tag = "f",
  show_ea_axis = TRUE,
  show_process_axis = TRUE
)

fig_estimates <-
  (estimate_lichen_gr + estimate_coral_gr + estimate_cell_gr) /
  (estimate_lichen_npp + estimate_coral_npp + estimate_cell_npp)

fig_estimates

ggsave(
  "figures/fig_estimates_GPP_R_NPP_with_matched_inset_silhouettes.png",
  fig_estimates,
  width = fig_w,
  height = fig_h,
  dpi = 600
)

# ============================================================
# 9) slope/curve panel builder
# ============================================================

make_curve_panel <- function(system_name,
                             processes = c("gpp", "r"),
                             tag,
                             show_y_axis = TRUE,
                             show_temp_axis = TRUE,
                             title = NULL,
                             show_legend = TRUE) {
  
  cdf <- curve_df_3x3 %>%
    filter(system == system_name, process %in% processes)
  
  panel_tref <- unique(cdf$Tref)[1]
  
  phylo_row <- phylopic_specs %>%
    filter(system == system_name)
  
  p <- ggplot() +
    geom_hline(yintercept = 1, linewidth = 0.5, linetype = 3) +
    geom_line(
      data = cdf,
      aes(x = temp, y = rate, colour = process),
      linewidth = 2.5,
      alpha = 1
    ) +
    coord_cartesian(xlim = c(-10, 42), ylim = c(0, 2.45)) +
    scale_x_continuous(breaks = c(-10, 0, 10, 20, 30, 40)) +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
    scale_colour_manual(
      values = proc_cols[processes],
      breaks = processes,
      labels = proc_labels[processes]
    ) +
    labs(
      title = title,
      tag = tag,
      x = if (show_temp_axis) "Temperature (°C)" else "",
      y = if (show_y_axis) paste0("Slope-aligned response\n(centred on ", panel_tref, " °C)") else "",
      colour = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = if (show_legend) c(0.17, 0.84) else "none",
      legend.background = element_blank(),
      legend.text = element_text(size = 18),
      axis.text = element_text(size = 20),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold"),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 25, face = "bold"),
      plot.tag.position = c(0.02, 0.98),
      plot.margin = margin(8, 14, 8, 8)
    )
  
  if (!show_temp_axis) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  if (!show_y_axis) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      )
  }
  
  if (nrow(phylo_row) == 1) {
    p <- p +
      add_phylopic(
        uuid = phylo_row$uuid,
        x = phylo_row$x,
        y = phylo_row$y,
        height = phylo_row$height,
        alpha = 1,
        fill = "black"
      )
  }
  
  p
}

# ============================================================
# 10) slope/curve figure
# ============================================================

curve_lichen_gr <- make_curve_panel(
  "lichen", c("gpp", "r"),
  tag = "a",
  title = "Lichen",
  show_y_axis = TRUE,
  show_temp_axis = FALSE,
  show_legend = TRUE
)

curve_coral_gr <- make_curve_panel(
  "coral", c("gpp", "r"),
  tag = "b",
  title = "Coral",
  show_y_axis = TRUE,
  show_temp_axis = FALSE,
  show_legend = TRUE
)

curve_cell_gr <- make_curve_panel(
  "cell", c("gpp", "r"),
  tag = "c",
  title = "Foraminifera",
  show_y_axis = TRUE,
  show_temp_axis = FALSE,
  show_legend = TRUE
)

curve_lichen_npp <- make_curve_panel(
  "lichen", c("npp"),
  tag = "d",
  title = NULL,
  show_y_axis = TRUE,
  show_temp_axis = TRUE,
  show_legend = FALSE
)

curve_coral_npp <- make_curve_panel(
  "coral", c("npp"),
  tag = "e",
  title = NULL,
  show_y_axis = TRUE,
  show_temp_axis = TRUE,
  show_legend = FALSE
)

curve_cell_npp <- make_curve_panel(
  "cell", c("npp"),
  tag = "f",
  title = NULL,
  show_y_axis = TRUE,
  show_temp_axis = TRUE,
  show_legend = FALSE
)

fig_slopes <-
  (curve_lichen_gr + curve_coral_gr + curve_cell_gr) /
  (curve_lichen_npp + curve_coral_npp + curve_cell_npp)

fig_slopes

ggsave(
  "figures/fig_slopes_GPP_R_NPP_with_silhouettes.png",
  fig_slopes,
  width = fig_w,
  height = fig_h,
  dpi = 600
)



#################
# thermal trait comparison between lichen and coral --> topt and tbreadth?

get_coefs_CI_breadth <- function(mod, metabolic_category, dataset) {
  sm <- summary(mod)
  
  if (inherits(mod, "averaging")) {
    # MuMIn::model.avg object
    coefs <- as.data.frame(sm$coefmat.full)
    
    # column names are usually "Estimate" and "Std. Error"
    est_col <- "Estimate"
    se_col  <- if ("Std. Error" %in% names(coefs)) "Std. Error" else "SE"
    
  } else if (inherits(mod, "lme")) {
    # nlme::lme object
    coefs <- as.data.frame(sm$tTable)
    
    # column names are "Value" and "Std.Error"
    est_col <- "Value"
    se_col  <- "Std.Error"
    
  } else {
    stop("Unsupported model class: ", paste(class(mod), collapse = ", "))
  }
  
  coefs %>%
    rownames_to_column("term") %>%
    transmute(
      term,
      estimate = .data[[est_col]],
      se       = .data[[se_col]],
      lower_95 = estimate - 1.96 * se,
      upper_95 = estimate + 1.96 * se,
      metabolic_category = metabolic_category,
      dataset  = dataset
    )
}


lichen_gpp_coefs_breadth <- get_coefs_CI_breadth(lichen_gpp_breadth_top, "gpp", "lichen")
lichen_npp_coefs_breadth <- get_coefs_CI_breadth(lichen_npp_breadth_top, "npp", "lichen")

coral_gpp_coefs_breadth  <- get_coefs_CI_breadth(coral_gpp_breadth_top_avg,  "gpp", "coral")
coral_npp_coefs_breadth  <- get_coefs_CI_breadth(coral_npp_breadth_top,  "npp", "coral")
coral_r_coefs_breadth    <- get_coefs_CI_breadth(coral_r_breadth_top,    "r",   "coral")

all_breadth_coefs <- bind_rows(lichen_gpp_coefs_breadth, lichen_npp_coefs_breadth, coral_gpp_coefs_breadth, coral_npp_coefs_breadth, coral_r_coefs_breadth)


get_coefs_CI_topt <- function(mod, metabolic_category, dataset) {
  sm <- summary(mod)
  
  if (inherits(mod, "averaging")) {
    # MuMIn::model.avg object
    coefs <- as.data.frame(sm$coefmat.full)
    
    # column names are usually "Estimate" and "Std. Error"
    est_col <- "Estimate"
    se_col  <- if ("Std. Error" %in% names(coefs)) "Std. Error" else "SE"
    
  } else if (inherits(mod, "lme")) {
    # nlme::lme object
    coefs <- as.data.frame(sm$tTable)
    
    # column names are "Value" and "Std.Error"
    est_col <- "Value"
    se_col  <- "Std.Error"
    
  } else {
    stop("Unsupported model class: ", paste(class(mod), collapse = ", "))
  }
  
  coefs %>%
    rownames_to_column("term") %>%
    transmute(
      term,
      estimate = .data[[est_col]],
      se       = .data[[se_col]],
      lower_95 = estimate - 1.96 * se,
      upper_95 = estimate + 1.96 * se,
      metabolic_category = metabolic_category,
      dataset  = dataset
    )
}

lichen_gpp_coefs_topt <- get_coefs_CI_topt(lichen_gpp_topt_top_avg, "gpp", "lichen")
lichen_npp_coefs_topt <- get_coefs_CI_topt(lichen_npp_topt, "npp", "lichen")

coral_gpp_coefs_topt <- get_coefs_CI_topt(coral_gpp_topt,  "gpp", "coral")
coral_npp_coefs_topt  <- get_coefs_CI_topt(coral_npp_topt,  "npp", "coral")
coral_r_coefs_topt    <- get_coefs_CI_topt(coral_r_topt,    "r",   "coral")

all_topt_coefs <- bind_rows(lichen_gpp_coefs_topt, lichen_npp_coefs_topt, coral_gpp_coefs_topt, coral_npp_coefs_topt, coral_r_coefs_topt)


lichen_raw_breadth <- full_Ea %>% filter(dataset == "lichen") %>% filter(metabolic_category != "r") %>% 
  mutate(metabolic_category = factor(metabolic_category, levels = c("gpp", "npp", "r")))
lichen_breadth_mod_ea <- all_breadth_coefs %>% filter(dataset == "lichen" & term == "(Intercept)") %>% 
  mutate(metabolic_category = factor(metabolic_category, levels = c("gpp", "npp", "r")))

labels_df_l_breadth <- data.frame(
  dataset = c("Lichen", "Lichen"),           
  metabolic_category = c("gpp", "npp"),
  y = c(37, 37),  
  x= c(1.3, 2.3),
  label = c("n = 43", "n = 53"))

lichen_breadth <- 
  ggplot()+
  geom_jitter(data=lichen_raw_breadth,aes(x=metabolic_category, y=as.numeric(breadth), colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=lichen_breadth_mod_ea, aes(x=metabolic_category, y=estimate, ymin=lower_95, ymax=upper_95, colour=metabolic_category), size=1, linewidth=1.7,lineend='round')+
  coord_flip()+
  scale_x_discrete(limits=c("gpp", "npp", "r"), labels=c("GPP", "NPP", "R"), drop=FALSE)+
  scale_y_continuous(limits = c(-3, 40), breaks = c(0, 10, 20, 30, 40)) +
  scale_colour_manual(values=c("yellowgreen","lightgoldenrod4"))+
  scale_fill_manual(values=c("yellowgreen","lightgoldenrod4"))+
  labs( x="",
        y = "", 
        title ="Thermal breadth") +
  theme_classic() +
  theme(legend.position="none",
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.15), 
        plot.title.position = "plot")+
  #add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=2.4, y=30, height=0.35, alpha=1)+
  geom_label(data = labels_df_l_breadth,
             aes(x = x, y = y, label = label, fill=metabolic_category),
             color = "white", hjust = -0.15, vjust=1, size = 3)+
  annotate("text", y=40, x=3.4, label= "a", size=6, fontface = "italic")+
  add_phylopic(
    uuid   = "a208bba4-f4bf-4810-bcc9-c5868836fc76",
    x      = 3.2,
    y      = -2.85,
    height = 1,
    alpha  = 1
  )

coral_raw_breadth <- full_Ea %>% filter(dataset == "coral")
coraL_breadth_mod_ea <- all_breadth_coefs %>% filter(dataset == "coral" & term == "(Intercept)") 

labels_df_c_breadth <- data.frame(
  dataset = c("Coral", "Coral", "Coral"),           
  metabolic_category = c("gpp", "npp", "r"),
  y = c(35, 35, 35),  
  x= c(1.3, 2.3, 3.3),
  label = c("n = 46", "n = 38", "n=56"))
  

coral_breadth <- 
  ggplot()+
  geom_jitter(data=coral_raw_breadth,aes(x=metabolic_category, y=as.numeric(breadth), colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=coraL_breadth_mod_ea, aes(x=metabolic_category, y=estimate, ymin=lower_95, ymax=upper_95, colour=metabolic_category), size=1, linewidth=1.7,lineend='round')+
  coord_flip()+
  scale_x_discrete(limits=c("gpp", "npp", "r"), labels=c("GPP", "NPP", "R"))+
  scale_y_continuous(limits = c(-3, 40), breaks = c(0, 10, 20, 30, 40)) +
  scale_colour_manual(values=c("yellowgreen","lightgoldenrod4","honeydew4"))+
  scale_fill_manual(values=c("yellowgreen","lightgoldenrod4","honeydew4"))+
  #annotate("text", x=3.4, y=5.5, label= "Thermal breadth", size=4.5, fontface="italic") + 
  labs( x="",
        y = "Temperature (°C)") +
  theme_classic() +
  theme(legend.position="none",
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.15),
        plot.title.position = "plot")+
 # add_phylopic(uuid = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104", x=3.4, y=30, height=0.4, alpha=1)+
  geom_label(data = labels_df_c_breadth,
             aes(x = x, y = y, label = label, fill=metabolic_category),
             color = "white", hjust = -0.15, vjust=1, size = 3)+
  annotate("text", y=40, x=3.4, label= "b", size=6, fontface = "italic")+
  add_phylopic(
    uuid   = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",
    x      = 3.2,
    y      = -2.5,
    height = 0.7,
    alpha  = 1
  )


lichen_raw_topt <- full_Ea %>% filter(dataset == "lichen") %>% filter(metabolic_category != "r")
lichen_topt_mod_ea <- all_topt_coefs %>% filter(dataset == "lichen" & term == "(Intercept)") 

labels_df_l_topt <- data.frame(
  dataset = c("Lichen", "Lichen"),           
  metabolic_category = c("gpp", "npp"),
  y = c(37, 37),  
  x= c(1.3, 2.3),
  label = c("n = 43", "n = 53"))


lichen_topt <- 
  ggplot()+
  geom_jitter(data=lichen_raw_topt,aes(x=metabolic_category, y=as.numeric(topt), colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=lichen_topt_mod_ea, aes(x=metabolic_category, y=estimate, ymin=lower_95, ymax=upper_95, colour=metabolic_category), size=1, linewidth=1.7,lineend='round')+
  coord_flip()+
  scale_x_discrete(limits=c("gpp", "npp", "r"), labels=c("GPP", "NPP", "R"), drop=FALSE)+
  scale_y_continuous(limits = c(-3, 40), breaks = c(0, 10, 20, 30, 40)) +
  scale_colour_manual(values=c("yellowgreen","lightgoldenrod4"))+
  scale_fill_manual(values=c("yellowgreen","lightgoldenrod4"))+
  #annotate("text", x=2.4, y=-2, label= "Thermal optimum", size=4.5, fontface="italic") + 
  labs(y = "Temperature (°C)", 
       x="", 
       title = "Thermal optimum") +
  theme_classic() +
  theme(legend.position="none",
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face = "italic", hjust = 0))+
  #add_phylopic(uuid = "a208bba4-f4bf-4810-bcc9-c5868836fc76", x=2.4, y=42, height=0.35, alpha=1)+
  geom_label(data = labels_df_l_topt,
             aes(x = x, y = y, label = label, fill=metabolic_category),
             color = "white", hjust = -0.15, vjust=1, size = 3)+
  annotate("text", y=40, x=3.4, label= "c", size=6, fontface = "italic")+
  add_phylopic(
    uuid   = "a208bba4-f4bf-4810-bcc9-c5868836fc76",
    x      = 3.2,
    y      = -2.85,
    height = 1,
    alpha  = 1
  )
coral_raw_topt <- full_Ea %>% filter(dataset == "coral") 
coral_topt_mod_ea <- all_topt_coefs %>% filter(dataset == "coral" & term == "(Intercept)") 

labels_df_c_topt <- data.frame(
  dataset = c("Coral", "Coral", "Coral"),           
  metabolic_category = c("gpp", "npp", "r"),
  y = c(35, 35, 35),  
  x= c(1.3, 2.3, 3.3),
  label = c("n = 46", "n = 38", "n=56"))


coral_topt <- 
  ggplot()+
  geom_jitter(data=coral_raw_topt,aes(x=metabolic_category, y=as.numeric(topt), colour=metabolic_category), alpha=0.15, width=0.05)+
  geom_pointrange(data=coral_topt_mod_ea, aes(x=metabolic_category, y=estimate, ymin=lower_95, ymax=upper_95, colour=metabolic_category), size=1, linewidth=1.7,lineend='round')+
  coord_flip()+
  scale_x_discrete(limits=c("gpp", "npp", "r"), labels=c("GPP", "NPP", "R"))+
  scale_colour_manual(values=c("yellowgreen","lightgoldenrod4", "honeydew4"))+
  scale_fill_manual(values=c("yellowgreen","lightgoldenrod4", "honeydew4"))+
  #annotate("text", x=3.4, y=5.5, label= "Thermal optimum", size=4.5, fontface="italic") + 
  scale_y_continuous(limits = c(-3, 40), breaks = c(0, 10, 20, 30, 40)) +
  labs( x="",
        y = "Temperature (° C)") +
  theme_classic() +
  theme(legend.position="none",
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  geom_label(data = labels_df_c_topt,
             aes(x = x, y = y, label = label, fill=metabolic_category),
             color = "white", hjust = -0.15, vjust=1, size = 3)+
  annotate("text", y=40, x=3.4, label= "d", size=6, fontface = "italic")+
  add_phylopic(
    uuid   = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",
    x      = 3.2,
    y      = -2.5,
    height = 0.7,
    alpha  = 1
  )

fig3_tempdep <- (fig3a + fig3b)

ggsave(fig3_tempdep, filename = "./figures/figure3.png", dpi=700, width=12, height=6)


fig4_breadth <- (lichen_breadth / coral_breadth) 
fig4_topt <- (lichen_topt / coral_topt) 
fig4 <- fig4_breadth | fig4_topt + plot_layout(heights = c(0.5, 0.5)) 

ggsave(fig4, filename = "./figures/figure4.png", dpi=700, width=12, height=4)


####### figure 4 -> additional covariates included in top models

####  general function to extract the slopes and generate CI ribbons ####

fe_ci_line <- function(model, data, xvar, n = 200) {
  if (!xvar %in% names(data)) {
    stop("xvar '", xvar, "' not found in data")
  }
  
  # --- get fixed-effect coef + vcov depending on model type ---
  if (inherits(model, "averaging")) {
    # MuMIn::model.avg object
    beta <- coef(model, full = TRUE)  # named vector
    V    <- vcov(model)               # matrix
  } else if (inherits(model, c("lme", "lmerMod", "lmerModLmerTest", "glmerMod"))) {
    beta <- fixef(model)
    V    <- vcov(model)
  } else {
    # fallback for lm/glm/etc.
    beta <- coef(model)
    V    <- vcov(model)
  }
  
  if (!all(c("(Intercept)", xvar) %in% names(beta))) {
    stop(
      "Model must have an intercept and a slope named '", xvar, "'.\n",
      "Available fixed effects: ", paste(names(beta), collapse = ", ")
    )
  }
  
  # --- prediction grid over xvar ---
  x_seq <- seq(
    min(data[[xvar]], na.rm = TRUE),
    max(data[[xvar]], na.rm = TRUE),
    length.out = n
  )
  
  newdat <- data.frame(x = x_seq)
  names(newdat)[1] <- xvar
  
  # --- design matrix for intercept + xvar ---
  X <- model.matrix(reformulate(xvar), newdat)  # ~ xvar
  
  idx   <- c("(Intercept)", xvar)
  beta2 <- beta[idx]
  V2    <- V[idx, idx, drop = FALSE]
  
  # --- fitted line + 95% CI ---
  fit <- as.numeric(X %*% beta2)
  se  <- sqrt(diag(X %*% V2 %*% t(X)))
  
  newdat$fit <- fit
  newdat$lwr <- fit - 1.96 * se
  newdat$upr <- fit + 1.96 * se
  
  newdat
}


ci_coral_gpp_br <- fe_ci_line(model = coral_gpp_model_top, data  = gpp_coral_Eas,xvar  = "scaled_breadth")

label_df1 <- data.frame(
  x = 3.7,    # choose where in x-space you want the text
  y = 2.8,     # choose where in y-space you want the text
  label = "slope=-0.04, *"
)

figure4.1 <- ggplot() +
  geom_point(data = gpp_coral_Eas, aes(x = as.numeric(scaled_breadth), y = as.numeric(e)), colour = "yellowgreen", alpha  = 0.6) +
  geom_ribbon(data = ci_coral_gpp_br, aes(x = scaled_breadth, ymin = lwr, ymax = upr), fill  = "yellowgreen",alpha = 0.15) +
  geom_line(data = ci_coral_gpp_br, aes(x = scaled_breadth, y = fit), colour = "yellowgreen") +
  labs(x = "Scaled thermal breadth", y = "Temperature 
       dependence (Ea)",title = "", tag = "a" ) +
  #xlim(-15,15)+
  ylim(-0.5,3.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",x = 6.7, y = 3.2,height = 0.7, alpha  = 1)+
  geom_text(
    data = label_df1,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = 4
  )

ci_coral_npp_br <- fe_ci_line(model = coral_npp_model_top_avg, data  = npp_coral_Eas,xvar  = "scaled_breadth")

label_df2 <- data.frame(
  x = 1.5,    # choose where in x-space you want the text
  y = 2.7,     # choose where in y-space you want the text
  label = "slope=-0.065, ***"
)


figure4.2 <- ggplot() +
  geom_point(data = npp_coral_Eas, aes(x = as.numeric(scaled_breadth), y = as.numeric(e)), colour = "lightgoldenrod4", alpha  = 0.6) +
  geom_ribbon(data = ci_coral_npp_br , aes(x = scaled_breadth, ymin = lwr, ymax = upr), fill  = "lightgoldenrod4",alpha = 0.15) +
  geom_line(data = ci_coral_npp_br , aes(x = scaled_breadth, y = fit), colour = "lightgoldenrod4") +
  labs(x = "Scaled thermal breadth", y = "Temperature 
       dependence (Ea)",title = "", tag = "c" ) +
  #xlim(-15,15)+
  ylim(-0.5,3.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",x = 4.5, y = 3.2,height = 0.7, alpha  = 1)+
  geom_text(
    data = label_df2,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = 4)


ci_coral_r_br <- fe_ci_line(model = coral_r_model_top, data  = r_coral_Eas,xvar  = "scaled_breadth")

label_df3 <- data.frame(
  x = 2.5,    # choose where in x-space you want the text
  y = 2.7,     # choose where in y-space you want the text
  label = "slope=-0.05, *"
)


figure4.3 <- ggplot() +
  geom_point(data = r_coral_Eas, aes(x = as.numeric(scaled_breadth), y = as.numeric(e)), colour = "honeydew4", alpha  = 0.6) +
  geom_ribbon(data = ci_coral_r_br , aes(x = scaled_breadth, ymin = lwr, ymax = upr), fill  = "honeydew4",alpha = 0.15) +
  geom_line(data = ci_coral_r_br , aes(x = scaled_breadth, y = fit), colour = "honeydew4") +
  labs(x = "Scaled thermal breadth", y = "Temperature 
       dependence (Ea)",title = "", tag = "e" ) +
  #xlim(-15,15)+
  ylim(-0.5,3.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",x = 4.5, y = 3.2, height = 0.7, alpha  = 1)+
  geom_text(
    data = label_df3,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = 4)

ci_lichen_gpp_br <- fe_ci_line(model = lichen_gpp_model_top_avg, data  = gpp_lichen_Eas,xvar  = "scaled_breadth")

label_df4 <- data.frame(
  x = 5.5,    # choose where in x-space you want the text
  y = 2.7,     # choose where in y-space you want the text
  label = "slope=-0.02, ***"
)

figure4.4 <- ggplot() +
  geom_point(data = gpp_lichen_Eas, aes(x = as.numeric(scaled_breadth), y = as.numeric(e)), colour = "yellowgreen", alpha  = 0.6) +
  geom_ribbon(data = ci_lichen_gpp_br  , aes(x = scaled_breadth, ymin = lwr, ymax = upr), fill  = "yellowgreen",alpha = 0.15) +
  geom_line(data = ci_lichen_gpp_br  , aes(x = scaled_breadth, y = fit), colour = "yellowgreen") +
  labs(x = "Scaled thermal breadth", y = "",title = "", tag = "b" ) +
  #xlim(-15,15)+
  ylim(-0.5,3.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "a208bba4-f4bf-4810-bcc9-c5868836fc76",x = 13.5, y = 3.2,height = 0.9, alpha  = 1)+
  geom_text(
    data = label_df4,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = 4)

ci_lichen_npp_br <- fe_ci_line(model = lichen_npp_model_top_avg, data  = npp_lichen_Eas,xvar  = "scaled_breadth")

label_df5 <- data.frame(
  x = 4.5,    # choose where in x-space you want the text
  y = 2.7,     # choose where in y-space you want the text
  label = "slope=-0.02, ***")

figure4.5 <- ggplot() +
  geom_point(data = npp_lichen_Eas, aes(x = as.numeric(scaled_breadth), y = as.numeric(e)), colour = "lightgoldenrod4", alpha  = 0.6) +
  geom_ribbon(data = ci_lichen_npp_br  , aes(x = scaled_breadth, ymin = lwr, ymax = upr), fill  = "lightgoldenrod4",alpha = 0.15) +
  geom_line(data = ci_lichen_npp_br  , aes(x = scaled_breadth, y = fit), colour = "lightgoldenrod4") +
  labs(x = "Scaled thermal breadth", y = "",title = "", tag = "d" ) +
  #xlim(-15,15)+
  ylim(-0.5,3.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "a208bba4-f4bf-4810-bcc9-c5868836fc76",x = 11, y = 3.2,height = 0.9, alpha  = 1)+
  geom_text(
    data = label_df5,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = 4)


figure4 <- (figure4.1 | figure4.4) / (figure4.2 | figure4.5) / (figure4.3 | plot_spacer())

ggsave(figure4, filename = "./figures/figure4.png", dpi=700, width=8, height=11)












##########################################
##################### SI #################
##########################################


# model selection

model_selection_plus_table <- function(dredge_obj,
                                       n_models = NULL,
                                       include_formula = FALSE) {
  
  df <- as.data.frame(dredge_obj)
  
  # Optionally keep only top N models
  if (!is.null(n_models)) {
    df <- df[seq_len(min(n_models, nrow(df))), , drop = FALSE]
  }
  
  # Identify predictor columns
  term_names <- attr(dredge_obj, "term.names")
  term_names <- setdiff(term_names, "(Intercept)")
  term_names <- term_names[term_names %in% names(df)]
  
  # Term block: "+" / "" indicators
  term_block <- df[, term_names, drop = FALSE] %>%
    mutate(across(everything(), ~ ifelse(is.na(.x) | .x == 0, "", "+")))
  
  # Info block
  info_block <- tibble(
    k      = df[["df"]],
    logLik = round(df[["logLik"]], 3),
    AICc   = round(df[["AICc"]], 3),
    delta  = round(df[["delta"]], 3),
    weight = round(df[["weight"]], 3)
  )
  
  # Optional formula column
  if (include_formula) {
    formulas <- lapply(get.models(dredge_obj, subset = TRUE), formula)
    formula_strings <- vapply(
      formulas,
      function(f) paste(deparse(f), collapse = " "),
      FUN.VALUE = character(1)
    )
    info_block <- info_block %>%
      mutate(Formula = formula_strings[seq_len(nrow(df))])
  }
  
  # Model index as first column
  model_col <- tibble(Model = seq_len(nrow(df)))
  
  bind_cols(model_col, term_block, info_block)
}


# --------------------------
# Helper: build system table
# --------------------------
build_system_table <- function(system_models, delta_cutoff = 2) {
  # system_models: named list, e.g. list(
  #   NPP = lichen_npp_model_set,
  #   GPP = lichen_gpp_model_set,
  #   R   = lichen_r_model_set
  # )
  
  bind_rows(
    lapply(names(system_models), function(proc) {
      dredge_obj <- system_models[[proc]]
      
      # subset models by ΔAIC
      top_models <- subset(dredge_obj, delta < delta_cutoff)
      
      if (nrow(top_models) == 0) {
        return(NULL)
      }
      
      tbl <- model_selection_plus_table(top_models)
      
      # Add Process column as "subbanner" grouping
      tbl %>%
        mutate(Process = proc, .before = Model)
    }),
    .id = NULL
  )
}

# --------------------------
# Define your system-specific lists
# --------------------------
lichen_models <- list(
  NPP = lichen_npp_model_set,
  GPP = lichen_gpp_model_set,
  R   = lichen_r_model_set
)

coral_models <- list(
  NPP = coral_npp_model_set,
  GPP = coral_gpp_model_set,
  R   = coral_r_model_set
)

# Build data frames
lichen_df <- build_system_table(lichen_models, delta_cutoff = 2)
coral_df  <- build_system_table(coral_models,  delta_cutoff = 2)

# --------------------------
# Turn into flextables with "subbanner" formatting
# --------------------------
make_grouped_ft <- function(df) {
  ft <- flextable(df)
  
  # Merge Process cells vertically to create visual groups
  if ("Process" %in% names(df)) {
    ft <- merge_v(ft, j = "Process")
    ft <- valign(ft, j = "Process", valign = "top")
    ft <- bold(ft, j = "Process", part = "body")
  }
  
  ft <- autofit(ft)
  ft
}

ft_lichen <- make_grouped_ft(lichen_df)
ft_coral  <- make_grouped_ft(coral_df)

# --------------------------
# Write to DOCX with two tables
# --------------------------
doc <- read_docx()

doc <- doc %>%
  body_add_par("Lichen model selection", style = "heading 1") %>%
  body_add_flextable(ft_lichen) %>%
  body_add_par("") %>%
  body_add_par("Coral model selection", style = "heading 1") %>%
  body_add_flextable(ft_coral)

print(doc, target = "model_selection_tables_grouped.docx")




















# SI figures latitude
# Get CI line for this specific panel
ci_coral_gpp <- fe_ci_line(model = coral_gpp_model_top, data  = gpp_coral_Eas,xvar  = "scaled_latitude")

figure5a <- ggplot() +
  geom_point(data = gpp_coral_Eas, aes(x = as.numeric(scaled_latitude), y = as.numeric(e)), colour = "yellowgreen", alpha  = 0.6) +
  geom_ribbon(data = ci_coral_gpp, aes(x = scaled_latitude, ymin = lwr, ymax = upr), fill  = "yellowgreen",alpha = 0.15) +
  geom_line(data = ci_coral_gpp, aes(x = scaled_latitude, y = fit), colour = "yellowgreen") +
  labs(x = "Scaled latitude", y = "Temperature dependence (Ea)",title = "", tag = "d" ) +
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "f6a243aa-5cb1-41a2-a52c-c8d4c4300104",x = 18, y = 2.3,height = 0.4, alpha  = 1)

#scaled latitude was not in the best model for coral npp and r 


lichen_gpp_coefs <- get_coefs_CI(, "gpp", "lichen")
lichen_npp_coefs <- get_coefs_CI(, "npp", "lichen")
lichen_r_coefs   <- get_coefs_CI(lichen_r_model_top_avg,   "r",   "lichen")


ci_lichen_gpp <- fe_ci_line(model = lichen_gpp_model_top_avg, data  = gpp_lichen_Eas,xvar  = "scaled_latitude")

figure5b <- ggplot() +
  geom_point(data = gpp_lichen_Eas, aes(x = as.numeric(scaled_latitude), y = as.numeric(e)), colour = "yellowgreen", alpha  = 0.6) +
  geom_ribbon(data = ci_lichen_gpp, aes(x = scaled_latitude, ymin = lwr, ymax = upr), fill  = "yellowgreen",alpha = 0.15) +
  geom_line(data = ci_lichen_gpp, aes(x = scaled_latitude, y = fit), colour = "yellowgreen") +
  labs(x = "Scaled latitude", y = "Temperature dependence (Ea)",title = "", tag = "d" ) +
  ylim(0, 2.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "a208bba4-f4bf-4810-bcc9-c5868836fc76",x = 31, y = 2.3,height = 0.5, alpha  = 1)

ci_lichen_npp <- fe_ci_line(model = lichen_npp_model_top_avg, data  = npp_lichen_Eas, xvar  = "scaled_latitude")


figure4c <- ggplot() +
  geom_point(data = npp_lichen_Eas, aes(x = as.numeric(scaled_latitude), y = as.numeric(e)), colour = "lightgoldenrod4", alpha  = 0.6) +
  geom_ribbon(data = ci_lichen_npp, aes(x = scaled_latitude, ymin = lwr, ymax = upr), fill  = "lightgoldenrod4",alpha = 0.15) +
  geom_line(data = ci_lichen_npp, aes(x = scaled_latitude, y = fit), colour = "lightgoldenrod4") +
  labs(x = "Scaled latitude", y = "Temperature dependence (Ea)",title = "", tag = "d" ) +
  ylim(0, 2.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "a208bba4-f4bf-4810-bcc9-c5868836fc76",x = 31, y = 2.3,height = 0.5, alpha  = 1)


ci_lichen_r <- fe_ci_line(model = lichen_r_model_top_avg, data  = r_lichen_Eas,xvar  = "scaled_latitude")

figure4d <- ggplot() +
  geom_point(data = r_lichen_Eas, aes(x = as.numeric(scaled_latitude), y = as.numeric(e)), colour = "lightgoldenrod4", alpha  = 0.6) +
  geom_ribbon(data = ci_lichen_npp, aes(x = scaled_latitude, ymin = lwr, ymax = upr), fill  = "lightgoldenrod4",alpha = 0.15) +
  geom_line(data = ci_lichen_npp, aes(x = scaled_latitude, y = fit), colour = "lightgoldenrod4") +
  labs(x = "Scaled latitude", y = "Temperature dependence (Ea)",title = "", tag = "d" ) +
  ylim(0, 2.5)+
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)) +
  add_phylopic( uuid   = "a208bba4-f4bf-4810-bcc9-c5868836fc76",x = 31, y = 2.3,height = 0.5, alpha  = 1)





#####







##################### Figures ######################
#plots
#install.packages("remotes")
pacman::p_load(rphylopic)
#code to get the uuid's for the phylopic pngs
uuid <- rphylopic::get_uuid(name = "Sorites orbiculus")







