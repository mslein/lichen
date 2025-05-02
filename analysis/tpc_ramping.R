####ramping tpc analysis 
pacman::p_load(rTPC, nls.multstart, broom, tidyverse)
ramping <- read_csv("extraction/coral dataset/ramping/ramping_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("%_O2") ~ response_value,
                         gas_units %in% c("mg_O2") ~ response_value*(1/32000),
                         gas_units %in% c("umol_O2") ~ response_value*(1/10^6)),
         gC = case_when(`mass/area` %in% c("g") ~ 1*0.5, 
                          `mass/area` %in% c("cm^2") ~  (5*44), 
                          `mass/area` %in% c("m^2") ~ (44*0.01)),
         min = case_when(time %in% c("hr") ~ 60),
         mol_gCmin= mol/(gC*min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, ramp_rate_per_day, cnidarian_type, study_id, metabolic_category)

#breaking up metabolic categories
npp <- ramping %>% filter(metabolic_category == "npp") %>%
filter(ramp_rate_per_day != 0) #omitting studies that have zeros (controls)
r <- ramping %>% filter(metabolic_category == "r") %>%  
  mutate(mol_gCmin=abs(mol_gCmin)) %>%
  select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, ramp_rate_per_day, cnidarian_type, study_id, metabolic_category) %>%
  filter(ramp_rate_per_day != 0) #omitting studies that have zeros (controls)
gpp <- ramping %>% filter(metabolic_category == "gpp") %>%
filter(ramp_rate_per_day != 0) #omitting studies that have zeros (controls)

# gpp and npp arithmetic
# gpp
gpp_calcs <-  read_csv("extraction/coral dataset/ramping/ramping_mrcalcs24mar25.csv") %>%
  filter(study_id != "aichelman2019") %>%
  mutate(gpp= npp_value + abs(r_value)) %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("%_O2") ~ gpp,
                         gas_units %in% c("mg_O2") ~ gpp*(1/32000),
                         gas_units %in% c("umol_O2") ~ gpp*(1/10^6)),
         gC = case_when(`mass/area` %in% c("g") ~ 1*0.5, 
                          `mass/area` %in% c("cm^2") ~ (5*44), 
                          `mass/area` %in% c("m^2") ~ (44*0.01)),
         min = case_when(time %in% c("hr") ~ 60),
         mol_gCmin= mol/(gC*min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove = FALSE) %>%
  select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, ramp_rate_per_day, cnidarian_type, study_id)

gpp_collect <- gpp %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, ramp_rate_per_day, cnidarian_type, study_id)
gpp_all <- rbind(gpp_collect, gpp_calcs)

#npp
npp_calcs <-  read_csv("extraction/coral dataset/ramping/ramping_mrcalcs24mar25.csv") %>%
  filter(study_id == "aichelman2019") %>%
  mutate(npp = gpp_value - abs(r_value)) %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("%_O2") ~ npp,
                         gas_units %in% c("mg_O2") ~ npp*(1/32000),
                         gas_units %in% c("umol_O2") ~ npp*(1/10^6)),
         gC = case_when(`mass/area` %in% c("g") ~ 1*0.5, 
                          `mass/area` %in% c("cm^2") ~ (5*44), 
                          `mass/area` %in% c("m^2") ~ (44*0.01)),
         min = case_when(time %in% c("hr") ~ 60),
         mol_gCmin= mol/(gC*min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, ramp_rate_per_day, cnidarian_type, study_id)

npp_collect <- npp %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, ramp_rate_per_day, cnidarian_type, study_id)
npp_all <- rbind(npp_collect, npp_calcs)


fit_and_predict <- function(data) {
  temp_data <- data %>%
    group_by(tpc_grp) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    mutate(Temperature = as.numeric(temp),
           rate = abs(as.numeric(mol_gCmin)))  # ensure positive rates
  
  # Skip if not enough data
  if (nrow(temp_data) < 4) return(NULL)
  
  # Get tailored start values and bounds for this curve
  start_vals <- get_start_vals(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  low_lims   <- get_lower_lims(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  upper_lims <- get_upper_lims(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  
  # Fit TPC model
  fit <- tryCatch({
    nls_multstart(rate ~ pawar_2018(temp = Temperature, r_tref, e, eh, th, tref = 15),
                  data = temp_data,
                  iter = 500,
                  start_lower = start_vals - 10,
                  start_upper = start_vals + 10,
                  lower = low_lims,
                  upper = upper_lims,
                  supp_errors = "Y")
  }, error = function(e) return(NULL))
  
  # Return NULL if model failed
  if (is.null(fit)) return(NULL)
  
  # Extract params and make predictions
  est <- calc_params(fit) %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  new_data <- data.frame(Temperature = seq(min(temp_data$Temperature),
                                           max(temp_data$Temperature),
                                           by = 0.5))
  
  preds <- broom::augment(fit, newdata = new_data) %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  return(list(data = temp_data, preds = preds, est = est))
}

########################################################
#######################   NPP   ########################
########################################################

#generating tpcs for NPP ####################
results_npp <- npp_all %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)


#visualizing TPCs for NPP ####################
combined_data_npp <- do.call(rbind, lapply(results_npp, function(x) x$data))
combined_preds_npp <- do.call(rbind, lapply(results_npp, function(x) x$preds))

plot_data_npp <- bind_rows(
  combined_data_npp %>% mutate(Type = "Observed"),
  combined_preds_npp %>% mutate(Type = "Fitted")
)
ggplot(plot_data_npp, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_npp, Type == "Observed")) +
  geom_line(data = subset(plot_data_npp, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "bottom")


#pulling the activation energies out for NPP####################
combined_est_npp <- do.call(rbind, lapply(results_npp, function(x) x$est)) %>%
  filter(!tpc_grp %in% c("7_jurriaans2020", "5_jurriaans", "5_jurriaans2020",  #these studies didn't have an Ea estimates
                           "2_puntin2023", "8_jurriaans2019", "1_jurriaans2019", 
                         "3_jurriaans2019", "3_jurriaans2021")) %>%
  select(tpc_grp, e, eh, topt)

#creating a dataframe of the TPC parameters joined w/ attributes from the full dataset ####################
model_comp <- left_join(combined_est_npp, npp_all, by="tpc_grp") %>%
  select(tpc_grp, e, eh, topt, temp, ramp_rate_per_day, cnidarian_type, latitude) %>%
  group_by(tpc_grp) %>%
  summarise(avg_temp = mean(temp), 
            e=mean(e), 
            eh=mean(eh), 
            topt=mean(topt), 
            ramp_rate_per_day=mean(ramp_rate_per_day), 
            latitude=mean(latitude)) 

pacman::p_load(moments)
#checking for normality
shapiro.test(model_comp$ramp_rate_per_day) #not normal


model_comp$RAMP <- log10(model_comp$ramp_rate_per_day)
skewness(model_comp$RAMP, na.rm = TRUE)#improves skewness slightly



#fitting simple linear models to assess if rate of ramping affects Eas ####################
model00npp_r<- lm(e ~ log10(ramp_rate_per_day), data=model_comp) 
model00npp_r2<- lm(eh ~ log10(ramp_rate_per_day), data=model_comp, na.action="na.exclude") 
model00npp_r3<- lm(topt ~ log10(ramp_rate_per_day), data=model_comp, na.action="na.exclude") 
confint(model00npp_r) #intervals span zero
confint(model00npp_r2) #intervals span zero
confint(model00npp_r3)#intervals span zero

#CONCLUSION: ramping rate does not affect Ea for NPP --> include in full analysis





########################################################
#######################    R    ########################
########################################################

#generating tpcs for R ####################
results_r <- r %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

#visualizing TPCs for R ####################
combined_data_r <- do.call(rbind, lapply(results_r, function(x) x$data))
combined_preds_r <- do.call(rbind, lapply(results_r, function(x) x$preds))

plot_data_r <- bind_rows(
  combined_data_r %>% mutate(Type = "Observed"),
  combined_preds_r %>% mutate(Type = "Fitted")
)
ggplot(plot_data_r, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_r, Type == "Observed")) +
  geom_line(data = subset(plot_data_r, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "bottom")

#pulling the activation energies out for R ####################
combined_est_r <- do.call(rbind, lapply(results_r, function(x) x$est))%>%
  select(tpc_grp, e, eh, topt)

#creating a dataframe of the TPC parameters joined w/ attributes from the full dataset ####################
model_comp2 <- left_join(combined_est_r, r, by="tpc_grp") %>%
  select(tpc_grp, e, eh, topt, temp, ramp_rate_per_day, cnidarian_type, latitude) %>%
  group_by(tpc_grp) %>%
  summarise(avg_temp = mean(temp), 
            e=mean(e), 
            eh=mean(eh), 
            topt=mean(topt), 
            ramp_rate_per_day=mean(ramp_rate_per_day), 
            latitude=mean(latitude)) 

#fitting simple linear models to assess if rate of ramping affects Eas ####################
model00r_r<- lm(e ~ log10(ramp_rate_per_day), data=model_comp2) 
model00r_r2<- lm(eh ~ log10(ramp_rate_per_day), data=model_comp2) 
model00r_r3<- lm(topt ~ log10(ramp_rate_per_day), data=model_comp2) 
confint(model00r_r) #intervals span zero
confint(model00r_r2) #intervals span zero
confint(model00r_r3) #intervals span zero



#CONCLUSION: ramping rate does not affect Ea for R --> include in full analysis




########################################################
#######################   GPP   ########################
########################################################


results_gpp <- gpp_all %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

# Combine results
combined_data_gpp <- do.call(rbind, lapply(results_gpp, function(x) x$data))
combined_preds_gpp <- do.call(rbind, lapply(results_gpp, function(x) x$preds))

plot_data_gpp <- bind_rows(
  combined_data_gpp %>% mutate(Type = "Observed"),
  combined_preds_gpp %>% mutate(Type = "Fitted")
)


ggplot(plot_data_gpp, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_gpp, Type == "Observed")) +
  geom_line(data = subset(plot_data_gpp, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "bottom")

#pulling the activation energies out for R ####################
combined_est_gpp <- do.call(rbind, lapply(results_gpp, function(x) x$est))%>%
  select(tpc_grp, e, eh, topt)

#creating a dataframe of the TPC parameters joined w/ attributes from the full dataset ####################
model_comp3 <- left_join(combined_est_gpp, gpp_all, by="tpc_grp") %>%
  select(tpc_grp, e, eh, topt, temp, ramp_rate_per_day, cnidarian_type, latitude) %>%
  group_by(tpc_grp) %>%
  summarise(e=mean(e), 
            eh=mean(eh), 
            topt=mean(topt), 
            ramp_rate_per_day=mean(ramp_rate_per_day), 
            latitude=mean(latitude)) 

#fitting simple linear models to assess if rate of ramping affects Eas ####################
model00gpp_r<- lm(e ~ log10(ramp_rate_per_day), data=model_comp3) 
model00gpp_r2<- lm(eh ~ log10(ramp_rate_per_day), data=model_comp3) 
model00gpp_r3<- lm(topt ~ log10(ramp_rate_per_day), data=model_comp3) 
confint(model00gpp_r) #intervals span zero
confint(model00gpp_r2) #intervals span zero
confint(model00gpp_r3) #intervals span zero


#CONCLUSION: ramping rate does not affect Ea for GPP --> include in full analysis
npp_export <- npp_all %>% mutate(metabolic_category = "npp")
gpp_export <- gpp_all %>% mutate(metabolic_category = "gpp")
r_export <- r %>% mutate(metabolic_category = "r")


final_dataset <- rbind(r_export, npp_export, gpp_export) %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, cnidarian_type, metabolic_category, study_id) 

write_csv(final_dataset, "analysis/tidy data/tidy_ramping_coral_24marc25.csv")


