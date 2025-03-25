####ramping tpc analysis 
pacman::p_load(rTPC, nls.multstart, broom, tidyverse)
ramping <- read_csv("extraction/coral dataset/ramping/ramping_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("%_O2") ~ response_value,
                         gas_units %in% c("mg_O2") ~ response_value*(1/32000),
                         gas_units %in% c("umol_O2") ~ response_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("cm^2") ~ 1, 
                          `mass/area` %in% c("m^2") ~ 1000),
         min = case_when(time %in% c("hr") ~ 60),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  unite("tpc_grp", c(paper_unit, study_id))

#breaking up metabolic categories
npp <- ramping %>% filter(metabolic_category == "npp")
r <- ramping %>% filter(metabolic_category == "r") %>%  select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_grammin, inv_T, temp, ramp_rate_per_day, cnidarian_type)
gpp <- ramping %>% filter(metabolic_category == "gpp")

# gpp and npp arithmetic
# gpp
gpp_calcs <-  read_csv("extraction/coral dataset/ramping/ramping_mrcalcs24mar25.csv") %>%
  filter(study_id != "aichelman2019") %>%
  mutate(r_value = abs(r_value)*-1,
         response_value= npp_value + r_value) %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("%_O2") ~ response_value,
                         gas_units %in% c("mg_O2") ~ response_value*(1/32000),
                         gas_units %in% c("umol_O2") ~ response_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("cm^2") ~ 1, 
                          `mass/area` %in% c("m^2") ~ 1000),
         min = case_when(time %in% c("hr") ~ 60),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  unite("tpc_grp", c(paper_unit, study_id)) %>%
  select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_grammin, inv_T, temp, ramp_rate_per_day, cnidarian_type)

gpp_collect <- gpp %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_grammin, inv_T, temp, ramp_rate_per_day, cnidarian_type)
gpp_all <- rbind(gpp_collect, gpp_calcs)

#npp
npp_calcs <-  read_csv("extraction/coral dataset/ramping/ramping_mrcalcs24mar25.csv") %>%
  filter(study_id == "aichelman2019") %>%
  mutate(r_value = abs(r_value)*-1,
         response_value= npp_value + r_value) %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("%_O2") ~ response_value,
                         gas_units %in% c("mg_O2") ~ response_value*(1/32000),
                         gas_units %in% c("umol_O2") ~ response_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("cm^2") ~ 1, 
                          `mass/area` %in% c("m^2") ~ 1000),
         min = case_when(time %in% c("hr") ~ 60),
         mol_grammin= mol/(gram*min),
         abs_molco2_g_min = abs(mol_grammin), 
         log=log(abs_molco2_g_min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral)) %>%
  unite("tpc_grp", c(paper_unit, study_id)) %>%
  select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_grammin, inv_T, temp, ramp_rate_per_day, cnidarian_type)

npp_collect <- npp %>% select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_grammin, inv_T, temp, ramp_rate_per_day, cnidarian_type)
npp_all <- rbind(npp_collect, npp_calcs)


fit_and_predict <- function(data) {
  
  temp_data <- data %>%
    mutate(Temperature = as.numeric(temp), 
           rate=abs(mol_grammin)) # just in case integer causes problems
  
  # Get start values and limits
  start_vals <- get_start_vals(temp_data$Temperature, temp_data$rate, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(temp_data$Temperature, temp_data$rate, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(temp_data$Temperature, temp_data$rate, model_name = 'sharpeschoolhigh_1981')
  
  # Fit the model
  fit <- nls_multstart(rate ~ sharpeschoolhigh_1981(temp = Temperature, r_tref, e, eh, th, tref = 15),
                       data = temp_data,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  
  est <- calc_params(fit) %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  # Generate new data for predictions using augment
  new_data <- data.frame(Temperature = seq(min(temp_data$Temperature), max(temp_data$Temperature), 0.5))
  preds <- broom::augment(fit, newdata = new_data)
  
  # Add species information to the predictions
  preds <- preds %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  
  return(list(data = temp_data, preds = preds, est=est))
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
  select(paper_id, tpc_grp, e, eh, topt, temp, ramp_rate_per_day, cnidarian_type, latitude) %>%
  group_by(tpc_grp, paper_id) %>%
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
model00npp_r2<- lm(eh ~ log10(ramp_rate_per_day), data=model_comp) 
model00npp_r3<- lm(topt ~ log10(ramp_rate_per_day), data=model_comp) 
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
model_comp2 <- left_join(combined_est_r, tpc_r, by="tpc_grp") %>%
  select(paper_id, tpc_grp, e, eh, topt, temp, ramp_rate_per_day, cnidarian_type, latitude) %>%
  group_by(tpc_grp, paper_id) %>%
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
  summarise(avg_temp = mean(temp), 
            e=mean(e), 
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

r_final <- r %>% mutate(metabolic_category = "r")
npp_final <- npp_all %>% mutate(metabolic_category = "npp")
gpp_final <- gpp_all %>% mutate(metabolic_category = "gpp")

final_dataset <- rbind(r_final, npp_final, gpp_final)

write_csv(final_dataset, "tidy_ramping_coral_24marc25.csv")


