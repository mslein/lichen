######lichen tpc 

#studies w/ less than 4 temps
# - colesie 2018
# - insarov: compiles one off temps from 3-4 studies (would not qualify for TPC)
# - coxson 2003 --> several of these have just 3 temps
# - lange 1980 --> only has 3 temps
# brown & kershaw 1984 --> a few have only 3 temps 
# larson 1980 --> only has 3 temps
# dodds 2007 --> only has 3 temps


# other study questions:
# in the coral analysis: becker 2021 has an enriched trt that I did not extract, should i?
#for hadjioannou 2019 i did include a high nutrient?
#jones 1998 --> exposure i included a 1 hour and 4 hour exposure...


#exclude nutrients
#fit arrhenius function for respiration (rather than sharpe-schoolfield)
#random effect of study 
#include average assay temperature 


#run a separate analysis for the ramping studies, get Ea, and include ramping exp trts as explanatory variables




## how should I be dealing with the fact that some of the older studies are just raw points and not averages?


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
         centre_celsius = I(temp-mean(temp))) %>%
  unite("tpc_grp", c(paper_unit, study_id))


lichen_tpc_npp <- lichen %>% filter(metabolic_category == "npp")
lichen_tpc_r <- lichen %>% filter(metabolic_category == "r") %>%
  mutate(mol_grammin=abs(mol_grammin))
lichen_tpc_gpp <- lichen %>% filter(metabolic_category == "gpp")


fit_and_predict <- function(data) {
  
  temp_data <- data %>%
    mutate(Temperature = as.numeric(temp), 
           rate=mol_grammin) # just in case integer causes problems
  
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


#########NPP#######
results_npp <- lichen_tpc_npp %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

# Combine results
combined_est_npp <- do.call(rbind, lapply(results_npp, function(x) x$est))
combined_data_npp <- do.call(rbind, lapply(results_npp, function(x) x$data))
combined_preds_npp <- do.call(rbind, lapply(results_npp, function(x) x$preds))

plot_data_npp <- bind_rows(
  combined_data_npp %>% mutate(Type = "Observed"),
  combined_preds_npp %>% mutate(Type = "Fitted")
)


lichen_npp_plot <- ggplot(plot_data_npp, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_npp, Type == "Observed")) +
  geom_line(data = subset(plot_data_npp, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("lichen_npp_plot.png", lichen_npp_plot, height=15, width=20)




results_r <- lichen_tpc_r %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

# Combine results
combined_est_r <- do.call(rbind, lapply(results_r, function(x) x$est))
combined_data_r <- do.call(rbind, lapply(results_r, function(x) x$data))
combined_preds_r <- do.call(rbind, lapply(results_r, function(x) x$preds))

plot_data_r <- bind_rows(
  combined_data_r %>% mutate(Type = "Observed"),
  combined_preds_r %>% mutate(Type = "Fitted")
)


lichen_r_plot<-ggplot(plot_data_r, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_r, Type == "Observed")) +
  geom_line(data = subset(plot_data_r, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("lichen_r_plot.png", lichen_r_plot, height=15, width=20)



results_gpp <- lichen_tpc_gpp %>%
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


lichen_gpp_plot<-ggplot(plot_data_gpp, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_gpp, Type == "Observed")) +
  geom_line(data = subset(plot_data_gpp, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")



ggsave("lichen_gpp_plot.png", lichen_gpp_plot, height=15, width=20)
