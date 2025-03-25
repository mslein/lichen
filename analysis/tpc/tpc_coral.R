#coral tpc
coral <- read_csv("extraction/coral dataset/coral/coral_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("log umol_o2") ~ (exp(response_value))*(1/10^6),
                         gas_units %in% c("mg_o2") ~ response_value*(1/32000),
                         gas_units %in% c("ug_o2") ~ response_value*(1/32000000),
                         gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                         gas_units %in% c("umol_o2") ~ response_value*(1/10^6),
                         gas_units %in% c("log umol_o2") ~ (response_value^10)/10^-6, 
                         gas_units %in% c("mmol_o2") ~ response_value*(1/10^3),
                         gas_units %in% c("pg_DIC") ~ response_value*(1/(1.2011*10^13)),
                         gas_units %in% c("uL_o2") ~ response_value*(4*10^-14), 
                         gas_units %in% c("umol_c") ~ response_value*(1/10^6)),
         gram = case_when(`mass/area` %in% c("g") ~ 1, 
                          `mass/area` %in% c("468 mm^2") ~ 10^3/(468), 
                          `mass/area` %in% c("l") ~ 1/1000, 
                          `mass/area` %in% c("cm^2") ~ 1, 
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
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral))%>%
  unite("tpc_grp", c(paper_unit, study_id))



coral_tpc_npp <- coral %>% filter(metabolic_category == "npp")
coral_tpc_r <- coral %>% filter(metabolic_category == "r") %>%
  mutate(mol_grammin=abs(mol_grammin))
coral_tpc_gpp <- coral %>% filter(metabolic_category == "gpp")


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
results_npp <- coral_tpc_npp %>%
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


coral_npp_plot <- ggplot(plot_data_npp, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_npp, Type == "Observed")) +
  geom_line(data = subset(plot_data_npp, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("coral_npp_plot.png", coral_npp_plot, height=15, width=20)




results_r <- coral_tpc_r %>%
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


coral_r_plot<-ggplot(plot_data_r, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_r, Type == "Observed")) +
  geom_line(data = subset(plot_data_r, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("coral_r_plot.png", coral_r_plot, height=15, width=20)



results_gpp <- coral_tpc_gpp %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

# Combine results
combined_est_gpp <- do.call(rbind, lapply(results_gpp, function(x) x$est))
combined_data_gpp <- do.call(rbind, lapply(results_gpp, function(x) x$data))
combined_preds_gpp <- do.call(rbind, lapply(results_gpp, function(x) x$preds))

plot_data_gpp <- bind_rows(
  combined_data_gpp %>% mutate(Type = "Observed"),
  combined_preds_gpp %>% mutate(Type = "Fitted")
)


coral_gpp_plot<-ggplot(plot_data_gpp, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_gpp, Type == "Observed")) +
  geom_line(data = subset(plot_data_gpp, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")



ggsave("coral_gpp_plot.png", coral_gpp_plot, height=15, width=20)


