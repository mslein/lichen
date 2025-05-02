#coral tpc
pacman::p_load(rTPC, nls.multstart, broom, tidyverse)

coral_npp_all <- read_csv("analysis/tidy data/coral_npp_final.csv")
coral_r_all <- read_csv("analysis/tidy data/coral_r_final.csv")
coral_gpp_all <- read_csv("analysis/tidy data/coral_gpp_final.csv")


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


#########NPP#######
results_npp <- coral_npp_all %>%
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




results_r <- coral_r_all %>%
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



results_gpp <- coral_gpp_all %>%
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

#exporting version of the estimates w/ parameters of interest
combined_est_r = combined_est_r %>% mutate(metabolic_category="r")
combined_est_npp = combined_est_npp %>% mutate(metabolic_category="npp")
combined_est_gpp = combined_est_gpp %>% mutate(metabolic_category="gpp")

tpc_params <- rbind(combined_est_r, combined_est_npp, combined_est_gpp)
write_csv(tpc_params, "analysis/tidy data/coral_tpc_Eas.csv")
