#coral tpc
pacman::p_load(rTPC, nls.multstart, broom, tidyverse)

coral_npp <- read_csv("analysis/tidy data/coral_npp_final.csv")  
coral_r <- read_csv("analysis/tidy data/coral_r_final.csv") 
coral_gpp <- read_csv("analysis/tidy data/coral_gpp_final.csv") 

coral_npp_unweighted <- coral_npp %>% filter(is.na(sd_molgCmin))
coral_npp_weighted <- coral_npp %>% filter(!is.na(sd_molgCmin))

coral_r_unweighted <- coral_r %>% filter(is.na(sd_molgCmin))
coral_r_weighted <- coral_r %>% filter(!is.na(sd_molgCmin)) 

coral_gpp_unweighted <- coral_gpp %>% filter(is.na(sd_molgCmin))
coral_gpp_weighted <- coral_gpp %>% filter(!is.na(sd_molgCmin)) 

#creating my thermal performance curve models for raw and mean data

fit_and_predict <- function(data) {
  temp_data <- data %>%
    group_by(tpc_grp) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    mutate(Temperature = as.numeric(temp),
           rate = abs(as.numeric(mol_gCmin)))  # ensure positive rates
  
  if (nrow(temp_data) < 4) return(NULL)
  
  start_vals <- get_start_vals(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  low_lims   <- get_lower_lims(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  upper_lims <- get_upper_lims(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  
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
  
  if (is.null(fit)) return(NULL)
  
  # Extract main params
  est <- calc_params(fit) %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  # Extract SE from model summary
  fit_summary <- summary(fit)
  e_se <- fit_summary$coefficients["e", "Std. Error"]
  
  # Add Ea_se to est
  est$e_se <- e_se
  
  # Predictions
  new_data <- data.frame(Temperature = seq(min(temp_data$Temperature),
                                           max(temp_data$Temperature),
                                           by = 0.5))
  preds <- broom::augment(fit, newdata = new_data) %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  return(list(data = temp_data, preds = preds, est = est))
}


fit_and_predict_w <- function(data) {
  temp_data <- data %>%
    group_by(tpc_grp) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    mutate(Temperature = as.numeric(temp),
           rate = abs(as.numeric(mol_gCmin)), 
           sd= sd_molgCmin)  # ensure positive rates
  
  if (nrow(temp_data) < 4) return(NULL)
  
  start_vals <- get_start_vals(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  low_lims   <- get_lower_lims(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  upper_lims <- get_upper_lims(temp_data$Temperature, temp_data$rate, model_name = 'pawar_2018')
  
  fit <- tryCatch({
    nls_multstart(rate ~ pawar_2018(temp = Temperature, r_tref, e, eh, th, tref = 15),
                  data = temp_data,
                  iter = 500,
                  start_lower = start_vals - 10,
                  start_upper = start_vals + 10,
                  lower = low_lims,
                  upper = upper_lims,
                  supp_errors = "Y", 
                  modelweights = 1/(sd^2)) # weighting by the variance (which is 1/sd^2)
  }, error = function(e) return(NULL))
  
  if (is.null(fit)) return(NULL)
  
  # Extract main params
  est <- calc_params(fit) %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  # Extract SE from model summary
  fit_summary <- summary(fit)
  e_se <- fit_summary$coefficients["e", "Std. Error"]
  
  # Add Ea_se to est
  est$e_se <- e_se
  
  # Predictions
  new_data <- data.frame(Temperature = seq(min(temp_data$Temperature),
                                           max(temp_data$Temperature),
                                           by = 0.5))
  preds <- broom::augment(fit, newdata = new_data) %>%
    mutate(tpc_grp = unique(temp_data$tpc_grp))
  
  return(list(data = temp_data, preds = preds, est = est))
}

#########NPP#######

#unweighted
results_npp_unweighted <- coral_npp_unweighted %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

combined_est_npp_unweighted <- do.call(rbind, lapply(results_npp_unweighted, function(x) x$est))
combined_data_npp_unweighted <- do.call(rbind, lapply(results_npp_unweighted, function(x) x$data))
combined_preds_npp_unweighted <- do.call(rbind, lapply(results_npp_unweighted, function(x) x$preds))

plot_data_npp_unweighted <- bind_rows(
  combined_data_npp_unweighted %>% mutate(Type = "Observed"),
  combined_preds_npp_unweighted %>% mutate(Type = "Fitted")
)


coral_npp_unweighted_plot <- ggplot(plot_data_npp_unweighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_npp_unweighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_npp_unweighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/coral_npp_unweighted_plot.png", coral_npp_unweighted_plot, height=15, width=20)



#weighted
results_npp_weighted <- coral_npp_weighted %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict_w)

combined_est_npp_weighted <- do.call(rbind, lapply(results_npp_weighted, function(x) x$est))
combined_data_npp_weighted <- do.call(rbind, lapply(results_npp_weighted, function(x) x$data))
combined_preds_npp_weighted <- do.call(rbind, lapply(results_npp_weighted, function(x) x$preds))

plot_data_npp_weighted <- bind_rows(
  combined_data_npp_weighted %>% mutate(Type = "Observed"),
  combined_preds_npp_weighted %>% mutate(Type = "Fitted")
)


coral_npp_weighted_plot <- ggplot(plot_data_npp_weighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_npp_weighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_npp_weighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/coral_npp_weighted_plot.png", coral_npp_weighted_plot, height=15, width=20)


#########R#######

#unweighted
results_r_unweighted <- coral_r_unweighted %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

combined_est_r_unweighted <- do.call(rbind, lapply(results_r_unweighted, function(x) x$est))
combined_data_r_unweighted <- do.call(rbind, lapply(results_r_unweighted, function(x) x$data))
combined_preds_r_unweighted <- do.call(rbind, lapply(results_r_unweighted, function(x) x$preds))

plot_data_r_unweighted <- bind_rows(
  combined_data_r_unweighted %>% mutate(Type = "Observed"),
  combined_preds_r_unweighted %>% mutate(Type = "Fitted")
)


coral_r_unweighted_plot <- ggplot(plot_data_r_unweighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_r_unweighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_r_unweighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/coral_r_unweighted_plot.png", coral_r_unweighted_plot, height=15, width=20)


#weighted
results_r_weighted <- coral_r_weighted %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict_w)

combined_est_r_weighted <- do.call(rbind, lapply(results_r_weighted, function(x) x$est))
combined_data_r_weighted <- do.call(rbind, lapply(results_r_weighted, function(x) x$data))
combined_preds_r_weighted <- do.call(rbind, lapply(results_r_weighted, function(x) x$preds))

plot_data_r_weighted <- bind_rows(
  combined_data_r_weighted %>% mutate(Type = "Observed"),
  combined_preds_r_weighted %>% mutate(Type = "Fitted")
)


coral_r_weighted_plot <- ggplot(plot_data_r_weighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_r_weighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_r_weighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/coral_r_weighted_plot.png", coral_r_weighted_plot, height=15, width=20)


#########GPP#######

#unweighted
results_gpp_unweighted <- coral_gpp_unweighted %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict)

combined_est_gpp_unweighted <- do.call(rbind, lapply(results_gpp_unweighted, function(x) x$est))
combined_data_gpp_unweighted <- do.call(rbind, lapply(results_gpp_unweighted, function(x) x$data))
combined_preds_gpp_unweighted <- do.call(rbind, lapply(results_gpp_unweighted, function(x) x$preds))

plot_data_gpp_unweighted <- bind_rows(
  combined_data_gpp_unweighted %>% mutate(Type = "Observed"),
  combined_preds_gpp_unweighted %>% mutate(Type = "Fitted")
)


coral_gpp_unweighted_plot <- ggplot(plot_data_gpp_unweighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_gpp_unweighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_gpp_unweighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/coral_gpp_unweighted_plot.png", coral_gpp_unweighted_plot, height=15, width=20)



#weighted
results_gpp_weighted <- coral_gpp_weighted %>%
  group_by(as.factor(tpc_grp)) %>%
  group_split() %>%
  lapply(fit_and_predict_w)

combined_est_gpp_weighted <- do.call(rbind, lapply(results_gpp_weighted, function(x) x$est))
combined_data_gpp_weighted <- do.call(rbind, lapply(results_gpp_weighted, function(x) x$data))
combined_preds_gpp_weighted <- do.call(rbind, lapply(results_gpp_weighted, function(x) x$preds))

plot_data_gpp_weighted <- bind_rows(
  combined_data_gpp_weighted %>% mutate(Type = "Observed"),
  combined_preds_gpp_weighted %>% mutate(Type = "Fitted")
)


coral_gpp_weighted_plot <- ggplot(plot_data_gpp_weighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_gpp_weighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_gpp_weighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/coral_gpp_weighted_plot.png", coral_gpp_weighted_plot, height=15, width=20)



#exporting version of the estimates w/ parameters of interest
combined_est_r = rbind(combined_est_r_weighted, combined_est_r_unweighted) %>% mutate(metabolic_category="r")
combined_est_npp = rbind(combined_est_npp_weighted, combined_est_npp_unweighted) %>% mutate(metabolic_category="npp")
combined_est_gpp = rbind(combined_est_gpp_weighted, combined_est_gpp_unweighted) %>% mutate(metabolic_category="gpp")

tpc_params <- rbind(combined_est_r, combined_est_npp, combined_est_gpp)
write_csv(tpc_params, "analysis/tidy data/coral_tpc_Eas.csv")
