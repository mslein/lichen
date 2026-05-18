######lichen tpc 
pacman::p_load(tidyverse, rTPC, nls.multstart)
lichen_npp <- read_csv("analysis/tidy data/lichen_npp_final.csv") 
lichen_r <- read_csv("analysis/tidy data/lichen_r_final.csv") 
lichen_gpp <- read_csv("analysis/tidy data/lichen_gpp_final.csv") 


lichen_npp_unweighted <- lichen_npp %>% filter(is.na(sd_molgCmin))
lichen_npp_weighted <- lichen_npp %>% filter(!is.na(sd_molgCmin))

lichen_r_unweighted <- lichen_r %>% filter(is.na(sd_molgCmin))
lichen_r_weighted <- lichen_r %>% filter(!is.na(sd_molgCmin)) 

lichen_gpp_unweighted <- lichen_gpp %>% filter(is.na(sd_molgCmin))
lichen_gpp_weighted <- lichen_gpp %>% filter(!is.na(sd_molgCmin)) 



lichen_gpp_unweighted %>%
  filter(study_id == "tretiach1997") %>%
  ggplot(aes(x=temp, y=mol_gCmin, colour=tpc_grp))+
  geom_point()+
  facet_wrap(~tpc_grp, scales = "free")+
  theme(legend.position = "none")



#creating my thermal performance curve models for raw and mean data

fit_and_predict <- function(data) {
  grp <- unique(as.character(data$tpc_grp))
  grp <- grp[!is.na(grp)][1]
  
  message("--------------------------------------------------")
  message("Processing: ", grp)
  
  temp_data <- data %>%
    dplyr::mutate(tpc_grp = as.character(tpc_grp)) %>%
    dplyr::group_by(tpc_grp) %>%
    dplyr::filter(dplyr::n() >= 5) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Temperature = as.numeric(temp),
      rate = abs(as.numeric(mol_gCmin))
    ) %>%
    dplyr::filter(!is.na(Temperature), !is.na(rate))
  
  message("Rows after cleaning: ", nrow(temp_data))
  message("Unique temperatures: ", dplyr::n_distinct(temp_data$Temperature))
  
  if (nrow(temp_data) < 4 || dplyr::n_distinct(temp_data$Temperature) < 4) {
    return(list(
      data  = temp_data,
      preds = tibble::tibble(),
      est   = tibble::tibble(tpc_grp = grp, e = NA_real_, e_se = NA_real_)
    ))
  }
  
  message("Getting start values...")
  start_vals <- tryCatch(
    get_start_vals(temp_data$Temperature, temp_data$rate, model_name = "pawar_2018"),
    error = function(e) NULL
  )
  message("Getting lower limits...")
  low_lims <- tryCatch(
    get_lower_lims(temp_data$Temperature, temp_data$rate, model_name = "pawar_2018"),
    error = function(e) NULL
  )
  message("Getting upper limits...")
  upper_lims <- tryCatch(
    get_upper_lims(temp_data$Temperature, temp_data$rate, model_name = "pawar_2018"),
    error = function(e) NULL
  )
  
  if (is.null(start_vals) || is.null(low_lims) || is.null(upper_lims)) {
    message("Failed starts/limits for ", grp)
    return(list(
      data  = temp_data,
      preds = tibble::tibble(),
      est   = tibble::tibble(tpc_grp = grp, e = NA_real_, e_se = NA_real_)
    ))
  }
  
  message("Attempting nls_multstart fit...")
  fit <- tryCatch(
    {
      nls_multstart(
        rate ~ pawar_2018(temp = Temperature, r_tref, e, eh, th, tref = 15),
        data = temp_data,
        iter = 500,
        start_lower = start_vals - 10,
        start_upper = start_vals + 10,
        lower = low_lims,
        upper = upper_lims,
        supp_errors = "Y"
      )
    },
    error = function(e) {
      message("Fit failed for ", grp, ": ", e$message)
      return(NULL)
    }
  )
  message("nls_multstart step finished")
  
  if (is.null(fit)) {
    return(list(
      data  = temp_data,
      preds = tibble::tibble(),
      est   = tibble::tibble(tpc_grp = grp, e = NA_real_, e_se = NA_real_)
    ))
  }
  
  # first try calc_params; if that fails, fall back to raw coefficient extraction
  message("Extracting params...")
  est <- tryCatch(
    {
      calc_params(fit) %>%
        dplyr::mutate(tpc_grp = grp)
    },
    error = function(e) {
      message("calc_params failed for ", grp, ": ", e$message)
      
      coefs <- tryCatch(stats::coef(fit), error = function(e) NULL)
      
      tibble::tibble(
        tpc_grp = grp,
        e = if (!is.null(coefs) && "e" %in% names(coefs)) unname(coefs["e"]) else NA_real_
      )
    }
  )
  message("Param extraction done")
  
  message("Extracting e SE...")
  e_se <- tryCatch(
    {
      fit_summary <- summary(fit)
      if ("e" %in% rownames(fit_summary$coefficients)) {
        fit_summary$coefficients["e", "Std. Error"]
      } else {
        NA_real_
      }
    },
    error = function(e) {
      message("SE extraction failed for ", grp, ": ", e$message)
      NA_real_
    }
  )
  message("e SE extraction done")
  
  if (!"e" %in% names(est)) {
    est$e <- NA_real_
  }
  est$e_se <- e_se
  
  message("Building predictions...")
  preds <- tryCatch(
    {
      new_data <- data.frame(
        Temperature = seq(
          min(temp_data$Temperature, na.rm = TRUE),
          max(temp_data$Temperature, na.rm = TRUE),
          by = 0.5
        )
      )
      
      pred_vals <- predict(fit, newdata = new_data)
      
      new_data %>%
        dplyr::mutate(
          .fitted = as.numeric(pred_vals),
          tpc_grp = grp
        )
    },
    error = function(e) {
      message("Prediction failed for ", grp, ": ", e$message)
      tibble::tibble()
    }
  )
  message("Prediction step done")
  
  message("Finished: ", grp,
          " | obs rows = ", nrow(temp_data),
          " | pred rows = ", nrow(preds),
          " | est rows = ", nrow(est))
  
  list(data = temp_data, preds = preds, est = est)
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
results_npp_unweighted <- lichen_npp_unweighted %>%
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


lichen_npp_unweighted_plot <- ggplot(plot_data_npp_unweighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_npp_unweighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_npp_unweighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/lichen_npp_unweighted_plot.png", lichen_npp_unweighted_plot, height=15, width=20)



#weighted
results_npp_weighted <- lichen_npp_weighted %>%
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


lichen_npp_weighted_plot <- ggplot(plot_data_npp_weighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_npp_weighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_npp_weighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/lichen_npp_weighted_plot.png", lichen_npp_weighted_plot, height=15, width=20)


#########R#######

#unweighted
results_r_unweighted <- lichen_r_unweighted %>%
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


lichen_r_unweighted_plot <- ggplot(plot_data_r_unweighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_r_unweighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_r_unweighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/lichen_r_unweighted_plot.png", lichen_r_unweighted_plot, height=15, width=20)


#weighted ---> none of the weights applied here, so leaving this out for lichen


#########GPP#######

#unweighted
results_gpp_unweighted <- lichen_gpp_unweighted %>%
  dplyr::mutate(tpc_grp = as.character(tpc_grp)) %>%
  dplyr::group_by(tpc_grp) %>%
  dplyr::group_split() %>%
  lapply(fit_and_predict)

combined_est_gpp_unweighted <- dplyr::bind_rows(lapply(results_gpp_unweighted, function(x) x$est))
combined_data_gpp_unweighted <- dplyr::bind_rows(lapply(results_gpp_unweighted, function(x) x$data))
combined_preds_gpp_unweighted <- dplyr::bind_rows(lapply(results_gpp_unweighted, function(x) x$preds))

plot_data_gpp_unweighted <- bind_rows(
  combined_data_gpp_unweighted %>% mutate(Type = "Observed"),
  combined_preds_gpp_unweighted %>% mutate(Type = "Fitted")
)


lichen_gpp_unweighted_plot <- ggplot(plot_data_gpp_unweighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_gpp_unweighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_gpp_unweighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/lichen_gpp_unweighted_plot.png", lichen_gpp_unweighted_plot, height=15, width=20)



#weighted
results_gpp_weighted <- lichen_gpp_weighted %>%
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


lichen_gpp_weighted_plot <- ggplot(plot_data_gpp_weighted, aes(x = Temperature, y = rate, color=tpc_grp)) +
  #scale_color_manual(values=c( "#a559aa", "#59a89c", "#f0c571", "#ae282c","#082a54" ))+
  geom_point(data = subset(plot_data_gpp_weighted, Type == "Observed")) +
  geom_line(data = subset(plot_data_gpp_weighted, Type == "Fitted"), aes(y=.fitted)) +
  facet_wrap(~tpc_grp, scales="free_y") +
  theme_classic() +
  #labs(x = 'Temperature (ºC)', y = 'Dispersal rate at leading edge', color = 'Species') +
  theme(legend.position = "none")

ggsave("figures/qc/lichen_gpp_weighted_plot.png", lichen_gpp_weighted_plot, height=15, width=20)



#exporting version of the estimates w/ parameters of interest
combined_est_r = combined_est_r_unweighted %>% mutate(metabolic_category="r")
combined_est_npp = rbind(combined_est_npp_weighted, combined_est_npp_unweighted) %>% mutate(metabolic_category="npp")
combined_est_gpp = rbind(combined_est_gpp_weighted, combined_est_gpp_unweighted) %>% mutate(metabolic_category="gpp")

tpc_params <- rbind(combined_est_r, combined_est_npp, combined_est_gpp)
write_csv(tpc_params, "analysis/tidy data/lichen_tpc_Eas.csv")





