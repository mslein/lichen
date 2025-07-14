isl_algae <- read_csv("extraction/coral dataset/isolated algae/isolated_algae.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  filter(!gas_units %in% c("pg_c", "ug_C")) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol_o2 = case_when(gas_units %in% c("ug_o2") ~ response_value*(1/32000000), 
                            gas_units %in% c("nmol_co2") ~ response_value*(1/10^9),
                            gas_units %in% c("mmol_o2") ~ response_value*(1/10^3),
                            gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                            gas_units %in% c("fmol_o2") ~ response_value*(1/10^15),
                            gas_units %in% c("umol_o2") ~ response_value*(1/10^6), 
                            gas_units %in% c("pmol_o2") ~ response_value*(1/10^12)),
         gC = case_when(response_correction %in% c("cell") ~ 152.9*10^-9, 
                        response_correction %in% c("cell x 10^9 ") ~ (152.9*10^-9)*10^9,
                        response_correction %in% c("cell x10^6") ~ (152.9*10^-9)*10^6, 
                        response_correction %in% c("cell x10^9") ~ (152.9*10^-9)*10^9,
                        response_correction %in% c("ug chla") ~ 44*10^-6, 
                        response_correction %in% c("ug chl a") ~ 44*10^-6,
                        response_correction %in% c("mg chla") ~ 44*10^-3),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1),
         mol_gCmin= mol_o2/(gC*min)) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE)  %>%
  select(tpc_grp, latitude, species, mol_gCmin, inv_T, temp, metabolic_category, study_id, paper_unit) 

isl_algae_r_all <- isl_algae %>% filter(metabolic_category == "r")
isl_algae_gpp <- isl_algae %>% filter(metabolic_category == "gpp")
isl_algae_npp <- isl_algae %>% filter(metabolic_category == "npp")



gpp_calcs <- read_csv("extraction/coral dataset/isolated algae/isolated_algae_arithmetic_jul2025.csv") %>%
filter(study_id %in% c("goulet2005", "gregoire2017", "iglesios-prieto1992", "leggat2004", "	
perez2001", "russnak2021", "yang2022")) %>%
  select(-mol_gCmin) %>%
  mutate(mol_gCmin= npp_value + abs(r_value), 
         metabolic_category = "gpp") %>%
  select(-r_value, -npp_value, -gpp_value) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE)  

npp_calcs <-  read_csv("extraction/coral dataset/isolated algae/isolated_algae_arithmetic_jul2025.csv") %>%
  filter(study_id %in% c("pierangelini2020")) %>%
  select(-mol_gCmin) %>%
  mutate(mol_gCmin= gpp_value - abs(r_value), 
         metabolic_category = "npp")%>%
  select(-r_value, -npp_value, -gpp_value) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE)  

isl_algae_npp_all <- rbind(isl_algae_npp, npp_calcs)
isl_algae_gpp_all <- rbind(isl_algae_gpp, gpp_calcs)

#saving a final version of the dataframes for the analysis
write_csv(isl_algae_npp_all, "analysis/tidy data/isl_algae_npp_final.csv")
write_csv(isl_algae_r_all, "analysis/tidy data/isl_algae_r_final.csv")
write_csv(isl_algae_gpp_all, "analysis/tidy data/isl_algae_gpp_final.csv")

