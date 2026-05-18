#loading libraries and packages
pacman::p_load(rTPC, nls.multstart, broom, tidyverse)

##############################
####### LICHEN ###############
##############################
#original dataset
lichen_0 <- read_csv("extraction/lichen dataset/lichen_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE)

gpp_calcs <-  read_csv("extraction/lichen dataset/lichen_gpp_arithmetic.csv") %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(-sample_size) %>%
  mutate(metabolic_category = "gpp",
         response_value = npp_value + abs(r_value),
         error = npp_error + r_error,
         sample_size = n) %>%
  select(-npp_value, -r_value, -npp_error, -r_error, -inv_T, -n, -mol, -gram, -min)

#joining for one full df
lichen_1 <- rbind(lichen_0, gpp_calcs)

#converting to 1 universal unit 
lichen <- lichen_1 %>%
  mutate(sd = case_when(error_type == "SD" ~ error,
                        error_type == "SEM" ~ error*sqrt(sample_size), 
                        error_type == "2 SEM" ~ (error/2)*sqrt(sample_size), TRUE ~ NA)) %>%
  mutate(mol_o2_conversion = case_when(gas_units %in% c("mg_co2") ~ (1/44010), # units in mol o2; molar mass of co2 ~ 44 g per mol, 1000 mg in 1 g, 1 to 1 mol conversion from co2 to o2
                            gas_units %in% c("mg_o2") ~ (1/32000), # units in mol o2; molar mass of o2 ~32 g per mol, 1000 mg in 1 g
                            gas_units %in% c("nmol_co2") ~ (1/10^9), # units in mol o2; 10^9 nanomoles in 1 mol, 1 to 1 mol conversion from co2 to o2
                            gas_units %in% c("nmol_o2") ~ (1/10^9), # units in mol o2; 10^9 nanomoles in 1 mol
                            gas_units %in% c("umol_co2") ~ (1/10^6)), # units in mol o2; 10^6 micromoles in 1 mol, 1 to 1 mol conversion from co2 to o2
         gC_conversion = case_when(response_correction %in% c("g") ~ 1*0.5, #units in g carbon; 
                        response_correction  %in% c("kg") ~ 1000*0.5, #units in g carbon; 1000 g in 1 kg
                        response_correction  %in% c("m^2") ~ 0.062*0.5, #units in g carbon; 0.062 g per m^2 in lichen
                        response_correction  %in% c("mg") ~ (1/1000)*0.5, #units in g carbon; 1000 mg per 1 g
                        response_correction  %in% c("mg chla") ~ 2.13*0.5, #units in g carbon; 2.13 mg Chl a per g biomass, 50% is Carbon
                        response_correction  %in% c("mg chl") ~ 2.13*0.5), #units in g carbon; 2.13 mg Chl a per g biomass, 50% is Carbon
         min_conversion = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60),
         mass_area = case_when(gC_conversion %in% c("g", "kg", "mg", "mg chla","mg chl") ~ "mass", 
                               gC_conversion %in% c("m") ~ "area"),
         mol_gCmin = (response_value*mol_o2_conversion)/(gC_conversion*min_conversion),
         sd_molgCmin = (sd*mol_o2_conversion)/(gC_conversion*min_conversion),
         latitude_broad=as.factor(latitude_broad), 
         latitude_broad=as.factor(latitude_broad), 
         elevation_broad=as.factor(elevation_broad), 
         lichen_type=as.factor(lichen_type))   %>%
  select(tpc_grp, lichen_type, elevation_broad, latitude, species, mol_gCmin, sd_molgCmin, temp, metabolic_category, study_id, paper_unit) 


#breaking up metabolic categories
lichen_npp_final <- lichen %>% filter(metabolic_category == "npp")
lichen_r_final <- lichen %>% filter(metabolic_category == "r")  %>% mutate(mol_grammin=abs(mol_gCmin))
lichen_gpp_final <- lichen %>% filter(metabolic_category == "gpp") 


#saving a final version of the dataframes for the lichen analysis
write_csv(lichen_npp_final, "analysis/tidy data/lichen_npp_final.csv")
write_csv(lichen_r_final, "analysis/tidy data/lichen_r_final.csv")
write_csv(lichen_gpp_final, "analysis/tidy data/lichen_gpp_final.csv")


##############################
####### CORAL ################
##############################

#original dataset
coral_0 <- read_csv("extraction/coral dataset/coral/coral_data.csv") %>%
  filter(!study_id %in% c("bahr2018", "samiei2015")) %>% #removing because of CaCO3 correction, and L corrction not being convertible
  unite("gas_units", gas_unit, gas) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  mutate(response_value = case_when(response_units == "mmoles.core ^-1 h^-1 (x10^-3)" ~ response_value*1000, #fixing the units for rodolfo metalpa 2006 being several orders off
                                    TRUE ~ response_value))
#ramping studies where ramp rate was not significant as an effect, see tpc_ramping.R for this analysis 
ramping_0 <- read_csv("extraction/coral dataset/ramping/ramping_data.csv") %>% 
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  unite("gas_units", gas_unit, gas) %>%
  select(-c(ramp_duration, ramp_duration_units, ramp_rate_per_day))


  
###back calculating different MR types for studies w/ at least 2/3 metabolic categories
#back_calcstudies <- count(coral, study_id, metabolic_category)
#write_csv(back_calcstudies, "coral_arithmetic_studies.csv")

###creating a duplicate sheet to do MR arithmetic for the studies from back_calcstudies df
#gpp_arithmetic_0 <- coral %>% filter(study_id %in% c("bahr2018", "banc-prandi2022", "godefroid2023", "higuchi2015", 
#                                                    "hill2014", "howe2001", "juillet-leclerc2014", "kemp2011", "rodolfo-metalpa2006", 
#                                                    "samiei2015"))
#write_csv(gpp_arithmetic_0, "coral_gpp_arithmetic.csv")

coral_gpp_calcs <-  read_csv("extraction/coral dataset/coral/coral_gpp_arithmetic.csv") %>%
  filter(study_id %in% c("godefroid2023", "hill2014", "howe2001", "juillet-leclerc2014", "rodolfo-metalpa2006")) %>%
  mutate(response_value = case_when(response_units == "mmoles.core ^-1 h^-1 (x10^-3)" ~ response_value*1000, #fixing the units for rodolfo metalpa 2006 being several orders off
                                    TRUE ~ response_value)) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(-sample_size, -response_value, -error) %>%
  mutate(metabolic_category = "gpp",
         response_value = npp_value + abs(r_value),
         error = npp_error + r_error,
         sample_size = n) %>%
  select(-npp_value, -r_value, -npp_error, -r_error,  -n, -gpp_value, -gpp_error, -broad_coral)

coral_npp_calcs <-  read_csv("extraction/coral dataset/coral/coral_gpp_arithmetic.csv") %>% 
  filter(study_id %in% c("higuchi2015", "kemp2011")) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(-sample_size, -response_value, -error) %>%
  mutate(metabolic_category = "gpp",
         response_value = gpp_value - abs(r_value),
         error = gpp_error - abs(r_error),
         sample_size = n) %>%
  select(-npp_value, -r_value, -npp_error, -r_error,  -n, -gpp_value, -gpp_error, -broad_coral)

ramping_gpp_calcs <- read_csv("extraction/coral dataset/ramping/ramping_mrcalcs24mar25.csv") %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  unite("gas_units", gas_unit, gas) %>%
  select(-sample_size, -error) %>%
  mutate(metabolic_category = "gpp",
         response_value = npp_value + abs(r_value),
         error = npp_error + r_error,
         sample_size = n) %>%
  select(-npp_value, -r_value, -npp_error, -r_error,  -n, -gpp_value, -gpp_error)%>%
  select(-c(ramp_duration, ramp_duration_units, ramp_rate_per_day))


ramping_npp_calcs <- read_csv("extraction/coral dataset/ramping/ramping_mrcalcs24mar25.csv") %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  unite("gas_units", gas_unit, gas) %>%
  select(-sample_size, -error) %>%
  mutate(metabolic_category = "gpp",
         response_value = gpp_value - abs(r_value),
         error = gpp_error - abs(r_error),
         sample_size = n) %>%
  select(-npp_value, -r_value, -npp_error, -r_error,  -n, -gpp_value, -gpp_error)%>%
  select(-c(ramp_duration, ramp_duration_units, ramp_rate_per_day))

#uniting all these into one df
coral_1 <-rbind(coral_0, ramping_0, coral_gpp_calcs, coral_npp_calcs, ramping_gpp_calcs, ramping_npp_calcs)

count(coral_1, `mass/area`)

coral_1 %>%
  filter(gas_units == "pg_DIC") %>%
  count(study_id)

#doing the conversion to a single unit
coral <-  coral_1 %>%
  mutate(sd = case_when(error_type == "SD" ~ error,
                        error_type == "SEM" ~ error*sqrt(sample_size), 
                        error_type %in% c("95 CI", "95% CI") ~ (error/1.96)*sqrt(sample_size), TRUE ~ NA)) %>%
    mutate(mol_o2_conversion = case_when(gas_units %in% c("mg_o2", "mg_O2") ~ (1/32000), # units in mol o2; 1000 mg in 1g, molar mass is 32 g / mol
                         gas_units %in% c("ug_o2") ~ (1/32000000), # units in mol o2; 1 million ug in 1g, molar mass is 32 g / mol
                         gas_units %in% c("nmol_o2") ~ (1/10^9), # units in mol o2; 10^9 nanomoles in 1 mole
                         gas_units %in% c("umol_o2", "umol_O2") ~ (1/10^6), # units in mol o2; 10^6 micromoles in 1 mol
                         gas_units %in% c("mmol_o2") ~ (1/10^3), # units in mol o2; 10^6 micromoles in 1 mol
                         gas_units %in% c("pg_DIC") ~ (1/(1.2011*10^13)),# units in mol o2; 10^12 pg in 1 g, molar mass of c ~12.011 g per mol, 1 mol C to 1 mol O2 
                         gas_units %in% c("umol_c") ~ (1/10^6)), #units in mol o2; 1 for 1 mol C to o2, 1 million micromoles in 1mol
         gC_conversion = case_when(`mass/area` %in% c("g") ~ 1*0.5, 
                        `mass/area` %in% c("mm^2") ~ (5*41)/100, # units in cm^2 = g; mm^2 (1cm^2 / 100mm^2), 5 ug Chl a per cm2 in coral, ratio is 44 C to 1 Chla  
                        `mass/area` %in% c("cm^2") ~ (5*41), # units in cm^2 = g
                        `mass/area` %in% c("m^2") ~ (5*41)*10000,  # units in cm^2 = g; m^2 (10000cm^2 / 1m^2), 5 ug Chl a per cm2 in coral, ratio is 44 C to 1 Chla
                        `mass/area` %in% c("larvae^-1", "cell") ~  152.9*(1/10^9),  #152.9 ng C per cell, covert to from nanograms to grams,
                        `mass/area` %in% c("ug chla") ~ 44*(1/10^6), # ratios 44 C to 1 Chl a, scaling from ug to g divide by 1 million
                        `mass/area` %in% c("g ADFW", "g") ~ 1*0.5, # units in g; half of weight in g is carbon
                        `mass/area` %in% c("mg") ~ 0.5*1000, # units in g; multiply by 1000, half of weight in g is carbon
                        `mass/area` %in% c("zoox x 10^6") ~ (152.9*(1/10^9)*10^6)), #152.9 ng C per cell, covert to from nanograms to grams, scale up to 1 million cells
         min_conversion = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60, 
                         time %in% c("day") ~ 1440),
         mol_gCmin= (response_value*mol_o2_conversion)/(gC_conversion*min_conversion),
         sd_molgCmin = (sd*mol_o2_conversion)/(gC_conversion*min_conversion),
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral))%>%
  dplyr::select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, sd_molgCmin, temp, cnidarian_type, metabolic_category, study_id)


#breaking up metabolic categories
coral_npp_all <- coral %>% filter(metabolic_category == "npp") 
coral_r_all <- coral %>% filter(metabolic_category == "r")  %>% mutate(mol_grammin=abs(mol_gCmin))
coral_gpp_all <- coral %>% filter(metabolic_category == "gpp") 


#saving a final version of the dataframes for the lichen analysis
write_csv(coral_npp_all, "analysis/tidy data/coral_npp_final.csv")
write_csv(coral_r_all, "analysis/tidy data/coral_r_final.csv")
write_csv(coral_gpp_all, "analysis/tidy data/coral_gpp_final.csv")




