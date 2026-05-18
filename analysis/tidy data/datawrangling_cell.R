#loading libraries and packages
pacman::p_load(rTPC, nls.multstart, broom, tidyverse)

##############################
####### SINGLE CELL ###############
##############################
#original dataset
cell_0 <- read_csv("extraction/single cell/cell_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) 

cell_gpp_calcs <-  read_csv("extraction/single cell/cell_arithmetic.csv") %>%
  filter(study_id != "deldicq2025") %>%
  unite("gas_units", gas_unit, gas) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(-response_value, -error) %>%
  mutate(metabolic_category = "gpp",
         response_value = npp_value + abs(r_value),
         error= NA) %>%
  select(-npp_value, -r_value, -npp_error, -r_error, -gpp_value, -gpp_error)

cell_npp_calcs <-  read_csv("extraction/single cell/cell_arithmetic.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  filter(study_id == "deldicq2025") %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(-response_value, -error) %>%
  mutate(metabolic_category = "npp",
         response_value = npp_value + abs(r_value),
         error = sqrt((npp_error*sqrt(sample_size))^2 + (r_error*sqrt(sample_size))^2)) %>%
  select(-npp_value, -r_value, -npp_error, -r_error, -gpp_value, -gpp_error)

#uniting all these into one df
cell_1 <-rbind(cell_0, cell_gpp_calcs, cell_npp_calcs)

#doing the conversion to a single unit
cell <-  cell_1 %>%
  mutate(sd = case_when(error_type == "SD" ~ error,
                        error_type == "SEM" ~ error*sqrt(sample_size),
                        TRUE ~ NA)) %>%
  mutate(mol_o2_conversion = case_when(gas_units %in% c("ug  mL^-1_O2") ~ (15/32000000), # units in mol o2; 15 mL flasks were used and included to cancel mL, 1 million ug in 1g, molar mass is 32 g / mol
                                       gas_units %in% c("ug L^-1_O2") ~ (0.06/32000000), # units in mol o2; 60 mL flasks were used and included to cancel the L, 1 million ug in 1g, molar mass is 32 g / mol
                                       gas_units %in% c("ug_O2") ~ (1/32000000), # units in mol o2; 1 million ug in 1g, molar mass is 32 g / mol
                                       gas_units %in% c("nmol_O2") ~ (1/10^9), # units in mol o2; 10^9 nanomoles in 1 mole, 
                                       gas_units %in% c("pmol_O2") ~ (1/10^12)), # units in mol o2; 10^12 picomoles in 1 mole, 
         gC_conversion = case_when(`mass/area` %in% c("mm^2") ~ (5*41)/100, # units in cm^2 = g; mm^2 (1cm^2 / 100mm^2), 5 ug Chl a per cm2 in coral, ratio is 44 C to 1 Chla --> assuming coral context holds for this
                                   `mass/area` %in% c("mg") ~ 0.5*1000, # units in g; multiply by 1000, half of weight in g is carbon
                                   `mass/area` %in% c("um^3", "indiv") ~ 0.113/10^6), # units in g; Foraminifera have coversion from biovolume (um) to g C of 0.113
         min_conversion = case_when(time %in% c("hr") ~ 60, 
                                    time %in% c("day") ~ 1440,
                                    time %in% c("week") ~ 10080),
         mol_gCmin= (response_value*mol_o2_conversion)/(gC_conversion*min_conversion),
         sd_molgCmin = (sd*mol_o2_conversion)/(gC_conversion*min_conversion)) %>%
  dplyr::select(tpc_grp, latitude, species, mol_gCmin, sd_molgCmin, temp, metabolic_category, study_id)


#breaking up metabolic categories
cell_npp_all <- cell %>% filter(metabolic_category == "npp") 
cell_r_all <- cell %>% filter(metabolic_category == "r")  %>% mutate(mol_grammin=abs(mol_gCmin))
cell_gpp_all <- cell %>% filter(metabolic_category == "gpp") 


#saving a final version of the dataframes for the lichen analysis
write_csv(cell_npp_all, "analysis/tidy data/cell_npp_final.csv")
write_csv(cell_r_all, "analysis/tidy data/cell_r_final.csv")
write_csv(cell_gpp_all, "analysis/tidy data/cell_gpp_final.csv")

min(cell_npp_all$temp)



