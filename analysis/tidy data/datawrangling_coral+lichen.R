#loading libraries and packages
pacman::p_load(rTPC, nls.multstart, broom, tidyverse)

##############################
####### LICHEN ###############
##############################
#original dataset
lichen_0 <- read_csv("extraction/lichen dataset/lichen_data.csv") %>%
  unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mol_o2 = case_when(gas_units %in% c("mg_co2") ~ response_value*(1/44010), 
                            gas_units %in% c("mg_o2") ~ response_value*(1/32000), 
                            gas_units %in% c("nmol_co2") ~ response_value*(1/10^9), 
                            gas_units %in% c("nmol_o2") ~ response_value*(1/10^9),
                            gas_units %in% c("umol_co2") ~ response_value*(1/10^6)),
         gC = case_when(response_correction %in% c("g") ~ 1*0.5, 
                        response_correction  %in% c("kg") ~ 1000*0.5, 
                        response_correction  %in% c("m") ~ 0.062*0.5,
                        response_correction  %in% c("mg") ~ (1/1000)*0.5, 
                        response_correction  %in% c("mg chla") ~ 2.13*0.5, 
                        response_correction  %in% c("mg chl") ~ 2.13*0.5),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60),
         mass_area = case_when(gC %in% c("g", "kg", "mg", "mg chla","mg chl") ~ "mass", 
                               gC %in% c("m") ~ "area"),
         mol_gCmin= mol_o2/(gC*min),
         latitude_broad=as.factor(latitude_broad), 
         latitude_broad=as.factor(latitude_broad), 
         elevation_broad=as.factor(elevation_broad), 
         lichen_type=as.factor(lichen_type), ) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE)  %>%
  select(tpc_grp, lichen_type, elevation_broad, latitude, species, mol_gCmin, inv_T, temp, metabolic_category, study_id, paper_unit) 


###back calculating different MR types for studies w/ at least 2/3 metabolic categories
#back_calcstudies <- count(lichen, study_id, metabolic_category)
#write_csv(back_calcstudies, "lichen_arithmetic_studies.csv")

###creating a duplicate sheet to do MR arithmetic for the studies from back_calcstudies df
#gpp_arithmetic_0 <- lichen_0 %>% filter(study_id %in% c("adams_1971", "aubert_2007", "brownkershaw1984", 
#         "colesie2018", "coxson_2003", "delprado_1995", 
#        "domaschke_2013", "eickmeieradams1973", "kappen1983", 
#        "kappen_2000", "kershaw1977", "kershaw1983", 
#        "lange2000", "lange2004", "lange_1980", "lange_1991", 
#        "larson1980", "lechowicz1973", "lechowiczadams1974", 
#        "macfarlane1983", "reiter2000", "sancho1997", 
#        "sanchokappen1989", "schroeter1995", "smith2001", 
#        "sonesson1989", "sundberg1997", "tretiach1997", 
#        "uchida2006", "zotz1998"))
#write_csv(gpp_arithmetic_0, "lichen_gpp_arithmetic.csv")

#gpp arithmetic
gpp_calcs <-  read_csv("lichen_gpp_arithmetic.csv") %>%
  mutate(gpp= npp_value + abs(r_value)) %>%
  #unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)), 
         metabolic_category = "gpp",
         mol_o2 = case_when(gas_units %in% c("mg_co2") ~ gpp*(1/44010), 
                            gas_units %in% c("nmol_co2") ~ gpp*(1/10^9),
                            gas_units %in% c("umol_co2") ~ gpp*(1/10^6)),
         gC = case_when(response_correction %in% c("g") ~ 1*0.5, 
                        response_correction  %in% c("kg") ~ 1000*0.5, 
                        response_correction  %in% c("m") ~ 0.062*0.5,
                        response_correction  %in% c("mg") ~ (1/1000)*0.5, 
                        response_correction  %in% c("mg chla") ~ 2.13*0.5, 
                        response_correction  %in% c("mg chl") ~ 2.13*0.5),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60),
         mass_area = case_when(gC %in% c("g", "kg", "mg", "mg chla","mg chl") ~ "mass", 
                               gC %in% c("m") ~ "area"),
         mol_gCmin= mol_o2/(gC*min),
         latitude_broad=as.factor(latitude_broad), 
         elevation_broad=as.factor(elevation_broad), 
         lichen_type=as.factor(lichen_type)) %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  select(tpc_grp, lichen_type, elevation_broad, latitude, species, mol_gCmin, inv_T, temp, metabolic_category, study_id, paper_unit)


#breaking up metabolic categories
lichen_npp_final <- lichen_0 %>% filter(metabolic_category == "npp")
lichen_r_final <- lichen_0 %>% filter(metabolic_category == "r")  
lichen_gpp_final <- lichen_0 %>% filter(metabolic_category == "gpp") %>% rbind(gpp_calcs)



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
         gC = case_when(`mass/area` %in% c("g") ~ 1*0.5, 
                        `mass/area` %in% c("mm^2") ~ (5*44)/100, 
                        `mass/area` %in% c("cm^2") ~ (5*44), 
                        `mass/area` %in% c("larvae^-1") ~  152.9*10^-9, 
                        `mass/area` %in% c("cell") ~ 152.9*10^-9, 
                        `mass/area` %in% c("ug chla") ~ 44*10^-6, 
                        `mass/area` %in% c("g ADFW") ~ 1*0.5,
                        `mass/area` %in% c("mg") ~ 0.5*0.001,
                        `mass/area` %in% c("zoox x 10^6") ~ ( 152.9*10^-9)*10^6),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1, 
                         time %in% c("sec") ~ 1/60, 
                         time %in% c("day") ~ 1440),
         mol_gCmin= mol/(gC*min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral))%>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  dplyr::select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, cnidarian_type, metabolic_category, study_id)
#adding in the ramping dataset
coral_ramping <- read_csv("analysis/tidy data/tidy_ramping_coral_24marc25.csv")
coral <- rbind(coral_0, coral_ramping)

###back calculating different MR types for studies w/ at least 2/3 metabolic categories
#back_calcstudies <- count(coral, study_id, metabolic_category)
#write_csv(back_calcstudies, "coral_arithmetic_studies.csv")

###creating a duplicate sheet to do MR arithmetic for the studies from back_calcstudies df
#gpp_arithmetic_0 <- coral %>% filter(study_id %in% c("bahr2018", "banc-prandi2022", "godefroid2023", "higuchi2015", 
#                                                    "hill2014", "howe2001", "juillet-leclerc2014", "kemp2011", "rodolfo-metalpa2006", 
#                                                    "samiei2015"))
#write_csv(gpp_arithmetic_0, "coral_gpp_arithmetic.csv")

#breaking up metabolic categories
npp <- coral %>% filter(metabolic_category == "npp") 
r <- coral %>% filter(metabolic_category == "r")  
gpp <- coral %>% filter(metabolic_category == "gpp") 


coral_r_all <- r %>% mutate(mol_grammin=abs(mol_gCmin))

#gpp and npp arithmetic
#coral
gpp_calcs <-  read_csv("coral_gpp_arithmetic.csv") %>%
  filter(study_id %in% c("godefroid2023", "hill2014", "howe2001", "juillet-leclerc2014", "rodolfo-metalpa2006")) %>%
  mutate(gpp= npp_value + abs(r_value))%>%
  #unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("mg_o2") ~ gpp*(1/32000),
                         gas_units %in% c("mmol_o2") ~ gpp*(1/10^3),
                         gas_units %in% c("ug_o2") ~ gpp*(1/32000000),
                         gas_units %in% c("umol_o2") ~ gpp*(1/10^6)),
         gC = case_when(`mass/area` %in% c("468 mm^2") ~ (5*44)/100,
                        `mass/area` %in% c("chla") ~ 44*10^-6, 
                        `mass/area` %in% c("cm^2") ~ (5*44)),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("sec") ~ 1/60),
         mol_gCmin= mol/(gC*min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral), 
         metabolic_category = "gpp") %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  dplyr::select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, cnidarian_type, metabolic_category, study_id)

coral_gpp_all <- rbind(gpp, gpp_calcs)

#npp
npp_calcs <-  read_csv("coral_gpp_arithmetic.csv") %>% 
  filter(study_id %in% c("higuchi2015", "kemp2011")) %>%
  mutate(npp = gpp_value - abs(r_value)) %>%
  #unite("gas_units", gas_unit, gas) %>%
  mutate(inv_T=(1/((8.617333262145*10^-5)*(temp+273.15))),
         mol = case_when(gas_units %in% c("mg_o2") ~ npp*(1/32000),
                         gas_units %in% c("umol_o2") ~ npp*(1/10^6)),
         gC = case_when(`mass/area` %in% c("cm^2") ~ (5*44)/100),
         min = case_when(time %in% c("hr") ~ 60, 
                         time %in% c("min") ~ 1),
         mol_gCmin= mol/(gC*min), 
         depth_broad = case_when(avg_depth_m < 20 ~ "shallow", 
                                 TRUE ~ "deep"), 
         depth_broad = as.factor(depth_broad), 
         cnidarian_type = as.factor(cnidarian_type), 
         broad_coral = case_when(cnidarian_type == "soft coral" ~ "soft", 
                                 TRUE ~ "hard"), 
         broad_coral =as.factor(broad_coral), 
         metabolic_category = "npp") %>%
  unite("tpc_grp", c(paper_unit, study_id), remove=FALSE) %>%
  dplyr::select(tpc_grp, depth_broad, broad_coral, latitude, species, mol_gCmin, inv_T, temp, cnidarian_type, metabolic_category, study_id) 


coral_npp_all <- rbind(npp, npp_calcs)


#saving a final version of the dataframes for the lichen analysis
write_csv(coral_npp_all, "analysis/tidy data/coral_npp_final.csv")
write_csv(coral_r_all, "analysis/tidy data/coral_r_final.csv")
write_csv(coral_gpp_all, "analysis/tidy data/coral_gpp_final.csv")




