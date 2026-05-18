#loading necessary packages
pacman::p_load(tidyverse, stringr, readr)

### lichen

gpp_lichen <- read_csv("analysis/tidy data/lichen_gpp_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, metabolic_category) %>%
  separate(tpc_grp, into =c("rep", "study_id"), remove=FALSE) %>%
  mutate(coupling_unit = interaction(study_id, species, latitude, elevation_broad)) %>% 
  distinct(coupling_unit, .keep_all = TRUE)
r_lichen  <- read_csv("analysis/tidy data/lichen_r_final.csv") %>% select(tpc_grp, elevation_broad, lichen_type, latitude, species, metabolic_category) %>%
  separate(tpc_grp, into =c("rep", "study_id"), remove=FALSE) %>%
  mutate(coupling_unit = interaction(study_id, species, latitude, elevation_broad)) %>% 
  distinct(coupling_unit, .keep_all = TRUE)
ea_lichen  <- read_csv("analysis/tidy data/lichen_tpc_Eas.csv") %>% filter(metabolic_category != "npp") %>%
  select(e, e_se, tpc_grp, metabolic_category)


lichen <- rbind(gpp_lichen, r_lichen) %>%
  left_join(ea_lichen, by = c("tpc_grp", "metabolic_category")) %>%
  pivot_wider(id_cols = coupling_unit,
    names_from = metabolic_category,
    values_from = c(e, e_se),
    names_sep = "_") %>%
  mutate(delta = e_r - e_gpp) %>%
  drop_na()

t.test(lichen$delta) #0.3301134

### coral 

gpp_coral <- read_csv("analysis/tidy data/coral_gpp_final.csv") %>% select(tpc_grp, latitude, species, metabolic_category) %>%
  separate(tpc_grp, into =c("rep", "study_id"), remove=FALSE) %>%
  mutate(coupling_unit = interaction(study_id, species, latitude)) %>% 
  distinct(coupling_unit, .keep_all = TRUE)
r_coral  <- read_csv("analysis/tidy data/coral_r_final.csv") %>% select(tpc_grp, latitude, species, metabolic_category) %>%
  separate(tpc_grp, into =c("rep", "study_id"), remove=FALSE) %>%
  mutate(coupling_unit = interaction(study_id, species, latitude)) %>% 
  distinct(coupling_unit, .keep_all = TRUE)

ea_coral  <- read_csv("analysis/tidy data/coral_tpc_Eas.csv") %>% filter(metabolic_category != "npp") %>%
  select(e, e_se, tpc_grp, metabolic_category)


coral <- rbind(gpp_coral, r_coral) %>%
  left_join(ea_coral, by=c("tpc_grp", "metabolic_category")) %>%
  select(-c(tpc_grp, rep)) %>%
  pivot_wider(id_cols = coupling_unit,
              names_from = metabolic_category,
              values_from = c(e, e_se),
              names_sep = "_") %>%
  mutate(delta = e_r - e_gpp) %>%
  drop_na()

hist(coral$delta, breaks=30)
  

t.test(coral$delta) #0.1379857 

### cell

gpp_cell <- read_csv("analysis/tidy data/cell_gpp_final.csv") %>% select(tpc_grp, latitude, species, metabolic_category) %>%
  separate(tpc_grp, into =c("rep", "study_id"), remove=FALSE) %>%
  mutate(coupling_unit = interaction(study_id, species, latitude)) %>% 
  distinct(coupling_unit, .keep_all = TRUE)
r_cell  <- read_csv("analysis/tidy data/cell_r_final.csv") %>% select(tpc_grp, latitude, species, metabolic_category) %>%
  separate(tpc_grp, into =c("rep", "study_id"), remove=FALSE) %>%
  mutate(coupling_unit = interaction(study_id, species, latitude)) %>% 
  distinct(coupling_unit, .keep_all = TRUE)
ea_cell  <- read_csv("analysis/tidy data/cell_tpc_Eas.csv") %>% filter(metabolic_category != "npp") %>%
  select(e, e_se, tpc_grp, metabolic_category)


cell <- rbind(gpp_cell, r_cell) %>%
  left_join(ea_cell, by=c("tpc_grp", "metabolic_category")) %>%
  select(-c(tpc_grp, rep)) %>%
  pivot_wider(id_cols = coupling_unit,
              names_from = metabolic_category,
              values_from = c(e, e_se),
              names_sep = "_") %>%
  mutate(delta = e_r - e_gpp) %>%
  drop_na()

t.test(cell$delta) #-0.5239927
