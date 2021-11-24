install.packages("revtools")
library(revtools)
library(tidyverse)
pacman::p_load(ggpmisc)

#
lichen1 <- read_bibliography("lichen24sept21.bib")
lichen2 <- read_bibliography("wos155.bib")
screen_abstracts(lichen2)

##ask mary about the difference in units between the spreadsheet 
## and the plos one figure 1....
lichen_data <- read_csv("extraction/lichen_data.csv") %>%
  mutate(temp_dep=1/((8.617333262145*10^-5)*(temp+273.15)),
         mmol_co2 = case_when(gas_unit %in% c("mg") ~ (1/44.01), 
                                  gas_unit %in% c("nmol") ~ (0.000001),
                                  gas_unit %in% c("umol") ~ (0.0001)),
         mg_sample = case_when(mass %in% c("g") ~ (1000), 
                               mass %in% c("kg") ~ (1000000),
                               mass %in% c("mg") ~ (1)),
         per_min = case_when(time %in% c("hr") ~ (1/60),
                             time %in% c("sec") ~ (60), 
                             TRUE ~ (1)),
         std_response = response_value*(mmol_co2*mg_sample*per_min))
       

subset_lichen<- lichen_data %>%
  filter(response_def %in% c("dark respiration", 
                             "respiration", 
                             "net photosynthesis",
                             "net photosynthesis rate")) %>%
  mutate(broad_responses = case_when(response_def %in% 
                    c("dark respiration","respiration") ~ "respiration",
TRUE ~ "photosynthesis"))
  
                                
ggplot(subset_lichen, aes(x=temp_dep, y=std_response, color=broad_responses))+
  geom_smooth(method="lm")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  geom_point()
  #facet_wrap(~broad_responses)
