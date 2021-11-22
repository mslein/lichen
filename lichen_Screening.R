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
lichen_data <- read_csv("lichen_data.csv") 
count(lichen_data, study_id, response_def, response_units)
count(lichen_data, response_units)
lichen_calcs <- lichen_data %>%
  mutate(temp_dep=1/((8.617333262145*10^-5)*(temp+273.15)),
         metric= case_when(response_units %in% c("nmol O2 *g of dried weight ^1 *min^-1") ~ "vo2", 
                           response_units %in% c("ratio") ~ "other", 
                           TRUE ~ "vco2")) %>%
  filter(!metric %in% c("vo2"),
         !response_def %in% c("max photosynthesis"))
test_1<- lichen_calcs %>%
  filter(response_def %in% c("dark respiration", "net photosynthesis"),
         !response_units %in% c("umol CO2 * m^-2 *s^-1")) %>%
  mutate(std_units = response_value/1000)
count(test_1, response_def, response_units)

ggplot(lichen_calcs,aes(x=temp_dep, y=response_value, group=response_def, color=response_def))+
  geom_smooth(method="lm")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  geom_point()



  mutate(value_converted_o2=value*1375) %>%
  mutate(temp_dep=1/((8.617333262145*10^-5)*(temp+273.15)))
view(lichen1)
screen_abstracts(lichen1)

my.formula <- y ~ x

ggplot(lichen_metabolic, aes(x=temp_dep, y=value_converted_o2, colour=response))+
  geom_smooth(method="lm")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  geom_point()
  #facet_wrap(~response)
