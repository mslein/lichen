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
lichen_metabolic <- read_csv("lichen temp metabolism.csv") %>%
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
