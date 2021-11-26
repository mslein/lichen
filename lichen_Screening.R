install.packages("revtools")
library(revtools)
library(tidyverse)
pacman::p_load(ggpmisc)

lichen1 <- read_bibliography("lichen24sept21.bib")
lichen2 <- read_bibliography("wos155.bib")
screen_abstracts(lichen2)

##ask mary about the difference in units between the spreadsheet 
## and the plos one figure 1....
lichen_data <- read_csv("extraction/lichen_data.csv") %>%
  mutate(inv_T=1/((8.617333262145*10^-5)*(temp+273.15)),
         mmol_co2 = case_when(gas_unit %in% c("mg") ~ (1/44.01), 
                                  gas_unit %in% c("nmol") ~ (0.000001),
                                  gas_unit %in% c("umol") ~ (0.0001)),
         mg_sample = case_when(mass %in% c("g") ~ (1000), 
                               mass %in% c("kg") ~ (1000000),
                               mass %in% c("mg") ~ (1)),
         per_min = case_when(time %in% c("hr") ~ (1/60),
                             time %in% c("sec") ~ (60), 
                             TRUE ~ (1)),
         mmol_mg_min = response_value*(mmol_co2*(mg_sample^-1)*(per_min^-1)), 
         mmol_mg_minplus1 = mmol_mg_min + 1)
       

subset_lichen<- lichen_data %>%
  filter(response_def %in% c("dark respiration", 
                             "respiration", 
                             "net photosynthesis",
                             "net photosynthesis rate")) %>%
  mutate(broad_responses = case_when(response_def %in% 
                    c("dark respiration","respiration") ~ "respiration",
TRUE ~ "photosynthesis"))

my.formula <- y ~ x
subset_lichen %>%
ggplot(aes(x=desc(inv_T), y=log(mmol_mg_min), color=broad_responses))+
  geom_smooth(method="lm")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  geom_point(alpha=0.5)
  #facet_wrap(~avg_ppfd)
ggplot(subset_lichen, aes(x=avg_ppfd, fill=broad_responses))+
  geom_histogram(binwidth=20)


library(nlme)
modPBF <-lme(log(mmol_mg_minplus1) ~ 1 + broad_responses*I(inv_T), 
             random = ~ 1 | study_id, data=subset_lichen, method="ML", na.action=na.omit) 
modPB <- modPBF

IPB1 <- fixef(modPB)[1] - fixef(modPB)[1]*mean(subset_lichen$inv_T)
SlPB1 <- fixef(modPB)[2]

IPB2 <- fixef(modPB)[1] + fixef(modPB)[2] - fixef(modPB)[1]*mean(subset_lichen$inv_T) - fixef(modPB)[1]*mean(subset_lichen$inv_T)
SlPB2 <- fixef(modPB)[2] + fixef(modPB)[4]

PB.funcP <- function(x) {IPB1 + SlPB1*x} # for trophic level 1
PBvalsP <- PB.funcP(subset_lichen[(subset_lichen$broad_responses=="respiration"),]$inv_T)

PB.funcPZ <- function(x) { IPB2 + SlPB2*x } # for trophic level 2
PBvalsPZ <- PB.funcPZ(subset_lichen[(subset_lichen$broad_responses=="photosynthesis"),]$inv_T)




k <- 8.617342*10^-5  # eV/K
  
Fig2A <- ggplot(data = subset_lichen, aes(x = -inv_T, y = log(mmol_mg_minplus1), ymin = 0, ymax = 3)) + #, xmin = -40.2, xmax = -38.2
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5), 
        legend.key = element_rect(fill = NA), 
        axis.text = element_text(size = 10),
        #strip.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        plot.title = element_text(hjust = -1)) +
  ggtitle("A. Phytoplankton Biomass [Chla]") +
  geom_line(data = subset(subset_lichen, broad_responses == "respiration"), aes(x = -inv_T, y = as.vector(PBvalsP)), color = 'black', size = 1) +
  geom_line(data = subset(subset_lichen, broad_responses == "photosynthesis"), aes(x = -inv_T, y = PBvalsPZ), color = 'black', linetype = 3, size = 1) 
  #geom_line(data = subset(data, trophic.level == "PZN"), aes(x = -invTT, y = PBvalsPZN), color = 'black', linetype = 2, size = 1) +
  #scale_linetype_identity() +
  #ylab("ln(ug Chl a / L)") +
  #labs(title = "Phytoplankton Biomass [Chla]", size = 3) +
  #annotate("text", x = c(-39.5,-39.5), y = c(1, 0.7), label = c("solid = respiration", "dotted = photosynthesis"), size = 4) +
  #scale_x_continuous("Temperature (1/kTi)", sec.axis = sec_axis(~(-(1/(k*+.))-273), name = xlab))







