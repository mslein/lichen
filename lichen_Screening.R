install.packages("revtools")
library(revtools)
library(tidyverse)
pacman::p_load(ggpmisc, viridis)

lichen1 <- read_bibliography("lichen24sept21.bib")
lichen2 <- read_bibliography("wos155.bib")
screen_abstracts(lichen2)

##ask mary about the difference in units between the spreadsheet 
## and the plos one figure 1....
library(tidyverse)
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
         abs_mmol_mg_min = abs(response_value*(mmol_co2*(mg_sample^-1)*(per_min^-1))), 
         mmol_mg_minplus1 = abs_mmol_mg_min+1)
       

subset_lichen<- lichen_data %>%
  filter(!response_def %in% c("gross photosynthesis")) %>%
  mutate(broad_responses = as.factor(case_when(ppfd %in% 
                    c(0) ~ "respiration",
TRUE ~ "photosynthesis")))

#looking at rates across ppfd values 
subset_lichen %>%
ggplot(aes(x=desc(inv_T), y=log(abs_mmol_mg_min), color=as.numeric(avg_ppfd)))+
  geom_point(alpha=0.5)+
  facet_wrap(~broad_responses)+
  scale_color_viridis(option="D")
#trying my hand at random slope + random intercepts?
library(lme4)
model <- lmer(log(abs_mmol_mg_min) ~ broad_responses * inv_T + (1 + inv_T|study_id), data=subset_lichen)
subset_lichen$fit <- predict(model) 
ggplot(subset_lichen,aes(inv_T, log(abs_mmol_mg_min), 
                         group=interaction(study_id, broad_responses), color = study_id, shape=broad_responses )) + 
  facet_grid(~broad_responses) +
  geom_line(aes(y=fit, lty=broad_responses), size=0.8) +
  geom_point(alpha = 0.3) + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()+
  scale_color_viridis(option= "D", discrete=TRUE)
  
#mary's model
library(nlme)
mod1 <-lme(log(abs_mmol_mg_min) ~ 1 + broad_responses*I(inv_T - mean(inv_T)), 
             random = ~ 1 | study_id, data=subset_lichen, method="REML", na.action=na.omit) 
mod <- mod1
slope_photo <- fixef(mod)[3]
slope_resp <- fixef(mod)[3] + fixef(mod)[4]
int_photo <- fixef(mod)[1] - fixef(mod)[3]*mean(subset_lichen$inv_T)
int_resp <- fixef(mod)[1] + fixef(mod)[2] - fixef(mod)[3]*mean(subset_lichen$inv_T) - fixef(mod)[4]*mean(subset_lichen$inv_T)


funcphoto <- function(x) {int_photo + slope_photo*x} # for trophic level 1
valsphoto <- funcP1(subset_lichen[(subset_lichen$broad_responses=="photosynthesis"),]$inv_T)

funcresp <- function(x) {int_resp + slope_resp*x } # for trophic level 2
valsresp <- funcP2(subset_lichen[(subset_lichen$broad_responses=="respiration"),]$inv_T)

k <- 8.617342*10^-5  # eV/K
  
Fig2A <- ggplot(data = subset_lichen, aes(x = -inv_T, y = log(abs_mmol_mg_min), color=avg_ppfd)) + #, xmin = -40.2, xmax = -38.2
  theme_bw() +
  geom_point()+
  scale_color_viridis()+
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
  geom_line(data = subset(subset_lichen, broad_responses == "photosynthesis"), aes(x = -inv_T, y = as.vector(valsphoto)), color = 'black', size = 1)+
  geom_line(data = subset(subset_lichen, broad_responses == "respiration"), aes(x = -inv_T, y = as.vector(valsresp)), color = 'black', linetype = 3, size = 1) 
  #geom_line(data = subset(data, trophic.level == "PZN"), aes(x = -invTT, y = PBvalsPZN), color = 'black', linetype = 2, size = 1) +
  #scale_linetype_identity() +
  #ylab("ln(ug Chl a / L)") +
  #labs(title = "Phytoplankton Biomass [Chla]", size = 3) +
  #annotate("text", x = c(-39.5,-39.5), y = c(1, 0.7), label = c("solid = respiration", "dotted = photosynthesis"), size = 4) +
  #scale_x_continuous("Temperature (1/kTi)", sec.axis = sec_axis(~(-(1/(k*+.))-273), name = xlab))




