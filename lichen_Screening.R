install.packages("revtools")
library(revtools)
library(tidyverse)
pacman::p_load(viridis, revtools, nlme, lme4, MuMIn, patchwork)


#######screening abstracts
#lichen1 <- read_bibliography("lichen24sept21.bib")
#lichen2 <- read_bibliography("wos155.bib")
#screen_abstracts(lichen2)

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
         mmol_mg_minplus1 = abs_mmol_mg_min+1, 
         log=log(abs_mmol_mg_min))
       

subset_lichen<- lichen_data %>%
  filter(!response_def %in% c("gross photosynthesis")) %>%
  mutate(broad_responses = as.factor(case_when(ppfd %in% 
                    c(0) ~ "respiration",
TRUE ~ "photosynthesis")), 
elevation_broad = as.factor(elevation_broad),
lichen_type = as.factor(lichen_type),
latitude_broad = as.factor(latitude_broad)) %>%
  mutate(centre_temp = I(inv_T-mean(inv_T)), 
         ppfd_group = case_when(avg_ppfd %in% c(0) ~ "dark", 
                                avg_ppfd %in% c(1:150) ~ "low",
                                avg_ppfd %in% c(200:500) ~ "medium", 
                                TRUE ~ "high")) 

photosynthesis <- subset_lichen %>%
  filter(broad_responses == "photosynthesis")

respiration <- subset_lichen %>%
  filter(broad_responses == "respiration")




#subset_lichen$broad_responses<- relevel(subset_lichen$broad_responses, "respiration")
photosynthesis$latitude_broad<- relevel(photosynthesis$latitude_broad, "subarctic")
photosynthesis$elevation_broad<- relevel(photosynthesis$elevation_broad, "neutral")
photosynthesis$lichen_type<- relevel(photosynthesis$lichen_type, "green algae")
respiration$latitude_broad<- relevel(respiration$latitude_broad, "tropical")
respiration$elevation_broad<- relevel(respiration$elevation_broad, "neutral")
respiration$lichen_type<- relevel(respiration$lichen_type, "green algae")

#random effects for photosynthesis 
model00p<- lmer(log(abs_mmol_mg_min) ~abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|species/ study_id) , data=photosynthesis, REML=FALSE)
model0p <- lmer(log(abs_mmol_mg_min) ~abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=photosynthesis, REML=FALSE)
model1p <- lmer(log(abs_mmol_mg_min) ~abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|study_id) , data=photosynthesis, REML=FALSE)
model2p <- gls(log(abs_mmol_mg_min) ~abs(latitude) +lichen_type*centre_temp + elevation_broad, data=photosynthesis)
#random effects for respiration
model00r<- lmer(log(abs_mmol_mg_min) ~abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|species/ study_id) , data=respiration, REML=FALSE)
model0r <- lmer(log(abs_mmol_mg_min) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=respiration, REML=FALSE)
model1r <- lmer(log(abs_mmol_mg_min) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|study_id) , data=respiration, REML=FALSE)
model2r <- gls(log(abs_mmol_mg_min) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad, data=respiration)

#model selection -- random effects --> all the data togetehr (not relevant anymore)
#model00 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|centre_temp) , data=subset_lichen, REML=FALSE)
#model0 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=subset_lichen, REML=FALSE)
#model1 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|centre_temp) + (1 |study_id) , data=subset_lichen, REML=FALSE)
#model1 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1|study_id) , data=subset_lichen, REML=FALSE)
#model2 <- gls(log(abs_mmol_mg_min) ~broad_responses*centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad, data=subset_lichen)

#model selection for photosynthesis random effects
model.sel(model00p, model0p, model1p, model2p)

#model selection for respiration random effects
model.sel(model00r, model0r, model1r, model2r)

#model0 prevails

#model selection for photosynthesis
#missing elevation
model3p <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=photosynthesis)
#missing latitude
model4p <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=photosynthesis)
#missing both latitude and elevation
model5p <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=photosynthesis)
#missing interaction between temp and response
model6p <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type*centre_temp + abs(latitude)+ (1+centre_temp|study_id), data=photosynthesis)
#missing temp
model7p <- lmer(log(abs_mmol_mg_min) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=photosynthesis)
#missing lichen type
model10p <- lmer(log(abs_mmol_mg_min) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=photosynthesis)
#missing interaction between lichen type + centre temp
model11p <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=photosynthesis)
full_modelp <- lmer(log(abs_mmol_mg_min) ~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=photosynthesis)
model.sel(model3p, model4p, model5p, model6p, model7p, model10p, model11p, full_modelp)

#model7p  is best by < 2 AIC
#summary(model7p)
#confint(model7p)

#picking the model with the same structure as respiration
summary(model11p)
confint(model11p)


#model selection for respiration
#missing elevation
model3r <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=respiration)
#missing latitude
model4r <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=respiration)
#missing both latitude and elevation
model5r <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=respiration)
#missing temp
model7r <- lmer(log(abs_mmol_mg_min) ~abs(latitude) +lichen_type + elevation_broad + (1+centre_temp|study_id), data=respiration)
#missing lichen type
model10r <- lmer(log(abs_mmol_mg_min) ~centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=respiration)
#missing interaction between lichen type + centre temp
model11r <- lmer(log(abs_mmol_mg_min) ~centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=respiration)
full_modelr <- lmer(log(abs_mmol_mg_min) ~lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=respiration)
model.sel(model3r, model4r, model5r, model6r, model7r, model10r, model11r, full_modelr)


#model11r is best by < 2 AIC
summary(model11r)
confint(model11r)




#model1 prevails
#missing elevation
#model3 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp +lichen_type*centre_temp +abs(latitude)+ (1+centre_temp|study_id), data=subset_lichen)
#missing latitude
#model4 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id) , data=subset_lichen)
#missing both latitude and elevation
#model5 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp +lichen_type*centre_temp + (1+centre_temp|study_id), data=subset_lichen)
#missing interaction between temp and response
#model6 <- lmer(log(abs_mmol_mg_min) ~broad_responses + centre_temp +lichen_type*centre_temp + abs(latitude)+ (1+centre_temp|study_id), data=subset_lichen)
#missing temp
#model7 <- lmer(log(abs_mmol_mg_min) ~broad_responses + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id), data=subset_lichen)
#missing response
#model8 <- lmer(log(abs_mmol_mg_min) ~centre_temp + abs(latitude) +lichen_type*centre_temp + elevation_broad + (1+centre_temp|study_id), data=subset_lichen)
#missing temp and response
#model9 <- lmer(log(abs_mmol_mg_min) ~abs(latitude) + elevation_broad +lichen_type*centre_temp + (1+centre_temp|study_id), data=subset_lichen)
#missing lichen type
#model10 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=subset_lichen)
#missing interaction between lichen type + centre temp
#model11 <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp +lichen_type +abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=subset_lichen)
#full_model <- lmer(log(abs_mmol_mg_min) ~broad_responses*centre_temp +lichen_type*centre_temp +  abs(latitude) + elevation_broad + (1+centre_temp|study_id), data=subset_lichen)
#model.sel(model3, model4, model5, model6, model7, model8, model9, model10, model11, full_model)

#model11 is best by < 2 AIC
#summary(model11)
#confint(model11)




photosynthesis_df<- count(photosynthesis, study_id, lichen_type, elevation_broad, latitude, centre_temp)


randslope_p<- random.effects(model11p) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

photosynthesis_joined<- left_join(randslope_p,photosynthesis_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "green algae" & elevation_broad == "neutral" ~-9.05784,
                                lichen_type == "green algae" & elevation_broad == "high" ~-9.05784+-9.07998, 
                                lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~-9.05784+0.18316,
                                lichen_type == "cyanobacteria" & elevation_broad == "high" ~-9.05784+0.18316-9.07998),
         latitudeii=abs(latitude)*-0.06160, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope+0.07049) %>%
  as.data.frame()

p_plot<- ggplot(data=photosynthesis, aes(x=centre_temp, y=log(abs_mmol_mg_min)), colour="olivedrab4")+
  geom_point(colour="olivedrab4", alpha=0.5)+
  geom_abline(data=photosynthesis_joined, aes(slope=f_slope, intercept=f_intercept), colour="olivedrab4")+
  geom_abline(aes(slope=0.07049, intercept=-9.05784), colour="black")+
  geom_hline(yintercept=-9.1055, size=1.5)+  
  xlab("Temperature (1/kT)")+
  ylab("Metabolic rate (mmol CO2 per mg per min)")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("Photosynthesis")+
  ylim(-31, -4)

respiration_df<- count(respiration, study_id, lichen_type, elevation_broad, latitude, centre_temp)


randslope_r<- random.effects(model11r) %>%
  as.data.frame() %>%
  select(-condsd) %>%
  pivot_wider(names_from = term, values_from = condval) %>%
  rename(r_intercept = `(Intercept)`, 
         study_id= grp, 
         r_slope=centre_temp)

respiration_joined<- left_join(randslope_r,respiration_df, by="study_id") %>%
  mutate(c_intercept=case_when(lichen_type == "cyanobacteria" & elevation_broad == "high" ~-18.771838,
                               lichen_type == "cyanobacteria" & elevation_broad == "neutral" ~-18.771838+0.122152, 
                               lichen_type == "green algae" & elevation_broad == "high" ~-18.771838+-0.692983,
                               lichen_type == "green algae" & elevation_broad == "neutral" ~-18.771838+-0.692983+0.122152,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "high" ~-18.771838+-0.333793,
                               lichen_type == "cyanobacteria + green algae" & elevation_broad == "neutral" ~-18.771838+-0.333793+0.122152),
         latitudeii=abs(latitude)*-0.008383, 
         f_intercept=r_intercept+c_intercept+latitudeii, 
         f_slope=r_slope+-0.602350) %>%
  as.data.frame()


r_plot <- ggplot(data=respiration, aes(x=centre_temp, y=log(abs_mmol_mg_min)), colour="grey")+
  geom_point(colour="grey", alpha=0.5)+
  geom_abline(data=respiration_joined, aes(slope=f_slope, intercept=f_intercept), colour="grey")+
  geom_abline(slope=-0.602350, intercept = -18.771838, size=1.5)+  
  xlab("Temperature (1/kT)")+
  ylab("")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("Respiration")+
  ylim(-31, -4)
  #scale_x_continuous(sec.axis = sec_axis(~.*(1/(8.61*10^-5))-273.15), name = "Temperature (°C)")
  #scale_x_continuous(sec.axis = sec_axis(~ . * 20, name = "Temperature (°C)"))

#how do i plot this if this is the centred term?



p_plot + r_plot




