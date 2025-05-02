



#lichen
ggplot()+
  geom_abline(slope=0.65, intercept=0, colour="grey30", linetype="dotted", size=2, alpha=0.5)+
  geom_abline(slope=0.32, intercept=0, colour="#abbb80", linetype="dotted", size=2, alpha=0.5)+
  geom_abline(slope=0.47, intercept=0, colour="grey30", size=2)+
  geom_abline(slope=0.3824633, intercept=0, colour="#abbb80", size=2)+
  xlim(0,5)+
  ylim(-1,5)+
  theme_classic()


#coral
ggplot()+
  geom_abline(slope=0.65, intercept=0, colour="grey30", linetype="dotted", size=2, alpha=0.5)+
  geom_abline(slope=0.32, intercept=0, colour="#abbb80", linetype="dotted", size=2, alpha=0.5)+
  geom_abline(slope=0.66, intercept=0, colour="grey30", size=2)+
  geom_abline(slope=0.617, intercept=0, colour="#abbb80", size=2)+
  xlim(0,5)+
  ylim(-1,5)+
  theme_classic()
