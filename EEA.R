require(dplyr)
require(tidyr)
require(ggplot2)

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

################################################
#import enzymatic activity results
################################################
EEA_results<- read.table("data/EEA.txt", dec=",", fill=TRUE, sep="\t", header = TRUE) %>% 
  mutate(Treatment  = stringr::str_replace_all(Treatment, "C|J", function(x) ifelse(x == "C", "J", "C")))


################################################
#Plot Fig. 3 - enzymatic activity 
################################################
EEA_results %>% 
  group_by(Experiment, Treatment, Time.hours., Enzyme) %>% 
  summarise_at(vars("Measurement"), list(mean = mean, se = se)) %>% 
  mutate(Time.hours.= factor(Time.hours., levels=c("0","9","24","43","67","72","96","120","144")),
         Enzyme=factor(Enzyme, levels=c("AMA","APA","AGAL","BGAL","BGLU","CHIT","OLE"))) %>% 
  filter(!Enzyme %in% c("AGAL","BGAL","BGLU")) %>% 
  ggplot(aes(x=Time.hours., y=mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean +se), width = 0.2, colour = "gray50")+
  geom_point(size =5, colour = "black")+geom_point(size =4)+
  geom_line()+
  facet_grid(cols=vars(Experiment),rows = vars(Enzyme), scales="free",space="free_x",switch="x")+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Fig_3-EEA.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 180, 
       scale = 2,
       dpi = 300)


################################################
#Plot Fig. S1 - enzymatic activity 
################################################
EEA_results %>% 
  group_by(Experiment, Treatment, Time.hours., Enzyme) %>% 
  summarise_at(vars("Measurement"), list(mean = mean, se = se)) %>% 
  mutate(Time.hours.= factor(Time.hours., levels=c("0","9","24","43","67","72","96","120","144")),
         Enzyme=factor(Enzyme, levels=c("AMA","APA","AGAL","BGAL","BGLU","CHIT","OLE"))) %>% 
  filter(Enzyme %in% c("AGAL","BGAL","BGLU")) %>% 
  ggplot(aes(x=Time.hours., y=mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean +se), width = 0.2, colour = "gray50")+
  geom_point(size =5, colour = "black")+geom_point(size =4)+
  geom_line()+
  facet_grid(cols=vars(Experiment),rows = vars(Enzyme), scales="free",space="free_x",switch="x")+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Fig_S1-EEA.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 180, 
       scale = 2,
       dpi = 300)