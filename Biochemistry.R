require(dplyr)
require(tidyr)
require(ggplot2)

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

################################################
#import raw table
################################################
#biochemical parameters and cell counts
mic_biochem_results<- read.table("data/Mic_biochem_results.txt", dec=",", fill=TRUE, sep="\t", header = TRUE)

################################################
#summarize the biological replicates for each data point
################################################
mic_biochem_summary<- mic_biochem_results %>% filter(!Time.hours.=="48") %>%  #remove timpoint with no measurments
  reshape2::melt(id.vars=c("Experiment", "Treatment", "Bottle", "Time.hours.")) %>% 
  group_by(Experiment, Treatment, Time.hours., variable) %>% 
  summarise_at(vars("value"),funs(mean= mean(., na.rm=TRUE), se = se(., na.rm=TRUE)))

################################################
#plot Fig. 2 - Biochemical data
################################################
mic_biochem_summary%>% 
  filter(variable %in% c("DOC.umol.L.1.","NH4..umol.L.1.", "PO43..umol.L.1.",
                         "DON..umol.L.1.","Bacterial.production.ug.C.L.1.h.1", 
                         "PP..ug.C.L.1.h.1.")) %>% 
  mutate(Time.hours.= factor(Time.hours., levels=c("0","9","24","43","67","72","96","120","144"))) %>% 
  ggplot(aes(x=Time.hours., y=mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean +se), width = 0.2, colour = "gray50")+
  geom_point(size =7, colour = "black")+geom_point(size =5)+
  geom_line()+
  facet_grid(cols=vars(Experiment),rows = vars(variable), scales="free",space="free_x",switch="x")+
  theme_EF+
  theme(legend.position = "none")

#save the plot
ggsave("./Figures/Fig_2.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 180, 
       scale = 2,
       dpi = 300)

################################################
#plot Fig. 6 - Cell abundances
################################################
mic_biochem_results %>% filter(!Time.hours.=="48", Experiment =="B") %>%  #remove timpoint with no measurments
  reshape2::melt(id.vars=c("Experiment", "Treatment", "Bottle", "Time.hours.")) %>% 
  filter(grepl("FCM|euk", variable, ignore.case=TRUE)) %>% 
  mutate(variable=case_when(variable=="Abundance.FCM.cells.mL.1" ~ "Microbial cells",
                            variable=="Abundance.Syn.FCM.cells.mL.1" ~ "Synechococcus",
                            variable=="Abundance.PicoEukaryotes.cells.ml.1" ~ "Picoeukaryotes",
                            variable=="Abundance.phototrophic.nanoeukaryotes.cells.ml.1" ~ "Phot. nanoeukaryotes")) %>% 
  mutate(variable=factor(variable, c("Microbial cells","Synechococcus","Picoeukaryotes","Phot. nanoeukaryotes"))) %>% 
  group_by(Experiment, Treatment, Time.hours., variable) %>% 
  summarise_at(vars("value"), list(mean = mean, se = se)) %>% 
  mutate(Time.hours.= factor(Time.hours., levels=c("0","9","24","43","67","72","96","120","144"))) %>% 
  ggplot(aes(x=Time.hours., y=mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean +se), width = 0.2, colour = "gray50")+
  geom_point(size =7, colour = "black")+geom_point(size =5)+
  geom_line()+
  ylab("Abundance (cells mL-1)")+
  facet_wrap(.~variable, scales = "free_y", ncol=2)+
  #facet_grid(cols=vars(Experiment),rows = vars(variable), scales="free",space="free_x",switch="x")+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Fig_6-cell_counts.pdf",
       plot = last_plot(),
       units = "mm",
       #width = 90,
       #height = 180, 
       scale = 2,
       dpi = 300)

################################################
#Generate Tab. S1 - biochemical data
################################################
mic_biochem_summary %>% 
  mutate(value=paste0(round(mean,2), "\u00B1",round(se,2)),
         Treatment = case_when(Treatment=="C"~"Control", 
                               Treatment=="J"~"Jelly-OM"),
         Experiment = case_when(Experiment =="A"~"Bac. deg.", 
                                Experiment=="B"~"Phyt. res.")) %>% 
  select(-c(starts_with("Abundance"),mean,se)) %>%
  tidyr::spread(variable, value) %>% 
  ungroup() %>% 
  openxlsx::write.xlsx(.,  file= 'Tables/Tab_S1-Mic_biochem_summary.xlsx')
