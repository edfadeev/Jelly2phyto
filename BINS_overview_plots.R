require(dplyr)
require(ggplot2)
require(gridExtra)

#########################
##    Preparations     ##
#########################
#define work directory
wd <- "D:/UniVie/Projects/Microcosm_Piran/Analysis/MetaG/08_BINS"


#taxonomy of bins 
BINS_summary <- read.csv(paste0(wd,"/combined_bins_summary.txt"),
                         sep = "\t") %>% 
                mutate(group= factor(group,
                                        levels =c ("T0","T3J","T4J",
                                                   "T3C","T4C", "SJ","SC",
                                                   "M","ML")))
#ggplot theme and colours
theme_EF <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        title = element_text(size=25, face ="bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 22, face ="bold"),
        axis.title = element_text(size = 20))

#colours
tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", "#DD7788")



##################################
## Plot overview on Class level ##
##################################
BINS_summary %>% 
  group_by(group, t_class, t_family) %>% 
  summarize(n_bins=n()) %>% 
  ggplot(aes(x=group, y=n_bins, fill = t_class))+
  geom_col()+
  scale_fill_manual(values = tol21rainbow)+
  guides(fill=guide_legend(title = "Class"))+
  theme_EF
#save
ggsave(paste0("Figures/BINS_summary.png"),
       plot = last_plot(),
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)


###############################
## Plot bins in main Classes ##
###############################
#summarize bins in main classes
tax_by_class<- lapply(c("Alphaproteobacteria", "Bacteroidia","Gammaproteobacteria"),
       function(x){
                  return(BINS_summary %>% 
                          filter(t_class==x) %>% 
                          group_by(group, t_genus) %>% 
                          summarize(n_bins=n()) %>% 
                          ggplot(aes(x=group, y=n_bins, fill = t_genus))+
                          geom_col()+
                          scale_fill_manual(values = c(tol21rainbow, "gray50","yellow","red","blue"))+
                          guides(fill=guide_legend(title = paste0(x), ncol =1))+
                          theme_EF)
                          })

#plot
grid.arrange(tax_by_class[[1]],tax_by_class[[2]],
                        tax_by_class[[3]],nrow = 1)

#save
tax_by_class.p <- arrangeGrob(tax_by_class[[1]],tax_by_class[[2]],
                              tax_by_class[[3]],nrow = 1)
ggsave(paste0("Figures/BINS_summary_Class.png"),
       plot = tax_by_class.p,
       units = "cm",
       width = 30, height = 10, 
       scale = 2,
       dpi = 300)
