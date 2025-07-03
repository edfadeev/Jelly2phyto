require(ggplot2)
require(dplyr)

#configuration of plots
#theme
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
colour_pallette<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", "#DD7788")

names(colour_pallette)<- c("Alphaproteobacteria_uncl",  "Puniceispirillaceae", "Rhodobacteraceae", "Sphingomonadaceae",
                           "Bacteroidia_uncl", "Crocinitomicaceae", "Flavobacteriaceae", "Salibacteraceae","Schleiferiaceae", 
                           "Alteromonadaceae", "Gammaproteobacteria_uncl", "Halieaceae", "Litoricolaceae", "Nitrincolaceae", 
                           "Pseudoalteromonadaceae",  "Spongiibacteraceae", "Vibrionaceae", "Balneolaceae", "Other taxa")
##############################################
#import enrichment results from anvio
##############################################

#taxonomy of bins 
BINS_summary <- read.csv("data/combined_bins_summary.txt",sep = "\t") %>% 
  mutate(group= factor(group, levels =c ("T0","T3J","T4J",
                                         "T3C","T4C", "SJ","SC",
                                         "M","ML")))
#KEGG modules in each bin
BINS_modules <- read.csv("data/BINS_modules.txt", sep = "\t") 

#import enrichment test results
J_C_enr_T4_mod <- read.csv("data/Jelly_Control_T4_enrichment.txt", sep = "\t") %>% mutate(Time="T4") %>% 
  left_join(read.csv("data/modules_info.txt", sep = "\t"), by = c("accession" = "module")) #module definitions

##############################################
## Enriched KEGG pathways Jelly-OM vs. Control
##############################################
#combined taxa and metabolism
J_C_enr_T4_bins<- J_C_enr_T4_mod %>% tidyr::separate_longer_delim(sample_ids, ",") %>% 
  filter(adjusted_q_value< 0.1) %>% #FDR threshold
  left_join(BINS_summary, by = c("sample_ids"="db_name")) %>% 
  left_join(BINS_modules, by = c("sample_ids"="db_name", "accession"="module") )


#calculate how many enriched pathways were present in each bin
BINS_enr_modules_sum<- J_C_enr_T4_bins %>% 
  mutate(t_family=case_when(grepl("\\d", t_family) ~paste0(t_class, "_uncl"),
                            t_order!="" & t_family==""~ paste0(t_class, "_uncl"),
                            t_class!="" & t_order=="" & t_family==""~ paste0(t_class, "_uncl"),
                            t_class=="" & t_order=="" & t_family==""~ paste0(t_domain, "_uncl"),
                            t_genus=="Pseudoalteromonas"~ "Pseudoalteromonadaceae",
                            TRUE ~ t_family)) %>% 
  group_by(associated_groups, t_domain, t_class, t_order, 
           t_family, category, subcategory) %>% 
  summarize(Num_pathways=n())
  
#plot taxonomy of bins
BINS_summary %>% filter(t_domain!="", grepl("T4",group)) %>% 
  mutate(t_family=case_when(grepl("\\d", t_family) ~paste0(t_class, "_uncl"),
                            t_order!="" & t_family==""~ paste0(t_class, "_uncl"),
                            t_class!="" & t_order=="" & t_family==""~ paste0(t_class, "_uncl"),
                            t_class=="" & t_order=="" & t_family==""~ paste0(t_domain, "_uncl"),
                            t_genus=="Pseudoalteromonas"~ "Pseudoalteromonadaceae",
                            TRUE ~ t_family)) %>% 
  mutate(t_family=ifelse(t_family %in% tax_fam, t_family, "Other taxa")) %>% 
  group_by(group, t_class, t_family) %>% 
  summarize(n_bins=n()) %>% 
  filter(!t_family=="Bacteria_uncl") %>% 
  mutate(t_family=factor(t_family, levels=c("Alphaproteobacteria_uncl",  "Puniceispirillaceae", "Rhodobacteraceae", "Sphingomonadaceae",
                                            "Bacteroidia_uncl", "Crocinitomicaceae", "Flavobacteriaceae", "Salibacteraceae","Schleiferiaceae", 
                                            "Alteromonadaceae", "Gammaproteobacteria_uncl", "Halieaceae", "Litoricolaceae", "Nitrincolaceae", 
                                            "Pseudoalteromonadaceae",  "Spongiibacteraceae", "Vibrionaceae", "Balneolaceae", "Other taxa"))) %>% 
  ggplot(aes(x=group, y=n_bins, fill = t_family))+
  geom_col()+
  scale_fill_manual(values = colour_pallette)+
  guides(fill=guide_legend(title = "Class"))+
  theme_EF+
  theme(legend.position = "right")

#save plot
ggsave("Figures/Fig_5A-Bins_taxa.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180, height = 180, 
       scale = 2,
       dpi = 300)


#plot enriched KEGG pathways
BINS_enr_modules_sum %>% 
  mutate(t_family=ifelse(t_family %in% tax_fam, t_family, "Other taxa")) %>% 
  mutate(t_family=factor(t_family, levels=c("Alphaproteobacteria_uncl",  "Puniceispirillaceae", "Rhodobacteraceae", "Sphingomonadaceae",
                                            "Bacteroidia_uncl", "Crocinitomicaceae", "Flavobacteriaceae", "Salibacteraceae","Schleiferiaceae", 
                                            "Alteromonadaceae", "Gammaproteobacteria_uncl", "Halieaceae", "Litoricolaceae", "Nitrincolaceae", 
                                            "Pseudoalteromonadaceae",  "Spongiibacteraceae", "Vibrionaceae", "Balneolaceae", "Other taxa"))) %>% 
  ggplot(aes(x=category, y=Num_pathways, fill= t_family))+
  geom_col(position = "stack")+
  scale_fill_manual(values = colour_pallette)+
  #guides(fill=guide_legend(title = paste0(x), ncol =1))+
  facet_wrap(.~associated_groups)+
  theme_EF+
  theme(legend.position = "none")+
  coord_flip()

#save plot
ggsave("Figures/Fig_5B-met_enr.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180, height = 180, 
       scale = 2,
       dpi = 300)