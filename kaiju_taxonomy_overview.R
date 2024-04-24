require(dplyr)
require(ggplot2)


lapply(c("Phylum", "Class", "Order","Family", "Genus"), function(tax) {
kaiju_table<- read.csv(paste0("D:/UniVie/Projects/Microcosm_Piran/Analysis/MetaG/Taxonomy/Tax_Kaiju/kaiju_nr_euk_merged_",tax,".out"),
                       sep = "\t") %>% 
  mutate(Sample= gsub("/lisc/scratch/oceanography/efadeev/Microcosm_Piran/MicPir_sep/07_TAXONOMY/|_nr_euk_tax.out","",file)) %>% 
  mutate(Sample= gsub("[A-H]0[0-9]_","",Sample),
         Domain = case_when(grepl("Eukaryota", taxon_name)~"Eukaryota",
                            grepl("Bacteria", taxon_name)~"Bacteria",
                            grepl("Archaea", taxon_name)~"Archaea",
                            grepl("Viruses", taxon_name)~"Viruses")) %>%
  filter(!grepl("merged", Sample)) %>% 
  group_by(Sample, Domain, taxon_name) %>% 
  summarize(Total = sum(percent)) %>% 
  mutate(Sample=case_when(Sample=="T0_S" ~ "S_0", TRUE~ Sample),
          Sample=factor(Sample, levels=c("M01","M02","M03","ML1","ML2","ML3","T0","T3C1","T3C2","T4C1","T4C2","T4C3","T3C3","T3J1",
                                        "T3J2","T3J3","T4J1","T4J2","T4J3","S_0","SC1_2","SC2_2","SC3_2","SJ1_2","SJ2_2","SJ3_2")),
           Group = case_when(grepl("M",Sample) ~ "Environment",
                             grepl("T0|T3|T4",Sample) ~ "Incubation",
                             grepl("S",Sample) ~ "Succession")) %>% 
  tidyr::separate(taxon_name, into = c("Domain", "Phylum","Class","Order", "Family","Genus"), extra = 'drop', remove = FALSE)
  
  
kaiju_table%>% 
  filter(Total> 1, Domain == "Bacteria") %>% 
  ggplot(aes(x= Sample, y= Total, fill = get(tax)))+
  geom_col()+
  facet_wrap(~Group, scales = "free_x")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Percent of classified reads (>1%)")+
  labs(fill = paste(tax))+
  theme_EF+
  theme(legend.position = "right",
        axis.text.x = element_text(angle=90))

ggsave(paste0("Figures/Kaiju_Bac_",tax,".png"),
       plot = last_plot(),
       units = "mm",
       width = 250, height = 250, 
       scale = 2,
       dpi = 300)

kaiju_table%>% 
  filter(Total> 0.05, Domain == "Eukaryota") %>% 
  ggplot(aes(x= Sample, y= Total, fill = get(tax)))+
  geom_col()+
  facet_wrap(~Group, scales = "free_x")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Percent of classified reads (>0.05%)")+
  labs(fill = paste(tax))+
  theme_EF+
  theme(legend.position = "right",
        axis.text.x = element_text(angle=90))

ggsave(paste0("Figures/Kaiju_Euk_",tax,".png"),
       plot = last_plot(),
       units = "mm",
       width = 250, height = 250, 
       scale = 2,
       dpi = 300)

})
