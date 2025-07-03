################################################
## This script visualizes the taxonomic profiling
## of the fastq libraries based on NCBI nr database
################################################

require(dplyr)
require(ggplot2)

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
tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", "#DD7788")

#####################################
#import Kaiju output table
#####################################
kaiju_table<- read.csv("data/kaiju_nr_euk_merged_Order.out", sep = "\t") %>% 
  mutate(Sample= gsub("/lisc/scratch/oceanography/efadeev/Microcosm_Piran/MicPir_sep/07_TAXONOMY/|_nr_euk_tax.out","",file)) %>% 
  mutate(Sample= gsub("[A-H]0[0-9]_","",Sample),
         #add Taxonomic domain
         Domain = case_when(grepl("Eukaryota", taxon_name)~"Eukaryota",
                            grepl("Bacteria", taxon_name)~"Bacteria",
                            grepl("Archaea", taxon_name)~"Archaea",
                            grepl("Viruses", taxon_name)~"Viruses")) %>%
  filter(!grepl("merged", Sample)) %>% 
  group_by(Sample, Domain, taxon_name) %>% 
  #summarize percentage
  summarize(Reads = sum(reads)) %>% 
  mutate(Sample=case_when(Sample=="T0_S" ~ "S_0", TRUE~ Sample),
         Sample=factor(Sample, levels=c("M01","M02","M03","ML1","ML2","ML3","T0","T3C1","T3C2","T4C1","T4C2","T4C3","T3C3","T3J1",
                                        "T3J2","T3J3","T4J1","T4J2","T4J3","S_0","SC1_2","SC2_2","SC3_2","SJ1_2","SJ2_2","SJ3_2")),
         Type = case_when(grepl("M",Sample) ~ "Environment",
                          grepl("T0|T3|T4",Sample) ~ "Experiment A",
                          grepl("S",Sample) ~ "Experiment B")) %>% 
  tidyr::separate(taxon_name, into = c("Domain", "Phylum","Class","Order", "Family","Genus"), extra = 'drop', remove = FALSE, sep =";")

#####################################
# Plot taxonomy of bacterial communities in microcosm
#####################################
#total reads per sample
tax_total_reads<- kaiju_table
mutate(Sample=case_when(Sample=="T0_S" ~ "S_0", TRUE~ Sample)) %>% 
  group_by(Sample, Domain) %>% 
  #summarize percentage
  summarize(Total = sum(reads))

#summarize on order level
summary_order<- kaiju_table %>% filter(Domain=="Bacteria") %>% 
  left_join(tax_total_reads, by = c("Sample", "Domain")) %>% 
  mutate(Proportion=Reads/Total) %>% 
  mutate(Order=case_when(Proportion < 0.01 ~ "Other taxa", TRUE~Order),
         Sample=factor(Sample, levels=c("M01","M02","M03","ML1","ML2","ML3","T0","T3C1","T3C2","T3C3","T4C1","T4C2","T4C3","T3J1",
                                        "T3J2","T3J3","T4J1","T4J2","T4J3","S_0","SC1_2","SC2_2","SC3_2","SJ1_2","SJ2_2","SJ3_2")))
  

#plot
summary_order %>% 
  filter(!Sample %in% c("M01","M02","M03","ML1","ML2","ML3", "T3C1","T3C2","T3C3","T3J1", "T3J2","T3J3")) %>% 
  mutate(Order= factor(Order, levels=c(summary_order %>% 
                                         filter(Proportion> 0.01) %>% 
                                         pull(Order) %>% unique() %>% 
                                         sort(),"Other taxa"))) %>% 
  ggplot(aes(x= Sample, y= Proportion*100, fill = Order))+
  geom_col()+
  facet_wrap(~Type, scales = "free_x")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Percent of classified reads (>1%)")+
  labs(fill = "Order")+
  theme_EF+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90))

#save the plot
ggsave("./Figures/Fig_4-metaG_taxa_Order.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       #height = 90, 
       scale = 2,
       dpi = 300)

#####################################
# Plot taxonomy of bacterial communities in situ
#####################################
#plot
summary_order %>% 
  filter(Sample %in% c("M01","M02","M03","ML1","ML2","ML3")) %>% 
  mutate(Order= factor(Order, levels=c(summary_order %>% 
                                         filter(Proportion> 0.01) %>% 
                                         pull(Order) %>% unique() %>% 
                                         sort(),"Other taxa"))) %>% 
  ggplot(aes(x= Sample, y= Proportion*100, fill = Order))+
  geom_col()+
  facet_wrap(~Type, scales = "free_x")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Percent of classified reads (>1%)")+
  labs(fill = "Order")+
  theme_EF+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90))

#save the plot
ggsave("./Figures/Fig_S2-metaG_taxa_Order.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       #height = 90, 
       scale = 2,
       dpi = 300)