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

#define work directory
wd <- "C:/Users/fadee/ucloud/Projects/Microcosm_Piran/Analysis/MetaG/"


#####################################
## Summarize selected taxonomic rank
#####################################
tax <- "Class" # can be changed to "Phylum", "Class", "Order","Family", "Genus"

#import Kaiju output table
kaiju_table<- read.csv(paste0(wd,"07_TAXONOMY/Tax_Kaiju/kaiju_nr_euk_merged_",tax,".out"),
                         sep = "\t") %>% 
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
    summarize(Total = sum(percent)) %>% 
    mutate(Sample=case_when(Sample=="T0_S" ~ "S_0", TRUE~ Sample),
           Sample=factor(Sample, levels=c("M01","M02","M03","ML1","ML2","ML3","T0","T3C1","T3C2","T4C1","T4C2","T4C3","T3C3","T3J1",
                                          "T3J2","T3J3","T4J1","T4J2","T4J3","S_0","SC1_2","SC2_2","SC3_2","SJ1_2","SJ2_2","SJ3_2")),
           Group = case_when(grepl("M",Sample) ~ "Environment",
                             grepl("T0|T3|T4",Sample) ~ "Incubation",
                             grepl("S",Sample) ~ "Succession")) %>% 
    tidyr::separate(taxon_name, into = c("Domain", "Phylum","Class","Order", "Family","Genus"), extra = 'drop', remove = FALSE, sep =";")
  
#plot
kaiju_table%>% 
  filter(Total> 1#, Domain %in% c("Bacteria", "Eukaryota")
         ) %>% #threshold for visualization purposes
  ggplot(aes(x= Sample, y= Total, fill = get(tax)))+
  geom_col()+
  facet_wrap(~Group, scales = "free_x")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Percent of classified reads (>1%)")+
  labs(fill = paste(tax))+
  theme_EF+
  theme(legend.position = "right",
        axis.text.x = element_text(angle=90))

#####################################
## Proportions whithin bacterial community
#####################################
tax <- "Family" # can be changed to "Phylum", "Class", "Order","Family", "Genus"

#import Kaiju output table
kaiju_table<- read.csv(paste0(wd,"07_TAXONOMY/Tax_Kaiju/kaiju_nr_euk_merged_",tax,".out"),
                       sep = "\t") %>% 
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


tax_total_reads<- read.csv(paste0(wd,"07_TAXONOMY/Tax_Kaiju/kaiju_nr_euk_merged_",tax,".out"),
                       sep = "\t") %>% 
  mutate(Sample= gsub("/lisc/scratch/oceanography/efadeev/Microcosm_Piran/MicPir_sep/07_TAXONOMY/|_nr_euk_tax.out","",file)) %>% 
  mutate(Sample= gsub("[A-H]0[0-9]_","",Sample),
         #add Taxonomic domain
         Domain = case_when(grepl("Eukaryota", taxon_name)~"Eukaryota",
                            grepl("Bacteria", taxon_name)~"Bacteria",
                            grepl("Archaea", taxon_name)~"Archaea",
                            grepl("Viruses", taxon_name)~"Viruses")) %>%
  filter(!grepl("merged", Sample)) %>% 
  mutate(Sample=case_when(Sample=="T0_S" ~ "S_0", TRUE~ Sample)) %>% 
  group_by(Sample, Domain) %>% 
  #summarize percentage
  summarize(Total = sum(reads))


summary_order<- kaiju_table %>% filter(Domain=="Bacteria") %>% 
  left_join(tax_total_reads, by = c("Sample", "Domain")) %>% 
  mutate(Proportion=Reads/Total) %>% 
  mutate(Order=case_when(Proportion < 0.01 ~ "Other taxa", TRUE~Order),
         Sample=factor(Sample, levels=c("M01","M02","M03","ML1","ML2","ML3","T0","T3C1","T3C2","T3C3","T4C1","T4C2","T4C3","T3J1",
                                        "T3J2","T3J3","T4J1","T4J2","T4J3","S_0","SC1_2","SC2_2","SC3_2","SJ1_2","SJ2_2","SJ3_2"))) %>% 
  filter(!Sample %in% c("M01","M02","M03","ML1","ML2","ML3", "T3C1","T3C2","T3C3","T3J1", "T3J2","T3J3"))
  

summary_order.p<- summary_order %>% 
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
ggsave("./Figures/metaG_taxa_Order.pdf",
       plot = summary_order.p,
       units = "mm",
       width = 180,
       #height = 90, 
       scale = 2,
       dpi = 300)

kaiju_table %>% filter(Domain=="Bacteria") %>% 
  left_join(tax_total_reads, by = c("Sample", "Domain")) %>% 
  group_by(Type, Sample, Class, Order) %>% 
  summarize(Reads=sum(Reads), Total=Total) %>% 
  mutate(Proportion=Reads/Total) %>% 
  mutate(#Order=case_when(Proportion < 0.01 ~ "Other taxa", TRUE~Order),
         Sample=factor(Sample, levels=c("M01","M02","M03","ML1","ML2","ML3","T0","T3C1","T3C2","T3C3","T4C1","T4C2","T4C3","T3J1",
                                        "T3J2","T3J3","T4J1","T4J2","T4J3","S_0","SC1_2","SC2_2","SC3_2","SJ1_2","SJ2_2","SJ3_2"))) %>% 
  filter(!Sample %in% c("M01","M02","M03","ML1","ML2","ML3")) %>% 
  mutate(Group= case_when(grepl("T3C",Sample) ~ "T3C",
                          grepl("T4C",Sample) ~ "T4C",      
                          grepl("T3J",Sample) ~ "T3J",
                          grepl("T4J",Sample) ~ "T4J",      
                          grepl("SC",Sample) ~ "SC",     
                          grepl("SJ",Sample) ~ "SJ", TRUE ~ Sample)) %>% 
  group_by(Type, Group, Class, Order) %>% 
  summarize(Order_mean=mean(Proportion), Order_se=se(Proportion), Order_min=min(Proportion), Order_max=max(Proportion)) %>% View()





%>% 
  ggplot(aes(x= Group, y= Order_mean*100, fill = Order))+
  geom_col()+
  #geom_text(aes(y=-1, label=paste(round(Total_mean/10E6,2), round(Total_se/10E6,2), sep="±")))+
  facet_wrap(~Type, scales = "free_x")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Percent of classified reads (>1%)")+
  labs(fill = "Order")+
  theme_EF+
  theme(legend.position = "right",
        axis.text.x = element_text(angle=90))


#####################################
## Dissimilarity based on taxonomic composition
#####################################

class_tax_matrix<- kaiju_table %>% mutate(Class=case_when(is.na(Class)~"Unclassified", TRUE~Class)) %>% 
                                    group_by(Sample, Class) %>% 
                                    summarize(Total=sum(Total)/100) %>% 
                                    select(Sample, Class, Total) %>% 
                                    tidyr::spread(key="Sample", value="Total") %>% 
                                    tibble::column_to_rownames("Class") %>% 
                                    mutate_all(as.numeric) %>% 
                                    mutate_all(funs(tidyr::replace_na(.,0)))

# Running NMDS in vegan (metaMDS)
NMDS <-  metaMDS(t(class_tax_matrix),
                 distance = "bray",
                 #k = 3,
                 maxit = 999, 
                 trymax = 999,
                 wascores = FALSE,
                 autotransform = FALSE,
                 tidy= "sites",na.rm = TRUE)
#plot
NMDS.scores <- as.data.frame(scores(NMDS)) %>% 
  tibble::rownames_to_column(var ="Microcosm") %>% 
  mutate(group=case_when(grepl("M",Microcosm)~"Environmental",
                         grepl("T0|S_0",Microcosm)~"Inoculum",
                         grepl("T3|T4",Microcosm)~"Experiment A",
                         grepl("S",Microcosm)~"Experiment B"),
         Type=case_when(grepl("C", Microcosm) ~ "Control",
                        grepl("J", Microcosm) ~ "Jelly-OM",
                        TRUE~"Other"))

NMDS.scores %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, shape = Type, color=group,  label = Microcosm))+
  geom_point(size =5, colour = "black")+
  geom_point(size =4)+
  geom_text(nudge_y = -0.05, size =4)+
  scale_color_manual(values = c("#009E73", "#F0E442", "#0072B2", 
                                "#D55E00"))+
  theme_EF


#save the plot
ggsave("./Figures/metaG_class_NMDS.pdf",
       plot = last_plot(),
       units = "mm",
       #width = 180,
       #height = 90, 
       scale = 2,
       dpi = 300)


#test whether the differences between the runs are significant
distmat <- 
  vegdist(t(class_tax_matrix), method = "bray",na.rm = TRUE)

df <- NMDS.scores %>% select(Microcosm, group, Type) 

adonis_all <- adonis2(distmat ~ group+Type, df,
                      permutations=999)

adonis_all

#posthoc to check which groups are different
pairwiseAdonis::pairwise.adonis(distmat, factors=df$Type, 
                                p.adjust.m='bonferroni', perm = 999)



#test whether the differences between the runs are significant
df_sub <- NMDS.scores %>% select(Microcosm, group, Type) %>% filter(group %in% c("Experiment A"#,
                                                                                 #"Experiment B"
                                                                                 ),
                                                                    #Type=="Jelly-OM"
                                                                    )

distmat_sub  <- 
  vegdist(t(class_tax_matrix[,df_sub$Microcosm]), method = "bray",na.rm = TRUE)

adonis_all <- adonis2(distmat_sub ~ Type, df_sub,
                      permutations=999)

adonis_all




#posthoc to check which ponds are different
groups <- df[["group"]]
mod <- vegan::betadisper(distmat, groups)
vegan::permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)

#####################################
## Automartically summarize and plot all taxonomic ranks
#####################################
#import Kaiju output tables for the different taxonomic ranks
lapply(c("Phylum", "Class", "Order","Family", "Genus"), function(tax) {
kaiju_table<- read.csv(paste0(wd,"07_TAXONOMY/Tax_Kaiju/kaiju_nr_euk_merged_",tax,".out"),
                       sep = "\t") %>% 
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

ggsave(paste0("Figures/Kaiju_Bac_",tax,".pdf"),
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

ggsave(paste0("Figures/Kaiju_Euk_",tax,".pdf"),
       plot = last_plot(),
       units = "mm",
       width = 250, height = 250, 
       scale = 2,
       dpi = 300)

})
