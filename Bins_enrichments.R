require(ggplot2)
require(dplyr)

##############################################
## Enriched KEGG pathways Jelly vs. Control ##
##############################################
#taxonomy of bins 
BINS_summary <- read.csv("data/combined_bins_summary.txt",sep = "\t") %>% 
  mutate(group= factor(group, levels =c ("T0","T3J","T4J",
                                  "T3C","T4C", "SJ","SC",
                                  "M","ML")))
#KEGG modules in each bin
BINS_modules <- read.csv("data/BINS_modules.txt",
                         sep = "\t") 

#import enrichment test results
J_C_enr_T4_mod <- read.csv("data/Jelly_Control_T4_enrichment.txt",
                           sep = "\t") %>% mutate(Time="T4") %>% 
  left_join(read.csv("data/modules_info.txt", sep = "\t"), by = c("accession" = "module")) #module definitions

J_C_enr_T4_bins<- J_C_enr_T4_mod %>% tidyr::separate_longer_delim(sample_ids, ",") %>% 
  filter(adjusted_q_value< 0.1) %>% #FDR threshold
  left_join(BINS_summary, by = c("sample_ids"="db_name")) %>% 
  left_join(BINS_modules, by = c("sample_ids"="db_name", "accession"="module") )


#Calculate how many enriched pathways were present in each bin
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
  
#families of bins that have also enriched modules
tax_fam<- unique(BINS_enr_modules_sum$t_family)

colour_pallette <- c(tol21rainbow)
names(colour_pallette)<- c("Alphaproteobacteria_uncl",  "Puniceispirillaceae", "Rhodobacteraceae", "Sphingomonadaceae",
                           "Bacteroidia_uncl", "Crocinitomicaceae", "Flavobacteriaceae", "Salibacteraceae","Schleiferiaceae", 
                           "Alteromonadaceae", "Gammaproteobacteria_uncl", "Halieaceae", "Litoricolaceae", "Nitrincolaceae", 
                           "Pseudoalteromonadaceae",  "Spongiibacteraceae", "Vibrionaceae", "Balneolaceae", "Other taxa")

BINS_enr_modules.p <-BINS_enr_modules_sum %>% 
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



BINS_summary.p <- BINS_summary %>% filter(t_domain!="", grepl("T4",group)) %>% 
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

cowplot::plot_grid(BINS_summary.p, BINS_enr_modules.p,  ncol = 1)

ggsave(paste0("Figures/J_vs_C_enr_legend.pdf"),
       plot = last_plot(),
       units = "mm",
       width = 180, height = 180, 
       scale = 2,
       dpi = 300)


##############################################
## Enriched KEGG pathways JS vs. JC ##
##############################################
#taxonomy of bins 
BINS_summary <- read.csv("data/combined_bins_summary.txt",sep = "\t") %>% 
  mutate(group= factor(group, levels =c ("T0","T3J","T4J",
                                         "T3C","T4C", "SJ","SC",
                                         "M","ML")))
#KEGG modules in each bin
BINS_modules <- read.csv("data/BINS_modules.txt",
                         sep = "\t") 

#import enrichment test results
J_C_enr_T4_mod <- read.csv("data/T3_T4_combined/SJ_SC_met_enrichment.txt",
                           sep = "\t") %>% mutate(Time="T4") %>% 
  left_join(read.csv("data/modules_info.txt", sep = "\t"), by = c("accession" = "module")) #module definitions

J_C_enr_T4_bins<- J_C_enr_T4_mod %>% tidyr::separate_longer_delim(sample_ids, ",") %>% 
  filter(adjusted_q_value< 0.1) %>% #FDR threshold
  left_join(BINS_summary, by = c("sample_ids"="db_name")) %>% 
  left_join(BINS_modules, by = c("sample_ids"="db_name", "accession"="module") )


#Calculate how many enriched pathways were present in each bin
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

#families of bins that have also enriched modules
tax_fam<- unique(BINS_enr_modules_sum$t_family)

colour_pallette <- c(tol21rainbow)
names(colour_pallette)<- c("Alphaproteobacteria_uncl",  "Puniceispirillaceae", "Rhodobacteraceae", "Sphingomonadaceae",
                           "Bacteroidia_uncl", "Crocinitomicaceae", "Flavobacteriaceae", "Salibacteraceae","Schleiferiaceae", 
                           "Alteromonadaceae", "Gammaproteobacteria_uncl", "Halieaceae", "Litoricolaceae", "Nitrincolaceae", 
                           "Pseudoalteromonadaceae",  "Spongiibacteraceae", "Vibrionaceae", "Balneolaceae", "Other taxa")

BINS_enr_modules.p <-BINS_enr_modules_sum %>% 
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



BINS_summary.p <- BINS_summary %>% filter(t_domain!=""#, grepl("T4",group)
                                          ) %>% 
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

cowplot::plot_grid(BINS_summary.p, BINS_enr_modules.p,  ncol = 1)

ggsave(paste0("Figures/J_vs_C_enr_legend.pdf"),
       plot = last_plot(),
       units = "mm",
       width = 180, height = 180, 
       scale = 2,
       dpi = 300)






##############################################
## Enriched CAZymes Jelly vs. Control       ##
##############################################
#import enrichment results
J_C_CAZy_enr <- read.csv("data/Jelly_control_T4_CAZy_enrichment.txt",
                         sep = "\t") %>% select(-c("key", "accession")) %>% 
  rename(accession = function.) %>% 
  mutate(Class = case_when(grepl("GH", accession) ~ "Glycoside_Hydrolase",
                           grepl("GT", accession) ~ "Glycosyl_Transferase",
                           grepl("PL", accession) ~ "Polysaccharide_Lyase",
                           grepl("CE", accession) ~ "Carbohydrate_Esterase",
                           grepl("AA", accession) ~ "Auxiliary_Activities",
                           grepl("CBM", accession) ~ "Carbohydrate_binding_module"))


#explore table
str(J_C_CAZy_enr)

#filter the table according to enrichment groups or CAZy class
View(
  J_C_CAZy_enr %>% 
    filter(associated_groups == "Jelly", #can be changed to "Control"
           Class == "Glycoside_Hydrolase", # any of the classes above
    ))

#plot
J_C_CAZy_enr %>% 
  mutate(Class=factor(Class)) %>%  
  group_by(associated_groups, Class) %>% 
  count(Class, name = "Num_CAZymes", .drop = FALSE) %>%
  ggplot(aes(x=Class, y=Num_CAZymes, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()

##############################################
## Enriched Proteases SJ vs. SC               ##
##############################################
#import enrichment results
SJ_SC_Prot_enr <- read.csv("data/Jelly_control_T4_MEROPS_enrichment.txt",
                           sep = "\t") %>% 
  tidyr::separate(function., into = c('Annotation',"function."), sep ="_\\(", remove = TRUE, extra = "merge") %>% 
  tidyr::separate(function., into = c('Tax',"MEROPS_ID"), sep ="_\\[", remove = TRUE, extra = "merge") %>% 
  mutate(Tax=gsub("\\(|\\)|\\[|\\]|\\{|\\}", "", Tax),
         MEROPS_ID= gsub("\\[|\\]", "", MEROPS_ID),
         Family = gsub("\\..*", "", MEROPS_ID),
         Type = case_when(grepl("non-", Annotation, ignore.case = TRUE)~ "Non-peptidase",
                          grepl("aminopeptidase", Annotation, ignore.case = TRUE)~ "Aminopeptidase",
                          grepl("peptidase", Annotation, ignore.case = TRUE)~ "Peptidase",
                          TRUE ~ "Other"))

#explore table
str(SJ_SC_Prot_enr)

#filter the table according to enrichment groups
View(
  SJ_SC_Prot_enr %>% 
    filter(associated_groups == "Jelly", #can be changed to "Control"
           Type == "Aminopeptidase"
    ))

#plot
SJ_SC_Prot_enr %>% 
  filter(Type %in% c("Aminopeptidase", "Peptidase"), 
         associated_groups %in% c("Jelly","Control"), adjusted_q_value < 0.1) %>% 
  mutate(Family=factor(Family)) %>% 
  group_by(associated_groups, Family) %>% 
  count(Family, name = "Num_enr.genes", .drop = FALSE) %>%
  #filter(Num_enr.genes>5) %>% 
  ggplot(aes(x=Family, y=Num_enr.genes, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()

##############################################
## Enriched COGs Jelly vs. Control          ##
##############################################
#import Anvio COG20 definitions
COG20_table <- read.csv("data/COG20_definitions_table.txt",
                        sep= " ") %>% 
  rename(accession = COG20_FUNCTION_accession,
         function.=COG20_FUNCTION_function.,
         category_accession=COG20_CATEGORY_accession,
         category=COG20_CATEGORY_function.)

#import enrichment results table from Anvio
J_C_COG20_enr <- read.csv("data/Jelly_control_T4_COG20_FUN_enrichment.txt",
                          sep = "\t") %>% select(-c("key")) %>% 
  left_join(COG20_table,  by = c("accession","function."))

#explore the output 
str(J_C_COG20_enr)

#filter the table according to enrichment groups or COG categories
View(
  J_C_COG20_enr %>% 
    filter(associated_groups == "Jelly", #can be changed to "Control"
           category_accession == "T", # any of "T" "E" "L" "I" "S" "R" "Q" "K" "J" "C" "P" "G" "D" "X" "V" NA  "U" "F" "O" "A" "H" "M" "B" "W"
    )
)

#summarize by category the number of enriched COGs  
J_C_COG20_enr_CAT_sum<- J_C_COG20_enr %>% 
  group_by(associated_groups, category, category_accession) %>% 
  summarize(COGs_n=n())

#plot
J_C_COG20_enr %>% 
  filter(adjusted_q_value< 0.05, !category_accession %in% c("S","R"), !is.na(category_accession)) %>% 
  mutate(category=factor(category)) %>% 
  group_by(associated_groups, category) %>% 
  count(category, name = "Num_funs", .drop = FALSE) %>%
  ggplot(aes(x=category, y=Num_funs, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()
