################################################
## This script summarizes the enrichment tests
## between SJ and SC in the second incubation
################################################

require(dplyr)
require(ggplot2)
require(gridExtra)

#########################
##    Preparations     ##
#########################
#define work directory
wd <- "C:/Users/fadee/ucloud/Projects/Microcosm_Piran/Analysis/MetaG/"

#taxonomy of bins 
BINS_summary <- read.csv(paste0(wd,"08_BINS/combined_bins_summary.txt"),
                         sep = "\t") %>% 
  select(group, db_name, starts_with("t_"))

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

##############################################
## Enriched COGs SJ vs. SC                  ##
##############################################
#import Anvio COG20 definitions
COG20_table <- read.csv(paste0(wd,"09_METABOLISM/COG20/COG20_definitions_table.txt"),
                        sep= " ") %>% 
  rename(accession = COG20_FUNCTION_accession,
         function.=COG20_FUNCTION_function.,
         category_accession=COG20_CATEGORY_accession,
         category=COG20_CATEGORY_function.)

#import enrichment results table from Anvio
SJ_SC_COG20_enr <- read.csv(paste0(wd,"09_METABOLISM/COG20/SJ_SC_COG20_FUN_enrichment.txt"),
                      sep = "\t") %>% select(-c("key")) %>% 
  left_join(COG20_table,  by = c("accession","function."))

#explore the output 
str(SJ_SC_COG20_enr)

#filter the table according to enrichment groups or COG categories
View(
  SJ_SC_COG20_enr %>% 
    filter(associated_groups == "Jelly", #can be changed to "Control"
           category_accession == "T", # any of "T" "E" "L" "I" "S" "R" "Q" "K" "J" "C" "P" "G" "D" "X" "V" NA  "U" "F" "O" "A" "H" "M" "B" "W"
           )
  )

#summarize by category the number of enriched COGs  
SJ_SC_COG20_enr_CAT_sum<- SJ_SC_COG20_enr %>% 
  group_by(associated_groups, category, category_accession) %>% 
  summarize(COGs_n=n())

#plot
SJ_SC_COG20_enr %>% 
  mutate(category=factor(category)) %>% 
  group_by(associated_groups, category) %>% 
  count(category, name = "Num_funs", .drop = FALSE) %>%
  ggplot(aes(x=category, y=Num_funs, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()
  
##############################################
## Enriched CAZymes SJ vs. SC               ##
##############################################
#import enrichment results
SJ_SC_CAZy_enr <- read.csv(paste0(wd,"09_METABOLISM/CAZy/SJ_SC_CAZy_enrichment.txt"),
                       sep = "\t") %>% select(-c("key", "accession")) %>% 
                rename(accession = function.) %>% 
                mutate(Class = case_when(grepl("GH", accession) ~ "Glycoside_Hydrolase",
                                         grepl("GT", accession) ~ "Glycosyl_Transferase",
                                         grepl("PL", accession) ~ "Polysaccharide_Lyase",
                                         grepl("CE", accession) ~ "Carbohydrate_Esterase",
                                         grepl("AA", accession) ~ "Auxiliary_Activities",
                                         grepl("CBM", accession) ~ "Carbohydrate_binding_module"))
                                               
                                               
#explore table
str(SJ_SC_CAZy_enr)

#filter the table according to enrichment groups or CAZy class
View(
  SJ_SC_CAZy_enr %>% 
    filter(associated_groups == "SJ", #can be changed to "Control"
           Class == "Glycoside_Hydrolase", # any of the classes above
    ))

#plot
SJ_SC_CAZy_enr %>% 
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
SJ_SC_Prot_enr <- read.csv(paste0(wd,"09_METABOLISM/MEROPS/SJ_SC_MEROPS_enrichment.txt"),
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
    filter(associated_groups == "SJ", #can be changed to "SC"
           Type == "Aminopeptidase"
    ))

#plot
SJ_SC_Prot_enr %>% 
  filter(Type %in% c("Aminopeptidase", "Peptidase")) %>% 
  mutate(Family=factor(Family)) %>% 
  group_by(associated_groups, Family) %>% 
  count(Family, name = "Num_Prot", .drop = FALSE) %>%
  filter(Num_Prot>10) %>% 
  ggplot(aes(x=Family, y=Num_Prot, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()


##############################################
## Enriched KEGG pathways SJ vs. SC         ##
##############################################
#KEGG modules in each bin
BINS_modules <- read.csv(paste0(wd,"09_METABOLISM/KEGG/BINS_modules.txt"),
                         sep = "\t")
#module definitions
modules_info <- read.csv(paste0(wd,"09_METABOLISM/KEGG/modules_info.txt"),
                         sep = "\t")

#import enrichment test results
SJ_SC_enr_mod <- read.csv(paste0(wd,"09_METABOLISM/KEGG/SJ_SC_met_enrichment.txt"),
                        sep = "\t") %>% 
  left_join(modules_info, by = c("accession" = "module"))

#plot
SJ_SC_enr_mod %>% 
  mutate(subcategory=factor(subcategory)) %>% 
  group_by(associated_groups, subcategory) %>% 
  count(subcategory, name = "Num_pathways", .drop = FALSE) %>%
  ggplot(aes(x=subcategory, y=Num_pathways, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()

#check in which bins the enriched pathways are present
enriched_pathway_IDs<- SJ_SC_enr_mod %>% 
  select(accession, associated_groups)

BINS_enr_modules<- BINS_modules %>% 
  left_join(BINS_summary, by = "db_name") %>% 
  filter(grepl("T3|T4", group)) %>%  #subset only bins from the Jelly-Control
  left_join(enriched_pathway_IDs, by = c("module" ="accession")) %>% 
  filter(!is.na(associated_groups))


#Calculate in how many bins enriched pathways were present
enr_modules_bins_sum<- BINS_enr_modules %>% 
  group_by(associated_groups, t_class, t_order, 
           t_family, t_genus, t_species,
           module_subcategory, module_name) %>% 
  count(module_name, name = "Num_of_bins", 
        .drop = FALSE)

#Calculate how many enriched pathways were present in each bin
BINS_enr_modules_bin_sum<- BINS_enr_modules %>% 
  group_by(associated_groups, db_name, t_class, t_order, 
           t_family, t_genus, t_species,
           module_subcategory) %>% 
  count(db_name, name = "Num_pathways", 
        .drop = FALSE)


#generate plots for each class
enr_mod_by_taxa<- lapply(c("Alphaproteobacteria", "Bacteroidia","Gammaproteobacteria"),
                         function(x){
                           return(BINS_enr_modules_sum %>% 
                                    filter(t_class==x ) %>% #filter main classes
                                    ggplot(aes(x=module_subcategory, y=Num_pathways, fill= t_genus))+
                                    geom_col(position = "stack")+
                                    scale_fill_manual(values = c(tol21rainbow, "gray50","yellow","red","blue"))+
                                    guides(fill=guide_legend(title = paste0(x), ncol =1))+
                                    facet_wrap(.~associated_groups)+
                                    theme_EF+
                                    coord_flip())
                         })

#save plot of each class separately
#Alphaproteobacteria
ggsave(paste0("Figures/SJ_vs_SC_enr_pathway_Alpha.png"),
       plot = enr_mod_by_taxa[[1]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)

#Bacteroidia
ggsave(paste0("Figures/SJ_vs_SC_enr_pathway_Bac.png"),
       plot = enr_mod_by_taxa[[2]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)

#Gammaproteobacteria
ggsave(paste0("Figures/SJ_vs_SC_enr_pathway_Gamma.png"),
       plot = enr_mod_by_taxa[[3]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)


