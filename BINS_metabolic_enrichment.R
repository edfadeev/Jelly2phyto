require(dplyr)
require(ggplot2)
require(gridExtra)

#########################
##    Preparations     ##
#########################
#define work directory
wd <- "D:/UniVie/Projects/Microcosm_Piran/Analysis/MetaG/"

#import Anvio output
#KEGG modules in each bin
BINS_modules <- read.csv(paste0(wd,"09_METABOLISM/KEGG/BINS_modules.txt"),
                         sep = "\t")
#module definitions
modules_info <- read.csv(paste0(wd,"09_METABOLISM/KEGG/modules_info.txt"),
                         sep = "\t")

#taxonomy of bins 
BINS_summary <- read.csv(paste0(wd,"08_BINS/combined_bins_summary.txt"),
                         sep = "\t") %>% 
                select(group, db_name, starts_with("t_"))

##############################################
## Enriched KEGG pathways Jelly vs. Control ##
##############################################
#import enrichment test results
J_C_enr_mod <- read.csv(paste0(wd,"09_METABOLISM/KEGG/Jelly_control_met_enrichment.txt"),
                             sep = "\t") %>% 
  left_join(modules_info, by = c("accession" = "module"))

#plot
J_C_enr_mod %>% 
  mutate(subcategory=factor(subcategory)) %>% 
  group_by(associated_groups, subcategory) %>% 
  count(subcategory, name = "Num_pathways", .drop = FALSE) %>%
  ggplot(aes(x=subcategory, y=Num_pathways, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()

#check in which bins the enriched pathways are present
enriched_pathway_IDs<- J_C_enr_mod %>% 
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
ggsave(paste0("Figures/BINS_enr_pathway_Alpha.png"),
       plot = enr_mod_by_taxa[[1]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)

#Bacteroidia
ggsave(paste0("Figures/BINS_enr_pathway_Bac.png"),
       plot = enr_mod_by_taxa[[2]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)

#Gammaproteobacteria
ggsave(paste0("Figures/BINS_enr_pathway_Gamma.png"),
       plot = enr_mod_by_taxa[[3]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)



##############################################
## Enriched KEGG pathways SJ vs. SC ##
##############################################
#import enrichment test results
SJ_SC_enr_mod <- read.csv(paste0(wd,"09_METABOLISM/KEGG/SJ_SC_met_enrichment.txt"),
                        sep = "\t") %>% 
  left_join(modules_info, by = c("accession" = "module"))

#plot
SJ_SC_enr_mod %>% 
  filter(associated_groups!="NA") %>%  
  mutate(subcategory=factor(subcategory)) %>% 
  group_by(associated_groups, subcategory) %>% 
  count(subcategory, name = "Num_pathways", .drop = FALSE) %>%
  ggplot(aes(x=subcategory, y=Num_pathways, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()

#check in which bins the enriched pathways are present
SJ_SC_enriched_pathway_IDs<- SJ_SC_enr_mod %>% 
  filter(associated_groups!="NA") %>%  
  select(accession, associated_groups)

SJ_SC_BINS_enr_modules<- BINS_modules %>% 
  left_join(BINS_summary, by = "db_name") %>% 
  filter(grepl("T3|T4", group)) %>%  #subset only bins from the Jelly-Control
  left_join(SJ_SC_enriched_pathway_IDs, by = c("module" ="accession")) %>% 
  filter(!is.na(associated_groups))


#Calculate in how many bins enriched pathways were present
SJ_SC_enr_modules_bins_sum<- SJ_SC_BINS_enr_modules %>% 
  group_by(associated_groups, t_class, t_order, 
           t_family, t_genus, t_species,
           module_subcategory, module_name) %>% 
  count(module_name, name = "Num_of_bins", 
        .drop = FALSE)

#Calculate how many enriched pathways were present in each bin
SJ_SC_BINS_enr_modules_bin_sum<- SJ_SC_BINS_enr_modules %>% 
  group_by(associated_groups, db_name, t_class, t_order, 
           t_family, t_genus, t_species,
           module_subcategory) %>% 
  count(db_name, name = "Num_pathways", 
        .drop = FALSE)


#generate plots for each class
SJ_SC_enr_mod_by_taxa<- lapply(c("Alphaproteobacteria", "Bacteroidia","Gammaproteobacteria"),
                         function(x){
                           return(SJ_SC_BINS_enr_modules_bin_sum %>% 
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
ggsave(paste0("Figures/SJ_SC_BINS_enr_pathway_Alpha.png"),
       plot = SJ_SC_enr_mod_by_taxa[[1]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)

#Bacteroidia
ggsave(paste0("Figures/SJ_SC_BINS_enr_pathway_Bac.png"),
       plot = SJ_SC_enr_mod_by_taxa[[2]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)

#Gammaproteobacteria
ggsave(paste0("Figures/SJ_SC_BINS_enr_pathway_Gamma.png"),
       plot = SJ_SC_enr_mod_by_taxa[[3]],
       units = "cm",
       width = 25, height = 25, 
       scale = 2,
       dpi = 300)








heme and colours
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
