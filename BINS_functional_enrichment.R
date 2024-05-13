require(dplyr)
require(ggplot2)
require(gridExtra)

#########################
##    Preparations     ##
#########################
#define work directory
wd <- "D:/UniVie/Projects/Microcosm_Piran/Analysis/MetaG/"

#import Anvio output
COG20_table <- read.csv(paste0(wd,"09_METABOLISM/COG20/COG20_definitions_table.txt"),
                        sep= " ") %>% 
          rename(accession = COG20_FUNCTION_accession,
         function.=COG20_FUNCTION_function.,
         category_accession=COG20_CATEGORY_accession,
         category=COG20_CATEGORY_function.)

#taxonomy of bins 
BINS_summary <- read.csv(paste0(wd,"08_BINS/combined_bins_summary.txt"),
                         sep = "\t") %>% 
  select(group, db_name, starts_with("t_"))

##############################################
## Enriched functions Jelly vs. Control     ##
##############################################
J_C_COG20_enr <- read.csv(paste0(wd,"COG20/Jelly_control_COG20_FUN_enrichment.txt"),
                      sep = "\t") %>% 
  left_join(COG20_table,  by = c("accession","function."))

J_C_COG20_enr_CAT_sum<- J_C_COG20_enr %>% 
  group_by(associated_groups, category, category_accession) %>% 
  summarize(COGs_n=n())

#plot
J_C_COG20_enr %>% 
  mutate(category=factor(category)) %>% 
  group_by(associated_groups, category) %>% 
  count(category, name = "Num_funs", .drop = FALSE) %>%
  ggplot(aes(x=category, y=Num_funs, fill= associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  coord_flip()

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#check in which bins the enriched pathways are present
enriched_funs_IDs<- J_C_COG20_enr %>% 
  select(accession, associated_groups)

J_C_BINS_enr_CAT<- BINS_modules %>% 
  left_join(BINS_summary, by = "db_name") %>% 
  filter(grepl("T3|T4", group)) %>%  #subset only bins from the Jelly-Control
  left_join(enriched_funs_IDs, by = c("module" ="accession")) %>% 
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

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
##############################################
## Enriched CAZymes Jelly vs. Control       ##
##############################################
J_C_CAZy_enr <- read.csv(paste0(wd,"/CAZy/Jelly_control_CAZy_enrichment.txt"),
                       sep = "\t") %>% 
  left_join(COG20_table,  by = c("accession","function."))
