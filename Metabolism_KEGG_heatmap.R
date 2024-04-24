require(dplyr)
require(pheatmap)

wd <- "D:/UniVie/Projects/Microcosm_Piran/Analysis/MetaG/MicPir_sep/Metabolism"

KO_names <- read.csv(paste0(wd,"/KEGG_hits.txt"),
                   sep = "\t") %>% 
            select(enzyme,enzyme_definition) %>% 
            unique()

K_hits_matrix <- read.csv(paste0(wd,"/KEGG-enzyme_hits-MATRIX.txt"),
                          sep = "\t")


K_modules_melted <- K_modules <- read.csv(paste0(wd,"/KEGG_modules.txt"),
                                          sep = "\t") %>% 
  tidyr::separate_rows(enzyme_hits_in_module, sep = ',') %>% 
  dplyr::rename(enzyme=enzyme_hits_in_module) %>% 
  select(enzyme, module_name, module_category,module_subcategory) %>% 
  unique()

K_hits_merged<- K_hits_matrix %>% 
  left_join(KO_names, by = "enzyme") %>% 
  left_join(K_modules_melted, by = "enzyme") 


AA_met <- K_hits_merged %>% 
  #filter(module_category == "Amino acid metabolism") %>% 
  #"Amino acid metabolism"
  #"Carbohydrate metabolism"
  #Metabolism of cofactors and vitamins
  #Energy metabolism
  select(1:10) %>% 
  mutate_at(c(2:10), as.numeric) %>% 
  unique() %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("enzyme") #%>% 
  #filter(rowSums(.) >= 200) 

KO_subcategory<- K_hits_merged %>% 
  select(enzyme, module_subcategory) %>% 
  unique() %>% 
  group_by(enzyme) %>% 
  summarize(subcategory = paste0(module_subcategory, collapse = ",")) %>%  
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("enzyme")

KO_category<- K_hits_merged %>% 
  select(enzyme, module_category) %>% 
  unique() %>% 
  group_by(enzyme) %>% 
  summarize(category = paste0(module_category, collapse = ",")) %>%  
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("enzyme")


pheatmap::pheatmap(as.matrix(AA_met),
                   scale = "row",
                   cluster_rows = T,
                   show_rownames = F,
                   annotation_row= KO_category, 
                   annotation_legend = FALSE, filename = "heatmap_all_corr.png")

