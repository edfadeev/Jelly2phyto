require(dplyr)
require(ggplot2)

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


#taxonomy of bins 
BINS_summary <- read.csv("data/combined_bins_summary.txt",
                         sep = "\t")

BINS_summary %>% 
  group_by(group, t_class, t_family) %>% 
  summarize(n_bins=n()) %>% 
  ggplot(aes(x=group, y=n_bins, fill = t_class))+
  geom_col()+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF


BINS_summary %>% 
  filter(t_class=="Gammaproteobacteria") %>% 
  group_by(group, t_genus) %>% 
  summarize(n_bins=n()) %>% 
  ggplot(aes(x=group, y=n_bins, fill = t_genus))+
  geom_col()+
  scale_fill_manual(values = c(tol21rainbow, "gray50","yellow","red","blue"))+
  theme_EF





BINS_modules <- read.csv("data/BINS_modules.txt",sep = "\t")

modules_info <- read.csv("data/modules_info.txt",sep = "\t")

enriched_modules_J_C <- read.csv("data/Jelly_control_met_enrichment.txt", sep = "\t") %>% 
                      left_join(modules_info, by = c("accession" = "module"))%>% 
  filter(adjusted_q_value < 0.05) %>% 
  group_by(associated_groups, subcategory) %>% 
  summarize(n_pathways=n()) %>% 
  mutate(Comp="J_C")

enriched_modules_SJ_SC <- read.csv("data/SJ_SC_met_enrichment.txt",sep = "\t") %>% 
                      left_join(modules_info, by = c("accession" = "module"))%>% 
  filter(adjusted_q_value < 0.05) %>% 
  group_by(associated_groups, subcategory) %>% 
  summarize(n_pathways=n()) %>% 
  mutate(Comp="SJ_SC")

enriched_modules_J_SJ <- read.csv("data/Jelly_SJ_met_enrichment.txt",sep = "\t") %>% 
                              left_join(modules_info, by = c("accession" = "module"))%>% 
  filter(adjusted_q_value < 0.05) %>% 
  group_by(associated_groups, subcategory) %>% 
  summarize(n_pathways=n()) %>% 
  mutate(Comp="J_SJ")

rbind(enriched_modules_J_SJ, enriched_modules_SJ_SC, enriched_modules_J_C) %>% 
  ggplot(aes(x=subcategory, y=n_pathways, fill= associated_groups, group = associated_groups))+
  geom_col(width=0.8, position = "dodge")+
  scale_fill_manual(values = tol21rainbow)+
  facet_grid(Comp~.)+
  theme_EF+
  theme(axis.text.x = element_text(angle = 90))




COG20_enr <- read.csv(paste0(wd,"/Jelly_control_COG20_FUN_enrichment.txt"),
                      sep = "\t")








#import tables
wd <- "D:/UniVie/Projects/Microcosm_Piran/Analysis/MetaG/MicPir_sep/Metabolism/"

COG20_table <- read.table(paste0(wd,"COG20/COG20_definitions_table.txt"),
                         sep= " ") %>% 
                rename(accession = COG20_FUNCTION_accession,
                       function.=COG20_FUNCTION_function.,
                       category_accession=COG20_CATEGORY_accession,
                       category=COG20_CATEGORY_function.)

COG20_enr <- read.csv(paste0(wd,"COG20/Jelly_control_COG20_FUN_enrichment.txt"),
                        sep = "\t") %>% 
  left_join(COG20_table,  by = c("accession","function."))



COG20_enr <- read.csv(paste0(wd,"COG20/SJ_SC_COG20_FUN_enrichment.txt"),
                      sep = "\t") %>% 
  left_join(COG20_table,  by = c("accession","function."))


  
COG20_enr %>% 
  group_by(associated_groups, category, category_accession) %>% 
  summarize(COGs_n=n())



CAZy_enr <- read.csv("data/Jelly_control_CAZy_enrichment.txt",
                      sep = "\t") %>% 
  left_join(COG20_table,  by = c("accession","function."))





read.csv("./03_CONTIGS/COG20.txt",
         sep = "\t") %>%
  mutate(across(everything(),~ gsub("!!!.*","", .))) %>% 
  group_by(gene_callers_id, source) %>% 
  pivot_wider(names_from = source,
              values_from = c(accession,function.,e_value),
              names_glue = "{source}_{.value}") %>% 
  ungroup() %>% 
  select(COG20_FUNCTION_accession,COG20_FUNCTION_function.,
         COG20_CATEGORY_accession, COG20_CATEGORY_function.) %>% 
  unique() %>% 
  write.table("./03_CONTIGS/COG20_definitions_table.txt")

anvi-interactive --manual-mode \
-d $WORKDIR/03_CONTIGS/KEGG-module_pathwise_completeness-MATRIX.txt \
-t $WORKDIR/03_CONTIGS/KEGG-module_pathwise_completeness-MATRIX.txt.newick \
-p $WORKDIR/03_CONTIGS/Metabolism.db \
--title "Metabolism Heatmap" \
--server-only -P 5678


#produce interactive visualization
anvi-matrix-to-newick $WORKDIR/03_CONTIGS/KEGG-module_pathwise_presence-MATRIX.txt

# dry run to get the profile db:
anvi-interactive -d $WORKDIR/03_CONTIGS/KEGG-module_pathwise_presence-MATRIX.txt \
-p $WORKDIR/03_CONTIGS/Metabolism_presence.db \
--manual-mode \
--dry-run

# import the state file (changes from the Anvio tutorial)
anvi-import-state -s $WORKDIR/03_CONTIGS/metabolism_state.json \
-p $WORKDIR/03_CONTIGS/Metabolism_presence.db \
-n default

#import modules to the profile
anvi-import-misc-data $WORKDIR/modules_info.txt \
-p $WORKDIR/03_CONTIGS/Metabolism_presence.db -t items

anvi-interactive --manual-mode \
-d $WORKDIR/03_CONTIGS/KEGG-module_pathwise_presence-MATRIX.txt \
-t $WORKDIR/03_CONTIGS/KEGG-module_pathwise_presence-MATRIX.txt.newick \
-p $WORKDIR/03_CONTIGS/Metabolism_presence.db \
--title "Metabolism Heatmap presence" \
--server-only -P 5678