require(dplyr)
require(tidyr)
require(ggplot2)

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#colours
greens <- c(
  "#006400",  # DarkGreen
  "#228B22",  # ForestGreen
  "#32CD32",  # LimeGreen
  "#7CFC00",  # LawnGreen
  "#ADFF2F",  # GreenYellow
  "#98FB98"   # PaleGreen
)

################################################
#import raw table
################################################
phyto_counts<- read.table("data/Phytoplankton_counts.txt", dec=",", fill=TRUE, sep="\t", header = TRUE) 

################################################
#plot Fig. 7 - Phytoplankton composition
################################################
phyto_counts %>% 
  reshape2::melt()%>% 
  mutate(Taxa=tolower(Taxa)) %>%
  mutate(Taxa= factor(Taxa, c("coccolithophorids","diatoms","dinoflagellates","nanoflagellates","silicoflagellates","algae non ident.")),
         variable=factor(variable,c("SC1","SC2","SC3","SJ1","SJ2","SJ3","Fito.Mix","MO2","ML2")),
         panel= case_when(variable %in% c("Fito.Mix","MO2","ML2")~ "B",
                          TRUE~"A")) %>% 
  ggplot(aes(x=variable, y=value, group = variable, fill = Taxa))+
  geom_col()+
  scale_fill_manual(values=greens)+
  facet_grid(cols=vars(panel),scales="free",space="free_x",switch="x")+
  ylab("Abundance (cells mL-1)")+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Fig_7-phyto_prop.pdf",
       plot = last_plot(),
       units = "mm",
       #width = 90,
       #height = 90, 
       scale = 2,
       dpi = 300)