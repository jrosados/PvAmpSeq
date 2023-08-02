
library(tidyverse)
library(RColorBrewer) 
library(ggpubr)
library(cowplot)
library(here)

rm(list = ls())
SI_controls <- readRDS(here("raw_data","seq_tbl_end.rds"))

idea<- read.csv(here("raw_data","Ideal_data.csv"))

#Samples used for mimic mixed infections
#"AR-079-0", "AR-093-0"
cbp2 <- c("#459CC1", "#FF9933","#D55E00","#E69F00","#56B4E9")


#Aqua dark
AR079 <- c("Chr05-2" = "#459CC1","Chr03-2" ="#459CC1",
           "Chr02-1" ="#459CC1", "Chr09-1" ="#459CC1",
           "Chr01-1" = "#459CC1","Chr14-2" = "#459CC1")


#dark yellow
AR093 <- c("Chr05-1" = "#FCE515", "Chr03-1" = "#FCE515",
           "Chr02-2" = "#FCE515","Chr09-2" ="#FCE515",
           "Chr01-2" ="#FCE515","Chr14-1" ="#FCE515")

others<- c("Chr05-3" = "#D55E00", "Chr05-4" = "#E69F00",
           "Chr02-4" = "#D55E00", "Chr03-3" = "#D55E00","Chr03-4" = "#E69F00",
           "Chr14-3" = "#D55E00")


SI_controls %>%
    select(marker_id,count,sample, Haplotype, Frequency, info)%>%
    filter(marker_id%in% c("Chr01", "Chr05", "Chr02","Chr03",
                           "Chr09","Chr13","Chr14"),
           sample %in% c("AR79_093_1_1","AR79_093_1_10","AR79_093_1_50","AR79_093_1_100","AR79_093_1_500", "AR79_093_1_1000", 
                         "AR79_093_10_1", "AR79_093_50_1","AR79_093_100_1", "AR79_093_100_1","AR79_093_500_1","AR79_093_1000_1")) %>%
    rbind(idea) %>%
    mutate(Frequency = Frequency*100) %>%
    mutate(infoD = factor(info, levels = c("1_1","1_10","1_50","1_100","1_500","1_1000","10_1","50_1","100_1", "500_1", "1000_1"), 
                          labels = c( "1:1","1:10","1:50","1:100",
                                      "1:500","1:1000","10:1","50:1","100:1",
                                      "500:1","1000:1"))) %>%
    group_by(marker_id,sample)%>%
    arrange(desc(count), .by_group = T) %>% 
    mutate(marker_id = as.factor(marker_id), 
    Haplotype = forcats::fct_reorder(Haplotype,Frequency)) %>%
    ggplot(aes(x = marker_id, y = Frequency, fill = Haplotype))+
    geom_bar(position = "stack", stat = "identity", colour ="black") + 
    ylab("Clone relative frequency (%)") + xlab(" ") + 
    scale_fill_manual(values = c(AR079, AR093,others)) +
    theme_bw(base_size = 10)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 7, angle=90, vjust=0.5), 
          legend.position ="none", 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 10, color = "black"),
          panel.border = element_rect(colour = "black", fill = NA)) +
    coord_cartesian(ylim=c(0, 100))+
    scale_y_continuous( breaks = c(0,50,100), labels = c("0" = "0", "50" = "50","100"= "100"))  +
    facet_wrap(facets = vars(infoD), ncol = 6)  

ggsave("figs/figure3.tiff", units="cm", width=18, height=6.6, dpi=300, compression = 'lzw')
  
dev.off()



