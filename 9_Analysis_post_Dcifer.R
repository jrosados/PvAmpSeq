library(tidyverse)
library(here)
library(cowplot)
library(ggpubr)
library(table1)
library(tableone)
library(ggplot2)




#Loading data from previous runs 

sols_data <- readRDS(here("intermediate_data", "sols_data_merged.rds"))
peru_data <- readRDS(here("intermediate_data", "peru_data_merged.rds"))

#Loading Dcifer 
sols_all_meta <- readRDS(here("intermediate_data", "sols_all_meta.rds"))
peru_all_meta <- readRDS(here("intermediate_data", "peru_all_meta.rds"))
sols_sig_meta <- readRDS(here("intermediate_data", "sols_sig_meta.rds")) 
peru_sig_meta <- readRDS(here("intermediate_data", "peru_sig_meta.rds"))


#Loading PNH classification data
sols_PNH_data <- readRDS(here("intermediate_data", "sols_data_classification.rds"))
peru_PNH_data <- readRDS(here("intermediate_data", "peru_data_classification.rds"))

#Loading amplification success data
sols_amp_success_data <- readRDS(here("intermediate_data", "sols_amp_success.rds"))
peru_amp_success_data <- readRDS(here("intermediate_data", "peru_amp_success.rds"))

#Merging Dcifer data and PNH data

sols_PNH_r_all<- 
  sols_all_meta %>%
  filter(comparison_type=="paired") %>%
  filter(sampleid2 %in% sols_data$sample[sols_data$follow.x=="baseline"]) %>%
  right_join(sols_PNH_data, by =c("sampleid1" ="sample"))


peru_PNH_r_all<- 
  peru_all_meta %>%
  filter(comparison_type=="paired") %>% #there are 
  #two comparisons from diff homes but with similar patient name. CAH vs LUP
  filter(sampleid2 %in% peru_data$sample_id[peru_data$day=="Day 0"]) %>%
  right_join(peru_PNH_data, by =c("sampleid1" ="sample"))

#Merging Dcifer data and PNH data

#########
#SI     #
#########

amp_success_SI_r_all<-
sols_PNH_r_all %>%
  left_join(sols_amp_success_data, by = c("sampleid1" ="sample")) %>%
  select(sampleid1,sig_est, amp_success)

#########
#Peru   #
#########

amp_success_PE_r_all<-
  peru_PNH_r_all %>%
  left_join(peru_amp_success_data, by = c("sampleid1" ="sample")) %>%
  select(sampleid1,sig_est, amp_success)


#Merging both SI and PE database
amp_success_sols_peru<- rbind(amp_success_SI_r_all,amp_success_PE_r_all)

#Table success amplification vs IBD pairwise significant comparison 

table_amp_success_sols_peru<- CreateTableOne(vars = c("amp_success"), 
                                      strata = "sig_est", 
                                      data = amp_success_sols_peru)

print(table_amp_success_sols_peru,showAllLevels = TRUE,
      nonnormal = "amp_success")

#######
#Plot #
#######

#Plot r Peru
peru_PNH_r_plot<-
peru_PNH_r_all %>% 
  #mutate(CI_range = CI_upper-CI_lower) %>% 
  ggplot(aes(x = estimate, y = Prop_het)) +#, color = factor(sig_est))) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("navajowhite3", "mediumorchid3")) +
  labs(x = "r^",
       y = bquote(1 - P [NH]),
       color = "") +
  ggtitle("Peru") +
  theme_minimal()+
  stat_cor(label.x = 0.10) #+
  #facet_wrap(~clones)

#Plot r total Peru
peru_PNH_rtotal_sig<- 
  peru_sig_meta %>%
  filter(comparison_type=="paired") %>% #there are 
  #two comparisons from diff homes but with similar patient name. CAH vs LUP
  filter(sampleid2 %in% peru_data$sample_id[peru_data$day=="Day 0"]) %>%
  right_join(peru_PNH_data, by =c("sampleid1" ="sample"))


peru_PNH_rtotal_plot<-
peru_PNH_rtotal_sig %>% 
  #mutate(CI_range = CI_upper-CI_lower) %>% 
  ggplot(aes(x = rtotal1, y = Prop_het #, color = factor(sig_est)
             )) +
  geom_point(alpha = 0.4) +
  geom_smooth(color = "indianred4", ) +
  stat_cor(label.y = 1.9)+
  labs(x = "r^total",
       y = bquote(1 - P [NH]),
       color = "M'") +
  ggtitle("Peru") +
  theme_minimal()



#Plot r Solomons
sols_PNH_r_plot<-
  sols_PNH_r_all %>% 
  #mutate(CI_range = CI_upper-CI_lower) %>% 
  ggplot(aes(x = estimate, y = Prop_het)) +# , color = factor(sig_est))) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("navajowhite3", "mediumorchid3")) +
  labs(x = "r^",
       y = bquote(1 - P [NH]),
       color = "") +
  ggtitle("Solomon Islands") +
  stat_cor(label.x = 0.10)+
  theme_minimal() #+
  #facet_wrap(~clones)

#Plot r total Solomons
sols_PNH_rtotal_sig<- 
  sols_sig_meta %>%
  filter(comparison_type=="paired") %>%
  filter(sampleid2 %in% sols_data$sample[sols_data$follow.x=="baseline"]) %>%
  right_join(sols_PNH_data, by =c("sampleid1" ="sample"))


sols_PNH_rtotal_plot<-
  sols_PNH_rtotal_sig %>% 
  #mutate(CI_range = CI_upper-CI_lower) %>% 
  ggplot(aes(x = rtotal1, y = Prop_het)) +
  geom_point(alpha = 0.4) +
  geom_smooth(color = "indianred4", ) +
  stat_cor(label.y = 1.9)+
    labs(x = "r^total",
       y = bquote(1 - P [NH]),
       color = "M'") +
  ggtitle("Solomon Islands") +
  theme_minimal() 


PNH_r_plot <-ggarrange(sols_PNH_r_plot, peru_PNH_r_plot,
                                      ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

PNH_rtotal_plot<-ggarrange(sols_PNH_rtotal_plot, peru_PNH_rtotal_plot,
                                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")


figure_ibd_new_hap <- 
ggdraw() +
  draw_plot(PNH_r_plot, x = 0.005,y = 0.5, width = .99, height = 0.49) +
  draw_plot(PNH_rtotal_plot, x = 0.005,y = 0, width = .99, height = 0.49) +
  draw_plot_label(label = c("A", "B"), size = 12,
                  x = c(0,0), y = c(1, 0.5))

print(figure_ibd_new_hap)

ggsave("figs/Supp_figure6.tiff", units="cm", width=25.14, height=20.12, dpi=300, compression = 'lzw')
  