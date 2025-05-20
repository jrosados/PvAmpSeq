library(tidyverse)
library(here)
library(cowplot)
library(ggpubr)
library(table1)
library(tableone)
library(ggplot2)




#Loading data from previous runs 

sols_data <- readRDS(here("intermediate_data", "sols_data_merged.rds"))
#write.csv(sols_data, "sols_data_raw.csv", row.names = F)

peru_data <- readRDS(here("intermediate_data", "peru_data_merged.rds"))
#write.csv(peru_data, "peru_data.csv", row.names = F)

#Loading Dcifer 
sols_all_meta <- readRDS(here("intermediate_data", "sols_all_meta.rds"))
peru_all_meta <- readRDS(here("intermediate_data", "peru_all_meta.rds"))
sols_sig_meta <- readRDS(here("intermediate_data", "sols_sig_meta.rds")) 
peru_sig_meta <- readRDS(here("intermediate_data", "peru_sig_meta.rds"))


#Loading PNH classification data
sols_PNH_data <- readRDS(here("intermediate_data", "sols_data_classification.rds"))
#write.csv(sols_PNH_data, "sols_PNH_data.csv", row.names = F)
peru_PNH_data <- readRDS(here("intermediate_data", "peru_data_classification.rds"))
#write.csv(peru_PNH_data, "peru_PNH_data.csv", row.names = F)

#Loading amplification success data
sols_amp_success_data <- readRDS(here("intermediate_data", "sols_amp_success.rds"))
peru_amp_success_data <- readRDS(here("intermediate_data", "peru_amp_success.rds"))

#Merging Dcifer data and PNH data

sols_PNH_r_all<- 
  sols_all_meta %>%
  filter(comparison_type=="paired") %>%
  filter(sampleid2 %in% sols_data$sample[sols_data$follow.x=="baseline"]) %>%
  right_join(sols_PNH_data, by =c("sampleid1" ="sample")) %>%
 mutate(IBD_classification= case_when(estimate >=0.5 ~ "Homologous",
                                      estimate <=0.25 ~"Heterologous", 
                                      estimate >0.25 |estimate <0.5  ~"Difficult to define")) %>%
  rename(IBD = estimate) %>%
  select(sampleid1,Patient, treatment, treatment2, recurrence, 
         delay_since_prev_ep, change_moi,max_moi,clones,IBD,IBD_classification)

write.csv(sols_PNH_r_all, "Text_S3_SolomonsIsland_IBD_data_06_02_25.csv", row.names = F)

peru_PNH_r_all<- 
  peru_all_meta %>%
  filter(comparison_type=="paired") %>% #there are 
  #two comparisons from diff homes but with similar patient name. CAH vs LUP
  filter(sampleid2 %in% peru_data$sample_id[peru_data$day=="Day 0"]) %>%
  right_join(peru_PNH_data, by =c("sampleid1" ="sample")) %>%
  mutate(IBD_classification= case_when(estimate >=0.5 ~ "Homologous",
                                       estimate <=0.25 ~"Heterologous", 
                                       estimate >0.25 |estimate <0.5  ~"Difficult to define")) %>%
  rename(IBD = estimate) %>%
  select(sampleid1,Patient, treatment, recurrence, 
         delay_since_prev_ep, change_moi,max_moi,clones,IBD,IBD_classification)

write.csv(peru_PNH_r_all, "Text_S4_Peru_IBD_data_06_02_25.csv", row.names = F)


#Classification vs treatment
table1(~ IBD_classification|treatment2, data =sols_PNH_r_all)


SI_table_treatment_classification<- CreateTableOne(vars = c("IBD_classification"), strata = "treatment2", data = sols_PNH_r_all, 
                                                factorVars =c("IBD_classification"))

print(SI_table_treatment_classification,showAllLevels = TRUE,
      nonnormal = "IBD_classification")

fisher.test(sols_PNH_r_all$treatment2, sols_PNH_r_all$IBD_classification)

#Clonality and classification in SI
table1(~ clones|IBD_classification, data =sols_PNH_r_all)

CreateTableOne(vars = c("clones"), strata = "IBD_classification", data = sols_PNH_r_all, 
               factorVars =c("clones"))

fisher.test(sols_PNH_r_all$clones, sols_PNH_r_all$IBD_classification)

#Clonality and classification in PE

table1(~ clones|IBD_classification, data =peru_PNH_r_all)

CreateTableOne(vars = c("clones"), strata = "IBD_classification", data = peru_PNH_r_all, 
               factorVars =c("clones"))

fisher.test(peru_PNH_r_all$clones, peru_PNH_r_all$IBD_classification)



