
# Load pkgs ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(janitor)
library(ggpubr)
library(reshape2)
library(cowplot)
library(table1)
library(tableone)


# Load data ---------------------------------------------------------------
load("intermediate_data/intermediate_data.RData")
sols_all_meta <- readRDS(here("intermediate_data", "sols_all_meta.rds"))
peru_all_meta <- readRDS(here("intermediate_data", "peru_all_meta.rds"))
sols_sig_meta <- readRDS(here("intermediate_data", "sols_sig_meta.rds")) 
peru_sig_meta <- readRDS(here("intermediate_data", "peru_sig_meta.rds"))

# Wrangle data ------------------------------------------------------------

sols_ibd <- sols_all_meta %>% 
  # filter to include only paired samples, and only pw comparisons to baseline
  filter(comparison_type == "paired", day_id2 == 0) %>% 
  mutate(subject_id = paste0("AR-", patientid1),
         day_id1 = as.integer(day_id1),
         day_id2 = as.integer(day_id2)) %>% 
  group_by(subject_id) %>%  
  arrange(day_id1) %>%  
  # create episode_number accounting for baseline as episode 1
  mutate(episode_number = row_number() + 1) %>% 
  ungroup()

peru_ibd <- peru_all_meta %>% 
  # filter to include only paired samples
  filter(comparison_type == "paired") %>% 
  # merge with analysis_data to get patient_name and 'day' info
  left_join(analysis_data %>% 
              select(sample_id, subject_id, 
                     day = days_since_enrolment) %>% 
              distinct(),
            by = c("sampleid1" = "sample_id")) %>% 
  rename("subject_id1" = "subject_id",
         "day_id1" = "day") %>% 
  # do the same for sample id 2 in pair - merge with analysis_data to get patient_name and 'day' info
  left_join(analysis_data %>% 
              select(sample_id, subject_id, 
                     day = days_since_enrolment) %>% 
              distinct(),
            by = c("sampleid2" = "sample_id")) %>% 
  rename("day_id2" = "day") %>% 
  select(-subject_id1) %>% 
  # filter to include only pw comparisons to first episode in study period
  filter(day_id2 == "0") %>% 
  group_by(subject_id) %>%  
  arrange(date1) %>%  
  # create episode_number accounting for baseline as episode 1
  mutate(episode_number = row_number() + 1) %>% 
  ungroup()

# Plots -------------------------------------------------------------------


# peru_ibd_modified <- peru_ibd %>%
#   mutate(ibd_range = cut(estimate,
#                          breaks = c(0, 0.25, 0.5, 0.75, 1),
#                          include.lowest = TRUE,
#                          labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")))

SI_all_samples_day_recurrence <- Epi_data_SI %>%
  mutate(
    SampleID = paste0("AR-", str_pad(studyid, width = 3, side = "left", pad = "0"), "-", vis_visit_day),
  studyid= as.character(studyid)) %>%
    #mutate(SampleName = paste0("AR-",studyid,"-",vis_visit_day))
  inner_join(AmpSeq_samples_data_SI,
             by =c("studyid" = "SampleName",
                    "SampleID")) %>%
  select(vis_date, studyid,SampleID) %>%
  mutate(date = as.Date(vis_date, format = "%m/%d/%Y")) %>%
  arrange(studyid, SampleID, date) %>% 
  distinct() %>% 
  # Calculate delay since previous infection, within each ID
  group_by(studyid) %>% 
  mutate(delay_since_prev_ep = as.numeric( date - first(date), 
                                           units = "days", na.rm=T),
         recurrence =rank(date)) %>%
  select(studyid, SampleID, delay_since_prev_ep,recurrence) %>%
  #filter( !delay_since_prev_ep ==0) %>%
  filter(recurrence >1) %>%
  select(-recurrence) %>%
  rename(patientid1 = studyid, 
         sampleid1 = SampleID,
         day_id1 = delay_since_prev_ep) %>%
  left_join(sols_ibd %>% select(sampleid1,estimate, patientid1),
            by=c("sampleid1", "patientid1")) 



PE_all_samples_day_recurrence <- Epi_data_PE %>% 
  filter(cod_samp %in% AmpSeq_samples_data_PE$SampleID)  %>%
  select(cod_per, cod_samp, date) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  arrange(cod_per, cod_samp, date) %>% 
  distinct() %>% 
  # Calculate delay since previous infection, within each ID
  group_by(cod_per) %>% 
  mutate(delay_since_prev_ep = as.numeric( date - first(date), 
                                           units = "days", na.rm=T),
           recurrence =rank(date)) %>%
  select(cod_per, cod_samp, delay_since_prev_ep,recurrence) %>%
  #filter( !delay_since_prev_ep ==0) %>%
  filter(recurrence >1) %>%
  select(-recurrence) %>%
  rename(patientid1 = cod_per, 
         sampleid1 = cod_samp,
         day_id1 = delay_since_prev_ep) %>%
  left_join(peru_ibd %>% select(sampleid1,estimate, patientid1, day_id1)) 
  

# Define your desired colors
my_colors <- c("0-0.25" = "#66C2A5",  # Color 1 (from Set2 palette)
               "0.25-0.5" = "#FC8D62", # Color 2
               "0.5-0.75" = "#8DA0CB", # Color 3
               "0.75-1" = "#E78AC3",   # Color 4
               "ND" = "gray50")        # Gray for Not Detected


## IBD distribution
sols_ibd_plot_distribution<- sols_ibd %>% 
  filter(comparison_type == "paired") %>% 
  ggplot(aes(x = estimate)) +
  geom_histogram(bins = 20) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "IBD estimate",
       y = "Frequency",
       title = "Solomon Islands") +
  theme_bw(base_size = 13) 

peru_ibd_plot_distribution<-
peru_ibd %>% 
  filter(comparison_type == "paired") %>% 
  ggplot(aes(x = estimate)) +
  geom_histogram(bins = 20) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "IBD estimate",
       y = "Frequency",
       title = "Peru") +
  theme_bw(base_size = 13) 

## Days since baseline and IBD
sols_ibd_days_plot<- sols_ibd %>% 
  #SI_all_samples_day_recurrence %>%
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  #mutate(ibd_range = fct_explicit_na(ibd_range, na_level = "ND"))  %>%
  ggplot(aes(x = day_id1, y = reorder(patientid1, estimate), color = ibd_range)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  #scale_color_manual(values = my_colors) + # Apply custom colors
      labs(x = "Days since baseline infection",
         y = "Patient ID",
         color = "IBD estimate",
         title = "") +
    theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(face = "plain", color="black", size = 5)) 

table_time_ibd_SI<- CreateTableOne(vars = c("day_id1"), strata = "ibd_range", data = sols_ibd_days_plot[["data"]], 
                                   factorVars =c("ibd_range"))

print(table_time_ibd_SI,showAllLevels = TRUE,
      nonnormal = "day_id1")

peru_ibd_days_plot<- peru_ibd %>% 
  #PE_all_samples_day_recurrence %>%
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>%
  #mutate(ibd_range = fct_explicit_na(ibd_range, na_level = "ND"))  %>%
  ggplot(aes(x = day_id1, y = reorder(patientid1, estimate), color = ibd_range)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  #scale_color_manual(values = my_colors) + # Apply custom colors
  labs(x = "Days since first infection during study period",
       y = "Patient ID",
       color = "IBD estimate",
       title = "") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(face = "plain", color="black", size = 5)) 
  
table1(~ day_id1|ibd_range, data =peru_ibd_days_plot[["data"]])

table_time_ibd_pe<- CreateTableOne(vars = c("day_id1"), strata = "ibd_range", data = peru_ibd_days_plot[["data"]], 
                                              factorVars =c("ibd_range"))

print(table_time_ibd_pe,showAllLevels = TRUE,
      nonnormal = "day_id1")

sols_ibd_boxplot_days_plot<- sols_ibd %>% 
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  ggplot() +
  geom_boxplot(aes(x = ibd_range, y = day_id1,color = ibd_range),
               width = 0.5) +
  stat_compare_means(aes(x = ibd_range, y = day_id1), 
                     label = "p.format", size = 3,
                     label.x = 3.5, label.y = 150)+
  scale_color_brewer(palette = "Set2") +
  labs(y = "Days since baseline infection",
       x = "IBD estimate",
       color = "IBD estimate",
       title = "") +
  theme_bw(base_size = 14) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip() 
  
peru_ibd_boxplot_days_plot<- peru_ibd %>% 
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  ggplot() +
  geom_boxplot(aes(x = ibd_range, y = day_id1,color = ibd_range),
               width = 0.5) +
  stat_compare_means(aes(x = ibd_range, y = day_id1), 
                     label = "p.format", size = 3,
                     label.x = 4, label.y = 250)+
  scale_color_brewer(palette = "Set2") +
  labs(y = "Days since first infection during study period",
       x = "IBD estimate",
       color = "IBD estimate",
       title = "") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip() 

ibd_plot_distribution<-ggarrange(sols_ibd_plot_distribution, peru_ibd_plot_distribution,
                                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

ibd_days_plot<-ggarrange(sols_ibd_days_plot, peru_ibd_days_plot,
                                 ncol=2, nrow=1, common.legend = TRUE, legend="right")

ibd_boxplot_days_plot<-ggarrange(sols_ibd_boxplot_days_plot, peru_ibd_boxplot_days_plot,
                         ncol=2, nrow=1, common.legend = TRUE, legend="right")
figure4 <- ggdraw() +
  draw_plot(ibd_plot_distribution, x = 0.005,y = 0.75, width = .99, height = 0.245) +
  draw_plot(ibd_days_plot, x = 0.005,y = 0.25, width = .99, height = 0.5) +
  draw_plot(ibd_boxplot_days_plot, x = 0.005,y = 0, width = .99, height = 0.25) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 14,
                  x = c(0,0.5,0,0.5,0,0.5), y = c(1,1,0.75,0.75,0.25,0.25))

print(figure4)

#ggsave("figs/figure4_IBD_v2.pdf", units="cm", width=20.14, height=25.12, dpi=300)
ggsave("figs/figure4_IBD_v2.png", units="cm", width=20.14, height=25.12, dpi=300)

# Supplementary plots


Supp_sols_ibd_days_plot<- #sols_ibd %>% 
  SI_all_samples_day_recurrence %>%
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  mutate(ibd_range = fct_explicit_na(ibd_range, na_level = "ND"))  %>%
  ggplot(aes(x = day_id1, y = reorder(patientid1, estimate), color = ibd_range)) +
  geom_point(size = 1) +
  #scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = my_colors) + # Apply custom colors
  labs(x = "Days since baseline infection",
       y = "Patient ID",
       color = "IBD estimate",
       title = "") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(face = "plain", color="black", size = 5)) 


Supp_peru_ibd_days_plot<- #peru_ibd %>% 
  PE_all_samples_day_recurrence %>%
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>%
  mutate(ibd_range = fct_explicit_na(ibd_range, na_level = "ND"))  %>%
  ggplot(aes(x = day_id1, y = reorder(patientid1, estimate), color = ibd_range)) +
  geom_point(size = 1) +
  #scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = my_colors) + # Apply custom colors
  labs(x = "Days since first infection during study period",
       y = "Patient ID",
       color = "IBD estimate",
       title = "") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(face = "plain", color="black", size = 5)) 

supp_ibd_days_plot<-ggarrange(Supp_sols_ibd_days_plot, Supp_peru_ibd_days_plot,
                         ncol=2, nrow=1, common.legend = TRUE, legend="right")

Supp_fig_ibd_days <- ggdraw() +
  draw_plot(supp_ibd_days_plot, x = 0.005,y = 0, width = .99, height = 0.95) +
  draw_plot_label(label = c("A", "B"), size = 14,
                  x = c(0,0.5), y = c(1,1))

print(Supp_fig_ibd_days)

#ggsave("figs/figure4_IBD_v2.pdf", units="cm", width=20.14, height=25.12, dpi=300)
ggsave("figs/Supp_fig_ibd_days.png", units="cm", width=20.14, height=25.12, dpi=300)


# Number of previous infections in Peru -----------------------------------

PE_all_samples <- Epi_data_PE %>% 
  mutate(PvPCR = if_else(PCR == 2 |PCR == 4, "Pv+","Pv-",
                         missing = "unknown")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  group_by(cod_per) %>% 
  summarise(
    PCR_until_Nov2014 = sum(PvPCR == "Pv+" & date <= as.Date("2014-11-30"), na.rm = TRUE),
    PCR_Jan2013_Dec2015 = sum(PvPCR == "Pv+" & 
                                date >= as.Date("2013-01-01") & 
                                date <= as.Date("2015-12-31"), na.rm = TRUE),
    
    Pv_clinic_until_Nov2014 = sum(clinic_Pvinf == "1" & date <= as.Date("2014-11-30"), na.rm = TRUE),
    Pv_clinic_Jan2013_Dec2015 = sum(clinic_Pvinf == "1" & date <= as.Date("2015-12-31"), na.rm = TRUE),
                                  
    Pv_treatment_until_Nov2014 = sum(medicamento2 == "P" |medicamento2 == "C" & date <= as.Date("2014-11-30"), na.rm = TRUE),
    Pv_treatment_Jan2013_Dec2015 = sum(medicamento2 == "P" |medicamento2 == "C" & date <= as.Date("2015-12-31"), na.rm = TRUE))

write.csv(PE_all_samples, "PE_all_samples_num_infec.csv", row.names = F)

PE_sample_subset <- Epi_data_PE %>% 
  filter(cod_per %in% AmpSeq_samples_data_PE$PatientName)  %>%
  mutate(PvPCR = if_else(PCR == 2 |PCR == 4, "Pv+","Pv-",
                         missing = "unknown")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  group_by(cod_per) %>% 
  summarise(
    PCR_until_Nov2014 = sum(PvPCR == "Pv+" & date <= as.Date("2014-11-30"), na.rm = TRUE),
    PCR_Jan2013_Dec2015 = sum(PvPCR == "Pv+" & 
                                date >= as.Date("2013-01-01") & 
                                date <= as.Date("2015-12-31"), na.rm = TRUE),
    
    Pv_clinic_until_Nov2014 = sum(clinic_Pvinf == "1" & date <= as.Date("2014-11-30"), na.rm = TRUE),
    Pv_clinic_Jan2013_Dec2015 = sum(clinic_Pvinf == "1" & date <= as.Date("2015-12-31"), na.rm = TRUE),
    
    Pv_treatment_until_Nov2014 = sum(medicamento2 == "P" |medicamento2 == "C" & date <= as.Date("2014-11-30"), na.rm = TRUE),
    Pv_treatment_Jan2013_Dec2015 = sum(medicamento2 == "P" |medicamento2 == "C" & date <= as.Date("2015-12-31"), na.rm = TRUE))


write.csv(PE_sample_subset, "PE_sample_subset_num_infec.csv", row.names = F)



# Just selected samples in this study
PE_all_samples <- Epi_data_PE %>% 
  filter(cod_samp %in% AmpSeq_samples_data_PE$SampleID)  %>%
  select(cod_per, cod_samp, date) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  arrange(cod_per, cod_samp, date) %>% 
  distinct() %>% 
  # Calculate delay since previous infection, within each ID
  group_by(cod_per) %>% 
  mutate(delay_since_prev_ep = as.numeric( date - first(date), 
                                           units = "days", na.rm=T),
         recurrence =rank(date)) %>%
  select(cod_per, cod_samp, delay_since_prev_ep,recurrence) %>%
  #filter( !delay_since_prev_ep ==0) %>%
  filter(recurrence >1) %>%
  select(-recurrence) %>%
  rename(patientid1 = cod_per, 
         sampleid1 = cod_samp,
         day_id1 = delay_since_prev_ep) %>%
  left_join(peru_ibd %>% select(sampleid1,estimate, patientid1, day_id1)) 


# Proportions -------------------------------------------------------------

sols_ibd %>% 
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  tabyl(ibd_range) %>% 
  adorn_totals("row")

peru_ibd %>% 
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  tabyl(ibd_range) %>% 
  adorn_totals("row")


# Pv3Rs results -----------------------------------------------------------

## Sols proportions 
load("outputs/Pv3Rs_sols_posteriors_20240827_181146.RData")

joint_summary <- indiv_posteriors_joint %>% 
  group_by(subject_id, state_pair, joint_probability) %>% 
  filter(episode_number == max(episode_number)) %>% 
  mutate(percentage = round(joint_probability*100, 3),
         total_recurrences = episode_number-1) %>%
  select(-episode_number) %>% 
  arrange(subject_id, total_recurrences, percentage) %>% 
  relocate(total_recurrences, .before = state_pair)

joint_top3_summary <- joint_summary %>% 
  select(-percentage) %>% 
  group_by(subject_id) %>% 
  # Get the highest prob for each subject and store ranking
  mutate(ranking = rank(-joint_probability, ties.method = "first")) %>% 
  arrange(subject_id, ranking) %>% 
  filter(ranking %in% c(1, 2, 3))

joint_recurrence_summary <- joint_top3_summary %>% 
  select(-total_recurrences) %>% 
  filter(ranking == 1) %>% 
  mutate(state_pair_split = strsplit(state_pair, "")) %>%  # Split state_pair into individual letters
  unnest(state_pair_split) %>%  # Unnest to create a new row for each letter
  group_by(subject_id) %>%
  mutate(
    episode_number = row_number() + 1,  # Create episode number starting from n+1 to reflect number of recurrence
    joint_probability = joint_probability,  # Keep the same joint_probability for all episodes\
    state_classification = case_when(state_pair_split == "L" ~ "Relapse",
                                     state_pair_split == "C" ~ "Recrudescence",
                                     state_pair_split == "I" ~ "Reinfection")) %>%
  select(subject_id, episode_number, state_classification, joint_probability)

classification_summary_sols <- indiv_posteriors_marginal %>% 
  rename("posterior_probability_recrudescence" = "Posterior_marginal_prob_C",
         "posterior_probability_relapse" = "Posterior_marginal_prob_L",
         "posterior_probability_reinfection" = "Posterior_marginal_prob_I") %>% 
  # merge with epi variables
  left_join(analysis_data %>% select(subject_id, sample_id, episode_number, treatment_arm, visit_date, days_since_enrolment, days_since_last_episode) %>% distinct(), 
            by = c("subject_id", "episode_number")) %>% 
  # merge with IBD r
  left_join(sols_ibd,
            by = c("subject_id", "episode_number")) %>% 
  # create bins for IBD r
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  # merge with joint probabilities
  left_join(joint_recurrence_summary %>%
              rename("joint_classification" = "state_classification"),
            by = c("subject_id", "episode_number"))

classification_summary_sols

save(classification_summary_sols, file = "outputs/sols_infection_classification.RData")

# PERU --------------------------------------------------------------------

load("outputs/Pv3Rs_peru_posteriors_20240904_010047.RData")

joint_summary <- indiv_posteriors_joint %>% 
  group_by(subject_id, state_pair, joint_probability) %>% 
  filter(episode_number == max(episode_number)) %>% 
  mutate(percentage = round(joint_probability*100, 3),
         total_recurrences = episode_number-1) %>%
  select(-episode_number) %>% 
  arrange(subject_id, total_recurrences, percentage) %>% 
  relocate(total_recurrences, .before = state_pair)

joint_top3_summary <- joint_summary %>% 
  select(-percentage) %>% 
  group_by(subject_id) %>% 
  # Get the highest prob for each subject and store ranking
  mutate(ranking = rank(-joint_probability, ties.method = "first")) %>% 
  arrange(subject_id, ranking) %>% 
  filter(ranking %in% c(1, 2, 3))

joint_recurrence_summary <- joint_top3_summary %>% 
  select(-total_recurrences) %>% 
  filter(ranking == 1) %>% 
  mutate(state_pair_split = strsplit(state_pair, "")) %>%  # Split state_pair into individual letters
  unnest(state_pair_split) %>%  # Unnest to create a new row for each letter
  group_by(subject_id) %>%
  mutate(
    episode_number = row_number() + 1,  # Create episode number starting from n+1 to reflect number of recurrence
    joint_probability = joint_probability,  # Keep the same joint_probability for all episodes\
    state_classification = case_when(state_pair_split == "L" ~ "Relapse",
                                     state_pair_split == "C" ~ "Recrudescence",
                                     state_pair_split == "I" ~ "Reinfection")) %>%
  select(subject_id, episode_number, state_classification, joint_probability)

classification_summary_peru <- indiv_posteriors_marginal %>% 
  rename("posterior_probability_recrudescence" = "Posterior_marginal_prob_C",
         "posterior_probability_relapse" = "Posterior_marginal_prob_L",
         "posterior_probability_reinfection" = "Posterior_marginal_prob_I") %>% 
  # merge with epi variables
  left_join(analysis_data %>% select(subject_id, sample_id, episode_number, community, visit_date, days_since_enrolment, days_since_last_episode) %>% distinct(), 
            by = c("subject_id", "episode_number")) %>% 
  # merge with IBD data
  left_join(peru_ibd,
            by = c("subject_id", "episode_number")) %>% 
  # create bins for IBD
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  # merge with joint probabilities
  left_join(joint_recurrence_summary %>%
              rename("joint_classification" = "state_classification"),
            by = c("subject_id", "episode_number"))

classification_summary_peru

save(classification_summary_peru, file = "outputs/peru_infection_classification.RData")

# Results classification --------------------------------------------------

## Sols
classification_summary_sols %>% 
  tabyl(episode_number, joint_classification) %>% 
  adorn_totals("col") %>% 
  adorn_percentages()

classification_summary_sols %>% 
  filter(episode_number == 2, joint_classification == "Recrudescence") %>% 
  summarise(min = min(posterior_probability_recrudescence),
            avg = mean(posterior_probability_recrudescence),
            max = mean(posterior_probability_recrudescence))

classification_summary_sols %>% 
  filter(episode_number == 2, joint_classification == "Reinfection") %>% 
  summarise(min = min(posterior_probability_reinfection),
            avg = mean(posterior_probability_reinfection),
            max = mean(posterior_probability_reinfection))

classification_summary_sols %>% 
  filter(episode_number == 2, joint_classification == "Relapse") %>% 
  summarise(min = min(posterior_probability_relapse),
            avg = mean(posterior_probability_relapse),
            max = mean(posterior_probability_relapse))

## Peru
classification_summary_peru %>% 
  mutate(ibd_classification = case_when(ibd_range == "0-0.25" ~ "Heterologous",
                                        ibd_range == "0.25-0.5" ~ "Difficult to define (0.25-0.5)",
                                        (ibd_range == "0.5-0.75" | ibd_range == "0.75-1") ~ "Homologous")) %>% 
  tabyl(ibd_classification) %>% 
  adorn_totals("col") %>% 
  adorn_percentages()

classification_summary_peru %>% 
  tabyl(episode_number, joint_classification) %>% 
  adorn_totals("col") %>% 
  adorn_percentages()

classification_summary_peru %>% 
  filter(episode_number == 2, joint_classification == "Reinfection") %>% 
  summarise(min = min(posterior_probability_reinfection),
            avg = mean(posterior_probability_reinfection),
            max = mean(posterior_probability_reinfection))

classification_summary_peru %>% 
  filter(episode_number == 2, joint_classification == "Relapse") %>% 
  summarise(min = min(posterior_probability_relapse),
            avg = mean(posterior_probability_relapse),
            max = mean(posterior_probability_relapse))

classification_summary_peru %>% 
  filter(episode_number >2, joint_classification == "Reinfection") %>% 
  summarise(min = min(posterior_probability_reinfection),
            avg = mean(posterior_probability_reinfection),
            max = mean(posterior_probability_reinfection))


# Treatment arm results ---------------------------------------------------

## Sols
classification_summary_sols %>% 
  mutate(ibd_classification = case_when(ibd_range == "0-0.25" ~ "Heterologous",
                                        ibd_range == "0.25-0.5" ~ "Difficult to define (0.25-0.5)",
                                        (ibd_range == "0.5-0.75" | ibd_range == "0.75-1") ~ "Homologous")) %>% 
  tabyl(ibd_classification, treatment_arm) %>% 
  adorn_totals("row") 

classification_summary_sols %>% 
  mutate(ibd_classification = case_when(ibd_range == "0-0.25" ~ "Heterologous",
                                         ibd_range == "0.25-0.5" ~ "Difficult to define (0.25-0.5)",
                                         (ibd_range == "0.5-0.75" | ibd_range == "0.75-1") ~ "Homologous"),
         tmt = case_when(treatment_arm == "AL" ~ "non-PQ",
                         TRUE ~ "PQ")) %>% 
  tabyl(ibd_classification, tmt) %>% 
  adorn_totals("row") %>% 
  adorn_percentages()

classification_summary_sols %>% 
  mutate(tmt = case_when(treatment_arm == "AL" ~ "non-PQ",
                         TRUE ~ "PQ")) %>% 
  tabyl(joint_classification, tmt) %>% 
  adorn_totals("row") %>% 
  adorn_percentages()



# Confusion matrix --------------------------------------------------------

## IBD classification Sols
confusion_data_ibdstate_sols <- 
  classification_summary_sols %>% 
  select(subject_id, sample_id, episode_number, visit_date, days_since_enrolment, days_since_last_episode, starts_with("posterior_"), starts_with("joint"), ibd_range) %>% 
  pivot_longer(cols = (starts_with("posterior")), 
               names_to = "marginal_classification", 
               values_to = "marginal_probability") %>% 
  mutate(marginal_classification = case_when(marginal_classification == "posterior_probability_recrudescence" ~ "Recrudescence",
                                             marginal_classification == "posterior_probability_relapse" ~ "Relapse",
                                             marginal_classification == "posterior_probability_reinfection" ~ "Reinfection"),
         marginal_classification = factor(marginal_classification, 
                                          levels = c("Relapse", "Recrudescence", "Reinfection")),
         marginal_probability_range = cut(marginal_probability, 
                                          breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                          include.lowest = TRUE, 
                                          labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")),
         joint_probability_range = cut(joint_probability, 
                                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                       include.lowest = TRUE, 
                                       labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")),
         ibd_classification = case_when(ibd_range == "0-0.25" ~ "Heterologous",
                                        ibd_range == "0.25-0.5" ~ "Difficult to define (0.25-0.5)",
                                        (ibd_range == "0.5-0.75" | ibd_range == "0.75-1") ~ "Homologous")) %>% 
  select(subject_id, episode_number, days_since_enrolment, ibd_range, ibd_classification, marginal_classification, marginal_probability, marginal_probability_range, joint_classification, joint_probability, joint_probability_range) %>% 
  # now only keep the highest 'ranking' marginal classification
  group_by(subject_id, episode_number) %>% 
  slice_max(order_by = marginal_probability, n = 1, with_ties = F) %>% 
  ungroup() %>% 
  group_by(ibd_classification, marginal_classification, joint_classification) %>%
  tally() %>%
  dcast(ibd_classification ~ marginal_classification, value.var = "n", fill = 0)

ggplot(melt(confusion_data_ibdstate_sols), aes(ibd_classification, variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = value)) +
  theme_minimal() +
  labs(title = "Solomon Islands",
       x = "IBD Classification",
       y = "Probabilistic classification",
       fill = "Count")

## IBD range Sols
confusion_data_ibdrange_sols <- classification_summary_sols %>% 
  select(subject_id, sample_id, episode_number, visit_date, days_since_enrolment, days_since_last_episode, starts_with("posterior_"), starts_with("joint"), ibd_range) %>% 
  pivot_longer(cols = (starts_with("posterior")), 
               names_to = "marginal_classification", 
               values_to = "marginal_probability") %>% 
  mutate(marginal_classification = case_when(marginal_classification == "posterior_probability_recrudescence" ~ "Recrudescence",
                                             marginal_classification == "posterior_probability_relapse" ~ "Relapse",
                                             marginal_classification == "posterior_probability_reinfection" ~ "Reinfection"),
         marginal_classification = factor(marginal_classification, 
                                          levels = c("Relapse", "Recrudescence", "Reinfection")),
         marginal_probability_range = cut(marginal_probability, 
                                          breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                          include.lowest = TRUE, 
                                          labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")),
         joint_probability_range = cut(joint_probability, 
                                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                       include.lowest = TRUE, 
                                       labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  select(subject_id, episode_number, days_since_enrolment, ibd_range, marginal_classification, marginal_probability, marginal_probability_range, joint_classification, joint_probability, joint_probability_range) %>% 
  # now only keep the highest 'ranking' marginal classification
  group_by(subject_id, episode_number) %>% 
  slice_max(order_by = marginal_probability, n = 1, with_ties = F) %>% 
  ungroup() %>% 
  group_by(ibd_range, marginal_classification, joint_classification) %>%
  tally() %>%
  dcast(ibd_range ~ marginal_classification, value.var = "n", fill = 0)


ggplot(melt(confusion_data_ibdrange_sols), aes(ibd_range, variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = value)) +
  theme_minimal() +
  labs(title = "Solomon Islands",
       x = "IBD range",
       y = "Probabilistic classification",
       fill = "Count")

## IBD classification Peru
confusion_data_ibdstate_peru <- 
  classification_summary_peru %>% 
  select(subject_id, sample_id, episode_number, visit_date, days_since_enrolment, days_since_last_episode, starts_with("posterior_"), starts_with("joint"), ibd_range) %>% 
  pivot_longer(cols = (starts_with("posterior")), 
               names_to = "marginal_classification", 
               values_to = "marginal_probability") %>% 
  mutate(marginal_classification = case_when(marginal_classification == "posterior_probability_recrudescence" ~ "Recrudescence",
                                             marginal_classification == "posterior_probability_relapse" ~ "Relapse",
                                             marginal_classification == "posterior_probability_reinfection" ~ "Reinfection"),
         marginal_classification = factor(marginal_classification, 
                                          levels = c("Relapse", "Recrudescence", "Reinfection")),
         marginal_probability_range = cut(marginal_probability, 
                                          breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                          include.lowest = TRUE, 
                                          labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")),
         joint_probability_range = cut(joint_probability, 
                                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                       include.lowest = TRUE, 
                                       labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")),
         ibd_classification = case_when(ibd_range == "0-0.25" ~ "Heterologous",
                                        ibd_range == "0.25-0.5" ~ "Difficult to classify (0.25-0.5)",
                                        (ibd_range == "0.5-0.75" | ibd_range == "0.75-1") ~ "Homologous")) %>% 
  select(subject_id, episode_number, days_since_enrolment, ibd_range, ibd_classification, marginal_classification, marginal_probability, marginal_probability_range, joint_classification, joint_probability, joint_probability_range) %>% 
  # now only keep the highest 'ranking' marginal classification
  group_by(subject_id, episode_number) %>% 
  slice_max(order_by = marginal_probability, n = 1, with_ties = F) %>% 
  ungroup() %>% 
  group_by(ibd_classification, marginal_classification, joint_classification) %>%
  tally() %>%
  dcast(ibd_classification ~ marginal_classification, value.var = "n", fill = 0)

ggplot(melt(confusion_data_ibdstate_peru), aes(ibd_classification, variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = value)) +
  theme_minimal() +
  labs(title = "Peru",
       x = "IBD Classification",
       y = "Probabilistic classification",
       fill = "Count")

## IBD range Peru
confusion_data_ibdrange_peru <- classification_summary_peru %>% 
  select(subject_id, sample_id, episode_number, visit_date, days_since_enrolment, days_since_last_episode, starts_with("posterior_"), starts_with("joint"), ibd_range) %>% 
  pivot_longer(cols = (starts_with("posterior")), 
               names_to = "marginal_classification", 
               values_to = "marginal_probability") %>% 
  mutate(marginal_classification = case_when(marginal_classification == "posterior_probability_recrudescence" ~ "Recrudescence",
                                             marginal_classification == "posterior_probability_relapse" ~ "Relapse",
                                             marginal_classification == "posterior_probability_reinfection" ~ "Reinfection"),
         marginal_classification = factor(marginal_classification, 
                                          levels = c("Relapse", "Recrudescence", "Reinfection")),
         marginal_probability_range = cut(marginal_probability, 
                                          breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                          include.lowest = TRUE, 
                                          labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")),
         joint_probability_range = cut(joint_probability, 
                                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                       include.lowest = TRUE, 
                                       labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  select(subject_id, episode_number, days_since_enrolment, ibd_range, marginal_classification, marginal_probability, marginal_probability_range, joint_classification, joint_probability, joint_probability_range) %>% 
  # now only keep the highest 'ranking' marginal classification
  group_by(subject_id, episode_number) %>% 
  slice_max(order_by = marginal_probability, n = 1, with_ties = F) %>% 
  ungroup() %>% 
  group_by(ibd_range, marginal_classification, joint_classification) %>%
  tally() %>%
  dcast(ibd_range ~ marginal_classification, value.var = "n", fill = 0)


ggplot(melt(confusion_data_ibdrange_peru), aes(ibd_range, variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = value)) +
  theme_minimal() +
  labs(title = "Peru",
       x = "IBD range",
       y = "Probabilistic classification",
       fill = "Count")


## final supp plot
ggplot(melt(confusion_data_ibdstate_sols), aes(ibd_classification, variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = value)) +
  theme_minimal() +
  labs(title = "Solomon Islands",
       x = "IBD Classification",
       y = "Probabilistic classification",
       fill = "Count") +
ggplot(melt(confusion_data_ibdstate_peru), aes(ibd_classification, variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = value)) +
  theme_minimal() +
  labs(title = "Peru",
       x = "IBD Classification",
       y = "Probabilistic classification",
       fill = "Count") +
  plot_annotation(tag_levels = "A", )

ggsave("outputs/plots/supp_confusion_matrix.tiff", width = 12, height = 5)

# Further analysis --------------------------------------------------------

classification_summary_sols %>% 
  ggplot(aes(x= joint_classification, y = days_since_last_episode, group = joint_classification, fill = joint_classification)) +
  geom_boxplot(alpha = 0.5) +
  stat_kruskal_test() +
  scale_fill_manual(values = c("Relapse" = "turquoise3",
                               "Recrudescence" = "skyblue4",
                               "Reinfection" = "magenta3")) +
  theme_minimal()


# For those with recrudescence treatment type?
# PQ vs no PQ – any trends? 

classification_summary_sols %>% 
  ggplot(aes(x= joint_classification, y = days_since_last_episode, group = joint_classification, fill = joint_classification)) +
  geom_boxplot(alpha = 0.5) +
  stat_kruskal_test() +
  scale_fill_manual(values = c("Relapse" = "turquoise3",
                               "Recrudescence" = "skyblue4",
                               "Reinfection" = "magenta3")) +
  facet_wrap(treatment_arm~.) +
  theme_minimal()

  
classification_summary_sols %>% 
  ggplot(aes(x = treatment_arm, y = days_since_last_episode, group = treatment_arm, fill = joint_classification)) +
  geom_boxplot(alpha = 0.5) +
  stat_kruskal_test() +
  scale_fill_manual(values = c("Relapse" = "turquoise3",
                               "Recrudescence" = "skyblue4",
                               "Reinfection" = "magenta3")) +
  facet_wrap(joint_classification~.) +
  theme_minimal()
