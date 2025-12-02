
# Load pkgs ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(janitor)
library(ggpubr)
library(reshape2)

# Load data ---------------------------------------------------------------

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
                     day = days_since_last_episode) %>% 
              distinct(),
            by = c("sampleid1" = "sample_id")) %>% 
  rename("subject_id1" = "subject_id",
         "day_id1" = "day") %>% 
  # do the same for sample id 2 in pair - merge with analysis_data to get patient_name and 'day' info
  left_join(analysis_data %>% 
              select(sample_id, subject_id, 
                     day = days_since_last_episode) %>% 
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

## IBD distribution
sols_ibd %>% 
  filter(comparison_type == "paired") %>% 
  ggplot(aes(x = estimate)) +
  geom_histogram(bins = 20) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "IBD",
       y = "Frequency",
       title = "Solomon Islands") +
  theme_bw(base_size = 14) +

peru_ibd %>% 
  filter(comparison_type == "paired") %>% 
  ggplot(aes(x = estimate)) +
  geom_histogram(bins = 20) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "IBD",
       y = "Frequency",
       title = "Peru") +
  theme_bw(base_size = 14) 

## Days since baseline and IBD
sols_ibd %>% 
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  ggplot(aes(x = day_id1, y = reorder(patientid1, estimate), color = ibd_range)) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Days since baseline infection",
         y = "Patient ID",
         color = "IBD estimate",
         title = "Solomon Islands") +
    theme_bw(base_size = 13) +

peru_ibd %>% 
  mutate(ibd_range = cut(estimate,  
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         include.lowest = TRUE, 
                         labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))) %>% 
  ggplot(aes(x = day_id1, y = reorder(patientid1, estimate), color = ibd_range)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Days since first infection during study period",
       y = "Patient ID",
       color = "IBD estimate",
       title = "Peru") +
  theme_bw(base_size = 13)  


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
