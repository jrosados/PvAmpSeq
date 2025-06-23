# Set global variables ----------------------------------------------------
set_global_vars <- function(
    RUN_EXAMPLE = TRUE,
    WITHIN_INDIVIDUAL_ALLELE_FREQ_THR = 0,
    BENCHMARK_MARKERS = c("Chr05","Chr07","Chr09","Chr10","Chr08","Chr13","Chr11","Chr03","Chr01","Chr02","Chr14"),
    MAX_MOI_TO_INCLUDE = 8,
    PRIOR_3RS = c(C = 1/3, L = 1/3, I = 1/3)
) {
  list(
    RUN_EXAMPLE = RUN_EXAMPLE,
    WITHIN_INDIVIDUAL_ALLELE_FREQ_THR = WITHIN_INDIVIDUAL_ALLELE_FREQ_THR,
    BENCHMARK_MARKERS = BENCHMARK_MARKERS,
    MAX_MOI_TO_INCLUDE = MAX_MOI_TO_INCLUDE,
    PRIOR_3RS = PRIOR_3RS
  )
}

# Haplotype frequencies ---------------------------------------------------

calculate_freqs <- function(
    analysis_data
){
  fs <- analysis_data %>% 
    # split data frame by marker
    group_by(marker_id) %>% 
    group_split() %>% 
    
    # Derive a within-marker list of frequencies, by individual
    lapply(function(x) {
      x <- x %>% 
        # Build a within-individual frequency table that 
        # always includes every haplotype (even the ones absent)
        mutate(haplotype = factor(haplotype)) %>% 
        select(sample_id, haplotype, frequency) %>% 
        pivot_wider(names_from = haplotype, 
                    values_from = frequency, 
                    values_fill = 0) %>% 
        pivot_longer(cols = -sample_id, 
                     names_to = "haplotype", 
                     values_to = "frequency") %>% 
        # Get population-level haplotype frequency, 
        # correcting for when within-individual sum is not equal 
        # to 1, as can happen when a minority clone is <2%
        group_by(sample_id) %>% 
        mutate(frequency = frequency / sum(frequency, na.rm = TRUE)) %>% 
        group_by(haplotype) %>% 
        summarise(frequency_pop_mean = mean(frequency, na.rm = TRUE))
      
      return(deframe(x))
    }) %>% 
    # get marker_id from group_keys from the group dfs 
    setNames(nm = analysis_data %>% group_by(marker_id) %>% group_keys() %>% pull(marker_id))
  
  # DATA FRAME: Here we can also save as dataframe for easier printing and table-ready for paper
  fs_df <-analysis_data %>%
    # Group by marker_id and sample_id for further calculations
    group_by(marker_id, sample_id) %>%

    # Build a within-individual frequency table that always includes every haplotype (even the ones absent)
    mutate(haplotype = factor(haplotype)) %>%
    select(marker_id, sample_id, haplotype, frequency) %>%
    pivot_wider(names_from = haplotype, values_from = frequency, values_fill = list(frequency = 0)) %>%
    pivot_longer(cols = -c(marker_id, sample_id), names_to = "haplotype", values_to = "frequency") %>%

    # Get population-level haplotype frequency, correcting for when within-individual sum is not equal to 1
    group_by(marker_id, sample_id) %>%
    mutate(frequency = frequency / sum(frequency, na.rm = TRUE)) %>%
    group_by(marker_id, haplotype) %>%
    summarise(frequency_pop_mean = mean(frequency, na.rm = TRUE), .groups = 'drop') %>%

    # Ensure haplotypes are correctly associated with their marker_id - note that this works for us because , in future would have to make this flexible to allow for haplotype names that are not reliant on having marker_id
    filter(str_detect(haplotype, marker_id))

  list(fs = fs,
       fs_df = fs_df)
}

# Run Pv3Rs ---------------------------------------------------------------

compute_all <- function(
    analysis_data,
    fs,
    global_vars,
    filename = "./outputs/Pv3Rs_posteriors_"
){

list2env(global_vars, envir = environment())

# Start timer
t_start <- Sys.time()

# Initialize an empty list to store the results
indiv_posteriors <- list()

# Loop through each unique individual name
for (indiv_name in unique(analysis_data$subject_id)) {
  # Verbose
  cat("ID : ", indiv_name, "...\n", sep = "")
  
  ## Prepare the data
  # 1- Subset haplotype data to specific individual and apply filters
  indiv_haplotype_data <- analysis_data %>% 
    # Restrict to a single patient
    filter(subject_id == indiv_name) %>% 
    # Ensure episodes are in order
    arrange(episode_number) %>%  
    # Restrict to a subset of markers
    filter(marker_id %in% BENCHMARK_MARKERS) %>% 
    # Restrict to summed MOI below threshold
    group_by(subject_id, episode_number, marker_id) %>% 
    mutate(MOI_per_marker = sum(n())) %>% 
    group_by(subject_id, episode_number) %>% 
    mutate(MOI_per_episode = max(MOI_per_marker, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(marker_id = factor(marker_id, levels = BENCHMARK_MARKERS))
  
  # 2- Calculate per-episode and per-participant MOI for PvR3S eligibility
  indiv_MOI <- indiv_haplotype_data %>% 
    select(subject_id, episode_number, 
           marker_id, starts_with("MOI_")) %>% 
    distinct() %>% 
    group_by(subject_id, episode_number) %>% 
    # Get highest per-marker MOI only for each episode
    # (drop marker_id in case of ties with highest per-marker MOI)
    select(-marker_id) %>% 
    distinct() %>% 
    filter(MOI_per_marker == max(MOI_per_marker, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(MOI_summed = sum(MOI_per_episode, na.rm = TRUE))
  
  ## Run Pv3Rs only if MOI below threshold         
  if (unique(indiv_MOI$MOI_summed) <= MAX_MOI_TO_INCLUDE) {
    # Verbose
    cat("Running Pv3Rs for: ", 
        indiv_name, 
        " (summed MOI = ", 
        unique(indiv_MOI$MOI_summed), 
        ").\n", 
        sep = "")
    
    # Preserve episode numbers for naming the output list
    indiv_episode <- unique(indiv_haplotype_data$episode_number)
    cat("The episode number is: ", indiv_episode) # printing to check episode order number is correct
    
    # Finish data preparation
    indiv_haplotype_data <- indiv_haplotype_data %>% 
      group_by(episode_number) %>% 
      group_split() %>% 
      lapply(function(x) {
        res <- x %>% 
          select(sample_id, episode_number, 
                 marker_id, haplotype, frequency) %>% 
          # For sensitivity analysis, allow to include/drop allele 
          # based on their within-individual frequency
          filter(frequency >= WITHIN_INDIVIDUAL_ALLELE_FREQ_THR) %>% 
          select(-sample_id, -episode_number, -frequency) %>% 
          distinct() %>% 
          # Prevent dropping of markers that are not characterized 
          # by setting .drop to FALSE
          group_by(marker_id, .drop = FALSE) %>% 
          group_split() %>% 
          lapply(function(y) {
            unique(y$haplotype)
          })
        
        # Returned a list named with each episode, 
        # setting marker allele to NA in case none are observed
        return(lapply(setNames(res, BENCHMARK_MARKERS), 
                      function(y) {
                        if (length(y) == 0) return(NA) else return(y)
                      }))
      })
    
    # Run Aimee's posterior estimation
    indiv_posterior <- compute_posterior(y     = indiv_haplotype_data,
                                         fs    = fs[BENCHMARK_MARKERS], 
                                         prior = matrix(PRIOR_3RS, 
                                                        nrow     = length(indiv_haplotype_data), 
                                                        ncol     = length(PRIOR_3RS), 
                                                        byrow    = TRUE, 
                                                        dimnames = list(c(1:length(indiv_haplotype_data)), 
                                                                        names(PRIOR_3RS))))
    
  } else {
    # Verbose
    cat("NOT Running Pv3Rs for: ", 
        indiv_name, 
        " because summed MOI exceeds threshold (observed = ", 
        unique(indiv_MOI$MOI_summed), 
        ", MAX_MOI_TO_INCLUDE = ", 
        MAX_MOI_TO_INCLUDE, 
        ").\n", 
        sep = "")
    
    # Return NULL
    indiv_episode    <- unique(indiv_haplotype_data$episode_number)
    indiv_posterior <- NULL
  }
  
  # Append the results to the list
  indiv_posteriors[[indiv_name]] <- list("subject_id" = indiv_name, 
                                         "episode_number"       = indiv_episode, 
                                         "Pv3Rs"       = indiv_posterior)
}

# End timer
t_end <- Sys.time()
cat("Pv3Rs for the whole dataset took : ", as.numeric(difftime(time1 = t_end, 
                                                               time2 = t_start, 
                                                               units = "secs"))/60, " mins", "\n", sep = "")

# Present the marginal data in a clearer format
indiv_posteriors_marginal <- do.call(rbind, 
                                     lapply(indiv_posteriors, function(x) {
                                       if (!is.null(x[["Pv3Rs"]])) {
                                         return(data.frame("subject_id"               = x[["subject_id"]], 
                                                           "episode_number"                     = x[["episode_number"]][-1], 
                                                           "Posterior_marginal_prob_C" = x[["Pv3Rs"]]$marg[, "C"], 
                                                           "Posterior_marginal_prob_L" = x[["Pv3Rs"]]$marg[, "L"], 
                                                           "Posterior_marginal_prob_I" = x[["Pv3Rs"]]$marg[, "I"]))
                                       } else {
                                         return(NULL)
                                       }
                                       
                                     }))
row.names(indiv_posteriors_marginal) <- 1:nrow(indiv_posteriors_marginal)

# Present the joint posterior estimates in a clearer format
indiv_posteriors_joint <- do.call(rbind, 
                                  lapply(indiv_posteriors, function(x) {
                                    if (!is.null(x[["Pv3Rs"]])) {
                                      joint_probs <- x[["Pv3Rs"]]$joint
                                      # Extract the state pairs and probabilities
                                      state_pairs <- names(joint_probs)
                                      prob_values <- as.numeric(joint_probs)
                                      
                                      # Create a data frame with subject_id, episode_number, and joint probabilities
                                      return(data.frame("subject_id"      = rep(x[["subject_id"]], length(state_pairs)), 
                                                        "episode_number"  = rep(x[["episode_number"]][-1], each = length(state_pairs)),
                                                        "state_pair"      = state_pairs, 
                                                        "joint_probability" = prob_values))
                                    } else {
                                      return(NULL)
                                    }
                                  }))
row.names(indiv_posteriors_joint) <- 1:nrow(indiv_posteriors_joint)

# Save list for use within 
output_list <- list(
  RUN_EXAMPLE                      = RUN_EXAMPLE,
  WITHIN_FREQ_THR                  = WITHIN_INDIVIDUAL_ALLELE_FREQ_THR,
  BENCHMARK_MARKERS                = BENCHMARK_MARKERS,
  MAX_MOI_TO_INCLUDE               = MAX_MOI_TO_INCLUDE,
  PRIOR_3RS                        = PRIOR_3RS,
  analysis_data                    = analysis_data,
  fs                               = fs,
  indiv_posteriors                 = indiv_posteriors,
  indiv_posteriors_marginal        = indiv_posteriors_marginal,
  indiv_posteriors_joint           = indiv_posteriors_joint
  )
# Also save Pv3Rs output because it's time consuming and 
# we don't want to re-run it every time.
save(output_list,
     file = paste0(filename,
                   strftime(Sys.time(), format = "%Y%m%d_%H%M%S"),
                   ".RData"))

return(output_list) 
}


# Run Pv3Rs null distribution iterations ----------------------------------

compute_null <- function(
    analysis_data,
    fs,
    global_vars
){
  
  list2env(global_vars, envir = environment())
  
  # Start timer
  t_start <- Sys.time()
  
  # Initialize an empty list to store the results
  indiv_posteriors <- list()
  
  # Loop through each unique individual name
  for (indiv_name in unique(analysis_data$subject_id)) {
    # Verbose
    cat("ID : ", indiv_name, "...\n", sep = "")
    
    ## Prepare the data
    # 1- Subset haplotype data to specific individual and apply filters
    indiv_haplotype_data <- analysis_data %>% 
      # Restrict to a single patient
      filter(subject_id == indiv_name) %>% 
      # Ensure episodes are in order
      arrange(episode_number) %>%  
      # Restrict to a subset of markers
      filter(marker_id %in% BENCHMARK_MARKERS) %>% 
      # Restrict to summed MOI below threshold
      group_by(subject_id, episode_number, marker_id) %>% 
      mutate(MOI_per_marker = sum(n())) %>% 
      group_by(subject_id, episode_number) %>% 
      mutate(MOI_per_episode = max(MOI_per_marker, na.rm = TRUE)) %>% 
      ungroup() %>% 
      mutate(marker_id = factor(marker_id, levels = BENCHMARK_MARKERS))
    
    # 2- Calculate per-episode and per-participant MOI for PvR3S eligibility
    indiv_MOI <- indiv_haplotype_data %>% 
      select(subject_id, episode_number, 
             marker_id, starts_with("MOI_")) %>% 
      distinct() %>% 
      group_by(subject_id, episode_number) %>% 
      # Get highest per-marker MOI only for each episode
      # (drop marker_id in case of ties with highest per-marker MOI)
      select(-marker_id) %>% 
      distinct() %>% 
      filter(MOI_per_marker == max(MOI_per_marker, na.rm = TRUE)) %>% 
      ungroup() %>% 
      mutate(MOI_summed = sum(MOI_per_episode, na.rm = TRUE))
    
    ## Run Pv3Rs only if MOI below threshold         
    if (unique(indiv_MOI$MOI_summed) <= MAX_MOI_TO_INCLUDE) {
      # Verbose
      cat("Running Pv3Rs for: ", 
          indiv_name, 
          " (summed MOI = ", 
          unique(indiv_MOI$MOI_summed), 
          ").\n", 
          sep = "")
      
      # Preserve episode numbers for naming the output list
      indiv_episode <- unique(indiv_haplotype_data$episode_number)
      cat("The episode number is: ", indiv_episode) # printing to check episode order number is correct
      
      # Finish data preparation
      indiv_haplotype_data <- indiv_haplotype_data %>% 
        group_by(episode_number) %>% 
        group_split() %>% 
        lapply(function(x) {
          res <- x %>% 
            select(sample_id, episode_number, 
                   marker_id, haplotype, frequency) %>% 
            # For sensitivity analysis, allow to include/drop allele 
            # based on their within-individual frequency
            filter(frequency >= WITHIN_INDIVIDUAL_ALLELE_FREQ_THR) %>% 
            select(-sample_id, -episode_number, -frequency) %>% 
            distinct() %>% 
            # Prevent dropping of markers that are not characterized 
            # by setting .drop to FALSE
            group_by(marker_id, .drop = FALSE) %>% 
            group_split() %>% 
            lapply(function(y) {
              unique(y$haplotype)
            })
          
          # Returned a list named with each episode, 
          # setting marker allele to NA in case none are observed
          return(lapply(setNames(res, BENCHMARK_MARKERS), 
                        function(y) {
                          if (length(y) == 0) return(NA) else return(y)
                        }))
        })
      
      # Run Aimee's posterior estimation
      indiv_posterior <- compute_posterior(y     = indiv_haplotype_data,
                                           fs    = fs[BENCHMARK_MARKERS], 
                                           prior = matrix(PRIOR_3RS, 
                                                          nrow     = length(indiv_haplotype_data), 
                                                          ncol     = length(PRIOR_3RS), 
                                                          byrow    = TRUE, 
                                                          dimnames = list(c(1:length(indiv_haplotype_data)), 
                                                                          names(PRIOR_3RS))))
      
    } else {
      # Verbose
      cat("NOT Running Pv3Rs for: ", 
          indiv_name, 
          " because summed MOI exceeds threshold (observed = ", 
          unique(indiv_MOI$MOI_summed), 
          ", MAX_MOI_TO_INCLUDE = ", 
          MAX_MOI_TO_INCLUDE, 
          ").\n", 
          sep = "")
      
      # Return NULL
      indiv_episode    <- unique(indiv_haplotype_data$episode_number)
      indiv_posterior <- NULL
    }
    
    # Append the results to the list
    indiv_posteriors[[indiv_name]] <- list("subject_id" = indiv_name, 
                                           "episode_number"       = indiv_episode, 
                                           "Pv3Rs"       = indiv_posterior)
  }
  
  # End timer
  t_end <- Sys.time()
  cat("Pv3Rs for the whole dataset took : ", as.numeric(difftime(time1 = t_end, 
                                                                 time2 = t_start, 
                                                                 units = "secs"))/60, " mins", "\n", sep = "")
  
  # Present the marginal data in a clearer format
  indiv_posteriors_marginal <- do.call(rbind, 
                                       lapply(indiv_posteriors, function(x) {
                                         if (!is.null(x[["Pv3Rs"]])) {
                                           return(data.frame("subject_id"               = x[["subject_id"]], 
                                                             "episode_number"                     = x[["episode_number"]][-1], 
                                                             "Posterior_marginal_prob_C" = x[["Pv3Rs"]]$marg[, "C"], 
                                                             "Posterior_marginal_prob_L" = x[["Pv3Rs"]]$marg[, "L"], 
                                                             "Posterior_marginal_prob_I" = x[["Pv3Rs"]]$marg[, "I"]))
                                         } else {
                                           return(NULL)
                                         }
                                         
                                       }))
  row.names(indiv_posteriors_marginal) <- 1:nrow(indiv_posteriors_marginal)
  
  # Present the joint posterior estimates in a clearer format
  indiv_posteriors_joint <- do.call(rbind, 
                                    lapply(indiv_posteriors, function(x) {
                                      if (!is.null(x[["Pv3Rs"]])) {
                                        joint_probs <- x[["Pv3Rs"]]$joint
                                        # Extract the state pairs and probabilities
                                        state_pairs <- names(joint_probs)
                                        prob_values <- as.numeric(joint_probs)
                                        
                                        # Create a data frame with subject_id, episode_number, and joint probabilities
                                        return(data.frame("subject_id"      = rep(x[["subject_id"]], length(state_pairs)), 
                                                          "episode_number"  = rep(x[["episode_number"]][-1], each = length(state_pairs)),
                                                          "state_pair"      = state_pairs, 
                                                          "joint_probability" = prob_values))
                                      } else {
                                        return(NULL)
                                      }
                                    }))
  row.names(indiv_posteriors_joint) <- 1:nrow(indiv_posteriors_joint)
  
  # Save list for use within 
  output_list <- list(
    RUN_EXAMPLE                      = RUN_EXAMPLE,
    WITHIN_FREQ_THR                  = WITHIN_INDIVIDUAL_ALLELE_FREQ_THR,
    BENCHMARK_MARKERS                = BENCHMARK_MARKERS,
    MAX_MOI_TO_INCLUDE               = MAX_MOI_TO_INCLUDE,
    PRIOR_3RS                        = PRIOR_3RS,
    analysis_data                    = analysis_data,
    # fs                               = fs,
    indiv_posteriors                 = indiv_posteriors,
    indiv_posteriors_marginal        = indiv_posteriors_marginal,
    indiv_posteriors_joint           = indiv_posteriors_joint
  )
  # Also save Pv3Rs output because it's time consuming and 
  # we don't want to re-run it every time.
  # save(output_list,
  #      file = paste0(filename,
  #                    strftime(Sys.time(), format = "%Y%m%d_%H%M%S"),
  #                    ".RData"))
  
  return(output_list) 
}

# Get/wrangle data --------------------------------------------------------

get_analysis_data_sols <- function(random_pairs,baseline_haps)
{
  # merge randomized pairs back with full dataset
  sols_recurrent_null <- random_pairs %>% 
    left_join(sols %>% 
                clean_names() %>% 
                mutate(vis_date = mdy(vis_date)), 
              by = c("sample", "info", "info_d", "episodes", "trtgrp", "trt_label", "vis_date")) 
  
  # get baseline haps
  baseline_haps <- sols_recurrent_null %>% 
    filter(episode_type == "baseline") %>% 
    pull(haplotype)
  
  # Data on episode number and time since last episode for new shuffled pairs
  # Note that the enrolment period was long for Solomon Islands, so we want to make sure we impose a time-order constraint for randomized pairs.
  episode_summary_null <- sols_recurrent_null %>%
    distinct(shuffled_pair_id, info_d, vis_date, sample) %>%
    arrange(shuffled_pair_id, vis_date) %>%
    group_by(shuffled_pair_id) %>%
    mutate(
      # get the episode number
      episode_number = row_number(),
      # format date
      first_date = first(vis_date),
      second_date = last(vis_date),
      days_since_last_episode = as.integer(last(vis_date) - first(vis_date))
    ) %>%
    ungroup()
  
  # get final data for inference
  analysis_data <- sols_recurrent_null %>%
    select(subject_id = shuffled_pair_id, # using "new ID" instead of original participant ID
           sample_id = sample, 
           treatment_arm = trt_label, 
           days_since_treatment = info_d,
           timepoint = follow_x,
           visit_date = vis_date,
           age_years = age_y, 
           sex = gender, 
           # episodes, # this is not accurate anymore
           marker_id, 
           haplotype, 
           type,
           frequency,
           mean_moi,
           max_moi) %>% 
    # add extra epi info on episode number and time since last episode
    left_join(episode_summary_null %>% select(subject_id = shuffled_pair_id, 
                                              sample_id = sample, 
                                              episode_number,
                                              days_since_last_episode),
              by = c("subject_id", "sample_id"))  %>% 
    # keep only haplotypes present at baseline
    filter(haplotype %in% baseline_haps)
  
  return(analysis_data)
}

get_analysis_data_peru <- function(random_pairs, all_haps)
{
  # merge randomized pairs back with full datset
  peru_recurrent_null <- random_pairs %>% 
    left_join(peru %>% 
                clean_names() %>% 
                mutate(date = mdy(date)), 
              by = c("sample", "info" = "patient_name", "date", "episodes")) 
  
  # get all haps
  all_haps <- peru_recurrent_null %>% 
    pull(haplotype)
  
  ### Data on episode number and time since last episode for new shuffled pairs
  # Note that in the Peru study this is not necessarily a person's "first" episode, rather the first episode in the evaluated study period
  episode_summary_null <- peru_recurrent_null %>%
    distinct(shuffled_pair_id, day, date, sample) %>%
    arrange(shuffled_pair_id, date) %>%
    group_by(shuffled_pair_id) %>%
    mutate(
      # get the episode number
      episode_number = row_number(),
      first_date = first(date),
      second_date = last(date),
      days_since_last_episode = as.integer(last(date) - first(date))
    ) %>%
    ungroup()
  
  # get final data for inference
  analysis_data <- peru_recurrent_null %>%
    select(subject_id = shuffled_pair_id, # using "new ID" instead of original participant ID
           sample_id = sample, 
           visit_date = date,
           age_years = age, 
           sex, 
           # episodes, # this is not accurate anymore
           marker_id, 
           haplotype, 
           frequency,
           mean_moi,
           max_moi) %>% 
    # add extra epi info on episode number and time since last episode
    left_join(episode_summary_null %>% select(subject_id = shuffled_pair_id, 
                                              sample_id = sample,
                                              episode_number, days_since_last_episode),
              by = c("subject_id", "sample_id")) 
  
  return(analysis_data)
}

# Get summaries -----------------------------------------------------------

get_joint_summary <- function(indiv_posteriors_joint)
{
  indiv_posteriors_joint %>% 
    group_by(subject_id, state_pair, joint_probability) %>% 
    filter(episode_number == max(episode_number)) %>% 
    mutate(percentage = round(joint_probability*100, 3),
           total_recurrences = episode_number-1) %>%
    select(-episode_number) %>% 
    arrange(subject_id, total_recurrences, percentage) %>% 
    relocate(total_recurrences, .before = state_pair)
  
}

get_marginal_summary <- function (indiv_posteriors_marginal)
{
  indiv_posteriors_marginal %>% 
    pivot_longer(cols = !subject_id & !episode_number, 
                 names_to = "posterior_type", 
                 values_to = "posterior_value") %>% 
    mutate(posterior_classification = case_when(posterior_type == "Posterior_marginal_prob_C" ~ "Recrudescence",
                                                posterior_type == "Posterior_marginal_prob_L" ~ "Relapse",
                                                posterior_type == "Posterior_marginal_prob_I" ~ "Reinfection"),
           posterior_classification = factor(posterior_classification, 
                                             levels = c("Relapse", "Recrudescence", "Reinfection"))) %>% 
    select(-posterior_type)
}

get_marginal_stats_summary <- function(marginal_summary)
{
  per_sample <- marginal_summary %>% 
    group_by(posterior_classification) %>% 
    summarise(min_post = min(posterior_value),
              mean_post = mean(posterior_value),
              max_post = max(posterior_value))
  
  per_study <- marginal_summary %>% 
    group_by(subject_id) %>%
    filter(posterior_value == max(posterior_value)) %>%
    ungroup() %>% 
    tabyl(posterior_classification) %>% 
    adorn_totals()
  
  per_sample %>% left_join(per_study, 
                           by = join_by(posterior_classification)) %>% 
                 return() 
}

get_fdr_summary <- function(fdr_df)
{
  fdr_df %>%
    summarise(
      mean_fdr   = mean(fdr),
      se_fdr     = sd(fdr) / sqrt(n()),
      lower    = mean_fdr - qnorm(0.975) * se_fdr,
      upper    = mean_fdr + qnorm(0.975) * se_fdr
    )
}

get_fdr_bs_summary <- function(fdr_df, times = 1000)
{
  set.seed(5)
  
  # bootstrapping 1000 replicates
  boots <- bootstraps(fdr_df, times = 1000)
  
  boot_stats <- boots %>%
    mutate(stats = map(splits,
                       ~ analysis(.x) %>%
                         summarise(
                           term     = "mean_fdr",        
                           estimate = mean(fdr)          
                         )))
  # get 95%CI by bootstrapping
  fdr_ci_bs <- int_pctl(boot_stats,
                        statistics = stats,
                        alpha      = 0.05)
  
  return(fdr_ci_bs)
}

# Plots -------------------------------------------------------------------

plot_episode_prob <- function(marginal_summary)
{
  marginal_summary %>% 
    group_by(subject_id) %>%
    arrange(desc(posterior_value), .by_group = TRUE) %>%
    mutate(subject_id = factor(subject_id, levels = unique(subject_id[order(posterior_value, decreasing = TRUE)]))) %>%
    
    ggplot(aes(x = factor(episode_number), y = posterior_value, group = episode_number, fill = posterior_classification)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c("Relapse" = "turquoise3",
                                 "Recrudescence" = "skyblue4",
                                 "Reinfection" = "magenta3")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x     = "Episode number", 
         y     = "Posterior probability", 
         fill = "") + 
    theme_bw() +
    facet_wrap(~subject_id)
}

plot_new_pairs_sols <- function(random_pairs)
{
  # merge back with full dataset
  sols_recurrent_null <- random_pairs %>% 
    left_join(sols %>% 
                clean_names() %>% 
                mutate(vis_date = mdy(vis_date)), 
              by = c("sample", "info", "info_d", "episodes", "trtgrp", "trt_label", "vis_date")) 
  
  # Data on episode number and time since last episode for new shuffled pairs
  # Note that the enrolment period was long for Solomon Islands, so we want to make sure we impose a time-order constraint for randomized pairs.
  episode_summary_null <- sols_recurrent_null %>%
    distinct(shuffled_pair_id, info_d, vis_date, sample) %>%
    arrange(shuffled_pair_id, vis_date) %>%
    group_by(shuffled_pair_id) %>%
    mutate(
      # get the episode number
      episode_number = row_number(),
      # format date
      first_date = first(vis_date),
      second_date = last(vis_date),
      days_since_last_episode = as.integer(last(vis_date) - first(vis_date))
    ) %>%
    ungroup()
  
  # print intermediate plot of 'new' shuffled pairs: 
  plot_random_pairs <- episode_summary_null %>% 
    ggplot(aes(x = days_since_last_episode, y = reorder(shuffled_pair_id, days_since_last_episode))) +
    geom_line(aes(group = shuffled_pair_id), color = "darkgrey") +
    geom_point() +
    labs(x = "Days since baseline episode",
         y = "Participant") +
    theme_bw() 
  
  return(plot_random_pairs)
}

plot_new_pairs_peru <- function(random_pairs)
{
  # merge back with full dataset
  peru_recurrent_null <- random_pairs %>% 
    left_join(peru %>% 
                clean_names() %>% 
                mutate(date = mdy(date)), 
              by = c("sample", "info" = "patient_name", "date", "episodes")) 
  
  ### Data on episode number and time since last episode for new shuffled pairs
  # Note that in the Peru study this is not necessarily a person's "first" episode, rather the first episode in the evaluated study period
  episode_summary_null <- peru_recurrent_null %>%
    distinct(shuffled_pair_id, day, date, sample) %>%
    arrange(shuffled_pair_id, date) %>%
    group_by(shuffled_pair_id) %>%
    mutate(
      # get the episode number
      episode_number = row_number(),
      first_date = first(date),
      second_date = last(date),
      days_since_last_episode = as.integer(last(date) - first(date))
    ) %>%
    ungroup()
  
  # print intermediate plot of 'new' shuffled pairs:
  plot_random_pairs <- episode_summary_null %>% 
    ggplot(aes(x = days_since_last_episode, y = reorder(shuffled_pair_id, days_since_last_episode))) +
    geom_line(aes(group = shuffled_pair_id), color = "darkgrey") +
    geom_point() +
    labs(x = "Days since first episode during study period",
         y = "Participant") +
    theme_bw() 
  
  return(plot_random_pairs)
}

# Randomize pairs ---------------------------------------------------------

randomize_pairs_sols <- function(data, max_attempts = 100) 
{
  attempt <- 1
  
  repeat {
    # Shuffle
    shuffled_data <- data %>%
      mutate(rand_order = sample(n())) %>%
      arrange(rand_order) %>%
      select(-rand_order) %>%
      mutate(pair_index = ceiling(row_number() / 2))
    
    # Keep only pairs that meet constraints
    valid_pairs <- shuffled_data %>%
      group_by(pair_index) %>%
      filter(n() == 2) %>%
      mutate(
        vis_date = mdy(vis_date),
        date1 = first(vis_date),
        date2 = last(vis_date),
        info1 = first(info),
        info2 = last(info)
      ) %>%
      filter(date1 < date2 & info1 != info2) %>%
      ungroup()
    
    if (nrow(valid_pairs) > 0) {
      # Found at least some pairs, return them
      valid_pairs <- valid_pairs %>%
        mutate(shuffled_pair_id = paste0("NewID-", pair_index))
      
      return(valid_pairs)
    } else {
      # If we found no valid pairs, try again
      attempt <- attempt + 1
      if (attempt > max_attempts) {
        stop("Couldn't form any valid pairs after multiple attempts.")
      }
    }
  }
}

randomize_pairs_peru <- function(data, max_attempts = 100) 
{
  attempt <- 1
  
  repeat {
    # Shuffle
    shuffled_data <- data %>%
      mutate(rand_order = sample(n())) %>%
      arrange(rand_order) %>%
      select(-rand_order) %>%
      mutate(pair_index = ceiling(row_number() / 2))
    
    # Keep only pairs that meet constraints
    valid_pairs <- shuffled_data %>%
      group_by(pair_index) %>%
      filter(n() == 2) %>%
      mutate(
        date1 = first(date),
        date2 = last(date),
        info1 = first(info),
        info2 = last(info)
      ) %>%
      filter(date1 < date2 & info1 != info2) %>%
      ungroup()
    
    if (nrow(valid_pairs) > 0) {
      # Found at least some pairs, return them
      valid_pairs <- valid_pairs %>%
        mutate(shuffled_pair_id = paste0("NewID-", pair_index))
      
      return(valid_pairs)
    } else {
      # If we found no valid pairs, try again
      attempt <- attempt + 1
      if (attempt > max_attempts) {
        stop("Couldn't form any valid pairs after multiple attempts.")
      }
    }
  }
}
