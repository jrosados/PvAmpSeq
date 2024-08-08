
########################################################
# COMPARING AMPSEQ MOI VS MICROSATELLITE MOI , Table S5 #
#######################################################

#Microsatellite data extracted from Manrique et al. 
#Microsatellite analysis reveals connectivity among geographically 
#distant transmission zones of Plasmodium vivax in the Peruvian Amazon:
#A critical barrier to regional malaria elimination. 
#PLoS Negl Trop Dis, 2019. 13(11): p. e0007876 

library(here)
library(tidyverse)

load("./intermediate_data/intermediate_data.RData")

MOI_microsat<-read.csv(here("raw_data", "supp_info_MOI_microsatellites.csv"))

#Counting the number of AmpSeq and MOI markers used and classification of clonality
AmpSeq_micro_MOI<- MOI_marker_sample_PE %>%
  group_by(sample) %>%
  summarise(n_marker_poly = sum(n_hap >= 2)) %>%
  inner_join (MOI_max_PE, by = c("sample"))   %>%
  inner_join(MOI_microsat, by = c("sample" = "Sample")) %>% inner_join(amp_success_samp_PE, by = c("sample")) %>%
  select(sample,max_moi, n_marker_poly, n, COI, n_poly_loci,profile) %>%
  rename(Sample = sample,
         AmpSeq_MOI = max_moi ,
         AmpSeq_poly_markers = n_marker_poly ,
         AmpSeq_markers = n,
         MS_MOI = COI ,
         MS_poly_markers = n_poly_loci ,
         MS_markers = profile)

print(AmpSeq_micro_MOI)

