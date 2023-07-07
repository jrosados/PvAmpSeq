
###############################################
# COMPARING AMPSEQ MOI VS MICROSATELLITE MOI  #
###############################################

Paulo_AmpSeq_samples<- AmpSeq_samples_data_PE %>%
  inner_join(MOI_microsat, by = c("SampleID" = "Sample")) 


MOI_microsat<-read.csv("supp_info_MOI_microsatellites.csv")

AmpSeq_micro_MOI<- MOI_marker_sample_PE %>%
  group_by(sample) %>%
  summarise(n_marker_poly = sum(n_hap >= 2)) %>%
  inner_join (MOI_max_PE, by = c("sample"))   %>%
  inner_join(MOI_microsat, by = c("sample" = "Sample")) %>% inner_join(amp_success_samp_PE, by = c("sample")) %>%
  #mutate(sample_id = c(1,2,3,4,5)) %>%
  select(sample,max_moi, n_marker_poly, n, COI, n_poly_loci,profile) %>%
  rename(Sample = sample,
         AmpSeq_MOI = max_moi ,
         AmpSeq_poly_markers = n_marker_poly ,
         AmpSeq_markers = n,
         MS_MOI = COI ,
         MS_poly_markers = n_poly_loci ,
         MS_markers = profile)

AmpSeq_vs_MS <- write.csv(AmpSeq_micro_MOI, "AmpSeq_vs_MS.csv", row.names = F)