
library(tidyverse)
library(here)

rm(list = ls())
#Load data from replicates

iSeq_Rep1 <- readRDS(here("raw_data", "iSeq_Replicate1.rds"))

iSeq_Rep2 <- readRDS(here("raw_data","iSeq_Replicate2.rds"))




#Selecting marker from haplotype data

Rep1_Chr07 <- iSeq_Rep1 %>%
  group_by(marker_id, sample_name,haplotype) %>% 
  summarise(total_reads = sum(n_reads)) %>% 
  mutate(freq = total_reads/sum(total_reads)) %>% group_by(sample_name) %>% 
  filter(marker_id == "Chr07") 
  
Rep2_Chr07 <- iSeq_Rep2 %>%
  group_by(marker_id, sample_name,haplotype) %>% 
  summarise(total_reads = sum(n_reads)) %>% 
  mutate(freq = total_reads/sum(total_reads)) %>% group_by(sample_name) %>% 
  filter(marker_id == "Chr07") 


# Keeping matched haplotypes in both duplicates
iSeq_replicates_match_dist<-
  inner_join(iSeq_Rep1,iSeq_Rep2, by = c( "marker_id", "sample_name",
                                          "haplotype","country"))%>% 
  group_by(marker_id, sample_name,haplotype) %>%  
  summarise(total_reads.x = sum(n_reads.x), total_reads.y = sum(n_reads.y)) %>% 
  mutate(freq.x = total_reads.x /sum(total_reads.x),freq.y = total_reads.y /sum(total_reads.y)) %>% 
  group_by(sample_name)


#Estimating distance metric from all matching replicates

HAP <- iSeq_replicates_match_dist

sam_ID <- unique(HAP$sample_name)
N_sam <- length(sam_ID)

marker_ID <- unique(HAP$marker_id)
N_marker <- length(marker_ID)


d0 <- matrix(NA, nrow=N_sam, ncol=N_marker)
rownames(d0) <- sam_ID
colnames(d0) <- marker_ID

for(i in 1:N_sam)
{
  for(j in 1:N_marker)
  {
    HAP_ij <- HAP[which(HAP$sample_name == sam_ID[i] & HAP$marker_id == marker_ID[j]),]		
    
    d0[i,j] <- sum( 0.5*abs(HAP_ij$freq.x - HAP_ij$freq.y) )	
  }
}

#Number of haplotypes per sample / Selecting the sample with more haplotypes
N_hap <- length(iSeq_replicates_match_dist$haplotype[iSeq_replicates_match_dist$marker_id == "Chr07" &
                                                      iSeq_replicates_match_dist$sample_name == "AAKC033"])

hap_col <- rainbow(N_hap)


tiff(file="Supp_figure1.tif", width=30, height=12, units="cm", res=300, compression = "lzw")

par(oma=c(4,4,4,4),mfrow=c(2,4),mar=c(4,6,4,6))

barplot(cbind(iSeq_replicates_match_dist$freq.x[iSeq_replicates_match_dist$marker_id == "Chr07" 
                                          & iSeq_replicates_match_dist$sample_name == "AAKC033"],
        iSeq_replicates_match_dist$freq.y[iSeq_replicates_match_dist$marker_id == "Chr07"
                                          & iSeq_replicates_match_dist$sample_name == "AAKC033"]), col=hap_col,
        xlab = NULL, ylab = "Frequency", main = "AAKC033 Chr07", cex.main = 2,cex.lab=2,cex.axis=2 )
axis(1, at=c(.5, 2), labels = c("R1", "R2"),cex.axis=2)

mtext('A',at=.01,side=3,outer=T,cex=1.2) 

R1 <- iSeq_replicates_match_dist$freq.x[iSeq_replicates_match_dist$marker_id == "Chr07" 
                                        & iSeq_replicates_match_dist$sample_name == "AAKC033"]
R2 <- iSeq_replicates_match_dist$freq.y[iSeq_replicates_match_dist$marker_id == "Chr07"
                                        & iSeq_replicates_match_dist$sample_name == "AAKC033"]

#AL020379 , example of Peruvian sample with 2 clones

#N_hap <- length(iSeq_replicates_match_dist$haplotype)

mtext('B',at=.25,side=3,outer=T,cex=1.2) 

plot( x=R1, y=R2,
      pch=19, col=hap_col,
      xlim=c(0,1), ylim=c(0,1), xlab = "R1 frequency", ylab = "R2 frequency", 
      main = "AAKC033 Chr07", cex.main = 2,cex.axis=2, cex.lab=2)
#0.01462789
points(x=c(0,1), y=c(0,1), type='l', lty="dashed")
mtext("d0 = 0.015", side = 1, outer = FALSE, line = -4)

for(j in 1:N_hap)
{
  points( x = c(R1[j], R2[j] + 0.5*(R1[j]-R2[j])),
          y = c(R2[j], R2[j] + 0.5*(R1[j]-R2[j])),
          type='l', col=hap_col[j])
}

mtext('C',at=.50,side=3,outer=T,cex=1.2) 

hist(d0[,1], breaks=5, ylab = "count", xlab = "d0", xlim = c(0,.1),ylim = c(0,20), 
     main = " d0 Chr 03", cex.main = 2, cex.axis=2, cex.lab=2)

mtext('D',at=.75,side=3,outer=T,cex=1.2) 

hist(d0[,2], breaks=5, ylab = "count", xlab = "d0", xlim = c(0,.1),ylim = c(0,20), 
     main = " d0 Chr 05", cex.main = 2, cex.axis=2, cex.lab=2)

mtext('E',at=.01,side=3,line= -15,outer=T,cex=1.2) 

hist(d0[,3], breaks=5, ylab = "count", xlab = "d0", xlim = c(0,.1),ylim = c(0,20), 
     main = " d0 Chr 06", cex.main = 2, cex.axis=2, cex.lab=2)

mtext('F',at=.25,side=3,line= -15,outer=T,cex=1.2) 

hist(d0[,4], breaks=5, ylab = "count", xlab = "d0", xlim = c(0,.1),ylim = c(0,20), 
     main = " d0 Chr 07", cex.main = 2, cex.axis=2, cex.lab=2)

mtext('G',at=.5,side=3,line= -15,outer=T,cex=1.2) 

hist(d0[,5], breaks=5, ylab = "count", xlab = "d0", xlim = c(0,.1),ylim = c(0,20), 
     main = " d0 Chr 08", cex.main = 2, cex.axis=2, cex.lab=2)

mtext('H',at=.75,side=3, line= -15, outer=T,cex=1.2) 

hist(d0[,6], breaks=5, ylab = "count", xlab = "d0", xlim = c(0,.1),ylim = c(0,20), 
     main = " d0 Chr 13", cex.main = 2, cex.axis=2, cex.lab=2)

dev.off()

# save.image("Hap_frequency_replicates.RDta")

