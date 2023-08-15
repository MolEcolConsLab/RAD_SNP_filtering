blat.filter <- function(blat.df, distance = 150, min.overlap = 25){
  
#Going to remove sequences that align in the first 150 bases.
#Seq 1 gets split into blocks that might be very small, so we want to filter those first
#Every sequence gets BLATed against every other, so every sequence ends up in Seq1, so can focus filtering there. 
#Goal at end of this is to have a vector of loci that aligned with other loci within the region of interest (i.e. the distance specified by "distance")
  repetitive_seqs <- blat.df %>% dplyr::filter(seq1 != seq2) %>% #Remove instances where a sequence is compared against itself
    dplyr::filter(seq1_end - seq1_start >= min.overlap) %>% #The part of seq1 that aligns must be longer than min.overlap.
    dplyr::filter(seq1_start < distance - min.overlap) %>% #Alignments further from the cut site than "distance" - "min.overlap" don't count.
    distinct(seq1) %>% 
    pull(seq1)
  
    
  
    #Count number of times each locus shows up in the first column
  # seq1.df <- blat.df %>% dplyr::filter(seq1 == seq2 | seq1_end < distance | seq1_start < distance - min.bases) %>% 
  #   dplyr::count(seq1) %>% 
  #   dplyr::rename(seq = seq1, n_seq1 = n)
  # 
  # #Count number of times each locus shows up in the second column
  # seq2.df <- blat.df %>% dplyr::filter(seq1 == seq2 | seq2_end < distance | seq2_start < distance - min.bases) %>% 
  #   dplyr::count(seq2) %>% 
  #   dplyr::rename(seq = seq2, n_seq2 = n)
  # 
  # #Combine dataframes together and filter for loci that only aligned with itself
  # blat.counts <- seq1.df %>% left_join(seq2.df, by = "seq") %>%
  #   dplyr::filter(n_seq1 == 1 & n_seq2 == 1)
  
  #blat.counts$CHROM <- gsub(".*_", "", blat.counts$seq)
  
  #whitelist.blat.loci <- blat.counts %>% dplyr::select(CHROM)
  
  repetitive_seqs <- gsub(".*_", "", repetitive_seqs)
  
}



filter.microhaps <- function(fix_tidy.snps, whitelist.blat.loci, snp.distance, min.snps, max.snps){
  
  final_SNPs_per_locus <- fix_tidy.snps %>% count(ChromKey, name = "SNPs")
  
  #Isolate microhaplotypes
  (microhaps <- final_SNPs_per_locus %>% dplyr::filter(SNPs >= min.SNPs & SNPs <= max.SNPs))
  
  microhaps.all <- fix_tidy.snps %>% inner_join(microhaps, by = "ChromKey") %>% 
    dplyr::select(CHROM, POS)
  
  #Convert table of SNPs to wide format to filter for distance between SNPs.
  microhaps_wide <- microhaps.all %>% inner_join(whitelist.blat.loci, by = "CHROM") %>% 
    dplyr::mutate(CHROM = as.integer(CHROM)) %>% 
    group_by(CHROM) %>% 
    mutate(POS_index = paste0("Pos_", seq(1:n()))) %>% 
    pivot_wider(values_from = POS, names_from = POS_index)
  
  
  microhaps_wide2 <- microhaps_wide %>% dplyr::filter(is.na(Pos_2) == FALSE) %>% #Set maximum and minimum number of SNPs, respectively. The first is the maximum SNPs + 1 and the second is minimum SNPs + 0.
    dplyr::filter(Pos_2 - Pos_1 <= snp.distance | Pos_3 - Pos_2 <= snp.distance | Pos_4 - Pos_3 <= snp.distance) #%>% #Include only microhaps with two SNPs that are within 100 bp of one another
  
  #Make long dataframe
  microhaps_long <- microhaps_wide2 %>% 
    pivot_longer(cols = starts_with("Pos_"), names_to = c(), values_to = "POS") %>% 
    drop_na()
  
  #Make dataframe of loci with more than max.snps
  too_many_SNPs <- microhaps_long %>% group_by(CHROM) %>% 
    summarize(n.SNPs = n()) %>% 
    dplyr::filter(n.SNPs > max.snps)
  
  microhaps.final <- microhaps_long %>% anti_join(too_many_SNPs, by = "CHROM")
}




filter.snps <- function(fix_tidy.snps, blacklist.blat.loci, snp.distance){
  
  snps.final <- fix_tidy.snps %>% dplyr::filter(!CHROM %in% blacklist.blat.loci) %>%
    dplyr::select(CHROM, ChromKey, POS) %>% 
    dplyr::filter(POS <= snp.distance)

  }