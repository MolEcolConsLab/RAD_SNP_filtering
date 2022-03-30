####-------------Part 1: Remove failed individuals from gt_tidy before continuing. This ensures that failed samples don't skew additional filtering.-----------####
#Remove failed individuals 

remove_failed <- function(gt_tidy){
  failed_indv <- gt_tidy %>% group_by(Indiv) %>% 
  summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>% 
  filter(percent_missing > Ipass) %>% 
  left_join(popmap, by = "Indiv")

#Store individuals that pass "failed" filter
pass_indv <- gt_tidy %>% anti_join(failed_indv, by = "Indiv") %>%
  group_by(Indiv) %>% 
  summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>% 
  left_join(popmap, by = "Indiv")

#Keep only individuals that didn't fail
gt_tidy_fltr <- gt_tidy %>%  semi_join(pass_indv, by = "Indiv")

cat("Summary: There are", nrow(failed_indv), "individuals that failed at a threshold of", Ipass, "% missing data and were removed before filtering.")

return(gt_tidy_fltr)
}
####-----------------Part 2: re-calculate missingness and perform second round of filtering for individuals and loci---------------####

# #Checked 05/06/2021
calc_missing <- function(gt_tidy_fltr){

  missing_indv <- gt_tidy_fltr %>% group_by(Indiv) %>%
  summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
  filter(percent_missing > Imiss) %>%
  left_join(popmap, by = "Indiv")

keep_indv <- gt_tidy_fltr %>% anti_join(missing_indv, by = "Indiv") %>%
  group_by(Indiv) %>%
  summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
  left_join(popmap, by = "Indiv")

missing_loci <- gt_tidy_fltr %>% group_by(ChromKey, POS) %>%
  summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
  filter(percent_missing > Lmiss) %>%
  ungroup()


#Calculate locus depth and remove outliers (mean depth + 2*sd)
mean_locus_depth <- mean(gt_tidy_fltr$gt_DP_recode, na.rm = TRUE)
locus_depth_sd <- sd(gt_tidy_fltr$gt_DP_recode, na.rm = TRUE)

#Checked 05/06/2021
ridonkulous_coverage_loci <- gt_tidy_fltr %>% group_by(ChromKey, POS) %>%
  summarize(mean_depth = mean(gt_DP_recode, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(mean_depth > (mean_locus_depth + (2*locus_depth_sd)))


#Dataframe of loci to keep + calculate missingness across individuals
#Checked 05/06/2021
keep_loci <- gt_tidy_fltr %>% anti_join(missing_loci, by = c("ChromKey", "POS")) %>%
  anti_join(ridonkulous_coverage_loci, by = c("ChromKey", "POS")) %>%
  group_by(ChromKey, POS) %>%
  summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100)

    cat("There are", nrow(missing_indv), "samples with >", Imiss, "% missing data that will be removed from the dataset;\n", nrow(keep_indv), "samples remain.\n\n")
  
  cat("There are ", nrow(missing_loci), "SNPs that are missing from >= ", Lmiss, "% of individuals that will be removed from the dataset and ", nrow(ridonkulous_coverage_loci), "SNPs with coverage higher than the overall mean + 2 standard deviations.", nrow(keep_loci), "SNPs remain.")  
#return(list(keep_indv, keep_loci)) #Don't need to save these anymore
}


#Create filtered genotype table keeping only loci and individuals that passed the above filters
##Also calculate allele frequency for each locus so we can remove sites that have no variation after filtering. Some sites may have been called as SNPs due to one or two minor alleles in individuals that have now been removed due to missing data.
#Checked 05/06/2021
filter_all <- function(gt_tidy_fltr, Imiss = Imiss, Lmiss = Lmiss, interval = 5){

  #Make sequence of values to loop over for individual missingness
  indvLoopVec <- seq(from=100-interval, to=Imiss, by = -interval)
    #Make sequence of thresholds for filtering loci based on the length of the sequence being used to filter individuals
  locLoopVec <- seq(from=100, to=Lmiss, length.out = length(indvLoopVec)+1)
  locLoopVec <- locLoopVec[-1]
  
  gt_tidy.temp <- gt_tidy_fltr
    
  #To speed code, loop through all but last iteration, then add more detailed filtering at last iteration (e.g. removing ridonkulous coverage loci and calculating heterozygosity)
  
  #Loop over values of individual missingness
  for(i in 1:(length(indvLoopVec)-1)){
    ind.t <- indvLoopVec[i]
    loc.t <- locLoopVec[i]
        
    #Calculate missingness
    missing_indv <- gt_tidy.temp %>% group_by(Indiv) %>%
      summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
      filter(percent_missing > ind.t) %>%
      left_join(popmap, by = "Indiv")

    #Create dataframe of individuals to keep
    keep_indv <- gt_tidy.temp %>% anti_join(missing_indv, by = "Indiv") %>%
      group_by(Indiv) %>%
      summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
      left_join(popmap, by = "Indiv")
    
      #Calculate missingness for loci
      missing_loci <- gt_tidy.temp %>% group_by(ChromKey, POS) %>%
        summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
        filter(percent_missing > loc.t) %>%
        ungroup()
    
      #Create dataframe of loci to keep
      keep_loci <- gt_tidy.temp %>% anti_join(missing_loci, by = c("ChromKey", "POS")) %>%
        group_by(ChromKey, POS) %>%
        summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100)

      #Update gt_tidy.temp dataframe for next loop
      gt_tidy.temp <- gt_tidy.temp %>% semi_join(keep_loci, by = c("ChromKey", "POS")) %>% #Keep only loci that passed filter thresholds
        semi_join(keep_indv, by = "Indiv")
      
      cat("There are", nrow(missing_indv), "samples with >", ind.t, "% missing data that will be removed from the dataset;\n", nrow(keep_indv), "samples remain.\n")
      
      cat("There are ", nrow(missing_loci), "SNPs that are missing from >= ", loc.t, "% of individuals that will be removed from the dataset; \n", nrow(keep_loci), "SNPs remain.\n\n")
  }

  #Run last iteration outside of loop
  #Calculate missingness for individuals
  missing_indv <- gt_tidy.temp %>% group_by(Indiv) %>%
    summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
    filter(percent_missing > Imiss) %>%
    left_join(popmap, by = "Indiv")
  
  #Create dataframe of individuals to keep
  keep_indv <- gt_tidy.temp %>% anti_join(missing_indv, by = "Indiv") %>%
    group_by(Indiv) %>%
    summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
    left_join(popmap, by = "Indiv")
  
  #Calculate missingness for loci
  missing_loci <- gt_tidy.temp %>% group_by(ChromKey, POS) %>%
    summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>%
    filter(percent_missing > Lmiss) %>%
    ungroup()

  #Calculate locus depth and remove outliers (mean depth + 2*sd)
  mean_locus_depth <- mean(gt_tidy_fltr$gt_DP_recode, na.rm = TRUE)
  locus_depth_sd <- sd(gt_tidy_fltr$gt_DP_recode, na.rm = TRUE)
  
  #Checked 05/06/2021
  ridonkulous_coverage_loci <- gt_tidy_fltr %>% group_by(ChromKey, POS) %>%
    summarize(mean_depth = mean(gt_DP_recode, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(mean_depth > (mean_locus_depth + (2*locus_depth_sd)))
  
  #Create dataframe of loci to keep
  keep_loci <- gt_tidy.temp %>% anti_join(missing_loci, by = c("ChromKey", "POS")) %>%
    anti_join(ridonkulous_coverage_loci, by = c("ChromKey", "POS")) %>%
    group_by(ChromKey, POS) %>%
    summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100)
  
  #Create final dataframe to export
  gt_tidy_filt1 <- gt_tidy.temp %>% semi_join(keep_loci, by = c("ChromKey", "POS")) %>% #Keep only loci that passed filter thresholds
    semi_join(keep_indv, by = "Indiv") %>% #Keep only individuals that passed the filter thresholds
    dplyr::select(!gt_DP) %>% 
    dplyr::rename(gt_DP = gt_DP_recode) %>% #Save gt_recode as gt_DP so defined depth threshold for missingness is perpetuated.
    mutate(allele_A = ifelse(gt_GT == "0/0", 2, ifelse(gt_GT == "0/1", 1, 0)), allele_B = ifelse(gt_GT == "0/0", 0, ifelse(gt_GT == "0/1", 1, 2))) %>%
    mutate(zygosity = ifelse(gt_GT == "0/1", "heterozygous", "homozygous")) #Calculate quantity of alleles because some loci have no variation after filtering individuals
  
    #Remove SNPs that lost variation when missing individuals were removed
    not_SNPs <- gt_tidy_filt1 %>% group_by(ChromKey, POS) %>%
      summarize(total_allele_A = sum(allele_A, na.rm = TRUE), total_allele_B = sum(allele_B, na.rm = TRUE)) %>%
      mutate(total_alleles = total_allele_A + total_allele_B) %>%
      filter(total_allele_A == 0 | total_allele_B == 0)
    #filter(total_allele_A <= minor_allele_indv | total_allele_B <= minor_allele_indv) #If wanting to filter out based on minor allele freq
    
    #Keep remaining variants
    SNPs <- fix_tidy %>% anti_join(not_SNPs, by = c("ChromKey", "POS"))
    
    #Filter fixed table for loci to keep
    fix_tidy_filt1 <- fix_tidy %>% semi_join(keep_loci, by = c("ChromKey", "POS")) %>%
      anti_join(not_SNPs, by = c("ChromKey", "POS"))
    

        cat("There are", nrow(missing_indv), "samples with >", Imiss, "% missing data that will be removed from the dataset;\n", nrow(keep_indv), "samples remain.\n\n")
        
        cat("There are ", nrow(missing_loci), "SNPs that are missing from >= ", Lmiss, "% of individuals that will be removed from the dataset and ", nrow(ridonkulous_coverage_loci), "SNPs with coverage higher than the overall mean + 2 standard deviations.", nrow(keep_loci), "SNPs remain.")  
        return(list(fix_tidy_filt1, gt_tidy_filt1, keep_indv, keep_loci))
      }

merge_tables <- function(gt_tidy_filt1, fix_tidy_filt1){
gt_tidy_filt2 <- gt_tidy_filt1 %>% semi_join(fix_tidy_filt1, by = c("ChromKey", "POS"))
return(gt_tidy_filt2)
}



####----------Part 3: Calculate read depth per individual and locus after filtering------------####

#Calculate mean read depth per individual (across loci) after filtering
calc_depth_indv <- function(gt_tidy_filt2){
  mean_read_depth_indv <- gt_tidy_filt2 %>% group_by(Indiv) %>%
  summarize(mean_depth = mean(gt_DP, na.rm = TRUE)) %>%
  arrange(desc(mean_depth)) %>%
  left_join(popmap, by = "Indiv")
  
  return(mean_read_depth_indv)
}

calc_depth_loci <- function(gt_tidy_filt2){
#Calculate mean read depth per locus (across individuals)
mean_read_depth_loci <- gt_tidy_filt2 %>% group_by(ChromKey, POS) %>%
  summarize(mean_depth = mean(gt_DP, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(mean_depth))

return(mean_read_depth_loci)
}

sumstats <- function(){

  #Remove SNPs that lost variation when missing individuals were removed
  not_SNPs <- gt_tidy_filt1 %>% group_by(ChromKey, POS) %>%
    summarize(total_allele_A = sum(allele_A, na.rm = TRUE), total_allele_B = sum(allele_B, na.rm = TRUE)) %>%
    mutate(total_alleles = total_allele_A + total_allele_B) %>%
    filter(total_allele_A == 0 | total_allele_B == 0)
  #filter(total_allele_A <= minor_allele_indv | total_allele_B <= minor_allele_indv) #If wanting to filter out based on minor allele freq
  
cat("After filtering, ", nrow(not_SNPs), " positions are newly invariant (not SNPs) and were filtered.\n Overall, kept ", nrow(keep_indv), " individuals and ", nrow(fix_tidy_filt1), " SNPs over ", nrow(SNPs_per_locus), " loci.\n (assuming a de novo approach to locus formation - the number of loci will be different for reference-based approaches).\n\nThe breakdown of remaining samples by population is:")

  
  #Number of remaining individuals from each population
  table(keep_indv$Pop)
  
}

##########################################################
##################HWE and heterozygosity##################
##########################################################

###################################################################################
####Part 6.1: calculate and visualize individual heterozygosity after filtering####
###################################################################################

# Calculate heterozygosity
# (H[exp] - H[obs]) / H[exp]
# H[exp] = 2pq
##http://www.uwyo.edu/dbmcd/molmark/practica/fst.html

#Count alleles by population
#Checked 05/06/2021
calc_het <- function(gt_tidy_filt2, popmap){

#Calculate individual level heterozygosity. Accounts for missing data by only including SNPs that were labeled as heterozygous or homozygous.
Obs_het_indv <- gt_tidy_filt2 %>%
  group_by(Indiv) %>% 
  summarize(Obs_het = sum(zygosity == "heterozygous", na.rm = TRUE)/sum(zygosity == "heterozygous" | zygosity == "homozygous", na.rm = TRUE)) %>%
  left_join(popmap, by = "Indiv")


return(Obs_het_indv)
}


############################################################################
####Part 6.2: Calculate deviation from HWE and remove super deviant loci####
############################################################################

#Calculate frequency of AA, Aa, and aa genotypes at each locus and format for HWE calculations (which only accepts a matrix of genotype frequencies)
#Checked 05/06/2021
calc_hwe <- function(gt_tidy_filt2){
  obs_allele_count <- gt_tidy_filt2 %>%
  group_by(ChromKey, POS) %>% 
  summarize(obs_AA_count = sum(allele_A == 2, na.rm = TRUE),
            obs_Aa_count = sum(allele_A == 1 & allele_B == 1, na.rm = TRUE),
            obs_aa_count = sum(allele_B == 2, na.rm = TRUE)) %>%
  mutate(ChromPos = paste0(ChromKey, "_", POS)) %>% 
  column_to_rownames(var = "ChromPos") %>% 
  dplyr::select(obs_AA_count, obs_Aa_count, obs_aa_count)

#Calculate concordance with/deviation from HWE
hwe_stats <- HWExactStats(as.matrix(obs_allele_count))

#Reformat dataframe with HWE calculations so we can merge with genotype dataframe later
locus_HWE <- cbind(obs_allele_count, hwe_stats) %>% 
  rownames_to_column(var = "ChromPos") %>%
  separate(col = ChromPos, into = c("ChromKey", "POS"), sep = "_") %>% 
  mutate(ChromKey = as.numeric(ChromKey), POS = as.numeric(POS))

cat("Breakdown of loci:")
#Calculate how many loci are in and out of HWE based on the threshold set above
print(locus_HWE %>% mutate(hwe_sig = ifelse(hwe_stats < HWE_inquire, "out of HWE", "in HWE")) %>% 
  group_by(hwe_sig) %>% 
  summarize(number_loci = n()) %>% 
  mutate(percent_loci = number_loci/(sum(number_loci))*100))

return(locus_HWE)
}

final_sumstats <- function(){
  cat("Filtering for individuals with excess heterozygosity (>", max_het, ") removed ", nrow(gt_tidy_hetero_remove), "individuals.\n", nrow(keep_hetero_indiv), "individuals remain.\n\nFiltering for loci out of HWE at p <", HWE_threshold, "results in removal of", nrow(locus_HWE_remove),"SNPs.\n\nAfter filtering for short-distance LD (if selected), keeping one random SNP per locus,", nrow(fix_tidy_filt_final), "SNPs remain.")
}