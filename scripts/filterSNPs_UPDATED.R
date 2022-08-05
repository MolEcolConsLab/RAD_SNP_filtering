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
filter_all <- function(gt_tidy_fltr, Imiss = Imiss, Lmiss = Lmiss, interval = 5, min.mean.depth = min.mean.depth, max.mean.depth = max.mean.depth){

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
    loc.t <- round(locLoopVec[i], 0)
        
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
  keep_loci.temp <- gt_tidy.temp %>% anti_join(missing_loci, by = c("ChromKey", "POS")) %>%
    anti_join(ridonkulous_coverage_loci, by = c("ChromKey", "POS")) %>%
    group_by(ChromKey, POS) %>%
    summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100)
  
  
  #Create final dataframe to export
  gt_tidy_filt.temp2 <- gt_tidy.temp %>% semi_join(keep_loci.temp, by = c("ChromKey", "POS")) %>% #Keep only loci that passed filter thresholds
    semi_join(keep_indv, by = "Indiv") %>% #Keep only individuals that passed the filter thresholds
    dplyr::select(!gt_DP) %>% 
    dplyr::rename(gt_DP = gt_DP_recode) %>% #Save gt_recode as gt_DP so defined depth threshold for missingness is perpetuated.
    mutate(allele_A = ifelse(gt_GT == "0/0", 2, ifelse(gt_GT == "0/1", 1, 0)), allele_B = ifelse(gt_GT == "0/0", 0, ifelse(gt_GT == "0/1", 1, 2))) %>%
    mutate(zygosity = ifelse(gt_GT == "0/1", "heterozygous", "homozygous")) #Calculate quantity of alleles because some loci have no variation after filtering individuals
  
  #Calcualte mean depth of loci to filter
  mean_read_depth_loci.temp <- gt_tidy_filt.temp2 %>% group_by(ChromKey, POS) %>%
    summarize(mean_depth = mean(gt_DP, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(mean_depth))
  
  #Filter loci with low mean read depth that are prone to higher miscall rates
  lo_hi_depth_loci <- mean_read_depth_loci.temp %>% dplyr::filter(mean_depth < min.mean.depth | mean_depth > max.mean.depth)
  
  
  
  gt_tidy_filt1 <- gt_tidy_filt.temp2 %>% anti_join(lo_hi_depth_loci, by = c("ChromKey", "POS"))
  
  
  keep_loci <- gt_tidy_filt1 %>% anti_join(lo_hi_depth_loci, by = c("ChromKey", "POS")) %>%
    group_by(ChromKey, POS) %>%
    summarize(percent_missing = sum(is.na(gt_DP))/n() * 100)
    
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
  not_SNPs <- gt_tidy_filt2 %>% group_by(ChromKey, POS) %>%
    summarize(total_allele_A = sum(allele_A, na.rm = TRUE), total_allele_B = sum(allele_B, na.rm = TRUE)) %>%
    mutate(total_alleles = total_allele_A + total_allele_B) %>%
    filter(total_allele_A == 0 | total_allele_B == 0)
  #filter(total_allele_A <= minor_allele_indv | total_allele_B <= minor_allele_indv) #If wanting to filter out based on minor allele freq
  
cat("After filtering, kept ", nrow(keep_indv), " individuals and ", nrow(fix_tidy_filt1), " SNPs over ", nrow(SNPs_per_locus), " loci.\n (assuming a de novo approach to locus formation - the number of loci will be different for reference-based approaches).\n\nThe breakdown of remaining samples by population is:")

  
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

#Adding argument to potentially loop over populations. Need to keep working on this 06/28/2022
calc_hwe <- function(gt_tidy_filt2, by_pop = T){
  if(by_pop == F){
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
locus_HWE.df <- cbind(obs_allele_count, hwe_stats) %>% 
  rownames_to_column(var = "ChromPos") %>%
  separate(col = ChromPos, into = c("ChromKey", "POS"), sep = "_") %>% 
  mutate(ChromKey = as.numeric(ChromKey), POS = as.numeric(POS))

cat("Breakdown of loci:")
#Calculate how many loci are in and out of HWE based on the threshold set above
print(locus_HWE.df %>% mutate(hwe_sig = ifelse(hwe_stats < HWE_inquire, "out of HWE", "in HWE")) %>% 
  group_by(hwe_sig) %>% 
  summarize(number_loci = n()) %>% 
  mutate(percent_loci = number_loci/(sum(number_loci))*100))

return(locus_HWE.df)
  } else if(by_pop == T){
    #Split gt_tidy by population
    hwe.pop.list <- gt_tidy_filt2 %>%
      group_by(Pop) %>% 
      group_split()
    
    #Add names to list elements, ordering the same way as above
    pop.names <- gt_tidy_filt2 %>%
      group_by(Pop) %>% 
      distinct(Pop) %>% 
      pull(Pop)

    names(hwe.pop.list) <- pop.names
        
    calc.hwe.pop <- function(.x){
      .x %>% group_by(ChromKey, POS, Pop) %>% 
        summarize(obs_AA_count = sum(allele_A == 2, na.rm = TRUE),
                  obs_Aa_count = sum(allele_A == 1 & allele_B == 1, na.rm = TRUE),
                  obs_aa_count = sum(allele_B == 2, na.rm = TRUE)) %>%
        mutate(ChromPos = paste0(ChromKey, "_", POS)) %>% 
        column_to_rownames(var = "ChromPos") %>% 
        dplyr::select(obs_AA_count, obs_Aa_count, obs_aa_count)
    }
    
    obs_allele_count <- hwe.pop.list %>% map(.f = calc.hwe.pop)
    
    #Calculate concordance with/deviation from HWE for each population (population = list element)
    hwe_stats <- map(.x = as.matrix(obs_allele_count), .f = HWExactStats)
    
    names(hwe_stats) <- pop.names
    
    locus_HWE.list <- list()
    locus_HWE.df <- NULL
    #Reformat dataframe with HWE calculations so we can merge with genotype dataframe later
    for(l in 1:length(obs_allele_count)){
      locus_HWE.list[[l]] <- cbind(obs_allele_count[[l]], hwe_stats = hwe_stats[[l]]) %>%
        rownames_to_column(var = "ChromPos") %>%
        separate(col = ChromPos, into = c("ChromKey", "POS"), sep = "_") %>%
        mutate(ChromKey = as.numeric(ChromKey), POS = as.numeric(POS)) %>% 
        mutate(pop = pop.names[l])
      
      locus_HWE.df <- rbind(locus_HWE.df, locus_HWE.list[[l]])
      
    }
    
    number.all.loci <- nrow(fix_tidy_filt1)
    
    cat("Breakdown of loci:")
    #Calculate how many loci are in and out of HWE based on the threshold set above
    print(locus_HWE.df %>% mutate(hwe_sig = ifelse(hwe_stats < HWE_inquire, "out of HWE", "in HWE")) %>% 
            group_by(hwe_sig, pop) %>% 
            summarize(number_loci = n()) %>% 
            mutate(percent_loci = round((number_loci/number.all.loci)*100, 2)))
  }
  return(locus_HWE.df)
}


############################################################################
####Part 6.3: Filter for LD by removing loci with R2 values above a threshold####
############################################################################
#An alternative if we want to us bcftools to calculate different R^2 values
#In this function; set a parameter for the number of R2 thresholds used
#Separate the FILTER column based on the number of thresholds used, and name the columns based on the R2 value.
#Have the filter thresholds show up in the same order, then name the columns based on the threshold.
#Have an if statement to evaluate which threshold we're filtering; if it's the smallest (which will be column 1), then keep the rows that do NOT have the filter label; if it's one of the smallest, then filter by that column.

filter.LD <- function(df = fix_tidy_filt1, threshold = 0.4, metric = "D"){

  df2 <- df %>% mutate(CHROMPOS = paste0(CHROM, "_", POS))
  
  if(metric == "D"){
    loci2remove <- df2 %>% mutate(LINKAGE_POS = paste0(CHROM, "_", POS_LD)) %>% 
      dplyr::filter(LD >= threshold) %>% 
      pull(LINKAGE_POS)
    
    df_final <- df2 %>% dplyr::filter(!CHROMPOS %in% loci2remove)
  } else if(metric == "R2"){
    loci2remove <- df2 %>% mutate(LINKAGE_POS = paste0(CHROM, "_", POS_R2)) %>% 
      dplyr::filter(R2 >= threshold) %>% 
      pull(LINKAGE_POS)
    
    df_final <- df2 %>% dplyr::filter(!CHROMPOS %in% loci2remove)
    
  }

  return(df_final)
  }


final_sumstats <- function(){
  cat("Filtering for individuals with excess heterozygosity (>", max_het, ") removed ", nrow(gt_tidy_hetero_remove), "individuals.\n", nrow(keep_hetero_indiv), "individuals remain.\n\nFiltering for loci out of HWE at p <", HWE_threshold, "results in removal of", nrow(locus_HWE_remove),"SNPs.\n\nAfter filtering for short-distance LD (if selected), keeping one random SNP per locus,", nrow(fix_tidy_filt_final), "SNPs remain.")
}