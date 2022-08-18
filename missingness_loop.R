rm(list=ls())

# make sure you have Rtools installed
#if (!require("devtools")) install.packages("devtools")

# install latest versions of required packages from GitHub, etc. (NB need XCode command line tools installed for some of these)
#devtools::install_github("thierrygosselin/radiator")
#devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
#devtools::install_github("thibautjombart/adegenet")
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("SNPRelate")

library(vcfR)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(Polychrome)
library(HardyWeinberg)
library(strataG)
library(radiator)
library(SeqArray)
library(genetics)
library(whoa)

#Source custom functions (change paths accordingly if needed)
source("scripts/ggplot.R")
source("scripts/VCFfilterstats.R")
source("scripts/xtrafunctions.R")
source("scripts/plotstats.R")
source("scripts/filterSNPs_UPDATED.R")
source("scripts/HDplot.R")

options(dplyr.summarise.inform = FALSE) #Shut down notices for summarise (they're annoying)

file_directory <- "../data/snails/"
file_prefix <- "UroFull_M1n6"
popmap_file <- "popmap_full"

popmap <- read_tsv(paste0(file_directory, popmap_file, ".txt"), col_names = FALSE) %>% rename(Indiv = X1, Pop = X2)

raw_vcf <- vcfR::read.vcfR(paste0(file_directory,file_prefix,".recode.vcf")) #This takes a minute


#identify putative paralogs and plot - this may take a long time
paralog_calcs <- HDplot(raw_vcf) %>%
  arrange(D)

#Create heterozygosity vs read ratio deviation plot
paralog_calcs %>% ggplot()+geom_point(aes(x=H,y=D), alpha = 0.1) + 
  scale_y_continuous(name = "D", breaks = seq(-20, 20, by = 2)) +
  scale_x_continuous(name = "H", breaks = seq(0, 1, by = 0.05)) +
  ggtitle("Heterozygosity (H) vs Read Ratio Deviation (D)") +
  geom_hline(yintercept = c(-4, 4), color = "blue", linetype = "dashed", size = 1.1)

#Set limits for H and D based on plot from above
H_lim <- 0.45
D_lim <- 4


#Create heterozygosity vs read ratio plot, if wanted
#paralog_calcs %>% ggplot() + geom_point(aes(x=H,y=ratio))
#Make genotype frequency scatter plot
gfreqs <- exp_and_obs_geno_freqs(raw_vcf) #Calculate expected and observed genotype frequencies
geno_freqs_scatter(gfreqs)

#Infer heterozygote miscall rate - this can take a long time
het.miscall.all <- infer_m(raw_vcf, minBin = 1e15) #Set the minBin argument very high so it runs over all read depths

het.miscall.all$m_posteriors #mean column shows percent heterozygote miscall rate as a decimal (e.g. 0.1 = 10% miscall rate)

#Infer miscall rate for different read depth bins
#Can set the bin size to whatever we want; Consider targeting ~ 100 bins, so dividing n_total by 100 and setting the bin size to that (but can vary the bin size b to be whatever you want; see help for infer_m function for more info)
b <- het.miscall.all$m_posteriors$total_n/100
het.miscall.bins <- infer_m(raw_vcf, minBin = b)

#Plot heterozygote miscall rate at different read depths to help with filtering.
posteriors_plot(het.miscall.bins$m_posteriors)

#Create dataframe of putative paralogs based on thresholds set above
paralogs <- paralog_calcs %>% filter(H > H_lim & abs(D) > D_lim) %>% 
  mutate(POS = as.integer(POS))

print(paste0("Identified ", nrow(paralogs), " likely paralogs. These will be removed from the dataset."))

##If interested in viewing or extracting specific pieces of information, use the below code.
#queryMETA(raw_vcf)
#queryMETA(raw_vcf, element = "DP") #Look at meta region. Can specify arguments that return more details

#info <- data.frame(getFIX(raw_vcf)) #Query the fixed region
#raw_vcf@gt[1:6, 1:4] #Look at gt portion of vcf file

#Meanings of the fields in the gt portion of the vcf file output from Stacks and VCFTools, which is colon delimited in each cell.
#meta_meanings <- cbind(c("GT", "DP", "AD", "GQ", "GL"), c("Genotype", "Read depth", "Allelic depth", "Genotype quality", "Genotype likelihood")) 


#Convert vcf file to tibble
options(tibble.print_max = 30)
raw_vcf_tidy <- vcfR2tidy(raw_vcf) #This takes awhile


#Extract fixed and gt dataframes and subset for useful columns
## ChromKey is a link between the fixed portion and the gt portion of the files. Need to merge the paralog df with the fixed df to grab the ChromKey column to filter paralogs at the next step.
paralogs <- paralogs %>% inner_join(raw_vcf_tidy$fix, by = c("CHROM", "POS")) %>% 
  dplyr::select(ChromKey, CHROM, POS)

#Convert fixed element to tidy format and remove paralogs
fix_tidy <- raw_vcf_tidy$fix %>% 
  anti_join(paralogs, by = c("CHROM", "POS"))

#Convert gt element to tidy format and remove paralogs 
gt_tidy <- raw_vcf_tidy$gt %>%
  left_join(popmap, by = "Indiv") %>% 
  anti_join(paralogs, by = c("ChromKey", "POS"))

#Store meta region
gt_meta <- raw_vcf_tidy$meta

#After extracting the relevant dataframes, remove the original imported vcf files from the environment to save working memory (you can always re-run parts of this chunk if you want to work directly with those files)
#rm(raw_vcf)
#rm(raw_vcf_tidy)

## Step 1: set criteria below which very low-depth loci will be called as missing. Consider setting this to 2m
min_depth_indv <- 6 #Minimum depth for a locus to be included in subsequent analyses

#Recode loci in each individual as missing based on the threshold above
gt_tidy <- gt_tidy %>% mutate(gt_DP_recode = ifelse(gt_DP < min_depth_indv, NA, gt_DP))

#Calculate missingness per individual
all_indv <- gt_tidy %>% group_by(Indiv) %>% 
  summarize(percent_missing = sum(is.na(gt_DP_recode))/n() * 100) %>% 
  left_join(popmap, by = "Indiv")

#Plot missingness after recoding low-depth loci
p0.1 <- ggplot(all_indv, aes(x = percent_missing)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey95") +
  labs(x = "% missing data per indv")
p0.1 #look at this to gauge and tweak threshold for failed vs. passed samples


Ipass.vec <- c(70, 75, 80, 85)
Imiss.vec <- c(20, 30, 40, 50, 60) #After removing failed individuals, filter individuals that still have this much missing data or more
Lmiss.vec <- c(10, 15, 20, 25, 30) #After removing failed individuals, remove loci that are missing in this percent of individuals
min.mean.depth.vec <- c(25, 50)
max.mean.depth.vec <- c(75, Inf)
int <- 2


#For quick test of script
# Ipass.vec <- c(70,85)
# Imiss.vec <- c(20,60) #After removing failed individuals, filter individuals that still have this much missing data or more
# Lmiss.vec <- c(10,30) #After removing failed individuals, remove loci that are missing in this percent of individuals
# min.mean.depth.vec <- c(25,50)
# max.mean.depth.vec <- c(75,Inf)

stats.df = stats.df.temp <- NULL


#----------------start loop------------------
for(min.m in 1:length(min.mean.depth.vec)){
  for(max.m in 1:length(max.mean.depth.vec)){
    for(Ip in 1:length(Ipass.vec)){
      for(Im in 1:length(Imiss.vec)){
        for(Lm in 1:length(Lmiss.vec)){

          Ipass <- Ipass.vec[Ip]
          Imiss <- Imiss.vec[Im]
          Lmiss <- Lmiss.vec[Lm]
          min.mean.depth <- min.mean.depth.vec[min.m]
          max.mean.depth <- max.mean.depth.vec[max.m]
          
#Remove failed individuals before filtering
gt_tidy_fltr <- remove_failed(gt_tidy) #Checked 05/23/2021

keep <- filter_all(gt_tidy_fltr, Imiss = Imiss, Lmiss = Lmiss, interval = int, min.mean.depth = min.mean.depth, max.mean.depth = max.mean.depth) #loops over values of missingness from 100% to Imiss in intervals specified by the interval argument. Iterative filtering improves the number of individuals retained.

#Split output from above function into separate dataframes
fix_tidy_filt1 <- keep[[1]]
gt_tidy_filt1 <- keep[[2]]
(keep_indv <- keep[[3]] %>% arrange(percent_missing))
keep_loci <- keep[[4]]

#if curious, calculate individual and locus missingness, and how much we will remove if we filtered immediately without an iterative approach, just to compare. Will likely see fewer individuals and/or loci retained here than in the loop above.
#miss.stats <- calc_missing(gt_tidy_fltr) 

#filter gt dataframe so only good loci remain
gt_tidy_filt2 <- merge_tables(gt_tidy_filt1, fix_tidy_filt1)

#Identify one random SNP per locus
##JDS note: Comment out below if reference-based approach was used upstream
LD_loci_to_keep <- fix_tidy_filt1 %>% group_by(CHROM) %>%
  slice_sample(n=1)

fix_tidy_filt_final <- fix_tidy_filt1 %>% semi_join(LD_loci_to_keep, by = c("ChromKey", "CHROM", "POS"))

fix_tidy_4_final_merge <- fix_tidy_filt_final %>% dplyr::select(ChromKey, CHROM, POS, REF, ALT)

gt_tidy_filt_final <- gt_tidy_filt2 %>% 
  inner_join(fix_tidy_4_final_merge, by = c("ChromKey", "POS")) #%>% 
#dplyr::select(ChromKey, CHROM, POS, Indiv, gt_AD, gt_DP, gt_HQ, gt_GL, gt_GQ, gt_GT, gt_GT_alleles, allele_A, allele_B, zygosity, REF, ALT)#LMK: When the LD_loci_to_keep line is commented out, this line gives an error that ChromKey isn't found; may need to be adjusted for ref based/no short LD filtering option

loci.1 <- nrow(fix_tidy_filt1)
indvs.1 <- nrow(keep_indv)

loci.final <- nrow(fix_tidy_filt_final)
indvs.final <- gt_tidy_filt_final %>% group_by(Indiv) %>% 
  summarize(n()) %>% 
  nrow()

stats.df.temp <- tibble(min_mean_depth = min.mean.depth,
                        max_mean_depth = max.mean.depth,
                        min_depth_indv = min_depth_indv,
                        Ipass = Ipass,
                        Imiss = Imiss,
                        Lmiss = Lmiss,
                        Loci_all_SNPs = loci.1,
                        Indvs_final = indvs.final,
                        Loci_one_SNP = loci.final)

stats.df <- bind_rows(stats.df, stats.df.temp)

write.xlsx(stats.df, file = "G://My Drive/Personal_Drive/Collaborations/Snails/Filter_test.xlsx")

#Free up memory from unused objects
gc()
  }
  }
    }
  }
}

stats.df <- stats.df %>% mutate(max_mean_depth = as.character(max_mean_depth)) %>% arrange(desc(Indvs_final, Loci_one_SNP))

write.xlsx(stats.df, file = "G://My Drive/Personal_Drive/Collaborations/Snails/Filter_test.xlsx")
