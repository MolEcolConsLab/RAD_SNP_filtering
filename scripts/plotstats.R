####---------Plot raw stats---------####
###Plot raw stats from vcf tools output
plot_stats_init <- function(ind, loc){
  # plot missing data per indv ----
  p1 <- ggplot(ind, aes(x = MISS)) +
    geom_histogram(binwidth = .01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MISS, na.rm = TRUE)), color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "missing data per indv")
  
  # plot Fis per indv ----
  p2 <- ggplot(ind, aes(x = Fis)) +
    geom_histogram(binwidth = .01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(Fis, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "Fis per indv")
  
  # plot read depth per indv ----
  p3 <- ggplot(ind, aes(x = MEAN_DEPTH)) +
    geom_histogram(binwidth = 10, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "mean read depth per indv")
  
  # plot depth vs missing ----
  p4 <- ggplot(ind, aes(x = MEAN_DEPTH, y = MISS)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = mean(MISS, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "mean depth per indv", y = "% missing data")
  
  
  # plot Fis vs mean depth per indv ----
  #May not need this graph. Not sure why it's helpful ...
  p5 <- ggplot(ind, aes(x = Fis, y = MEAN_DEPTH)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(Fis, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "Fis per indv", y = "mean depth per indv")


  # plot distribution missing data per locus ----
  p6 <- ggplot(loc, aes(x = F_MISS)) +
    geom_histogram(binwidth = 0.01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(F_MISS, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "% missing data per locus")
  
  # plot distribution mean read depth ----
  p7 <- ggplot(loc, aes(x = MEAN_DEPTH)) +
    geom_histogram(binwidth = 5, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "mean read depth per locus")
  
  
  # plot no of SNPs per locus ----
  p8 <- loc %>%
    count(CHR) %>%
    ggplot(aes(x = n)) +
    geom_histogram(binwidth = 1, color = "black", fill = "grey95") + 
    labs(x = "number of SNPs per locus")

  m1 <- multiplot(p1, p2, p3, p4, p5, p6, p7, p8, cols=2)
}


####---------------Plot stats after filtering-----------------####
#### Individual level graphs first
plot_ind_stats_filt <- function(keep_indv, keep_loci, mean_read_depth_indv, mean_read_depth_loci){
  # plot missing data per indv ----
  p1.1 <- ggplot(keep_indv, aes(x = percent_missing)) +
    geom_histogram(binwidth = 1, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(percent_missing, na.rm = TRUE), color = "darkblue"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = Imiss, color = "red"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = target_miss_per_indv, color = "lightblue"), #arbitrary target of 10
               linetype = "dashed", size = 1) +
    labs(x = "% missing data per indv") +
    scale_color_identity(name = "",
                         breaks = c("darkblue", "lightblue", "red"),
                         labels = c("mean", "target", "threshold"),
                         guide = "legend") +
    theme_standard
  
  # plot missing data by locus ----
  # plot distribution missing data per locus ----
  p1.2 <- ggplot(keep_loci, aes(x = percent_missing)) +
    geom_histogram(binwidth = 1, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(percent_missing, na.rm = TRUE), color = "darkblue"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = Lmiss, color = "red"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = target_miss_per_locus, color = "lightblue"),
               linetype = "dashed", size = 1) +
    labs(x = "% missing data per locus") +
    scale_color_identity(name = "",
                         breaks = c("darkblue", "lightblue", "red"),
                         labels = c("mean", "target", "threshold"),
                         guide = "legend") +
    theme_standard
  
  # plot read depth per indv ----
  p1.3 <- ggplot(mean_read_depth_indv, aes(x = mean_depth)) +
    geom_histogram(binwidth = 5, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(mean_depth, na.rm = TRUE), color = "darkblue"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = target_depth, color = "lightblue"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = min_mean_depth_lim, color = "red"),
               linetype = "dashed", size = 1) +
    labs(x = "mean read depth per indv") +
    scale_color_identity(name = "",
                         breaks = c("darkblue", "lightblue", "red"),
                         labels = c("mean", "target", "threshold"),
                         guide = "legend") +
    theme_standard
  
  # plot read depth per locus ----
  p1.4 <- ggplot(mean_read_depth_loci, aes(x = mean_depth)) +
    geom_histogram(binwidth = 1, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(mean_depth, na.rm = TRUE), color = "darkblue"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = target_depth, color = "lightblue"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = min_mean_depth_lim, color = "red"),
               linetype = "dashed", size = 1) +
    labs(x = "mean read depth per locus") +
    scale_color_identity(name = "",
                         breaks = c("darkblue", "lightblue", "red"),
                         labels = c("mean", "target", "threshold"),
                         guide = "legend") +
    theme_standard
  
  #Print individual level graphs
  ggarrange(p1.1, p1.2, p1.3, p1.4, common.legend = TRUE)
}


###Population-level graphs second
plot_pop_stats_filt <- function(keep_indv, mean_read_depth_indv){

  #Set color scheme for graphs
num_pops <- length(levels(factor(keep_indv$Pop)))
pop_cols <- brewer.pal(n = num_pops, name = "Spectral")

#Missing data by population
p1.1_pop <- ggplot(keep_indv, aes(x = factor(Pop))) +
  geom_boxplot(aes(y = percent_missing), fill = pop_cols) + 
  geom_hline(aes(yintercept = mean(percent_missing, na.rm = TRUE), color = "darkblue"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept = target_miss_per_indv, color = "lightblue"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept = Imiss, color = "red"), linetype = "dashed", size = 1) +
  labs(x = "%missing data per population") +
  scale_color_identity(name = "",
                       breaks = c
                       ("darkblue", "lightblue", "red"),
                       labels = c("mean", "target", "threshold"),
                       guide = "legend") +
  theme_standard

#Read coverage by population
p1.2_pop <- ggplot(mean_read_depth_indv, aes(x = factor(Pop))) +
  geom_boxplot(aes(y = mean_depth), fill = pop_cols) + 
  geom_hline(aes(yintercept = mean(mean_depth, na.rm = TRUE), color = "darkblue"), linetype = "dashed", size = 1) + #red line is mean
  geom_hline(aes(yintercept = target_depth, color = "lightblue"), linetype = "dashed", size = 1) + 
  geom_hline(aes(yintercept = min_mean_depth_lim, color = "red"), linetype = "dashed", size = 1) +
  #blue line is target coverage
  labs(x = "mean read depth by population") +
  scale_color_identity(name = "",
                       breaks = c("darkblue", "lightblue", "red"),
                       labels = c("mean", "target", "threshold"),
                       guide = "legend") +
  theme_standard

# Print population level graphs
ggarrange(p1.1_pop, p1.2_pop, common.legend = TRUE)
}


####--------------Plot final stats-----------------####

###Individual level first
# plot missing data per indv ----
plot_final_indv <- function(final_indv_stats, final_loc_stats, final_mean_depth_indv, final_mean_depth_loc){
  p2.1 <- ggplot(final_indv_stats, aes(x = percent_missing)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey95") +
  geom_vline(aes(xintercept = mean(percent_missing, na.rm = TRUE), color = "darkblue"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = Imiss, color = "red"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = target_miss_per_indv, color = "lightblue"), #arbitrary target of 10
             linetype = "dashed", size = 1) +
  labs(x = "% missing data per indv") +
  scale_color_identity(name = "",
                       breaks = c("darkblue", "lightblue", "red"),
                       labels = c("mean", "target", "threshold"),
                       guide = "legend")

# plot missing data by locus ----
# plot distribution missing data per locus ----
p2.2 <- ggplot(final_loc_stats, aes(x = percent_missing)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey95") +
  geom_vline(aes(xintercept = mean(percent_missing, na.rm = TRUE), color = "darkblue"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = Lmiss, color = "red"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = target_miss_per_locus, color = "lightblue"),
             linetype = "dashed", size = 1) +
  labs(x = "% missing data per locus") +
  scale_color_identity(name = "",
                       breaks = c("darkblue", "lightblue", "red"),
                       labels = c("mean", "target", "threshold"),
                       guide = "legend")

# plot read depth per indv ----
p2.3 <- ggplot(final_mean_depth_indv, aes(x = mean_depth)) +
  geom_histogram(binwidth = 5, color = "black", fill = "grey95") +
  geom_vline(aes(xintercept = mean(mean_depth, na.rm = TRUE), color = "darkblue"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = target_depth, color = "lightblue"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = min_mean_depth_lim, color = "red"),
             linetype = "dashed", size = 1) +
  labs(x = "mean read depth per indv") +
  scale_color_identity(name = "",
                       breaks = c("darkblue", "lightblue", "red"),
                       labels = c("mean", "target", "threshold"),
                       guide = "legend")

# plot read depth per locus ----
p2.4 <- ggplot(final_mean_depth_loc, aes(x = mean_depth)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey95") +
  geom_vline(aes(xintercept = mean(mean_depth, na.rm = TRUE), color = "darkblue"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = target_depth, color = "lightblue"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = min_mean_depth_lim, color = "red"),
             linetype = "dashed", size = 1) +
  labs(x = "mean read depth per locus") +
  scale_color_identity(name = "",
                       breaks = c("darkblue", "lightblue", "red"),
                       labels = c("mean", "target", "threshold"),
                       guide = "legend")

#Print individual level graphs
ggarrange(p2.1, p2.2, p2.3, p2.4, common.legend = TRUE)
}


#### Population level second
plot_final_pop <- function(final_indv_stats, final_mean_depth_indv){
  #Determine color scheme for graph
  num_pops2 <- length(levels(factor(final_indv_stats$Pop)))
  pop_cols2 <- brewer.pal(n = num_pops2, name = "Spectral")
  
  #Missing data by population
  p2.1_pop <- ggplot(final_indv_stats, aes(x = factor(Pop))) +
    geom_boxplot(aes(y = percent_missing), fill = pop_cols2) + 
    geom_hline(aes(yintercept = mean(percent_missing, na.rm = TRUE), color = "darkblue"), linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = target_miss_per_indv, color = "lightblue"), linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = Imiss, color = "red"), linetype = "dashed", size = 1) +
    labs(x = "%missing data per population") +
    scale_color_identity(name = "",
                         breaks = c("darkblue", "lightblue", "red"),
                         labels = c("mean", "target", "threshold"),
                         guide = "legend") +
    theme_standard
  
  #Read coverage by population
  p2.2_pop <- ggplot(final_mean_depth_indv, aes(x = factor(Pop))) +
    geom_boxplot(aes(y = mean_depth), fill = pop_cols2) + 
    geom_hline(aes(yintercept = mean(mean_depth, na.rm = TRUE), color = "darkblue"), linetype = "dashed", size = 1) + #red line is mean
    geom_hline(aes(yintercept = target_depth, color = "lightblue"), linetype = "dashed", size = 1) + 
    geom_hline(aes(yintercept = min_mean_depth_lim, color = "red"), linetype = "dashed", size = 1) +
    #blue line is target of 20x coverage
    labs(x = "mean read depth by population") +
    scale_color_identity(name = "",
                         breaks = c("darkblue", "lightblue", "red"),
                         labels = c("mean", "target", "threshold"),
                         guide = "legend") +
    #scale_fill_manual(values = pop_cols) +
    theme_standard
  
  # Print population level information
  ggarrange(p2.1_pop, p2.2_pop, common.legend = TRUE)
}

#### Final heterozygosity and SNPs per locus
plot_final_het_denovo <- function(final_SNPs_per_locus, final_obs_het_indv){
  # plot # SNPs per locus
  SNPs2 <- final_SNPs_per_locus %>% 
    ggplot(aes(x = SNPs)) +
    geom_histogram(binwidth = 1, color = "black", fill = "grey95") + 
    labs(x = "number of SNPs per locus") +
    theme_standard
  
  #Vizualize heterozygosity
  pop_cols2 <- brewer.pal(n = 5, name = "Dark2")
  
  hetero2 <- ggplot(final_obs_het_indv, aes(x = factor(Pop), y = Obs_het, color = Pop)) + 
    geom_jitter() + 
    scale_color_manual(values = pop_cols2) + 
    labs(title = "Individual heterozygosity", x = "Population", y = "Obs het") +
    theme_facet
  
  ggarrange(SNPs2, hetero2)
}


plot_final_het_reference <- function(final_obs_het_indv){
  #Vizualize heterozygosity
  pop_cols2 <- brewer.pal(n = 5, name = "Dark2")
  
  ggplot(final_obs_het_indv, aes(x = factor(Pop), y = Obs_het, color = Pop)) + 
    geom_jitter() + 
    scale_color_manual(values = pop_cols2) + 
    labs(title = "Individual heterozygosity", x = "Population", y = "Obs het") +
    theme_facet
}