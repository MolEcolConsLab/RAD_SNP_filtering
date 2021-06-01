## ------------------------------------------------------------------------- ##
# --------------------  read VCFtools ind stats files ----------------------- #

library(dplyr)

# read individual stats files generated using vcftools
# specify data path (dir) and name of file to be loaded (vcf, no file ending)
# will load all het, idepth and imiss files and join into one data frame

#dir <- c("results")
#vcf <- c("SNAPPER_FIL-4")

read.ind.stats <- function(dir, vcf) {
  # read depth stats
  filename <- paste(vcf, ".idepth", sep = "")
  path <- file.path(dir, filename)
  idepth <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  
  # read missing stats
  filename <- paste(vcf, ".imiss", sep = "")
  path <- file.path(dir, filename)
  imiss <- read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(INDV, F_MISS)
  names(imiss) <- c("INDV", "MISS")
  # join stats
  temp <- left_join(imiss, idepth)
  
  # read missing stats
  filename <- paste(vcf, ".het", sep = "")
  path <- file.path(dir, filename)
  het <- read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(INDV, `F`)
  names(het) <- c("INDV", "Fis")
  # join stats
  final <- left_join(temp, het)
}

## ------------------------------------------------------------------------- ##



## ------------------------------------------------------------------------- ##
# --------------------  read VCFtools locus stats files ----------------------- #

library(dplyr)

# read locus stats files generated using vcftools
# specify data path (dir) and name of file to be loaded (vcf, no file ending)
# will load all het, idepth and imiss files and join into one data frame

read.loc.stats <- function(dir, vcf) {
  # read depth stats
  filename <- paste(vcf, ".ldepth.mean", sep = "")
  path <- file.path(dir, filename)
  ldepth <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  names(ldepth) <- c("CHR", "POS", "MEAN_DEPTH", "VAR_DEPTH")
  
  # read missing stats
  filename <- paste(vcf, ".lmiss", sep = "")
  path <- file.path(dir, filename)
  lmiss <- read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(CHR, POS, F_MISS)
  # join stats
  temp <- left_join(lmiss, ldepth, by = c("CHR", "POS"))
}

## ------------------------------------------------------------------------- ##