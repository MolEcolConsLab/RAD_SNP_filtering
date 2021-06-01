#!/bin/bash
#BSUB -J 09_filter_VCF_SNP
#BSUB -o ./logs/09_filter_VCF_SNP.log
#BSUB -e ./logs/09_filter_VCF_SNP.err
#BSUB -n 8
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=5000]
#BSUB -W 05:00
#BSUB -q interactive
#BSUB -Is singularity

module load vcftools/0.1.16

VCFTOOLSIMG=/share/pkg/vcftools/0.1.16/vcftools-0.1.16.sif

input_vcf=./08_FINAL_populations_out_m3_M3/populations.snps.vcf
fileprefix=CNR_raw

##############################################
### Filter VCF for use in pop gen analyses ###
##############################################

#--min-meanDP and --max-meanDP filter by the mean minimum and maximum depth values across all included individuals.
#--maf is minor allele frequency
#--hwe - check for hwe and remove sites with a p value below this threshold
#--max-missing removes sites on the proportion of missing data, where 0 means sites that are totally missing are fine, and 1 means no missing data allowed.
#--recode outputs a new vcf file after filtering

#First filtering step
singularity exec $VCFTOOLSIMG vcftools --vcf $input_vcf --max-missing 0.5 --min-meanDP 10 --max-meanDP 1000 --recode --mac 3 --out ./09_vcftools_out/${fileprefix}

#Second step calculates data missingness 
## Per individual ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --missing-indv

## Per site ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --missing-site


#Fourth step calculates deviation from HWE
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --hardy


#Fifth step calculates depth per individual and locus
## Mean depth per individual ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --depth

## mean depth per site averaged across individuals ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --site-mean-depth

## Depth for each genotype in the VCF file ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --geno-depth


# Sixth step calculates heterozygosity per individual
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --het


# Seventh step grabs SNP call quality
singularity exec $VCFTOOLSIMG vcftools --vcf ./09_vcftools_out/${fileprefix}.recode.vcf --out ./09_vcftools_out/${fileprefix} --site-quality