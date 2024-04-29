#!/bin/bash

vcf=$1
prefix=$2

vcftools="/120t/caix/anaconda3/bin/vcftools"

$vcftools  \
         --maxDP 500 \
         --minGQ  10 \
         --minQ 30 \
         --min-meanDP 3 \
         --min-alleles 2 --max-alleles 2   \
         --keep-filtered PASS \
         --out   ${prefix}.snps.maf0.05.vcf  \
         --gzvcf ${vcf}  \
         --recode --recode-INFO-all \
         --maf 0.05
bgzip ${prefix}.snps.maf0.05.vcf.recode.vcf
tabix -p vcf ${prefix}.snps.maf0.05.vcf.recode.vcf.gz


$vcftools  \
         --maxDP 500 \
         --minGQ  10 \
         --minQ 30 \
         --min-meanDP 3 \
         --min-alleles 2 --max-alleles 2   \
         --keep-filtered PASS \
         --out   ${prefix}.snps.miss0.5.maf0.05.vcf.recode.vcf  \
         --gzvcf ${vcf}  \
         --recode --recode-INFO-all \
         --max-missing 0.5 \
         --maf 0.05
bgzip  ${prefix}.snps.miss0.5.maf0.05.vcf.recode.vcf
tabix -p vcf ${prefix}.snps.miss0.5.maf0.05.vcf.recode.vcf.gz
