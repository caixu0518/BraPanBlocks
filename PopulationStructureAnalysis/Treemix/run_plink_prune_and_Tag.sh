#!/bin/bash 

vcfIn=$1
outprefix=$2

plink  --vcf  ${vcfIn} --recode    --allow-extra-chr  --out  ${outprefix}  --vcf-half-call  haploid 

plink  --file  ${outprefix}  --make-bed  --allow-extra-chr  --out    ${outprefix}

plink  --file  ${outprefix} --indep-pairwise  50 10 0.2  --out  ${outprefix}.ftLD  --geno  0.1  --allow-extra-chr

plink  --file  ${outprefix}  --blocks  no-pheno-req   --allow-extra-chr     

