#!/bin/bash

vcfIn=$1
popsam=$2
popindex=$3

vcftools   --gzvcf  ${vcfIn}   --recode  --recode-INFO-all  --stdout  --keep  ${popsam}    | bgzip  >  ${popindex}.snp.vcf.gz

tabix -p vcf  ${popindex}.snp.vcf.gz

plink  --vcf   ${popindex}.snp.vcf.gz   --recode  --allow-extra-chr  --out  ${popindex}.snp.tmp

plink  --file  ${popindex}.snp.tmp      --make-bed  --allow-extra-chr  --out  ${popindex}.snp.tmp

perl /120t/caix/src/PopLDdecay/bin/mis/plink2genotype.pl  -inPED  ${popindex}.snp.tmp.ped  -inMAP  ${popindex}.snp.tmp.map  -outGenotype  ${popindex}.snp.genotype

/120t/caix/src/PopLDdecay/PopLDdecay  -InGenotype   ${popindex}.snp.genotype    -OutStat  ${popindex}_LDdecay


