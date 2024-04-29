#!/usr/bin/bash 

vcfIn=$1
prefix=$2

plink  --vcf  ${vcfIn}  --recode --allow-extra-chr --out ${prefix}

plink  --file ${prefix}  --make-bed  --allow-extra-chr  --out  ${prefix}

plink --bfile  ${prefix}   --recode vcf-iid  --out  ${prefix}   --allow-extra-chr 


bgzip ${prefix}.vcf
tabix -p vcf ${prefix}.vcf.gz

##- clean
rm ${prefix}.bed  ${prefix}.bim  ${prefix}.fam  ${prefix}.log  ${prefix}.map  ${prefix}.nosex  ${prefix}.ped
