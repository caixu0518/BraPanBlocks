#!/bin/bash

vcf=$1
prefix=$2

gatk="/data/mg1/caix/src/NGSpipe/gatk"
vcftools="/data/mg1/caix/miniconda3/bin/vcftools"


##- select SNPs

snp=${prefix}.snp.gz
indel=${prefix}.indel.gz


$gatk   SelectVariants  \
        -V  ${vcf}  \
        -O  ${snp}  \
        --select-type-to-include SNP

$gatk   SelectVariants  \
        -V  ${vcf}  \
        -O  ${indel} \
        --select-type-to-include INDEL

$gatk   VariantFiltration \
        -V  ${snp}  \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${prefix}.snps_ft.vcf.gz


$gatk   VariantFiltration \
        -V  ${indel}   \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O ${prefix}.indels_ft.vcf.gz


$vcftools  \
        --maxDP 500 \
        --minGQ  10 \
        --minQ 30 \
        --min-meanDP 3 \
        --min-alleles 2 --max-alleles 2   \
        --keep-filtered PASS \
        --out   ${prefix}.snps_PASS_miss0.5.maf0.05.vcf  \
        --gzvcf ${prefix}.snps_ft.vcf.gz  \
        --recode --recode-INFO-all \
        --max-missing 0.5 \
        --maf 0.05


$vcftools  \
        --maxDP 500 \
        --minGQ  10 \
        --minQ 30 \
        --min-meanDP 3 \
        --min-alleles 2 --max-alleles 2  \
        --keep-filtered PASS \
        --out   ${prefix}.indels_PASS_miss0.5.maf0.05.vcf  \
        --gzvcf ${prefix}.indels_ft.vcf.gz   \
        --recode --recode-INFO-all \
        --max-missing 0.5 \
        --maf 0.05


bgzip   ${prefix}.snps_PASS_miss0.5.maf0.05.vcf.recode.vcf
tabix -p vcf  ${prefix}.snps_PASS_miss0.5.maf0.05.vcf.recode.vcf.gz

bgzip   ${prefix}.indels_PASS_miss0.5.maf0.05.vcf.recode.vcf
tabix -p vcf   ${prefix}.indels_PASS_miss0.5.maf0.05.vcf.recode.vcf.gz

rm ${prefix}.snps_ft.vcf.gz*  ${prefix}.indels_ft.vcf.gz*  ${snp}  ${indel}

