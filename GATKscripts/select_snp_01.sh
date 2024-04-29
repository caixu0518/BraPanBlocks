#!/bin/bash


genome=$1
vcf=$2
prefix=$3

gatk="/120t/caix/src/gatk-4.2.4.1/gatk"


if [ ! -d "snptmp" ];then
    mkdir snptmp
else
    rm -rf  snptmp
fi


$gatk  --java-options "-Xmx200g -Djava.io.tmpdir=./snptmp"  \
       SelectVariants \
       -R ${genome}  \
       -V ${vcf}  \
       --select-type SNP    \
       -O ${prefix}.raw.snp.vcf.gz  

$gatk  --java-options "-Xmx200g -Djava.io.tmpdir=./snptmp"   \
      VariantFiltration    \
      -R ${genome}     \
      -V ${prefix}.raw.snp.vcf.gz   \
      --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter'  \
      -O ${prefix}.filter.snp.vcf.gz 

rm ${prefix}.raw.snp.vcf.gz

$gatk  --java-options "-Xmx100g -Djava.io.tmpdir=./snptmp"  \
      SelectVariants   \
      -R  ${genome}     \
      -V  ${prefix}.filter.snp.vcf.gz --exclude-filtered   \
      -O  ${prefix}.filtered.snp.vcf.gz  

rm    ${prefix}.filter.snp.vcf.gz    snptmp   -rf
