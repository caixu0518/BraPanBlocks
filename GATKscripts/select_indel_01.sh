#!/bin/bash


genome=$1
vcf=$2
prefix=$3

gatk="/120t/caix/src/gatk-4.2.4.1/gatk"

if [ ! -d "indeltmp" ];then
    mkdir indeltmp
else
    rm -rf  indeltmp
fi


$gatk  --java-options "-Xmx200g -Djava.io.tmpdir=./indeltmp"   \
       SelectVariants  -R  ${genome}   \
       -V ${vcf}  \
       --select-type INDEL  \
       -O  ${prefix}.raw.indel.vcf.gz  

$gatk  --java-options "-Xmx200g -Djava.io.tmpdir=./indeltmp"   \
       VariantFiltration -R  ${genome}  \
       -V ${prefix}.raw.indel.vcf.gz  \
       --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 ||  ReadPosRankSum < -8.0"  \
       --filter-name 'INDEL_filter'  \
       -O  ${prefix}.filter.indel.vcf.gz

rm     ${prefix}.raw.indel.vcf.gz 

$gatk  --java-options "-Xmx200g -Djava.io.tmpdir=./indeltmp"  \
       SelectVariants  -R   ${genome}   \
       -V ${prefix}.filter.indel.vcf.gz  \
       --exclude-filtered   \
       -O ${prefix}.filtered.indel.vcf.gz

rm   ${prefix}.filter.indel.vcf.gz   indeltmp    -rf


