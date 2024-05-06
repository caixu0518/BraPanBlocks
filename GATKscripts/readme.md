# _A brief overview_:

First, we employed fastp (version 0.12.3) with the parameters "-z 4 -q 20 -u 30 -n 5" to filter out low-quality reads present in the sequencing data. Next, we aligned the clean reads from 1,812 samples separately to the reference genomes CFv4, Z1v2, ECD04, A03, PCE, LongYou, and OIB, generating gvcf files based on the clara-parabricks pipeline (version: 4.2.0-1, default parameters) (https://www.nvidia.com/en-us/clara/genomics/). Subsequently, we merged the gvcf files using the GenomicsDBImport module of the GATK software (version gatk-4.1.3.0, https://gatk.broadinstitute.org/hc/en-us) and then utilized the GenotypeGVCFs module within the GATK pipeline to generate vcf files. Furthermore, we extracted SNPs and InDels using the SelectVariants module of GATK. The final SNP and indel files were filtered for high quality and polymorphic SNPs and InDels across the 1,812 different B. rapa genomes using vcftools software (Danecek, et al., 2011) (version: 0.1.16; parameters: --maxDP 500, --minGQ 10, --minQ 30, --min-meanDP 3, --min-alleles 2, --max-alleles 2, --maf 0.05).

### The script below is used to convert reads from individual accessions into gVCF files using the clara-parabricks:4.2.1-1 package. The version used in this manuscript is 4.2.0-1 with default parameters.
```
#!/bin/bash
 
readId=$1

echo $(date)
echo "Start: ${readId} ..."

left=${readId}"_1.fq.ft.gz"
rigth=${readId}"_2.fq.ft.gz"

if [ -e ${left} ]; then
   echo "${left} exists "
else
  sshpass -p xxxxxx scp -P 9222  -r  TNAS@10.122.68.71:/Volume1/backup/00.WangLabData/all_Bra_resequencingData/oib_393/${left}  .   
fi

if [ -e ${rigth} ]; then
   echo "${rigth} exists "
else
  sshpass -p xxxxxx scp -P 9222 -r TNAS@10.122.68.71:/Volume1/backup/00.WangLabData/all_Bra_resequencingData/oib_393/${rigth}  . 
fi

##- 
docker run  \
       --gpus=all  \
       --rm \
       --volume $(pwd):/workdir \
       --volume $(pwd):/outputdir \
       --volume $(pwd):/tmpDir \
       nvcr.io/nvidia/clara/clara-parabricks:4.2.1-1   \
       pbrun fq2bam \
       --ref /workdir/OIB.fasta   \
       --in-fq /workdir/${left}  /workdir/${rigth} \
       --tmp-dir  /tmpDir  \
       --out-bam /workdir/${readId}.bam
      
docker run  \
       --gpus=all  \
       --rm \
       --volume $(pwd):/workdir \
       --volume $(pwd):/outputdir \
       --volume $(pwd):/tmpDir \
       nvcr.io/nvidia/clara/clara-parabricks:4.2.1-1  \
       pbrun haplotypecaller \
       --ref /workdir/OIB.fasta  \
       --in-bam /workdir/${readId}.bam \
       --tmp-dir /tmpDir \
        --gvcf \
       --out-variants /outputdir/${readId}.g.vcf.gz

rm -rf  /home/caix/OIB/${readId}.bam    /home/caix/OIB/${readId}.bam.bai  /home/caix/OIB/${readId}_chrs.txt /home/caix/OIB/${left}  /home/caix/OIB/${rigth} 

echo $(date)
echo "Finished: ${readId} ... "
```

#### _example: 1812 different _B. rapa_ assemblies were used along with 7 reference genomes to construct multiple reference variation maps_
<div align=center>
<img src="https://github.com/caixu0518/BraPanBlocks/blob/main/pngs/Figure-1.gif">
</div>
