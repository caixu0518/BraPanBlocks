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
  sshpass -p Admin1111 scp -P 9222  -r  TNAS@10.122.68.71:/Volume1/backup/00.WangLabData/all_Bra_resequencingData/oib_393/${left}  .   
fi

if [ -e ${rigth} ]; then
   echo "${rigth} exists "
else
  sshpass -p Admin1111 scp -P 9222 -r TNAS@10.122.68.71:/Volume1/backup/00.WangLabData/all_Bra_resequencingData/oib_393/${rigth}  . 
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
