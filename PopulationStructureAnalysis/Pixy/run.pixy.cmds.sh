#!/bin/bash

groups=$1
prefix=$2

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/10t/caix/.aspera/connect/lib/ 


pixy  --stats pi fst dxy --vcf  CFv4.snps.maf0.05.vcf.gz   --populations ${groups}  --bed_file   CFv4.chr.w200k.s5k.bed  --n_cores 40  --bypass_invariant_check yes  --output_prefix  ${prefix}    &>   ${prefix}.pixy.log
