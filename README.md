# Scripts related to the Brassica rapa Pan-Blocks manuscript.

The current directory contains scripts corresponding to the analysis pipelines used in the Brassica rapa Pan-Blocks manuscript. This includes scripts for constructing the _B. rapa_ Pan-Blocks (stored in the [ConstructPanBlocks](https://github.com/caixu0518/BraPanBlocks/tree/main/ConstructPanBlocks) directory), domestication signal identification pipeline (stored in the [DomesticationSignals](https://github.com/caixu0518/BraPanBlocks/tree/main/DomesticationSignals) directory), population analysis workflows (stored in the [PopulationStructureAnalysis](https://github.com/caixu0518/BraPanBlocks/tree/main/PopulationStructureAnalysis) directory), and variant detection scripts (stored in the [GATKscripts](https://github.com/caixu0518/BraPanBlocks/tree/main/GATKscripts) directory). Most of the scripts are written in Perl and incorporate calls to previously published software. These details will be referenced in our forthcoming manuscript.


## ConstructPanBlocks
The Pan-Blocks constructed in the current study are based on 21 published _B. rapa_ assemblies. To construct the _B. rapa_ Pan-Blocks, we used each genome as a reference genome and aligned the remaining genomes to it (the order used in the current study is CFv4, Z1v2, ECD04, A03, PCE, LongYou, OIB, TUE, BRO, CXB, CXA, PCB, PCA, PCD, PCC, TUB, MIZ, TUA, TCA, CCB). Those syntenic blocks present in the reference genome were extracted as part of the Pan-Blocks. Specifically, we first used the CFv4 genome as a reference and aligned the sequences of other genomes to CFv4 using the nucmer (version: 4.0.0beta2, -t 20 --mum). Next, we utilized delta-filter (parameters: -l 10000 -r -q) to extract "one-to-one alignments" on the same chromosome (defined as syntenic blocks in the current study). Following this, potential inversions were extracted from these one-to-one alignments based on the alignment coordinates between any two genomes. Finally, we used show-coords (parameters: -TrHcl) to identify and merge all syntenic blocks present on the reference genome. These syntenic blocks were then added to the Pan-Blocks, and the sequences contributed by CFv4 were labeled accordingly. Once the alignment with CFv4 as the reference genome was completed, CFv4 was removed from the list of genomes for comparison. In the second round of alignment, we used Z1v2 as the reference genome. Using the same method, we identified syntenic blocks present in the Z1v2 genome and added them to the _B. rapa_ Pan-Blocks, labeling them as contributed by Z1v2. It's important to note that blocks in Z1v2 that were syntenic with CFv4 were excluded, as these sequences are already present in the _B. rapa_ Pan-Blocks. Continuing in the same manner, we used each of the remaining genomes as references in succession. This iterative process enabled us to complete the construction of the _B. rapa_ Pan-Blocks.

#### example: Display of Pan-Blocks on chromosome A03

<div align=center>
<img src="https://github.com/caixu0518/BraPanBlocks/blob/main/pngs/A03_Pan-Blocks.gif" width="180" height="105"> width="180" height="105"/>
</div>

```
```
## DomesticationSignals
After obtaining population variations using seven different _B. rapa_ genomes as reference genomes, we further integrated the Pan-Blocks and genetic variation maps constructed with these reference genomes to identify potential selection signals. The steps were conducted as follows: Initially, we partitioned each Pan-Block into 10 kb bins and computed the enrichment levels of reference alleles and missing alleles within both the derived and control groups. A bin was designated as a potential selected region if all SNPs within the bin showed significant enrichment of both alleles in either the derived or control populations. Specifically, we assessed the frequencies of both alleles in these populations and discarded sites where both dominant and recessive allele frequencies were either greater than 0.9 or less than 0.1. Subsequently, we determined the mean (μdominant, μrecessive) and standard deviation (σdominant, σrecessive) of allele frequencies for each window in these populations. Following this, we performed _Fisher's exact test_ and adjusted the p-values using the Benjamini-Hochberg method. If (μdominant-3σdominant) > (μrecessive+3σrecessive) and Padjusted < 1e-5, the window was further considered as a putative selective signal.


#### example: Domestication Signals in Southern and Northern East Asian B. rapa Populations Using Pan-Blocks.
![image](https://github.com/caixu0518/BraPanBlocks/blob/main/pngs/signal_example.gif)
```
```
## PopulationStructureAnalysis
PCA analysis was performed using the plink software (version: v1.90b6.24, https://www.cog-genomics.org/plink2/) with the parameters "plink --noweb --bfile --pca 20 --allow-extra-chr". Fst and nucleotide diversity were calculated using pixy (version: 1.2.7.beta1), with a window size set to 200 kb, step size set to 5 kb, and parameters "pixy --stats pi fst dxy --vcf --populations --bed_file --n_cores 40 --bypass_invariant_check yes". Linkage disequilibrium decay analysis was conducted using PopLDdecay (https://github.com/BGI-shenzhen/PopLDdecay) with default parameters. Population structure analysis was conducted using the fastStructure algorithm. Estimating the size history of populations from whole-genome sequence data in B. rapa was conducted using the SMC++ program using the recommended pipeline (https://github.com/popgenmethods/smcpp). F4-ratio statistics were conducted using the Dsuite package (Malinsky, et al., 2021). TreeMix diagram for different B. rapa populations was conducted using TreeMix v. 1.12 with parameters "-root Rapini -k 500 -m 10 -bootstrap" . Finally, the Newick tree was generated using the PanKmer program (https://gitlab.com/salk-tm/pankmer) with the parameter "--newick --metric".

```

```
## GATKscripts
First, we employed fastp (version 0.12.3) with the parameters "-z 4 -q 20 -u 30 -n 5" to filter out low-quality reads present in the sequencing data. Next, we aligned the clean reads from 1,812 samples separately to the reference genomes CFv4, Z1v2, ECD04, A03, PCE, LongYou, and OIB, generating gvcf files based on the clara-parabricks pipeline (version: 4.2.0-1, default parameters) (https://www.nvidia.com/en-us/clara/genomics/). Subsequently, we merged the gvcf files using the GenomicsDBImport module of the GATK software (version gatk-4.1.3.0, https://gatk.broadinstitute.org/hc/en-us) and then utilized the GenotypeGVCFs module within the GATK pipeline to generate vcf files. Furthermore, we extracted SNPs and InDels using the SelectVariants module of GATK. The final SNP and indel files were filtered for high quality and polymorphic SNPs and InDels across the 1,812 different B. rapa genomes using vcftools software (Danecek, et al., 2011) (version: 0.1.16; parameters: --maxDP 500, --minGQ 10, --minQ 30, --min-meanDP 3, --min-alleles 2, --max-alleles 2, --maf 0.05). 


```


```
