# A brief overview:

PCA analysis was performed using the plink software (version: v1.90b6.24, https://www.cog-genomics.org/plink2/) with the parameters "plink --noweb --bfile --pca 20 --allow-extra-chr". Fst and nucleotide diversity were calculated using pixy (version: 1.2.7.beta1), with a window size set to 200 kb, step size set to 5 kb, and parameters "pixy --stats pi fst dxy --vcf --populations --bed_file --n_cores 40 --bypass_invariant_check yes". Population structure analysis was conducted using the fastStructure algorithm. F4-ratio statistics were conducted using the Dsuite package (Malinsky, et al., 2021). The Newick tree was generated using the PanKmer program (https://gitlab.com/salk-tm/pankmer) with the parameter "--newick --metric".

#### _example: population analysis results conducted in this study_
<div align=center>
<img src="https://github.com/caixu0518/BraPanBlocks/blob/main/pngs/population_analysis.png">
</div>
