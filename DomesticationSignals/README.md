# A brief overview:

After obtaining population variations using seven different _B. rapa_ genomes as reference genomes, we further integrated the Pan-Blocks and genetic variation maps constructed with these reference genomes to identify potential selection signals. The steps were conducted as follows: Initially, we partitioned each Pan-Block into 10 kb bins and computed the enrichment levels of reference alleles and missing alleles within both the derived and control groups. A bin was designated as a potential selected region if all SNPs within the bin showed significant enrichment of both alleles in either the derived or control populations. Specifically, we assessed the frequencies of both alleles in these populations and discarded sites where both dominant and recessive allele frequencies were either greater than 0.9 or less than 0.1. Subsequently, we determined the mean (μ<sub>dominant</sub>, μ<sub>recessive</sub>) and standard deviation (σ<sub>dominant</sub>, σ<sub>recessive</sub>) of allele frequencies for each window in these populations. Following this, we performed _Fisher's exact test_ and adjusted the p-values using the Benjamini-Hochberg method. If (μ<sub>dominant</sub>-3σ<sub>dominant</sub>) > (μ<sub>recessive</sub>+3σ<sub>recessive</sub>) and Padjusted < 1e-5, the window was further considered as a putative selective signal.


## The main script functions are briefly introduced as follows:

### _generate_interval_vcfs.pl_
We used the script generate_interval_vcfs.pl to segment VCF files based on the coordinates of Pan-Blocks, enabling parallel computation to reduce processing time. This script splits VCF files into different parts according to the positions defined by Pan-Blocks, facilitating efficient parallelized analysis.

### _generate_allele_frequency_and_Pi_v1.1.pl_
The script generate_allele_frequency_and_Pi_v1.1.pl is used to calculate the frequency of the reference allele and missing allele at each SNP site within each VCF file, as well as the original π (Pi) value.

### _generate_allele_frequency_and_Pi_v1.1.pl_
The script generate_allele_frequency_and_Pi_v1.1.pl performs a _Fisher's exact test_ and adjusts the p-values using the Benjamini-Hochberg method. If the condition (mean of dominant allele frequency - 3 * standard deviation of dominant allele frequency) > (mean of recessive allele frequency + 3 * standard deviation of recessive allele frequency) is met and the adjusted p-value <sub>(Padjusted)</sub> is less than 1e-5, then the window is further considered as a putative selective signal.

### _CalculateDomesticatedSignalsBasedOnPanBlocks_v1.6.pl_
The script CalculateDomesticatedSignalsBasedOnPanBlocks_v1.6.pl calls generate_allele_frequency_and_Pi_v1.1.pl and generate_allele_frequency_and_Pi_v1.1.pl, integrating these two steps into an automated workflow. It also utilizes the Perl threads module to enable parallel computation for efficient processing. 
The resulting files generated from this step are [Merged.group.lst.RefAllele.adjustPvalue.txt](https://github.com/caixu0518/BraPanBlocks/blob/main/DomesticationSignals/Data/Merged.group.lst.RefAllele.adjustPvalue.txt) and [Merged.group.lst.MissAllele.adjustPvalue.txt](https://github.com/caixu0518/BraPanBlocks/blob/main/DomesticationSignals/Data/Merged.group.lst.MissAllele.adjustPvalue.txt).

### _extract_candiate_genes_v1.3_3xigama.pl_
The script extract_candidate_genes_v1.3_3xigama.pl extracts candidate genes present in the selection signals and adds potential functional descriptions to these candidate genes using a provided list of syntenic genes with _Arabidopsis_. The resulting file generated from this step is [merged.outliers.signals.txt](https://github.com/caixu0518/BraPanBlocks/blob/main/DomesticationSignals/Data/merged.outliers.signals.txt).

#### _example: Domestication Signals in Southern and Northern East Asian B. rapa Populations Using Pan-Blocks_
<div align=center>
<img src="https://github.com/caixu0518/BraPanBlocks/blob/main/pngs/signal_example.gif">
</div>
