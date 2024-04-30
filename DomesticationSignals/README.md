
## The main script functions are briefly introduced as follows:

### generate_interval_vcfs.pl
We used the script generate_interval_vcfs.pl to segment VCF files based on the coordinates of Pan-Blocks, enabling parallel computation to reduce processing time. This script splits VCF files into different parts according to the positions defined by Pan-Blocks, facilitating efficient parallelized analysis.

### generate_allele_frequency_and_Pi_v1.1.pl
The script generate_allele_frequency_and_Pi_v1.1.pl is used to calculate the frequency of the reference allele and missing allele at each SNP site within each VCF file, as well as the original π (Pi) value.

### generate_allele_frequency_and_Pi_v1.1.pl
The script generate_allele_frequency_and_Pi_v1.1.pl performs a Fisher's exact test and adjusts the p-values using the Benjamini-Hochberg method. If the condition (mean of dominant allele frequency - 3 * standard deviation of dominant allele frequency) > (mean of recessive allele frequency + 3 * standard deviation of recessive allele frequency) is met and the adjusted p-value <sub>(Padjusted)</sub> is less than 1e-5, then the window is further considered as a putative selective signal.

### CalculateDomesticatedSignalsBasedOnPanBlocks_v1.6.pl
The script CalculateDomesticatedSignalsBasedOnPanBlocks_v1.6.pl calls generate_allele_frequency_and_Pi_v1.1.pl and generate_allele_frequency_and_Pi_v1.1.pl, integrating these two steps into an automated workflow. It also utilizes the Perl threads module to enable parallel computation for efficient processing.

### extract_candiate_genes_v1.3_3xigama.pl
The script extract_candidate_genes_v1.3_3xigama.pl extracts candidate genes present in the selection signals and adds potential functional descriptions to these candidate genes using a provided list of syntenic genes with _Arabidopsis_.
