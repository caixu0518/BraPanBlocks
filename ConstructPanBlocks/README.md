# A brief overview:

The Pan-Blocks constructed in the current study are based on 21 published B. rapa assemblies. To construct the B. rapa Pan-Blocks, we used each genome as a reference genome and aligned the remaining genomes to it (the order used in the current study is CFv4, Z1v2, ECD04, A03, PCE, LongYou, OIB, TUE, BRO, CXB, CXA, PCB, PCA, PCD, PCC, TUB, MIZ, TUA, TCA, CCB). Those syntenic blocks present in the reference genome were extracted as part of the Pan-Blocks. Specifically, we first used the CFv4 genome as a reference and aligned the sequences of other genomes to CFv4 using the nucmer (version: 4.0.0beta2, -t 20 --mum). Next, we utilized delta-filter (parameters: -l 10000 -r -q) to extract "one-to-one alignments" on the same chromosome (defined as syntenic blocks in the current study). Following this, potential inversions were extracted from these one-to-one alignments based on the alignment coordinates between any two genomes. Finally, we used show-coords (parameters: -TrHcl) to identify and merge all syntenic blocks present on the reference genome. These syntenic blocks were then added to the Pan-Blocks, and the sequences contributed by CFv4 were labeled accordingly. Once the alignment with CFv4 as the reference genome was completed, CFv4 was removed from the list of genomes for comparison. In the second round of alignment, we used Z1v2 as the reference genome. Using the same method, we identified syntenic blocks present in the Z1v2 genome and added them to the B. rapa Pan-Blocks, labeling them as contributed by Z1v2. It's important to note that blocks in Z1v2 that were syntenic with CFv4 were excluded, as these sequences are already present in the B. rapa Pan-Blocks. Continuing in the same manner, we used each of the remaining genomes as references in succession. This iterative process enabled us to complete the construction of the B. rapa Pan-Blocks.

### _run_nucmers.pl_
The script run_nucmers.pl is used to compare any two genomes and obtain one-to-one alignments.

### calculate_shared_haps.v1.pl
The script calculate_shared_haps.v1.pl is used to obtain the coordinates of _B. rapa_ Pan-Blocks.

### 03.plot_multi_syn.v1.1.pl
The script 03.plot_multi_syn.v1.1.pl uses the SVG module in Perl to generate plots of Pan-Blocks and their one-to-one alignments between different genomes. The results can be visualized as shown in the following example.

#### example: Display of Pan-Blocks on chromosome A03
<div align=center>
<img src="https://github.com/caixu0518/BraPanBlocks/blob/main/pngs/A03_Pan-Blocks.gif">
</div>
