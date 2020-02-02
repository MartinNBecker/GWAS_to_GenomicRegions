# GWAS_to_GenomicRegions
Permutation test underlying "Early developmental gene enhancers affect subcortical volumes in the adult human brain." publication. [Hum Brain Mapp. 2016 May;37(5):1788-800. doi: 10.1002/hbm.23136. Epub 2016 Feb 18.]

1. Arrange GWAS data in a format that the first column gives the chromosome and nucleotide position, separated by a colon ("chromosome:nucleotide"). The second column should contain the p-value.

2. standard bed file format should be used for genomic regions

3. update genome_enrichment_Py3_300120.py with your input file names

4. Set number of permutations in genome_enrichment_Py3_300120.py and run

5. USe Histogram_empirical_pValue to visualize results and obtain final p-value
