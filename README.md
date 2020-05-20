# methylMnM
methyMnM was developed for jointly analyzing MeDIP-seq and MRE-seq data, which were derived from methylated DNA immunoprecipitation (MeDIP) experiments [Weber et al., 2005] followed by sequencing (MeDIP-seq) and methyl-sensitive restriction enzymes experiments for unmethylated CpGs (MRE-seq). We have implemented the methylMnM method via a set of R functions with the computational intensive parts written in C. The method consists three steps:

Data Pre-processing:
Calculate the CpG count of each window.
Calculate the MRE CpG count of each window.
Calculate MeDIP-seq tag count of each window of control and treatment samples.
Calculate MRE-seq tag count of each window of control and treatment samples.
Calculating p-values of each window by the methylMnM test.
FDR control.

Link to MnM homepage: https://epigenome.wustl.edu/MnM/
Link to biocoductor: https://www.bioconductor.org/packages/release/bioc/html/methylMnM.html 

Citation: Functional DNA methylation differences between tissues, cell types, and across individuals discovered using the M&M algorithm. B Zhang, Y Zhou, N Lin, RF Lowdon, C Hong… - Genome research, 2013
