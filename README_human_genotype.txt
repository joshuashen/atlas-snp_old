Documentation 

Author: Zhengzheng Wan (zwan@bcm.edu)


1. The component specific for HapMap SNP calling


Introduction
This component provides a new optional functionality to evaluate the accuracy of candidate SNPs in the human genomic sequences for sites that have been extensively studied by the International HapMap Project (HapMap). Once enabled, comparing with traditional "atlas-snp-evaluate.rb", the prior probability of being a real SNP [P(SNP)] is determined based on the population frequencies of the particular site from HapMap, instead of using the default value. The purpose of doing this is to improve the sensitivity in calling the HapMap SNPs, without significant compromising on the specificity. 

Users can enable this component by using options '-f', and '-d' (followed by the path where the allele frequency information of the HapMap SNPs is stored).   

Command 
ruby atlas-snp-evaluate-hapmapspecific.rb -i SNP.list  [-e estimated_substition_error_rate][-s estimated prior of P(SNP)][-f option to enable the component][-d directory for hapmap][-c comparing Pr(e|c) produced by using and without using this component][other options] > output.

Input and options
-f	Used to turn on this entire component.

-d	Must be appended once '-f' is selected. It should be followed by the path where the allele frequency information for HapMap SNPs is stored. 

-c	An optional argument. It must work together with '-f' and '-d', otherwise won't be effective.  When '-c' is appended, P(SNP) and Pr(e|c) produced by using the default (i.e. without '-f') will be appended as the last two columns of the output for comparison purpose.  

-e	Estimated substitution error rate

-s	Estimated SNP rate

Outputs
The general output of this program has the same format as that produced by using the version without '-f', as the previous version of "Atlas-snp".  The overall Pr(e|c) of a SNP site that summarizes the Pr(e|c) values from each individual read (shown in the 13th column), and the Pr(e|c) for each SNP read (shown in the12th column), are different from the output produced by using the version without '-f'. 

Particularly, if '-c' is used, there will be two extra columns appended, which will be shown as "--" otherwise.   

For example:


Header line
refName coordinate refBase homopolymer refEnv coverage SNPBase adjustedQual oriQual numVariant numAlterReads reads_info less_strand Pr_e_c strandP snpRate hapmap_status default_SNP_rate old Pr_e_c

When using '-c'
chr12 38829960 c 4 tgagctcattttc 62 T 26 26 1 1 T(26)189958_3861_1871(138)(223.0/256)+tgagctTattttc(1.27/0.42/21)snp(5.34)(0.00477414616232097); 0 0.00477414616232097 0.5 0.542 known 0.001 0.722222222222222      

When not using '-c'
chr12 38829960 c 4 tgagctcattttc 62 T 26 26 1 1 T(26)189958_3861_1871(138)(223.0/256)+tgagctTattttc(1.27/0.42/21)snp(5.34)(0.00477414616232097); 0 0.00477414616232097 0.5 0.542 known -- --

2. Determine genotypes using p-values based on binomial exact test


Introduction
This is a stand-alone program to determine the genotypes based on the p-values from Binomial exact test. 

Command
ruby genotype-determination.rb -i SNP.list [-p binomial p value cutoff]

Input and options  
-I	Followed by the SNP.list and the full path for the file.

-p	Once '-p' is selected, it  must be appended with the cutoff for binomial p-value. Otherwise, the default p-value cutoff is 0.001.

Outputs 
This command creates the output named by the input file name with a suffix ".pvalue.gt". The output contains two extra columns--the binomial p-value for a SNP site, and the genotype of this SNP site, as the last two columns.

