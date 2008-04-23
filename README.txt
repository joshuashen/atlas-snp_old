*Quick start*

1.  ruby atlas-mapper-format-ref.rb  -r reference.fasta    ## format reference
2.1 ruby split-fasta-to-batches.rb -s input.fasta -p prefix_for_output -l length_for_each_batch -q  ## split the reads into smaller batches
2.2 ruby atlas-mapper.rb -q reads.fasta -r reference.fasta  ## map and align the reads (a batch) onto the reference
3.  ruby atlas-snp.rb -x cross_match_output -r reference.fasta -o prefix_of_outputs ## calling candidate SNPs
4.  ruby atlas-snp-evaluate.rb -i SNP.list -e overall_error_rate -s estimated_SNP_rate > SNP.list.eva  ## evaluate the error probability of each candidate SNP

*Requirement*
 - ruby
 - blat
 - cross_match
 - servers with large RAM

You can find blat binaries on Jim Kent's website: http://hgwdev.cse.ucsc.edu/~kent/exe/
And cross_match at: http://www.phrap.org/consed/consed.html#howToGet
After installing BLAT and cross_match, please make sure to put the path to these two programs in $PATH of your shell environment. 

-----------------------

*What's new*

Version 0.9.7, March 17th, 2008
 - Replace legacy perl scripts with ruby. 
 - Add README file in the package
 - Change the header rows of atlas-snp.rb and atlas-snp-evaluate.rb output

Version 0.9.6, March 08th, 2008
A few Incremental changes, including:

1. the -t option of atlas-mapper.rb. discarding reads with large tails or large insertions, since they are very likely to be from paralogs.
2. change default values of -s and -e of atlas-snp-evaluate.rb

Version 0.9.5, March 03rd, 2008
atlas-snp-evaluate.rb: SNP evaluation program based on logistic regression and Bayesian inference. The supervised learning procedure was trained on True/False SNP sets from Drosophila melanogaster 454-FLX runs. The logistic regression procedure on each read takes account of:
 quality scores
 longest nearby homopolymer size
 distance to 3' end
 "swapped" bases

P(error | all 454 substitutions on this base) is appended in the last column. For each read showing the "SNP", P(error | substitution) is appended in the last () field.

Version 0.9.4, February 29th, 2008
More information in the raw SNP output;
Retooled coverage calculation.

Version 0.9.3
Fix a bug in atlas-snp.rb output format

Version 0.9.2
1. atlas-mapper-format-ref.rb

Changes for limiting the memory usage of a single BLAT job to less than 2G. Record some meta information in a flat file. Note: re-running this step is a prerequisite of downstream steps.

2. atlas-mapper-do.rb is renamed to atlas-mapper.rb
Improvement on BLAT memory usage and cross_match output IO. Significant changes on how reads in a batch are distributed into smaller batches for doing cross_match against reference.

------------

*Protocol*

Code distribution: open access at http://code.google.com/p/atlas-snp/

Step 1. Creating the reference environment.

Command
    
    ruby atlas-mapper-format-ref.rb  -r reference.fasta [options]

Input and options
-r 	reference.fasta  required

This specify the file of the reference sequencing in fasta format, for example: 
>DDB0232429 |Chromosomal Sequence| on chromosome: 2 position 1 to 8470628
TTTTTTTTTTTTTTTTTTTTTTTTTATGTATGACACAATCATTAAATCATTACACATACC
AATTAGATTTTCTTTTTTTTTCTGATTTTAAAAACAAAAAAAAAACAAAAATTTATAAAT

-l 	length of pieces  optional, default value 100000
The program breaks the reference sequence into smaller pieces to avoid performance issues with cross_match alignment. 

-f 	frequency cutoff of 11-mers, optional, default 1024
The most important performance boost of BLAT comes from -ooc option, which enables the program to ignore over-represented kmers in the genome in the alignment seeding stage. 1024 is optimized for mammalian genomes. 100~200 is best for smaller or less complex genomes.
 
-b 	UNIX path to the BLAT program, optional

If BLAT is not accessible by default in the user's shell, the path to the blat program can be provided by this option. Note: it is best practice to add this path into one's PATH shell environment.

Output
This command creates a directory named by reference.fasta with a suffix "Env4mapping". It contains relevant information and data for subsequent blat and cross_match steps. The input reference fasta file is split into less than 900Mbps fragments to ensure running BLAT smoothly on servers with smaller than 2G RAM. The -ooc file is created according to command-line provided parameter. The reference sequence is also split into much smaller pieces, default at 100Kbps, to serve local alignment through cross_match.  Some meta-information was stored in a flat file inside the directory.

Computation requirement and performance
This step is fairly quick. The RAM usage depends on the size of the reference genome. It takes about 3.0G for the human genome. 

Step 2. Mapping and aligning the reads onto reference

2.1	(Optional) Split the reads into batches so that each batch can be executed on a cluster node in a typical parallel computing environment. The typical size of a batch is 5-10 Mbps.

ruby split-fasta-to-batches.rb -s input.fasta -p prefix_for_output -l length_for_each_batch -q 

Input and options
-s	original fasta file of the reads
-p 	prefix for the output read batches
-l 	the maximal length of a batch
-q 	optional with no argument; if provided, the program will also split the quality file correspondingly. 

Output
Batches of read fasta files, with prefix provided by -p option. 

2.2	 Run Atlas-mapper

ruby atlas-mapper.rb -q reads.fasta -r reference.fasta [options]

Input and options
Required:
-q 	reads.fasta, typically it is a fasta file containing a batch of reads.
-r 	reference.fasta it should be the same path of the reference as in Step 1.

Optional: 
-n 1 	default: none
This corresponds to the -oneOff option in BLAT, which allows one mismatch in the kmer "seeds", is useful to map reads that are from a different species or a more diverged strain to the genome of the reference organism
-i 	corresponds to the -minIdentity options in BLAT
-t 	min match ratio, default 0.85
After BLAT alignment, each read is assessed by dividing the matching length with read size, if the ratio is smaller than the cutoff value, the read is regarded as from a similar repeat region (or paralog) and thus discarded from following steps. 
-l 	masklevel option in cross_match, default 20
-s 	optimized for short reads, such as the ones from Solexa

Output
The output is a directory with prefix "Mapping_of" before name of reads.fasta. It contains the blat result in psl format and gzipped, blat best hits annotated as unique or repeat, and cross_match -discrep_lists output.

Note: This program assumes the user has access to BLAT and cross_match program by default. If BLAT or cross_match is not in the user's default shell $PATH, the path to them can be provided by -b and -x options respectively. 

Computation
(1) Here is a Bash shell command example to submit all batches of mapping jobs to an LSF cluster:

for f in batch_*.fa; do bsub -o lsf.o -e lsf.e -J $f "ruby atlas-mapper.rb -q $f -r ref.fa [ options ]"; done

After the computation on all batches are finished, all cross_match result from all batches should be concatenated into one file for Step 3.
(2) Memory usage: typically a batch takes less than 2G RAM.
(3) Performance: a batch of 10Mbps 454-FLX reads from human genome can be mapped to human reference in 3.3 hours on a Mac Pro Xeon 2.4GHz.  


2.2a (Alternative) For some project, the reference sequence is short enough for direct cross_match comparison without doing blat first. This command is:
    cross_match batch.fa reference.fa -minscore 30 -raw -discrep_lists -masklevel 50 > batch.fa.xm.ref

The result can be fed to Step 3 in the same way as the one from Step 2.2. 

Step 3. Calling SNPs

3.1 Calling candidate SNPs from cross_match result.

Note: if in Step 2 multiple batches of reads were mapped in parallel, the results from all batches should be concatenated into one file for following analysis. 

This program calls SNPs from the cross_match -discrep_list output produced from step 2:

ruby atlas-snp.rb -x cross_match_output -r reference.fasta -o prefix_of_outputs [options]


Input and options
-x 	cross_match result 
-r 	reference fasta file
-o 	prefix of output file name

-s 	max allowed substitution rate (in %), default 5
-g 	max allowed indel rate (in %), default 5

Output
It outputs a list of SNPs in the format:

refName coordinate refBase homopolymer refEnv coverage SNPBase adjustedQual oriQual numSNPReads numAlterReads reads_info

DDB0169550 317 A 5 ATTTATATTTTTA 92 T 16.48 19 1 1 T(19)089207_1601_3197(16)(134.0/140)-taaaatAataaat(1.43/0.0/1)swap;

U 6849 C 2 GGCAGACTTCAAG 12 T 37 37 1 1 T(37)086878_2543_1772(211)(237.0/243)+ggcagaTttcaag(0.82/0.0/1)snp;

refName		the name of the reference sequence, for example, chr12 
coordinate	the position of the SNP site on the reference sequence
refBase		the base of the reference on that SNP site
homopolymer	the size of the longest homopolymer within a 13-bp window centered on the SNP base on the reference
refEnv		the reference sequence of a 13-bp window centered on the SNP site
coverage	the number of reads covering the SNP base
SNPBase		the consensus SNP base
adjustedQual	(depreciated) the adjusted quality score of the SNP base
oriQual		the sum of the phred quality scores of all reads showing the consensus SNP base
numSNPReads	the number of reads that show the consensus SNP base
numAlterReads	the number of reads that differ from the reference on the SNP site
reads_info	information about the reads that differ from the reference on the SNP site

If more than one read showing the SNP base, the information of each read is printed in the reads_info field and separated by ';'.  For more help information and other options, type -h as the argument.

3.2 (Optional) Calculating depth-coverage on reference

Command
    ruby atlas-mapper-coverage.rb -x cross_match.output -o prefix_of_output [-r ref.fasta] [ -t targeted_genomic_regions ] [ options ]

Input and options
-r 	specifies the reference fasta file. If provided, the program will compute the depth-coverage for each base on the reference. NOTE this is time-consuming for large genomes. 
-t 	specifies a list of targeted genomic regions, with format:
    target_name	reference_name start_on_ref end_on_ref direction(1 or -1)

If -t is provided, the program will compute the depth-coverage only for bases in the targeted regions.  NOTE: one of -r and -t options is required.

Step 4. Evaluating the accuracy of candidate SNPs.

Command
    ruby atlas-snp-evaluate.rb -i SNP.list [ -e estimated_substitution_error_rate ] [ -s estimated_SNP_rate ] > output

Input and options
Required:
-i 	candidate SNP list from Step 3.1

Optional
-e	estimated overall substitution error rate. By default it is 0.0008, which is suitable for 454-FLX reads.
-s	estimated overall SNP rate. By default it is 0.001, which is reasonable for a typical human genome. For other species, the best way is to make an intelligent guess first, then apply the procedure and derive a better number iteratively. 

The most important output of this program is appended in the 14th column. It is an estimation of Pr(e|c), the probability of being substitution errors of the candidate SNP site. To collect a set of high-confidence SNPs, the user can set a cutoff value for the Pr(e|c). For an example, a cutoff of Pr(e|c) at 0.05 would put an upper bound of false positive rate at about 0.05. 



