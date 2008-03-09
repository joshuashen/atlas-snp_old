*What's new*

Version 0.9.6
A few Incremental changes, including:

1. the -t option of atlas-mapper.rb. discarding reads with large tails or large insertions, since they are very likely to be from paralogs.
2. change default values of -s and -e of atlas-snp-evaluate.rb

Version 0.9.5
atlas-snp-evaluate.rb: SNP evaluation program based on logistic regression and Bayesian inference. The supervised learning procedure was trained on True/False SNP sets from Drosophila melanogaster 454-FLX runs. The logistic regression procedure on each read takes account of:

 quality scores
 longest nearby homopolymer size
 distance to 3' end
 "swapped" bases

P(error | all 454 substitutions on this base) is appended in the last column. For each read showing the "SNP", P(error | substitution) is appended in the last () field.

Version 0.9.4
More information in the raw SNP output;
Retooled coverage calculation.

Version 0.9.3
Fix a bug in atlas-snp.rb output format

Version 0.9.2
1. atlas-mapper-format-ref.rb

Changes for limiting the memory usage of a single BLAT job to less than 2G. Record some meta information in a flat file. Note: re-running this step is a prerequisite of downstream steps.
2. atlas-mapper-do.rb is renamed to atlas-mapper.rb

Improvement on BLAT memory usage and cross_match output IO. Significant changes on how reads in a batch are distributed into smaller batches for doing cross_match against reference.


----


*Requirement*

ruby
blat
cross_match
large RAM

You can find blat binaries on Jim Kent's website: http://hgwdev.cse.ucsc.edu/~kent/exe/

For cross_match, http://www.phrap.org/consed/consed.html#howToGet

*Protocol*

The procedure has three steps:

1. Mapping and aligning 454 reads onto the reference genome.

2. Calling raw (a) SNPs and (b) indels based on results from step 1.

3. SNP/indel filtering based on quality scores and biology.

I have written a set of streamlined programs to do step 1 and step 2a. The step 2b (indels) is in the working. Step 3 is straightforward in terms of programming but requires some discussions on the biology.

1. Atlas-mapper

Map the reads (454 or Sanger) to reference sequence via BLAT and do refined sequence local alignment via cross_match.

1) Format the reference sequence (including making blat ooc file, splitting reference into smaller pieces for cross_match):

ruby /PATH/TO/atlas-mapper-format-ref.rb -r reference.fasta
Note: this program requires about 3.0G RAM for the human genome. The reference.fasta is the reference sequences in one fasta file. The output is a directory named by reference.fasta with a suffix "Env4mapping".

2) Map and align: (including doing blat, picking best hits and uniquely mapped reads, doing cross_match etc)

ruby /PATH/TO/atlas-mapper.rb -q query.fasta -r reference.fasta [ options ]
query.fasta is a fasta file containing a batch of 454 reads.

The output is a directory with prefix "Mapping_of". It contains the blat result, blat best hits annotated as unique or repeats, and cross_match -discrep_lists output.

Note the query.fasta should be in reasonable size. It's best to be no larger than 10Mbp. A typical 454 run produce 100+Mbp sequence. So it is necessary to split the run into smaller batches and do each batch on one cluster node.

Here is a Bash shell command example to submit all batches of mapping jobs to a cluster:

for f in batch_*.fa; do bsub -o lsf.o -e lsf.e -J $f "ruby ..../atlas-mapper.rb -q $f -r ref.fa [ options ]"; done
2. Atlas-SNP

This program calls SNPs from the cross_match -discrep_list output produced from step 1:

ruby /PATH/TO/atlas-snp.rb -x cross_match_output -r reference.fasta -o prefix_of_outputs [-a adjust_quality_slope]
( The -x specifies the cross_match result concatenated from the batch results from step 1. )

It outputs a list of SNPs.

The format of SNP list is like this:

refName coordinate refBase homopolymer refEnv coverage SNPBase adjustedQual oriQual numVariant numAlter reads_info..

DDB0169550 317 A 5 ATTTATATTTTTA 92 T 16.48 19 1 1 T(19)089207_1601_3197(16)(134.0/140)-taaaatAataaat(1.43/0.0/1)swap;

U 6849 C 2 GGCAGACTTCAAG 12 T 37 37 1 1 T(37)086878_2543_1772(211)(237.0/243)+ggcagaTttcaag(0.82/0.0/1)snp;

This file contains all the information required to filter SNPs. We can use a straightforward ruby/perl/shell script to do filtering and calling homozygous/heterozygous based on:

number of variant reads vs coverage
adjusted quality
the size of longest homopolymer run nearby
For more help information and other options, just type -h as the argument.

2.1 Compute depth-coverage on reference

ruby atlas-mapper-coverage.rb -x cross_match.output -o prefix_of_output [-r ref.fasta] [ -t targeted_genomic_regions ] [ options ]
-r or -t is required. -r specifies the reference fasta file. If provided, the program will compute the depth-coverage for each base on the reference. NOTE this is time-consuming for large genomes. -t specifies a list of targeted genomic regions, with format like this:
target_name reference_name	   start_on_ref	 end_on_ref  direction(1 or -1)
If provided, the program will compute the depth-coverage only for bases in the targeted regions.