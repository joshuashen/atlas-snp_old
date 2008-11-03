Quick procedure:
 # convert file format, split into 1Million-reads batches, and then do Maq map
1. ruby pipeline_solexa_step1.rb -i InputDir -r reference.fasta -d 0 -n 3 -a 400

 # merge batch results, then do Maq assemble, SNP calling, pileup, and mapview etc. 
2. ruby PREFIX_of_outPut reference.fasta.bfa Mapping_of_InputDir/s*.map 


Note: -d 0 means run Maq locally in a precedural loop; to submit all batches into a LSF cluster (bsub ..), use -d 2; to submit to a SUN Grid cluster (qsub ..), use -d 1.


