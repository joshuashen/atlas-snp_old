# input format: 
# read_ID    chromosome    start    end    direction










## In 454 sequencing, some shotgun fragments are occasionally duplicated in the emulation PCR step. This creates skewed data for diploid genomes at those heterozygous loci. To avoid this one needs to remove the duplicated copies from SNP calling. 
## The cause of duplicated reads (from George): In production, When you make the sheared fragments in the first place, you only end up with a limited number of fragments. Then when you do emulsion PCR (for 454), or the amplification to make clusters for Solexa, sometimes fragments will escape from the micelle and populate new beads. But the new bead now has a duplicate fragment. So when you sequence, if you sequence "too much" from a library (imagine sequencing everything) you pick up these duplicates. Our solution has been to make multiple libraries and sequence less from each one. That seems to help a lot.

# So what one would like is to QC libraries for duplicates to make sure you haven't sequenced too much or that there hasn't been too much escape from micelles to form duplicates.







