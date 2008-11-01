#!usr/bin/perl
##currently not in use. but if we could get this to work, it might be more efficient than calling on maq's
##sol2sanger script. It's less computation here if we construct the Fastq_scores array.
    open(SOLEXA,"$ARGV[0]") or die "Can't open Solexa file";
    open(FASTQ, ">$ARGV[1]") or die "Can't open Fastq file";

    for($i=0; $i<=127; $i++) {
        $Q=(10 * log(1 + 10 ** ($i - 64) / 10.0)) / log(10);
               
        #Phred2Fastq
        $Fastq_scores[$i]=chr(($Q<=93? $Q : 93) + 33);
    }

    print FASTQ "@Fastq_scores\n";

    while($line=<SOLEXA>)
    {
        chomp($line);
        @nodes=split(/\s+/,$line);
	@chars=split(//, @nodes[7]);
	print FASTQ "@chars\n";
	$chars_length=@chars;

        print FASTQ "\@$nodes[0]\:$nodes[1]\:$nodes[2]\:$nodes[3]\:$nodes[4]\:$nodes[5]\n$nodes[6]\n\+\n";

        for($i=0; $i<$chars_length; $i++) {             
            print FASTQ "$Fastq_scores[$chars[$i]]";
        }
        print FASTQ "\n";
    }
    close(SOLEXA);
    close(FASTQ);

