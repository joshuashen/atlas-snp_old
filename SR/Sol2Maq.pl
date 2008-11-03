#!/usr/bin/env perl

# Author: Ye Liu, and Yufeng "Joshua" Shen
#
# Takes all the solexa export files in a given directory that can be run through maq, (filters them by chromosome and/or
# position,) breaks them into 2,000,000-line files, and runs them through maq. Intermediate files like .fastq and 
# 2,000,000-line files, along with intermediate maps, are saved in a Run directory with each pair of file's basic name 
# (s_n_Run).
#
# Directories entered as command-line arguments should have the last slash. Reference files can be in either .fasta or .bfa
# format; if .fasta is entered, the script will check for a corresponding .bfa file named {name}.fasta.bfq. If it does not exist,
# a file with this name will be created in the same directory as the .fasta file. If it does exist, it will be used.
#
# Chromosome number, start/end position, and maximum quality score in files can be entered as optional arguments with -c, 
# -s/-e, and -q, respectively. The default filter is none (all files are run) and the default maximum quality score is 40.
#
# To run, enter 1. the full path name of the directory containing standard-named Solexa files, 2. the full path name of the directory
# that you wish to output maq files, and 3. the full path name leading to the reference file in either .bfa or .fasta. Optional filtering
# commands include -c, which should be followed by the capital letter Roman numeral of a chromosome (for example, V for chromosome V); and -s
# and -e, which must coexist, followed by integer values of a position range from low to high. -q gives the option of setting the highest 
# quality score occuring in the files; if this is set too low, maq may be unable to recognize some reads and refuse to map files. The 
# default value is 40.

#used by BreakFile, Solexa2Fastq and RunMaq
use File::Basename;
## use Switch;
use Getopt::Long;
use vars qw /$opt_chr $opt_start $opt_end $opt_maxQual $opt_help $opt_maxMis $opt_maxDist/;

#Takes as inputs: SolexaOutputDirectory, MaqOutputDirectory, ReferenceFile, (optional: -c ChromNum, -s Start, -e End, -q maxQual)
@FileList = &GetFileList();
my $MaqOutputDirectory=();
my $ref=();
my @fullmap=();
mkdir "$ARGV[1]";
my $runfile = "$ARGV[1]Run.txt";
open (RUNFILE, ">>$runfile");

$opt_maxQual = 40;
$opt_maxMis = 3;
$opt_maxDist = 600;

GetOptions("help|h", "chr|c=s", "start|s=i", "end|e=i", "maxQual|q=i", "maxMis|n=i", "maxDist|a=i");

if ($opt_help) {
    &usage();
    exit;
}

if (($opt_start || $opt_end) && !($opt_start && $opt_end))  {
    die "Please input BOTH a start and an end value, or neither."
}


    


if ($ARGV[2] =~ /bfa$/)  {
    $ref=$ARGV[2]; 
}
else {
    if (!(-e "$ARGV[2].bfa")) {
    	system "maq fasta2bfa $ARGV[2] $ARGV[2].bfa 2>>$runfile";
    	}
    $ref = "$ARGV[2].bfa";
}

foreach (@FileList)  {
    chomp;
    my $line = $_;
    ($n, $p) = fileparse("$_");
    $MaqOutputDirectory = "${ARGV[1]}${n}_Run/";
    mkdir "$MaqOutputDirectory";
    my $nChunks = &ProcessFilePair ("$ARGV[0]$line", $MaqOutputDirectory);
    &RunMaq ($line, $ref, $nChunks, $MaqOutputDirectory);
}

#pulls a consensus file from all sets of data combined
my $fullmaps = "@fullmap";
system "maq mapmerge ${ARGV[1]}all.map $fullmaps 2>>$runfile";
system "maq pileup $ref ${ARGV[1]}all.map > ${ARGV[1]}pileup.txt";
system "maq assemble ${ARGV[1]}all.cns $ref ${ARGV[1]}all.map 2>>$runfile";
system "maq cns2snp ${ARGV[1]}all.cns > ${ARGV[1]}all.cns.snp 2>>$runfile";
system "maq.pl SNPfilter ${ARGV[1]}all.cns.snp > ${ARGV[1]}all.cns.filtered.snp 2>>$runfile";
system "maq cns2fq ${ARGV[1]}all.cns > ${ARGV[1]}all.cns.fq";
system "maq indelpe $ref ${ARGV[1]}all.map > ${ARGV[1]}all.indelpe";

print STDERR "Processing outputs...\n";
print RUNFILE "Processing outputs...\n";
system "sh ProcessOutput.sh ${ARGV[1]} $opt_start $opt_end"; 

close (RUNFILE);

#Extracts unmapped regions
#system "cat ${ARGV[1]}pileup.txt | awk '((\$2>-s) && (\$2<-e) && (\$4==0)){print $2}' | awk '{if (\$1 != stop+1) { if ( start) { print start"\t-\t"stop; } start = \$1;  }  stop = \$1; }' > ${ARGV[1]}unmapped_regions.txt";

#Opens a directory (input) and places all Solexa output files into an array @FileList.
sub GetFileList() {
    opendir(DIR, "$ARGV[0]") or die "Can't open directory";

    my @dir = readdir(DIR);
    my $in= "@dir";
    my @FileList=();
	
    for ($i=1; $i<10; $i++) {
      if($in=~/s_$i.1/ && $in=~/s_$i.2/) {
	push (@FileList, "$_[0]"."s_$i\n");
      }
    }			

    close(DIR);
    return @FileList;
}

#Processes each file pair based on the basic file name (s_n..._export.txt) and output directory
sub ProcessFilePair($BasicFileName $OutputDir) {
    my $BasicFileName = $_[0];
    print STDERR "Processing files ${BasicFileName}_1_export.txt and ${BasicFileName}_2_export.txt... \n";
    print RUNFILE "Processing files ${BasicFileName}_1_export.txt and ${BasicFileName}_2_export.txt... \n";
    my $nChunks = &BreakFile("${BasicFileName}_1_export.txt", "${BasicFileName}_2_export.txt", "$_[1]");
    return $nChunks;
}

#Breaks input files into smaller chunks of 2million lines each... takes in a pair of files
sub BreakFile($File1 $File2 $OutputDir)  {
    my $length=0;
    my $nChunks=1;
    my $switch=4;
    my $return=0;
    ($name1, $path) = fileparse("$_[0]", ".txt");
    ($name2) = fileparse("$_[1]", ".txt");

    open (IN1, $_[0]) or die "Can't open file $_[0]";
    open (IN2, $_[1]) or die "Can't open file $_[1]";
    open (OUT1, ">$_[2]$name1.$nChunks.txt") or die "Can't open output file $_[2]$name1.$nChunks.txt";
    open (OUT2, ">$_[2]$name2.$nChunks.txt") or die "Can't open output file $_[2]$name2.$nChunks.txt";

    if ($opt_chr && $opt_start)  {
	$switch = 1;
    }
    elsif ($opt_chr)  {
	$switch = 2;
    }
    elsif ($opt_start)  {
	$switch = 3;
    }

    while ($l1=<IN1> and $l2=<IN2>)  {
            if ($length >=2000000) { 
	        $length = 0;
	        $nChunks += 1;
	        close (OUT1);
		close (OUT2);
	        open (OUT1, ">$_[2]$name1.$nChunks.txt") or die "Can't open output file $_[2]$name1.$nChunks.txt";
		open (OUT2, ">$_[2]$name2.$nChunks.txt") or die "Can't open output file $_[2]$name2.$nChunks.txt";
            }

	    if ($switch == 1 )  {
		$return = (&FilterByChrom($l1, $l2) && &FilterByPosition($l1, $l2, $opt_start, $opt_end));
	    } elsif ($switch == 2) {
		$return = &FilterByChrom($l1, $l2);
	    } elsif ($switch == 3) {
		$return = &FilterByPosition($l1, $l2, $opt_start, $opt_end);
	    } elsif ($switch == 4) {
		$return = 1;
	    }

	    if ($return)  {
		$length += 1;
		print OUT1 "$l1";
		print OUT2 "$l2";
	    }
    }

    close (IN1);
    close (IN2);
    close (OUT1);
    close (OUT2);
    return $nChunks;
}
    
#Filters a pair of lines by chromosome number. At the time this is written, the format of chromosomes in the solexa files are
#"chr{NUMBER}.fa" where NUMBER is the chromosome number in capitalized Roman numerals. arguments: Line1, Line2.
sub FilterByChrom ($line1, $line2)  {
    my $chrom = "chr$opt_chr.fa";
    if ($_[0] =~ /($chrom)/ or $_[1] =~ /($chrom)/) {
        return "1";
    }
    else  {
	return "0";
    }
}

#Filters a pair of lines by read position. arguments: Line1, Line2.
sub FilterByPosition ($line1, $line2)  {
    my @line1 = split (/\t/, "$_[0]");
    my @line2 = split (/\t/, "$_[1]");

    if (defined $line1[12] or defined $line2[12])  {
        if (($line1[12] >= $opt_start) && ($line1[12] <= $opt_end))  { 
	    return "1";
        }
        elsif (($line2[12] >= $opt_start) && ($line2[12] <= $opt_end))  {
	    return "1";
        }
        else  {
	    return "0";
        }
    }
}	

#input solexa export file, output solexa-fastq file (that needs to be converted!)
sub Solexa2Fastq ($Input, $Output)  {
    open(SOLEXA,"$_[0]") or die "Can't open Solexa file";
    open(FASTQ, ">$_[1]") or die "Can't open Fastq file";

    print STDERR "Converting from Solexa to Fastq, $_[0] ...\n";
    print RUNFILE "Converting from Solexa to Fastq, $_[0] ...\n";

    
    while($line=<SOLEXA>)
    {
    	chomp($line);
	@nodes=split(/\s+/,$line);
	print FASTQ "\@$nodes[0]\:$nodes[1]\:$nodes[2]\:$nodes[3]\:$nodes[4]\/$nodes[5]\n$nodes[6]\n\+\n$nodes[7]\n";
    }

    close(SOLEXA);
    close(FASTQ);
}

#runs a pair of files through Maq. arguments: BasicFileName, Ref, nChunks, OutputDirectory
sub RunMaq ($filename, $ref, $nChunks, $outputdir) {
    print STDERR "Processing $_[0] in maq ...\n";
    print RUNFILE "Processing $_[0] in maq ...\n";
    ($name, $path) = fileparse("$_[0]");
    my $outdir = $_[3];
    my @maps = ();

    for ($m=1; $m<=$_[2]; $m++)  {
        my @files=();
        for ($j=1; $j<=2; $j++)  { 
	    $filetemp = "${outdir}temp.txt";
	    $file0 = "$outdir$name"."_$j"."_export.$m.txt";
	    $file1 = "$outdir$name"."_$j"."_export.$m.fastq";
	    $file2 = "$outdir$name"."_$j"."_export.$m.bfq";
	    &Solexa2Fastq ("$file0", "$filetemp");
#${ARGV[1]}all.map");
	    system "maq sol2sanger $filetemp $file1";
	    system "maq fastq2bfq $file1 $file2 2>>$runfile";
	    push (@files, "$file2");
        }
        system "maq map -n $opt_maxMis -a $opt_maxDist ${outdir}part$m.map $_[1] $files[0] $files[1] 2>>$runfile";
        push (@maps, "${outdir}part$m.map");
	system "rm $filetemp";
    }

    my $maps = "@maps";

    system "maq mapmerge $outdir$name.all.map $maps 2>>$runfile";
    push (@fullmap, "$outdir$name.all.map")
}

sub usage {

 print "Usage: perl ...pl inputdir outputdir ref [-n maxMismatch -a maxDistance] [-c chromosome -s start -e end -q maxquality] \n\n Remember to end the directories with the backslash (\\)\n";

}
