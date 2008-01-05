# divide fasta file with quality

use strict;
use Getopt::Long;
use vars qw /$opt_help $opt_fa $opt_prefix $opt_qual $opt_length/ ;

my ($count, $read, $num, $out, %start, $line);

GetOptions("help|h", "fa|f=s", "prefix|p=s", "qual|q", "length|l=i");
#print "$opt_length\n";

if ($opt_help  || !$opt_length) {
    print "Usage: perl __.pl -f foo.fasta -p prefix_for_output -l length_for_part -q\n";
    exit;
}

if(!$opt_prefix) {
    $opt_prefix=$opt_fa.'_part';
}

$count=0;
$num=1;
$out=$opt_prefix.'_'.$num.'.fa';
open (OUT, ">$out") or die ("Cannot write to $out\n");

if ($opt_fa) {
    open(STDIN, $opt_fa) or die("no $opt_fa\n");
}

while(<STDIN>) {
    $line=$_;
    if($line=~ /^>(\S+)/) {
	$read=$1;
	$count+= length($line);
#	print "$count\n";
	if($count > $opt_length) {
	    $start{$read}++;
	    close(OUT);
	    $count=0;
	    $num++;
	    $out=$opt_prefix.'_'.$num.'.fa';
	    open(OUT, ">$out") or die ("bad $out\n");
	}
	print OUT $line;
    } else {
	$count+=length($line);
	print OUT $line;
    }
}
close(IN);
close(OUT);

if(!$opt_qual) {
    exit;
}

my $qual=$opt_fa.'.qual';

open (IN, $qual) or die ("no $qual\n");

$num=1;
$out=$opt_prefix.'_'.$num.'.fa.qual';
open (OUT, ">$out") or die ("Bad $out\n");
while(<IN>) {
    $line=$_;
    if($line=~ /^>(\S+)/) {
	$read=$1;
	if($start{$read} > 0) {
	    print OUT "\n";
	    close(OUT);
	    $num++;
	    $out=$opt_prefix.'_'.$num.'.fa.qual';
	    open (OUT, ">$out") or die ("bad $out\n");
	}
	print OUT $line;
    } else {
	print OUT $line;
    }
}
close (IN);
close (OUT);

