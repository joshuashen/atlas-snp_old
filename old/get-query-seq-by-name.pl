use Getopt::Long;
use vars qw/$opt_help $opt_source $opt_list $opt_name $opt_output $opt_qual/;
use strict;

GetOptions("help|h", "source|s=s", "list|l=s", "name|n=s", "output|o=s", "qual|q");

if($opt_help || !$opt_source || (!$opt_list && !$opt_name)) {
    print "Usage 1: perl __.pl -s source.fasta -l names.list -o output.fa\n";
    print "Usage 2: perl __.pl -s source.fasta -n name.string -o output.fa\n";
    print "   note: if there is a quality file related to the source, the program will write a quality file related to output fasta as well\n";
    exit;
}

if(!$opt_output) {
    if($opt_list) {
	$opt_output=$opt_list.'.fa';
    } elsif ($opt_name) {
	$opt_output=$opt_name.'.fa';
    }
}

my (%in, $flag, %done);

if($opt_list) {
    open (IN, $opt_list) or die ("Bad file $opt_list\n");
    while(<IN>) {
	if(/^(\S+)/) {
	    $in{$1}++;
	}
    }
    close(IN);
} elsif ($opt_name) {
    $in{$opt_name}++;
}

$flag=0;
%done=();
open (IN, $opt_source) or die ("Bad file $opt_source\n");
open (OUT, ">$opt_output") or die ("Bad file $opt_output\n");
while(<IN>) {
    if (/^>(\S+)/) {
	my $name = $1;
	my $alter = '';
	if ($name =~ /^(\S+).scf/) {
	    $alter= $1;
	}
	
	if (($in{$name} > 0 || $in{$alter} > 0 ) && !($done{$name} > 0)) {
	    $flag=1;
	    $done{$name}++;
	    print OUT $_;
	} else {
	    $flag=0;
	}
    } elsif ($flag == 1) {
	print OUT $_;
    }

}
close(IN);
close(OUT);
	  

%done=();

if ($opt_qual) {
    my $qual=$opt_source.'.qual';

    if(-e $qual) {
	$flag=0;
	my $outqual = $opt_output.'.qual';
	open (QIN, "$qual") or die ("Cannot open file $qual\n");
	open (QOUT, ">$outqual")  or die ("Bad file $outqual\n");
	
	
	while (<QIN>) {
	    if (/^>(\S+)/) {
		if ($in{$1} > 0 && !($done{$1} >0)) {
		    $flag=1;
		    $done{$1}++;
		    print QOUT $_;
		} else {
		    $flag=0;
		}
	    } elsif ($flag==1) {
		print QOUT $_;
	    }
	}
	
	close(QIN);
	close(QOUT);
    }
}
exit;
