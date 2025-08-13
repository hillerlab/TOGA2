#!/usr/bin/env perl

## Virag Sharma, December 2014
## This script extracts a strectch of 10 nts upstream of start codon from a geneAnnotation file

use strict;
use warnings;

my $usage = "perl $0 longest_transcript_file output_file twoBitFile_Reference_Species\n";
die $usage if (scalar(@ARGV) != 3);

my $in_file     = $ARGV[0];
my $out         = $ARGV[1];
my $twoBitFile  = $ARGV[2];

open(FI,$in_file) || die "Cannot open the inputFile '$in_file'\n";
open(FO,">$out") || die "Cannot write to the outFile '$out'\n";
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($chr,$strand,$trans_start,$trans_stop) = (split /\t/,$line)[1,2,5,6];
        $chr = "chr".$chr if ($chr !~/chr/);
	
	my $start_1 = my $start_2 = my $stop_1 = my $stop_2 = "";
	
	if ($strand eq "+") {
		$start_1 = $trans_start - 10;
		$start_2 = $trans_start;
	} else {
		$start_1 = $trans_stop;
		$start_2 = $trans_stop + 10;
	}
	
	my $seq_start = `twoBitToFa $twoBitFile -seq=$chr -start=$start_1 -end=$start_2 stdout|grep -v ">"|tr -d "\n"`;
	$seq_start = revComp($seq_start) if ($strand eq "-");
	$seq_start = uc($seq_start);
	print FO "$seq_start\n";
}
close FI;
close FO;


sub revComp {  ## returns the recoverse complement of a given sequence
	my $seq = shift;
	$seq =~ tr/ATGCatgc-/TACGtacg-/;
	$seq = reverse($seq);
	return $seq;
}
