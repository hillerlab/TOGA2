#!/usr/bin/env perl
use strict;
use warnings;

my $usage = "usage: perl $0 input_file (produced by extract_splice_sites.pl or get_start_context.pl) splice_site_profile\n";
die $usage if (scalar(@ARGV) != 2);

my $in_file  = $ARGV[0];
my $out_file = $ARGV[1];

## Determine the length from the first line in the input file
my $line = `head -n 1 $in_file|tr -d "\n"`;
my $length = length($line);

my $total_count = `wc -l < $in_file`; chomp $total_count;

open(FO,">$out_file") || die "Error writing to outfile '$out_file'\n";
print FO "A\tT\tC\tG\n";

my $i = 1;
while ($i <= $length) {
	my @characters = qw(A T C G);
	my @prob_list = ();
	my $field = "c".$i;
	
	foreach my $nt(@characters) {
		my $count = `cut -$field $in_file|grep -w $nt|wc -l|tr -d "\n"`;
		my $prob = $count/$total_count;
		$prob = sprintf("%0.3f",$prob);
		push(@prob_list,$prob);
	}
	
	my $string = join("\t",@prob_list);
	print FO "$string\n";
	$i++;
}
close FO;
