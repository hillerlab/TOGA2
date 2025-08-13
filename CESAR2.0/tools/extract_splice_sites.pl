#!/usr/bin/env perl
## Virag Sharma, October 2014

use strict;
use warnings;

my $usage = "usage: perl $0 longest_transcript_file acceptors_sequences donor_sequences two_bit_file_ref\n";
die $usage if (scalar(@ARGV) != 4);

my $refseq_file   = $ARGV[0];
my $acc_file      = $ARGV[1];
my $donor_file    = $ARGV[2];
my $two_bit_file  = $ARGV[3];

my $tmp_bed    = `mktemp`; chomp $tmp_bed;
my $tmp_fasta  = `mktemp`; chomp $tmp_fasta;

open(FI,"$refseq_file") || die "Error opening reference sequence file '$refseq_file'\n";
open(FO,">$tmp_bed") || die "Error writing to the temporary bed file '$tmp_bed'\n"; 
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my($gene,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,7];
	$chr = "chr".$chr if ($chr !~/chr/);

	my @exons_list = split(/,/,$exons);
	my $noe = scalar(@exons_list);
	if (scalar(@exons_list) > 1) {
		open(FO,">$tmp_bed") || die "Error writing to the bed_file '$tmp_bed'\n";
		if ($strand eq "+") {
			for(my $i = 0; $i < scalar(@exons_list); $i++) {
				my @exon_cds = split(/-/,$exons_list[$i]);
				
				my $acc_cds_1 = $exon_cds[0]-22; ## Acceptor business
				my $acc_cds_2 = $acc_cds_1 + 22;			
				print FO "$chr\t$acc_cds_1\t$acc_cds_2\t$gene#$i#Acc#$strand\n" if ($i != 0);
				
				my $donor_cds_1 = $exon_cds[1]; ## Donor business
				my $donor_cds_2 = $donor_cds_1 + 6;
				
				print FO "$chr\t$donor_cds_1\t$donor_cds_2\t$gene#$i#Donor#$strand\n" if ($i != ($noe - 1));
			}
		} else {
			@exons_list = reverse(@exons_list);
			
			for(my $i = 0; $i < scalar(@exons_list); $i++) {
				my @exon_cds = split(/-/,$exons_list[$i]);
				
				my $acc_cds_1 = $exon_cds[1] + 22; ## Acceptor business
				my $acc_cds_2 = $acc_cds_1 - 22;
				print FO "$chr\t$acc_cds_1\t$acc_cds_2\t$gene#$i#Acc#$strand\n" if ($i != 0);
				
				my $donor_cds_1 = $exon_cds[0]; ## Donor business
				my $donor_cds_2 = $donor_cds_1 - 6;
				print FO "$chr\t$donor_cds_2\t$donor_cds_1\t$gene#$i#Donor#$strand\n" if ($i != ($noe - 1));
			}
		}
	}
}
close FI;
close FO;
				
my $call = "twoBitToFa $two_bit_file -bed=$tmp_bed $tmp_fasta";	
system($call) == 0 || die "Error running call '$call'\n";
	
open(FOA,">$acc_file") || die "Error writing to acceptor_sequences file '$acc_file'\n";
open(FOD,">$donor_file") || die "Error writing to donor sequence file '$donor_file'\n";

my $ref2Hash = fasta_hash($tmp_fasta);
my %seqHash = %$ref2Hash;
foreach my $keys(keys(%seqHash)) {
	my $seq = $seqHash{$keys};
	my $strand = (split /#/,$keys)[3];
	$seq = revComp($seq) if ($strand eq "-");
	die "The strand value is not defined for '$keys'\n" if ($strand ne "-" && $strand ne "+");
	$seq = uc($seq);
			
	if ($keys =~/Acc/) {
		print FOA "$seq\n";
	} else {
		print FOD "$seq\n";
	}
}	
close FOA;
close FOD;
`rm -rf $tmp_bed $tmp_fasta`;

sub fasta_hash {
	my $file_input = shift;

    	my $id = my $seq = "";
    	my %seq_hash = ();

    	open(FIS,"$file_input") || die "Error opening input file '$file_input' in fasta_hash function\n";
    	while (my $line=<FIS>) {
	chomp $line;

	if ($line =~/>/) {
		if ($seq ne "") {
			$seq_hash{$id} = $seq;
			$id = $seq = "";
		}
			$id =$line;
			$id =~s/>//;
		} else { 
			$seq = $seq.$line;
		}
 	}
    	close FIS;
    
	$seq_hash{$id} = $seq;
	return(\%seq_hash);
}

sub revComp {  ## returns the recoverse complement of a given sequence
	my $seq = shift;
	$seq =~ tr/ATGCatgc-/TACGtacg-/;
	$seq = reverse($seq);
	return $seq;
}
