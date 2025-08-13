#!/usr/bin/env perl

### Virag Sharma, 2017 @MPI-CBG and MPI-PKS, Dresden.
### This script formats genePred files in such a way that only the coordinates of coding exons are printed to an output file
### Also, the script excludes those genes whose CDS length is not a multiple of 3 or where there is an intron shorter than 30 bp

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
$| = 1;		# == fflush(stdout)

# input options/parameters: two optional parameters: A comma separated list of species and a transcriptID
my $longest = 0;  ## two optional parameters:
GetOptions ("longest" => \$longest);

my $usage = "wrongs # args\n $0 - Format a genePred file such that only the coordinates of coding exons are printed.
Each line should contain gene predictions in the format prescribed at https://genome.ucsc.edu/FAQ/FAQformat.html#format9

$0 excludes those genes whose CDS length is not a multiple of 3 or where there is an intron shorter than 30 bp or where the CDS is incomplete
usage:
$0 inputFile outputFile excludeFile

Options:
  -longest - Instead of printing all isoforms for a gene, only print the longest isoform\n";
  
die $usage if (scalar(@ARGV) != 3);

my $fileIn      = $ARGV[0];
my $fileOut     = $ARGV[1];
my $excludeFile = $ARGV[2];

##### Step 0: If -longest is specified, then just run genePredSingleCover and get the longest transcript per gene
my $fileRead = $fileIn;
if ($longest) {
	my $tmpFile = `mktemp`; chomp $tmpFile;
	my $call = "genePredSingleCover $fileIn $tmpFile";
	system($call) == 0 || die "Error running genePredSingleCover\n";
	$fileRead = $tmpFile;
	
	my $ct = 0;
	open(FI, $fileRead) || die "Error opening the file '$fileRead'\n";
	while (my $line = <FI>) {
		$ct++;	
	}	
	close FI;
	
	if ($ct == 0) {
		open(FOE,">$excludeFile") || die "Error writing to exclude file '$excludeFile'\n";
		open(FI, "$fileIn") || die "Error opening input file '$fileIn'\n";
		while (my $line = <FI>) {
			$line =~s/\s+$//;
			next if ($line =~/#/);
			
			my ($isoform,$chr,$strand,$transStart,$transStop,$cdsStart,$cdsStop,$noe,$exon_start,$exon_end) = (split /\t/,$line)[0,1,2,3,4,5,6,7,8,9];
			if ($cdsStart == $cdsStop) { ## These are not protein-coding genes. For RNA coding genes, cds_start = cds_stop. And we do not process them.
				print FOE "Non coding RNA gene is excluded (CDS start is the same as CDS stop) --> $line\n";
			}
		}
		close FI;
		close FOE;
		exit;		
	}
}

#### Step 1 --> read the file into array and have only some fields of interest
open(FI,$fileRead) || die "Error opening input file '$fileIn'\n";
open(FOE,">$excludeFile") || die "Error writing to exclude file '$excludeFile'\n";

my @fileArray = ();
while (my $line=<FI>)
{		
	$line =~s/\s+$//;
	next if ($line =~/#/);
	
	if ($line =~/incmpl/) {
		print FOE "Incomplete CDS --> $line\n";
		next;
	}
	
	my @tmp = split(/\t/,$line);
	my $nf = scalar(@tmp);
	die "Incorrent format\n'$line' should have minimum 10 fields, found only $nf\n" if ($nf < 10);
	
	my ($isoform,$chr,$strand,$transStart,$transStop,$cdsStart,$cdsStop,$noe,$exon_start,$exon_end) = (split /\t/,$line)[0,1,2,3,4,5,6,7,8,9];
	if ($cdsStart == $cdsStop) { ## These are not protein-coding genes. For RNA coding genes, cds_start = cds_stop. And we do not process them.
		print FOE "Non coding RNA gene is excluded (CDS start is the same as CDS stop) --> $line\n";
		next;	
	}	
	$exon_start =~s/,$//;
	$exon_end =~s/,$//;		
	my @exon_start = split(",",$exon_start);
	my @exon_end   = split(",",$exon_end);
		
	my @exon_boundaries = ();
	for(my $j = 0; $j <@exon_start; $j++)			### One of the operations that is performed as the file is read- instead of two fields for exon start coordinates and exon end coordinates
	{											### this FOR loop concatenates the two fields- as a result, each exon is represented as $exon_start-$exon_end.
		$exon_boundaries[$j] = "$exon_start[$j]-$exon_end[$j]";
	}
	my $exon_boundaries_string = join(",",@exon_boundaries);
			
	push(@fileArray,"$isoform\t$chr\t$strand\t$transStart\t$transStop\t$cdsStart\t$cdsStop\t$noe\t$exon_boundaries_string");
}
close FI;
close FOE;

### Step 2 ###
### To format every transcript such that only the coordinates of the coding sequence (CDS) are written to the output file
my %lengthHash = my %transcriptsHash =  ();

foreach my $line(@fileArray)
{
	$line =~s/\s+$//;
	my ($isoform,$chr,$strand,$transStart,$transStop,$cdsStart,$cdsStop,$noe,$exonsList) = (split /\t/,$line)[0,1,2,3,4,5,6,7,8];
	
	my @exons = split(/,/,$exonsList);
	my $exonsPrint = "";
	
	if ($noe == 1) {
		$exonsPrint = "$cdsStart-$cdsStop";  ## If there is only one exon, the coordinates are simply the value of the cds_start and cds_stop.
	} else {
		if ($strand eq "+") {				## Situation 1: if the gene is on the plus strand	
			my $j = 0;
			my $first_exon = my $fe = my $le = my $last_exon = "";
			my $exon_cds = "";
		
			foreach $exon_cds(@exons)
			{
				my @exon_each = split(/-/,$exon_cds);
				
				if ($cdsStart >= $exon_each[0] && $cdsStart <= $exon_each[1])  {
					$first_exon = "$cdsStart-$exon_each[1]"; ## Find the exon or rather exon boundary where $cds_start lies
					$fe = $exon_cds;
					last;
				}
				$j++;
			}
			
			splice(@exons,0,$j);		## Splice the exon array based on the index of the first exon, this value is obtained from $j.
			$j = 0;						## Reset the counter to 0.
			
			foreach $exon_cds(@exons)
			{
				my @exon_each = split(/-/,$exon_cds);
				if ($cdsStop >= $exon_each[0] && $cdsStop <= $exon_each[1]) {
					$last_exon = "$exon_each[0]-$cdsStop";   ## Find the exon where cds_stop lies
					$le = $exon_cds;
					last;
				}
				$j++;
			}
			
			if ($fe eq $le) {
				$exonsPrint = "$cdsStart-$cdsStop";				## In this case, there is only one exon whose coordinates are $cds_start-$cds_stop
			} else {											## Rest of the exons are all non-coding exons.		
				splice(@exons,$j+1);							## Splice the exons array based on the index on the last exon, this value is obtained from $j.
				$exons[0] = $first_exon;						## The first exon is $first_exon , the last is $last_exon.
				$exons[$#exons] = $last_exon;
				$exonsPrint = join(",",@exons);
			}
		}
					
		else {						### If the gene is on the negative strand
			@exons = reverse(@exons);
			
			my $j = 0;
			my $first_exon = my $fe = my $le = my $last_exon = "";
			my $exon_cds = "";
			
			foreach $exon_cds(@exons) {
				my @exon_each = split(/-/,$exon_cds);			
				
				if ($cdsStop <= $exon_each[1] && $cdsStop >= $exon_each[0]) {
					$first_exon = "$exon_each[0]-$cdsStop";			## Find the exon or rather exon boundary where $cds_stop lies ($cds_stop because this is negative strand).
					$fe = $exon_cds;
					last;
				}
				$j++;
			}
			splice(@exons,0,$j);								## Splice the exon array based on the index of the first exon, this value is obtained from $j.	
			$j = 0;												## Reset the counter to 0
			
			foreach $exon_cds(@exons) {
				my @exon_each = split(/-/,$exon_cds);
		
				if ($cdsStart >= $exon_each[0] && $cdsStart <= $exon_each[1]) {
					$last_exon = "$cdsStart-$exon_each[1]";			## Find the exon or rather exon boundary where $cds_start lies ($cds_start because this is negative strand).
					$le = $exon_cds;
					last;
				}
				$j++;
			}
			
			if ($fe eq $le) {
				$exonsPrint = "$cdsStart-$cdsStop";					## In this case, there is only one exon whose coordinates are $cds_start-$cds_stop
			} else {
				splice(@exons,$j+1);							## Splice the exons array based on the index on the last exon, this value is obtained from $j.
				$exons[0] = $first_exon;
				$exons[$#exons] = $last_exon;
		
				@exons = reverse(@exons);							## Reverse the array since the gene is on the negative strand.
				$exonsPrint = join(",",@exons);
			}
		}
	}
	
	$transcriptsHash{$isoform} = "$isoform\t$chr\t$strand\t$transStart\t$transStop\t$cdsStart\t$cdsStop\t$exonsPrint";
	$lengthHash{$isoform} = getCdsLength($exonsPrint);
}

my ($ref2ExcludedGenes) = getExcludedGenes(\%transcriptsHash);
my %excludedGenes = %$ref2ExcludedGenes;

### Final printing

open(FO,">$fileOut") || die "Error writing to the outFile '$fileOut'\n";
foreach my $transcript(keys(%transcriptsHash))
{
	print FO "$transcriptsHash{$transcript}\n" if (! exists $excludedGenes{$transcript});
}
close FO;
	
## Writing to the exclude file:	
open(FOE,">>$excludeFile");
foreach my $transcript (keys (%excludedGenes)) {
	print FOE "$transcript\t$excludedGenes{$transcript}\n";
}
close FOE;
`rm $fileRead` if ($longest);


##################################################################
############################ functions follow ####################
##################################################################

sub getExcludedGenes		### This function returns a list of accessions that should be excluded. 
{							### Critera for exclusion: i) the length of the CDS is not a multiple of 3 ii) An intron in the gene is shorter than 30 bp
	my $ref2InputHash = shift;	
	
	my %hashExclude = ();
	my %transcriptsHash = %$ref2InputHash;
	
	
		foreach my $acc (keys (%transcriptsHash)) {

			my $line = $transcriptsHash{$acc};
			my($acc,$strand,$exons_list) = (split /\t/,$line)[0,2,7];
				
			my $flag = "";
			my @statementList = ();
			my $s = 0;
			
			my @exonic_cds = split(/,/,$exons_list);
			@exonic_cds = reverse(@exonic_cds) if ($strand eq "-");
			
			if (scalar(@exonic_cds) > 1) {  ## If the number of exons is greater than 1, check for the intron length.
				for(my $j = 0; $j < @exonic_cds -1; $j++) {
					my @exon_cds1 = split(/-/,$exonic_cds[$j]);
					my @exon_cds2 = split(/-/,$exonic_cds[$j+1]);
					
					my $intron_length = "";
					
					$intron_length = $exon_cds2[0] - $exon_cds1[1] if ($strand eq "+");
					$intron_length = $exon_cds1[0] - $exon_cds2[1] if ($strand eq "-");
					
					if ($intron_length <= 30) {  	### Condition 1
						$flag = "T";
						my $exonNumber = $j + 1;
						push(@statementList,"Intron downstream of exon $exonNumber is <= 30bp");
					}
				}
			}
			
			my $length_cds = getCdsLength($exons_list);
			if ($length_cds %3 != 0) {	### Condition 2
				$flag = "T";
				push(@statementList,"The reference transcript has a frameshift in this gene");
			}
			
			my $statementString = join(", ",@statementList);
			$hashExclude{$acc} = $statementString if ($flag eq "T"); 
		}
	
	return \%hashExclude;
}

sub getCdsLength
{
	my $exon_cds = shift;
	my @tmp = split(/,/,$exon_cds);
	
	my $length = 0;
	foreach my $exon(@tmp) {
		my @tmpExon = split(/-/,$exon);
		$length = $length + ($tmpExon[1] - $tmpExon[0]);
	}
	
	return $length;
}

