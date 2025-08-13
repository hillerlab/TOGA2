#!/usr/bin/env perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS.
## Tool that converts the CESAR2.0 annotated bed file to a genePred file.

use strict;
use warnings;

my $usage = "wrongs # args\n $0 - Convert CESAR2.0 annotated exons (as elements in a bed file) to a genePred file

Mandatory positional arguments:
1) species
2) name of CESAR 2.0 output directory (containing the exon coordinates)
3) output genePred file
";

$| = 1;		# == fflush(stdout)
die $usage if (scalar(@ARGV) != 3);

my $species  = $ARGV[0];
my $dir      = $ARGV[1];
my $genePred = $ARGV[2];

####
my %accListHash = my %hashFromFile = ();
## Collect all exons of that species and convert it to a 2 dimensional hash, the first dimension is the accession, the second is just a number, the value is the line.
my $call = "find $dir -type f -path '*/${species}/*' -exec cat {} \\;";
print STDERR "executing  '$call'  ...\n";
my @results = `$call`;
die "ERROR: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
print STDERR "gives ", $#results+1, " exons\n";
for (my $i = 0; $i<=$#results; $i++) {
	my $line = $results[$i];	
	$line =~s/\s+$//;
	my $accExon = (split /\t/,$line)[4];
	my @tmp = split(/_/,$accExon);
	pop @tmp;
	
	my $acc = join("_",@tmp);
	$accListHash{$acc} = "T";
	$hashFromFile{$acc}{$i} = $line;
}
my @accListAll = keys(%accListHash);

my $tmpBed     = `mktemp /dev/shm/XXXXX.bed`; chomp $tmpBed;
my $tmpBedFile = `mktemp /dev/shm/XXXXX.bed`; chomp $tmpBedFile;
my $tmpFileB2L = `mktemp /dev/shm/XXXXX.bed`; chomp $tmpFileB2L;

# create output dir
$call = "mkdir -p `dirname $genePred`";
`$call`;
open(FO,">$genePred") || die "Error printing to the genePred file 'genePred'\n";
foreach my $acc(@accListAll)
{
	my @linesList = ();  ## This array contains all the lines associated with a particular gene/accession
	foreach my $l(keys %{$hashFromFile{$acc}}) {
		push(@linesList,$hashFromFile{$acc}{$l});
	}
	
	my $strandOut = my $chrOut = "";
	my %chrHash = my %strHash = ();
	
	### Write everything to a bedFile, sort this bedFile later
	open(FOB,">$tmpBed");
	foreach my $l(@linesList) {
		my $cds = (split /\t/,$l)[3];  ## The field#4 contains the useful information for exon annotation
		my($chr,$start,$stop,$strand) = (split ",",$cds)[1,2,3,4];
		$strandOut = $strand;
		$chrOut    = $chr;
		
		$chrHash{$chr} = "T";
		$strHash{$strand} = "T";
		
		print FOB "$chr\t$start\t$stop\n";
	}
	close FOB;
	
	## A gene could come from more than one chromosome, sometimes more than one strand
	my @chrUnique = keys(%chrHash);
	my @strUnique = keys(%strHash);

	if (scalar(@chrUnique) == 1 && scalar(@strUnique) == 1) {
		my $listExons = bedToLine($tmpBed,$acc,$tmpFileB2L);	
		
		if ($listExons ne "") {
			my $line = "$acc\t$chrOut\t$strandOut\t$listExons";
			my $lineGP = lineToGenePredLine($line);
			print FO "$lineGP\n";
		}
	} elsif (scalar(@chrUnique) > 1 && scalar(@strUnique) == 1) {
		foreach my $chrU(@chrUnique) {
			my $ct = 0;
			open(FOB,">$tmpBedFile");
			foreach my $l(@linesList) {
				my $cds = (split /\t/,$l)[3];
				my($chr,$start,$stop,$strand) = (split ",",$cds)[1,2,3,4];

				if ($chr eq $chrU) {
					print FOB "$chr\t$start\t$stop\n";
					$ct++;
				} 
			}
			close FOB;

			if ($ct > 0) {
				my $listExons = bedToLine($tmpBedFile,$acc,$tmpFileB2L);
				if ($listExons ne "") {
					my $line = "$acc\t$chrU\t$strandOut\t$listExons";
					my $lineGP = lineToGenePredLine($line);
					print FO "$lineGP\n";
				}
			}	
		}
	} elsif (scalar(@strUnique) > 1 && scalar(@chrUnique) == 1) {
		foreach my $strU(@strUnique) {
			my $ct = 0;
			open(FOB,">$tmpBedFile");
			foreach my $l(@linesList) {
				my $cds = (split /\t/,$l)[3];
				my($chr,$start,$stop,$strand) = (split ",",$cds)[1,2,3,4];

				if ($strand eq $strU) {
					print FOB "$chr\t$start\t$stop\n";
					$ct++;
				}
			}
			close FOB;

			if ($ct > 0) {
				my $listExons = bedToLine($tmpBedFile,$acc,$tmpFileB2L);
				if ($listExons ne "") {
					my $line = "$acc\t$chrOut\t$strU\t$listExons";
					my $lineGP = lineToGenePredLine($line);
					print FO "$lineGP\n";
				}
			}	
		}
	} elsif (scalar(@strUnique) > 1 && scalar(@chrUnique) > 1) {
		foreach my $strU(@strUnique) {
			foreach my $chrU(@chrUnique) {
				my $ct = 0;
				
				open(FOB,">$tmpBedFile");		
				foreach my $l(@linesList) {
					my $cds = (split /\t/,$l)[3];
					my($chr,$start,$stop,$strand) = (split ",",$cds)[1,2,3,4];
					
					if ($strand eq $strU && $chr eq $chrU) {
						print FOB "$chrU\t$start\t$stop\n";
						$ct++;
					}
				}
				close FOB;
				
				if ($ct > 0) {
					my $listExons = bedToLine($tmpBedFile,$acc,$tmpFileB2L);
					if ($listExons ne "") {
						my $line = "$acc\t$chrU\t$strU\t$listExons";
						my $lineGP = lineToGenePredLine($line);
						print FO "$lineGP\n";
					}
				}
			}
		}
	}
	
}
close FO;
`rm $tmpBed $tmpBedFile $tmpFileB2L`;
print STDERR "Results for $species are in $genePred\n";

### Sub-routines

sub lineToGenePredLine  ## This function produces a line that is written to the final genePred file. Basically a format conversion exercise
{
	my $line = shift;
	my($name,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,3];
	
	my @exonsList = split(",",$exons);
	my $noe = scalar(@exonsList);

	my $firstExon = $exonsList[0];
	my $lastExon  = $exonsList[$#exonsList];

	my $cdsStart = (split /-/,$firstExon)[0];
	my $cdsStop  = (split /-/,$lastExon)[1];

	my @exonStart = my @exonEnd = ();
	foreach my $ex(@exonsList) {
		my($start,$stop) = (split "-",$ex);
		push(@exonStart,$start);
		push(@exonEnd,$stop);
	}

	my $startList = join(",",@exonStart);  ## Prepare two comma-separated strings: one contains exonStarts, the other contain exonStops
	my $stopList  = join(",",@exonEnd);
	my $lineReturn = "$name\t$chr\t$strand\t$cdsStart\t$cdsStop\t$cdsStart\t$cdsStop\t$noe\t$startList\t$stopList";
	return $lineReturn;
}

## Produce a line for genePred file from a bedFile
sub bedToLine
{
	my ($tmpBed,$acc,$tmpFile) = @_;
	
	my $sortBedCall = "bedSort $tmpBed stdout";
	my $sortedBedLines = `$sortBedCall`;
	die "ERROR!! '$sortBedCall' failed in bedToLine function\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	
	my @sortedBedLinesList = split(/\n/,$sortedBedLines);
	## Now read sortedBedFile
	my $chrPrint = "";
	my @cdsList = ();
		
	foreach my $l(@sortedBedLinesList) {
		$l =~s/\s+$//;
		my($chr,$start,$stop) = (split /\t/,$l);
		my $cds = "$start-$stop";
		push(@cdsList,$cds);
		$chrPrint = $chr;
	}
		
	my $exonsCdsList = join(",",@cdsList);
	### Filter out exons that overlap with each other, we do not print them to the final genePred file:
	my $listExons = filterOverlappingExons($exonsCdsList,$chrPrint,$acc,$tmpFile);
	return $listExons;
}


sub filterOverlappingExons
{
	my($exonsList,$chr,$acc,$tmpFileIn) = @_;
	
	my @exonsListArray = split(",",$exonsList);
	my $i = 1;
	open(FOB,">$tmpFileIn") || die "Error opening tmpFile in filterOverlappingExons function\n";
	foreach my $ex(@exonsListArray) {
		my($start,$stop) = (split "-",$ex)[0,1];
		print FOB "$chr\t$start\t$stop\tExon$i\n";
		$i++;
	}
	close FOB;
	
	## Run overlapSelect with mergeOutput option
	my $callOverlapSelect = "overlapSelect $tmpFileIn $tmpFileIn stdout -mergeOutput";
	my $overlapSelectOut = `$callOverlapSelect`;
	die "ERROR!! '$callOverlapSelect' failed in filterOverlappingExons function\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	
	my @overlapSelectOutList = split(/\n/,$overlapSelectOut);
	my %overlappingExons = ();

	foreach my $line(@overlapSelectOutList) {
		$line =~s/\s+$//;
		my($exon1,$exon2) = (split /\t/,$line)[3,7];
		if ($exon1 ne $exon2) {
			$overlappingExons{$exon1} = "T";
			$overlappingExons{$exon2} = "T";
		}
	}
	
	my @overlappingExonsList = keys(%overlappingExons);
	if (scalar(@overlappingExonsList) > 0) {
		my $i = 1;
		my @exonsListNew = ();
		foreach my $ex(@exonsListArray) {
			my $exonString = "Exon$i";
			push(@exonsListNew,$ex) if (! exists $overlappingExons{$exonString});
			$i++;
		}
		@exonsListArray = @exonsListNew;
	}
	
	my $listReturn = join(",",@exonsListArray);
	return $listReturn;
}
