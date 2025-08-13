#!/usr/bin/env perl

## Virag Sharma, 2017. MPI-CBG and MPI-PKS.
## The CESAR2.0 annotation engine

use strict;
use warnings;
use Scalar::Util::Numeric qw(isint);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
# this perl module must be located in the same dir as this script
use FindBin qw($Bin);
use lib "$Bin/";
use geneAnnotationFunctions;

# input options/parameters:
my $verbose = my $list = 0;
my $maxMem = 30;
GetOptions ("v|verbose" => \$verbose, "maxMemory=i" => \$maxMem, "list" =>\$list);


my $usage = "wrongs # args\n $0 - The core annotation engine that runs CESAR2.0 on a given genePred file and outputs 
coordinates of intact exons in a bed format.
Mandatory positional arguments:
1) input [transript identifier| a file containing transcripts where every line is one transcript]
2) the indexed alignment (.bb) file
3) the genepred file which contains information about all transcripts (output of formatGenePred.pl)
4) reference species --> for example: hg38
5) speciesList --> a comma separated string of species whose exons are being annotated. For example mm10,bosTau7 for mouse and cow
6) output directory --> a file with the same name as the query species is generated in this directory and coordinates of all intact exons are appended to this file
7) the twoBit directory that contains the genomes of the reference and all query species
8) the path to the CESAR binary. This path is expected to contain the extra/ subdirectory that contains donor and acceptor profiles and the substitution matrix
Options:
  -v|verbose
  -maxMemory int    allow CESAR 2.0 to allocate at most that much memory in Gb. Default is 30 (Gb).
  -list [when the input is a file containing one transcript per line] 
";

die $usage if (scalar(@ARGV) < 8);

my $input       = $ARGV[0];
my $mafIndex    = $ARGV[1];
my $genePred    = $ARGV[2];
my $reference   = $ARGV[3];
my $speciesList = $ARGV[4];
my $outDir      = $ARGV[5];
my $twoBitDir   = $ARGV[6];
my $cesarDir    = $ARGV[7];

## Check arguments
die "The mafIndex file '$mafIndex' does not exist\n" if (! -e $mafIndex);
die "The genePred file '$genePred' does not exist\n" if (! -e $genePred);
die "The genes list file '$input' does not exist, you specified that the input is a list\n" if (! -e $input && $list != 0);

my $clade = "human";

my @speciesArray = split(",",$speciesList);
my %speciesOfInterest = map{$_ => 1}@speciesArray;

## Step0: check if two bit files are present for all species of interest
foreach my $species(@speciesArray) {
	my $twoBitFile = "$twoBitDir/$species/$species.2bit";
	die "The two bit file ('$twoBitFile') for species '$species' is not present. Aborting\n" if (! -e $twoBitFile);
}
my $twoBitFileRef = "$twoBitDir/$reference/$reference.2bit";

my %genesOfInterest = ();
if ($list) {
	open(FI,$input) || die "Error opening genesList file '$input'\n";
	while (my $gene = <FI>) {
		chomp $gene;
		$genesOfInterest{$gene} = "T";
	}
	close FI;
} else {
	$genesOfInterest{$input} = "T";
}

## Step1: Get the information about the input transcripts/genes from the genePred file
my %geneInfoHash = ();
open(FI,$genePred) || die "Error opening genePred file\n";
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($gene,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,7];
	if (exists $genesOfInterest{$gene}) {
		$geneInfoHash{$gene} = "$chr#$strand#$exons";
	}
}
close FI;

## check if there is nothing in the geneInfoHash
my @genesInt = keys(%geneInfoHash);
if (scalar(@genesInt == 0)) {
	print "############## ERROR ##############\n";
	print "No information found for the input gene '$input' from the genePred file '$genePred'\n";
	print "Could it be that the input is a file and you forgot to specify the '-list' parameter\nCheck and re-run\n";
	print "###################################\n";
	exit;
}

my @tempFiles = ();  ## create all temp files upfront and store them in @tempFiles array
my $mafFile      = `mktemp /dev/shm/exon.maf.XXXXXXXXX`; chomp $mafFile;
my $speciesBed   = `mktemp /dev/shm/species.bed.XXXXXX`; chomp $speciesBed;
my $speciesFasta = `mktemp /dev/shm/species.fa.XXXXXX`; chomp $speciesFasta;
my $cesarInput   = `mktemp /dev/shm/cesarIn.fa.XXXXXX`; chomp $cesarInput;
my $cesarOutput  = `mktemp /dev/shm/cesarOut.aln.XXXXXX`; chomp $cesarOutput;
push(@tempFiles,$mafFile,$speciesBed,$speciesFasta,$cesarInput,$cesarOutput);
	
foreach my $gene(keys(%geneInfoHash)) {
	my($chrRef,$strandRef,$allExons) = (split /#/,$geneInfoHash{$gene});

	print "The gene '$gene' comes on chromosome '$chrRef', strand '$strandRef'. Exons are '$allExons'\n" if ($verbose);
	my @exonsList = split(",",$allExons);
	@exonsList = reverse(@exonsList) if ($strandRef eq "-");

	## Step2: Get exon phases: i.e the split codon bases for each exons on both the 5' and the 3' end
	my $noe = scalar(@exonsList);
	my $ref2PhaseList = getExonPhaseList(\@exonsList);
	my @exonPhasesList = @$ref2PhaseList;

	## Step3: Run mafExtract on every exon and store the cdsStart, cdsStop, the chr and the strand for every query species present in the maf
	## At this stage also identify species-exon pairs which should  be excluded because:
	## 1) exon comes on more than one scaffold|strand, 
	## 2) nothing  but e-Lines 
	## 3) large insertions in the exon sequence of the query species

	my %ignoreExonsHash = ();
	my %cdsStartHash = my %cdsStopHash = my %strandHash = my %chrHash = ();
	my %exonsLengthHash = my %exonNumberHash = ();
	my %refSeqHash = my %referenceSpliceSites = ();

	my %hashInfo = ();
	## Initialize every entry in the hashInfo to NA at the start
	foreach my $species(@speciesArray) {
		my $k = 0;
		while ($k < $noe) {
			$hashInfo{$species}{$k}     = "NA#NA#NA#0";
			$k++;
		}
	}
	print "Processing gene '$gene'\n";

	my %exonHasStop = ();
	print "Fetching the alignment coordinates for each species and for each exon from '$gene' from the indexed maf\n###\n" if ($verbose); 
	my $exonCt = 1;
	foreach my $e(@exonsList) {
		print "Processing exon '$e'\n" if ($verbose);
		$exonNumberHash{$e} = $exonCt;
		my($refStart,$refStop) = (split /-/,$e)[0,1];
		my $lengthExon = abs($refStart - $refStop);
		$exonsLengthHash{$exonCt} = $lengthExon;
		
		## process reference sequence
		my ($phase5Prime,$phase3Prime) = (split /-/,$exonPhasesList[$exonCt-1])[0,1];
		my $refStartSS = $refStart - 2;  # the coordinates of the splice site
		my $refStopSS  = $refStop + 2;
		
		my $seqRefPrint = getSequenceFrom2Bit($twoBitFileRef,$chrRef,$refStartSS,$refStopSS);
		$seqRefPrint = revComp($seqRefPrint) if ($strandRef eq "-");
		
		## Get splice site nucleotides and find out if the acceptor is AC or the donor is AT
		print "The reference sequence is '$seqRefPrint'\n" if ($verbose);
		
		my $accSite = lc(substr($seqRefPrint,0,2));
		my $donorSite = lc(substr($seqRefPrint,length($seqRefPrint)-2,2));
		$referenceSpliceSites{$exonCt}{"acc"} = $accSite;
		$referenceSpliceSites{$exonCt}{"do"}  = $donorSite;
		
		## Get reference sequence that will be used as input to CESAR
		$seqRefPrint = substr($seqRefPrint,2,length($seqRefPrint)-4);  ## Get rid of the splice sites from the refSequence
		$seqRefPrint = fixPhase($seqRefPrint,$phase5Prime,$phase3Prime);
		
		## check if exon has stop codon:
		my $exonSeqInRef = $seqRefPrint;
		$exonSeqInRef = substr($exonSeqInRef,$phase5Prime,length($exonSeqInRef)-$phase5Prime-$phase3Prime);
		$exonSeqInRef =~s/\|//g;
		
		my $s= 0;
		while ($s < length($exonSeqInRef)) {
			my $codon = substr($exonSeqInRef,$s,1).substr($exonSeqInRef,$s+1,1).substr($exonSeqInRef,$s+2,1);

			if ($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA") {	 
				if ( ! ($s == length($exonSeqInRef) - 3) ) { ## This should not be the last codon at the end of the sequence though:
					$exonHasStop{$exonCt} = "T";
					print "'$exonSeqInRef' has stop codon at position '$s'\n" if ($verbose);
				}
			}
			$s = $s+3;
		}
		
		$refSeqHash{$exonCt} = $seqRefPrint;
		
		my $mafExtractCall = "/bin/bash -c 'set -o pipefail; mafExtract -region=$chrRef:$refStart-$refStop $mafIndex stdout|mafSpeciesSubset stdin NULL $mafFile -speciesList=$speciesList,$reference'";
		print "Running $mafExtractCall\n" if ($verbose);
		system($mafExtractCall) == 0 || die "Error running '$mafExtractCall'\n";
		
		my $ref_species_list = get_maf_species_list($mafFile); ## check if the maf file contains some s-lines. If a maf is generated from nets without running mafAddIRows
		my @species_list_maf = @$ref_species_list;
		if (scalar(@species_list_maf) < 1) {
			print "Omitting exon '$e', since there is nothing in the maf for this exon\n";
			$exonCt++;
			next;
		}
		
		my($ref2cdsStartQuery,$ref2cdsStopQuery,$ref2strandQuery,$ref2chrQuery,$ref2moreThanOneStrand,$ref2moreThanOneChr,$ref2SpeciesList,$ref2cdsPrint,$ref2seqLength) = getCoordinatesFromMaf($mafFile,$reference);
		
		my @speciesPresentInMaf = @$ref2SpeciesList;
		foreach my $species(@speciesPresentInMaf) {
			
			my $cdsStart = $ref2cdsStartQuery->{$species};
			my $cdsStop  = $ref2cdsStopQuery->{$species};
			
			if ($ref2moreThanOneStrand->{$species} eq "T" || $ref2moreThanOneChr->{$species} eq "T" || ($cdsStart eq "NA" && $cdsStop eq "NA") || ($cdsStart > $cdsStop) || (abs($cdsStart-$cdsStop) - $lengthExon > 10000) ) {
				$ignoreExonsHash{$exonCt}{$species} = "T";
				
				if ($strandRef eq "-") {  ## The function (getIntronTable) that processes hashInfo uses a zero based index and does not need the exons to be reversed in case the gene is on minus strand
					$hashInfo{$species}{$noe-$exonCt}     = "$ref2cdsPrint->{$species}#$ref2chrQuery->{$species}#$ref2strandQuery->{$species}#$ref2seqLength->{$species}"; ## Store everything useful from the maf in the "hashInfo" hash
				} else {
					$hashInfo{$species}{$exonCt-1}     = "$ref2cdsPrint->{$species}#$ref2chrQuery->{$species}#$ref2strandQuery->{$species}#$ref2seqLength->{$species}";
				}
				
				print "This exon is ignored --> $exonCt ... $species\n" if ($verbose);
				print "Ignored exon attributes: Strand --> '$ref2moreThanOneStrand->{$species}', chromosome --> '$ref2moreThanOneChr->{$species}', cdsStart --> '$cdsStart, cdsStop --> $cdsStop, lengthRef exon --> $lengthExon\n'" if ($verbose);
			} else {
				$cdsStartHash{$species}{$exonCt} = $ref2cdsStartQuery->{$species};
				$cdsStopHash{$species}{$exonCt}  = $ref2cdsStopQuery->{$species};
				$strandHash{$species}{$exonCt}   = $ref2strandQuery->{$species};
				$chrHash{$species}{$exonCt}      = $ref2chrQuery->{$species};
				
				if ($strandRef eq "-") {
					$hashInfo{$species}{$noe-$exonCt}     = "$ref2cdsPrint->{$species}#$ref2chrQuery->{$species}#$ref2strandQuery->{$species}#$ref2seqLength->{$species}"; ## Store everything useful from the maf in the "hashInfo" hash
				} else {
					$hashInfo{$species}{$exonCt-1}     = "$ref2cdsPrint->{$species}#$ref2chrQuery->{$species}#$ref2strandQuery->{$species}#$ref2seqLength->{$species}";
				}
				
				print "This exon is considered by CESAR --> $exonCt .. $species\n" if ($verbose);
			}
		}
		$exonCt++;
	}

	print "Now fetching intron lengths for each exon of interest\n" if ($verbose);
	### get intron lengths

	my %intronLengthHash = ();
	my $ref2IntronLengthTable =  getIntronTable($gene,$chrRef,$strandRef,\@exonsList,\@speciesArray,\%hashInfo,$verbose);
	my @intronLengthInfo = @$ref2IntronLengthTable;

	foreach my $l(@intronLengthInfo) {
		my($gene,$end,$chr,$exon,$species,$length) = (split /\t/,$l);  ## This is how the output of getIntronTables.pl look like:
		my $exonNumber = $exonNumberHash{$exon};		       ## ABCB4	fivePrime	chr7	87475385-87475465	mm10	1212
		
		$intronLengthHash{$species}{$exonNumber}{$end} = $length;
	}

	## Now I know the coordinates, time to run twoBitToFa to extract the sequence for every exon for every species:
	my %seqHash = my %headerHash = (); ## A 2 dimensional hash where the first dimension is the exon number, the second is the species and the value is the sequence|header

	foreach my $species(@speciesArray) {
		
		my $twoBitFile = "$twoBitDir/$species/$species.2bit";
		my %ntsAtEitherEnd = ();
		open(FO,">$speciesBed") || die "Error writing to speciesBed file '$speciesBed' for species '$species'\n";
		
		my $i = 1;
		my $seqCt = 0;
		while ($i <= $noe) {
			if (exists $ignoreExonsHash{$species}{$i} || ! exists $cdsStartHash{$species}{$i}) {  ## If there were no s-Lines for the species of interest in maf, then it would neither make it to ignore exons
				$i++;									       ## nor to the cdsStart,cdsStop hash.
				next;
			}
			
			my $lengthExon      = $exonsLengthHash{$i};
			my $cdsStartQuery   = $cdsStartHash{$species}{$i};
			my $cdsStopQuery    = $cdsStopHash{$species}{$i};
			my $lengthExonQuery = $cdsStopQuery - $cdsStartQuery;
			my $strandQuery     = $strandHash{$species}{$i};
			my $chrQuery        = $chrHash{$species}{$i};

			if ($verbose) {	
				print "Fetching coordinates for exon$i for species '$species'\n";
				print "The length is '$lengthExon', the cdsStart is '$cdsStartQuery', the cdsEnd is '$cdsStopQuery', the lengthExonQuery is '$lengthExonQuery', the strand is '$strandQuery' and the chromosome is '$chrQuery'\n";
			}
			
			my $intronLength5Prime = my $intronLength3Prime = 0;
				
			$intronLength5Prime = $intronLengthHash{$species}{$i}{"fivePrime"} if (exists $intronLengthHash{$species}{$i}{"fivePrime"});
			$intronLength3Prime = $intronLengthHash{$species}{$i}{"threePrime"} if (exists $intronLengthHash{$species}{$i}{"threePrime"});
		
			print "The intron length in 5Prime direction is '$intronLength5Prime' while the intron length in 3Prime direction is '$intronLength3Prime'\n" if ($verbose);
			
			if ($intronLength5Prime < 0 || $intronLength3Prime < 0) {
				$i++;
				next;
			}
			
			$intronLength5Prime = 50 if ($i == 1);  ## For the first exon, you can not have a length in 5' direction
			$intronLength3Prime = 50 if ($i == $noe);  ## Same for the last exon
			$intronLength5Prime = $intronLength3Prime = 50 if ($noe == 1);  ## And for a single exon gene
			$intronLength5Prime = 50 if (! exists $intronLengthHash{$species}{$i}{"fivePrime"}); ## if the exon in question is exon#2 and the first exon is deleted, we cannot get any value for intronLength here.
			$intronLength3Prime = 50 if (! exists $intronLengthHash{$species}{$i}{"threePrime"}); ## if the exon in question is exon#12 and the last exon (13th) is deleted, we cannot get any value for intronLength here.

			### Now set the default flankingSequence length = 50
			my $flankDef = 50;
			$flankDef = $flankDef + ($lengthExon - $lengthExonQuery) if ($lengthExon > $lengthExonQuery);
			print "The flank length is '$flankDef'\n" if ($verbose);
			my $flank5Prime = my $flank3Prime = $flankDef; ## The flankDef is what we have computed above
			$flank5Prime = $intronLength5Prime  if ($flankDef > $intronLength5Prime);  ## If the flankDef is longer than the intronLength that is available, we just use the intronLength
			$flank3Prime = $intronLength3Prime  if ($flankDef > $intronLength3Prime);
			
			my $cdsStartFinal = my $cdsStopFinal = "";
			print "The exon is $i --> length5Prime is $flank5Prime -->  length3Prime is $flank3Prime\n" if ($verbose);

			if ($strandQuery eq "+") {  ## These tell you the actual coordinates from where the sequence should be extracted
				$cdsStartFinal = $cdsStartQuery - $flank5Prime;	
				$cdsStopFinal  = $cdsStopQuery + $flank3Prime;
			} else {
				$cdsStartFinal = $cdsStartQuery - $flank3Prime;
				$cdsStopFinal  = $cdsStopQuery + $flank5Prime;
			}
			
			my $chromSize = getChromSize($species,$chrQuery);
			my $startString = my $endString = "";
			
			if ($cdsStartFinal < 0) {
				$startString = createString(abs($cdsStartFinal),"N");
				$cdsStartFinal = 1;
			}
			
			if ($cdsStopFinal > $chromSize) {
				$endString = createString(($cdsStopFinal - $chromSize),"N");
				$cdsStopFinal = $chromSize;
			}
			
			$ntsAtEitherEnd{$i}{"start"} = $startString;
			$ntsAtEitherEnd{$i}{"end"}   = $endString;
			my $queryID = "$i#$species#$cdsStartFinal#$cdsStopFinal#$chrQuery#$strandQuery";
			
			print FO "$chrQuery\t$cdsStartFinal\t$cdsStopFinal\t$queryID\n";
			print "###\n\n bedEntry is $chrQuery\t$cdsStartFinal\t$cdsStopFinal\t$queryID\n####\n\n" if ($verbose);
			$i++;
			$seqCt++;
		}
		close FO;

		next if ($seqCt == 0);
			
		## runTwoBitToFa on this bedFile
		my $twoBitToFaCall = "twoBitToFa $twoBitFile $speciesFasta -bed=$speciesBed";
		system($twoBitToFaCall) == 0 || die "Error running twoBitToFa call for species '$species'\n";
		
		## store the fasta returned by twoBitToFa in a hash
		my $ref2seqFastaHash = fasta2Hash($speciesFasta);
		
		foreach my $keys(keys(%$ref2seqFastaHash)) {
			my $sequence = $ref2seqFastaHash->{$keys};
			my @keyArray = split(/#/,$keys);
			my $exonNumber = shift(@keyArray);
			my $header = join("#",@keyArray);
			
			my $strandQuery = $strandHash{$species}{$exonNumber};
			my $startString = $ntsAtEitherEnd{$exonNumber}{"start"};
			my $endString   = $ntsAtEitherEnd{$exonNumber}{"end"};
			
			my($sequenceFinal,$headerFinal) = adjustForAssemblyGaps($sequence,$header,$strandRef,$startString,$endString);
			$seqHash{$exonNumber}{$species}     =  $sequenceFinal;
			$headerHash{$exonNumber}{$species}  =  $headerFinal;
		}
	}

	## Now run CESAR exon-wise by using the sequences stored in seqHash
	

	# create the output dir for all species and remove a potentially existing output file
	foreach my $species(@speciesArray) {
		my $createDirCall = "mkdir -p $outDir/$chrRef/$species";
		system($createDirCall) == 0 || die "Error: '$createDirCall' failed\n";
		system("rm -f $outDir/$chrRef/$species/$gene");
	}

	my %bedDirHash = ();
	my $k = 1;
	while($k <= $noe) {
		
		if (exists $exonHasStop{$k}) {
			print "Omitting exon '$k' because it has a stop codon\n" if ($verbose);
			$k++;
			safeDie("NA");
			next;
		}
		
		my ($exonStart,$exonEnd) = (split /-/,$exonsList[$k-1]);
		open(FO,">$cesarInput") || die "Error writing to CESAR input file '$cesarInput'\n";
		print FO ">referenceExon\n$refSeqHash{$k}\n#####\n";
		
		my $seqCt = 0;
		foreach my $species(@speciesArray) {
			print FO ">$headerHash{$k}{$species}\n$seqHash{$k}{$species}\n" if(exists $seqHash{$k}{$species});
			$seqCt++;
		}
		close FO;
		next if ($seqCt == 0);
		
		my $accProfile = "extra/tables/$clade/acc_profile.txt";
		my $doProfile = "extra/tables/$clade/do_profile.txt";
		$accProfile = "extra/tables/$clade/firstCodon_profile.txt" if ($k == 1);
		$doProfile = "extra/tables/$clade/lastCodon_profile.txt" if ($k == $noe);			
		$doProfile = "extra/tables/$clade/u12_donor_profile.txt" if ($referenceSpliceSites{$k}{"do"} eq "AT");
		$accProfile = "extra/tables/$clade/u12_acceptor_profile.txt" if ($referenceSpliceSites{$k}{"acc"} eq "AC");
		
		$doProfile  = "$cesarDir/$doProfile";
		$accProfile = "$cesarDir/$accProfile";
		my $matrix = "$cesarDir/extra/tables/human/eth_codon_sub.txt";	

		## runCESAR
		my $cesarCall = "$cesarDir/cesar $cesarInput --matrix $matrix -p $accProfile $doProfile --max-memory $maxMem > $cesarOutput";
		if ($noe == 1) {
			$cesarCall = "$cesarDir/cesar $cesarInput -f -l --matrix $matrix -p $accProfile $doProfile --max-memory $maxMem > $cesarOutput";
		} else {
			if ($k == 1) {
				$cesarCall = "$cesarDir/cesar $cesarInput -f --matrix $matrix -p $accProfile $doProfile --max-memory $maxMem > $cesarOutput";
			} elsif ($k == $noe) {
				$cesarCall = "$cesarDir/cesar $cesarInput -l --matrix $matrix -p $accProfile $doProfile --max-memory $maxMem > $cesarOutput";
			}
		}
		system($cesarCall) == 0 || die "Error running cesar on exon, '$cesarCall' failed\n";
		print "Running cesar call '$cesarCall'\n" if ($verbose);

		## Split the alignment into pairwise alignments and check if the exon is intact for each species
		open(FIA,"$cesarOutput") || die "Error opening cesarOutput file for exon number '$k'\n";
		my @alignmentArray = <FIA>;
		close FIA;
		
		for(my $i = 0; $i < scalar(@alignmentArray); $i++) {
			if ($alignmentArray[$i] =~/referenceExon/) {
				my $refSequence   = $alignmentArray[$i+1]; chomp $refSequence;
				my $queryHeader   = $alignmentArray[$i+2]; chomp $queryHeader;
				my $querySequence = $alignmentArray[$i+3]; chomp $querySequence;
				$querySequence = chopNAtStartOrEnd($querySequence);
				my $species = (split /#/,$queryHeader)[0];
				$species =~s/>//;
				
				print "Alignment for exon $k. gene '$gene'\n$refSequence\n$querySequence\n####\n" if ($verbose);
				
				# create output subdirectory
				open(FOS,">>$outDir/$chrRef/$species/$gene") || die "Error: cannot write to $outDir/$chrRef/$species/$gene file\n";
				my $accRef = $referenceSpliceSites{$k}{"acc"};
				my $doRef  = $referenceSpliceSites{$k}{"do"};
			
				my $exonCDS = checkExonIntactness($refSequence,$querySequence,$queryHeader,$strandRef,$gene,$k,$exonPhasesList[$k-1],$accRef,$doRef,$noe,$verbose);
			
				if ($exonCDS ne "NI") {  ## for exons that are not intact - getIntactExons.pl returns "NI";
					
					my($chrQuery,$exonStartQuery,$exonEndQuery,$strandQuery,$exonN) = (split /\t/,$exonCDS);
					my $queryKey = "$species,$chrQuery,$exonStartQuery,$exonEndQuery,$strandQuery";
					my $exonKey = "$gene\_$exonN";
					print FOS "$chrRef\t$exonStart\t$exonEnd\t$queryKey\t$exonKey\n";
					print "The exon '$k' from gene '$gene' is INTACT in the species '$species'\n" if ($verbose);
				} else {
					print "The exon '$k' from gene '$gene' is NON-INTACT in species '$species'\n" if ($verbose);
				}
				close FOS;
				print "DONE\n\n#####\n" if ($verbose);
				
				$i = $i+3;
			}
		}
		
		$k++;
	}
}
safeDie("NA"); ## this call simply deletes all temp files


########################################################
####        End of the main body of the code        ####
####               Subroutines begin                ####
########################################################

sub get_maf_species_list {  ## given a maf, this function runs mafSpeciersList and returns reference to an array that contains the species in the maf:
	my $maf = shift;	
	## Get the list of species present in maf
	my @speciesList = ();			### create a species_array which contains the list of all the species present in the maf block.
	@speciesList = `mafSpeciesList $maf stdout`; 
	chomp(@speciesList);
	$" = " ";	# separate array elements by a space in the output and remove the reference from the list
	return \@speciesList;
}

sub getCoordinatesFromMaf { ## This function gives me the coordinates (the start/stop/chr/strand) for the different species in the maf

	my ($maf,$reference,$ref2TempFilesList) = @_;

	open(FIM,"$maf") || die "Error !! Cannot open mafFile '$maf' in getCoordinatesFromMaf function\n";
	my @mafExtractResult = <FIM>;
	close FIM;

	my $ref_species_list = get_maf_species_list($maf);
	my @speciesList = @$ref_species_list;
	
	@speciesList = grep{$_ ne $reference} @speciesList;

	my %cdsStartHash = my %cdsStopHash = ();
	my %chrQuery = my %strandQuery = ();
	my %chrCountHash = my %strandCountHash = ();
	my %cdsPrint = my %seqLength = my %ct = ();

	my %cdsListHash = ();
	foreach my $species(@speciesList) {
		my @cdsList = ();
		$cdsListHash{$species} 	= \@cdsList;
		$strandQuery{$species}	= "NA";
		$chrQuery{$species}		= "NA";
		
		$seqLength{$species} = 0;
		$ct{$species} = 0; ## A species specific counter
		$cdsPrint{$species} = "NA";
	}
	
	foreach my $line(@mafExtractResult) {
		$line =~s/\s+$//;
		next if ($line =~/$reference/);
		if ($line =~/^s/) {	## For s lines
			my ($species,$chr,$start,$stop,$strand,$srcSize,$seq) = getSLineData($line);
			
			if ($species ne $reference) {
				$chrQuery{$species} = $chr;
				$strandQuery{$species} = $strand;
				
				my $cdsStartQuery = my $cdsStopQuery = "";
																	## An example is
																	## ##maf version=1 scoring=mafExtract
				if ($strand eq "-") {								## a score=0.000000
																	## s hg19-galGal4.chr19  3261974 109 + 61431566 TTCAGGCAGATGGTGGCTGAGGCAGTAGCGGTGGCTACAGTGCATGCAGAACTGGCCCAGAGTGGTGGTGCTGGCCGAGCACTTGGAGAAGCTACAGGTGTTGTCAGCC	
					$cdsStartQuery = $srcSize - $start - $stop;		## s galGal4.chr19      49577681 109 - 51267647 TTCAGGCAGATGGTGGCTGAGACAGTAGCGGTGGCTACAGTGCATGCAGAACTGCCCTAGGGTGGTGGTACTGGCTGAACACTTGGAGAAGCTACAGGTCTTGTCAGCC
					$cdsStopQuery  = $srcSize - $start;				##	
				} else {												## a score=0.000000
					$cdsStartQuery = $start;						## s hg19-galGal4.chr19  3262083 61 + 61431566 ---TTCACCACAGCAGACACCAGGGCGTCGAAGTCCTCCTCACAGGGCAGAGCCATAACTGGAC
					$cdsStopQuery  = $start + $stop;				## s galGal4.chr19      49577790 64 - 51267647 TTCTTCACCACAGCAGACACCAGGGCGTCAAAATCCTCCTCACAGGGCAGTGCCATAACTGGAC
				}													##
																	## startFirstBlock  = 51267647 - 49577681 - 109		==>	1689857		
																	## stopFirstBlock   = 51267647 - 49577681			==> 	1689966	-- We want this
																	##
																	## startLastBlock  = 51267647 - 49577790 -64		==> 	1689793	-- We want this
																	## stopLastBlock   = 51267647 - 49577790			==> 	1689857
				my $string = "$cdsStartQuery-$cdsStopQuery";
				
				my $ref2CdsList = $cdsListHash{$species};  ## get the relevant cdsList
				my @cdsList = @$ref2CdsList;
				push(@cdsList,$string);  ## push the string to the cdsList
				$cdsListHash{$species} = \@cdsList;  ## and update the value of the cdsList in the cdsListHash
				
				$chrCountHash{$species}{$strand}++;
				$strandCountHash{$species}{$chr}++;
				
				$cdsPrint{$species} = $start if ($ct{$species} == 0); ## Get the beginning of the CDS. i.e. the first block from the maf
				$seqLength{$species} = $seqLength{$species} + $stop; ## Get the length of the aligning sequence through all blocks in the maf
				$ct{$species}++;
			}
		}
	}
	close FIM;
	
## + strand is straight-forward as always. Minus is harder, as always
##maf version=1 scoring=mafExtract
# a score=0.000000
# s hg19.chr11 118532101 14 + 135006516 TGTAGAAAGGCGGT
# s mm9.chr9    79584924 14 - 124076172 TGAAGGAAGGCGAC
# i mm9.chr9   N 0 C 0

# a score=2010932.000000
# s hg19.chr11 118532115 53 + 135006516 GTCATTGGTGTGAGTCAAGTAGCAATCCATCATGAGGGTCAAGAGTGGGGGCT
# s mm9.chr9    79584938 53 - 124076172 ATCCTTGGTATGAGCTACATATCGATCCATCATGAGAGTCAGGAGTGGGGGCT
# i mm9.chr9   C 0 C 0

# a score=0.000000
# s hg19.chr11 118532168 50 + 135006516 GGCTCCGCTGCAGGTAGTACACGCGCCCACCATTGGGGACATGCCCATAG
# s mm9.chr9    79584991 50 - 124076172 GGCTCCGTTGCAGGTAATATATGCGTCCACCGTTGGGGATATGTCCGTAG
# i mm9.chr9   C 0 N 0
	
## With the code above, the firstElement (see below) is '44491234-44491248' and the last element is  '44491131-44491181'. The cdsStart is 44491234 and stop is 44491181. 
## This needs to be flipped. So I reverse the array for minus strand
	
	my %cdsStartQuery = my %cdsStopQuery = ();
	foreach my $species(@speciesList) {
		my $ref2CdsList = $cdsListHash{$species};
		my @cdsList = @$ref2CdsList;
		
		if (scalar(@cdsList) > 0) {
			@cdsList = reverse(@cdsList) if ($strandQuery{$species} eq "-");

			my $firstElement = $cdsList[0];
			my $lastElement  = $cdsList[$#cdsList];
		
			$cdsStartQuery{$species} = (split /-/,$firstElement)[0];
			$cdsStopQuery{$species}  = (split /-/,$lastElement)[1];
		} else {
			$cdsStartQuery{$species} = "NA";
			$cdsStopQuery{$species} = "NA";
		}
	}
	
	my %moreThanOneStrand = my %moreThanOneChr = ();
	foreach my $species(@speciesList) {
		$moreThanOneStrand{$species} = "F";
		$moreThanOneChr{$species} 	 = "F";
	}
	
	## check which species come from more than one strand, species which more come from than one chromosome/scaffold and species which have nothing but eLines
	foreach my $species(@speciesList) {
		my @chrList 	= keys(%{$chrCountHash{$species}});
		my @strandList  = keys(%{$strandCountHash{$species}});
		
		$moreThanOneStrand{$species} = "T"  if (scalar(@chrList) > 1);
		$moreThanOneChr{$species} = "T"     if (scalar(@strandList) > 1);
	}
	
	return(\%cdsStartQuery,\%cdsStopQuery,\%strandQuery,\%chrQuery,\%moreThanOneStrand,\%moreThanOneChr,\@speciesList,\%cdsPrint,\%seqLength);
}


sub safeDie {  ## cleans all temp files when the argument is "NA", otherwise cleans all temp files and dies after printing the argument(the error message)

	my $message = shift;
	my $allTempFiles = join(" ",@tempFiles);
	`rm -rf $allTempFiles`;
	die $message if ($message ne "NA");
}

sub adjustForAssemblyGaps {
	my ($sequence,$queryHeader,$strandRef,$startString,$endString) = @_;
	## if the sequence has Ns at the beginning or the end (assembly gaps), strip them right here and adjust the coordinates accordingly
	my $seqQueryCopy1 = $sequence;
	my $seqQueryCopy2 = $sequence;
	
	my($species,$cdsQueryStart,$cdsQueryStop,$chrQuery,$strandQuery) = (split /#/,$queryHeader);
	
	$seqQueryCopy1 =~s/^N+//;
	$cdsQueryStart = $cdsQueryStart + (length($sequence)  - length($seqQueryCopy1)); ## adjusting coordinates
	$seqQueryCopy2 =~s/N+$//;
	$cdsQueryStop = $cdsQueryStop - (length($sequence)  - length($seqQueryCopy2)); ## adjusting coordinates
	
	$sequence =~s/^N+//;
	$sequence =~s/N+$//;
	$sequence = revComp($sequence) if ($strandRef ne $strandQuery); 
	
	## Fix the strand orientation
	my $queryID = "$species#$cdsQueryStart#$cdsQueryStop#$chrQuery#$strandQuery";

	$sequence = $startString.$sequence.$endString;
	$sequence = uc($sequence);
	return($sequence,$queryID);
}


sub fasta2Hash {
    my $fastaInput = shift;
    
    my $id = my $seq = "";
    my %seqHash = ();

    open(FIS,"$fastaInput") || die "Error opening input file '$fastaInput'\n";
    while (my $line=<FIS>) {
		$line =~s/\s+$//;

		if ($line =~/>/) {
			if ($seq ne "") {
				$seqHash{$id} = $seq;
				$id = $seq = "";
			}

			$id =$line;
			$id =~s/>//;
		} else { 
			$seq = $seq.$line;
		}
    }
    close FIS;

    $seqHash{$id} = $seq;  ## For the last sequence
    return \%seqHash
}


sub getExonPhaseList {		## This function gives the phases (5' and 3') of each exon.
	my $ref2ExonsList = shift;
	my @exonsList =  @$ref2ExonsList;
	
	my $i = my $lengthCum = 0;
	my @exonLength = my @exonPhaseList = ();
	
	foreach my $exon(@exonsList) {
		my ($start,$stop) = (split /-/,$exon)[0,1];
		my $length = abs($start-$stop);
		$lengthCum = $lengthCum + $length;
		
		$exonLength[$i] = $lengthCum;
		$i++;
	}
	
	for (my $i = 0;$i < scalar(@exonLength); $i++) {
		my $exonLengthCum = $exonLength[$i];	
		my $phase5Prime = my $phase3Prime = "";
		
		if ($i != 0) {
			$phase5Prime = 3 - ($exonLength[$i-1]%3);
			$phase3Prime = ($exonLength[$i]%3);
		} else {
			$phase5Prime = 0;
			$phase3Prime = ($exonLength[$i]%3);
		}
		
		$phase5Prime = 0 if ($phase5Prime == 3);
		$phase3Prime = 0 if ($phase3Prime == 3);
		
		my $phaseList = "$phase5Prime-$phase3Prime";
		$exonPhaseList[$i] = $phaseList;
	}
	return \@exonPhaseList;
}


sub createString {  ## This function creates a string given the character that makes this string and the length
	my ($length,$char) = @_;
	my $string = "";
	
	my $i = 1;
	while ($i <= $length) {
		$string = $string.$char;
		$i++;
	}
	return $string;
}

sub getSequenceFrom2Bit {  ## This function gets sequence from2  bit file given the twoBit file, the chrom, the start and the stop coordinate
	my ($twoBitFile,$chr,$start,$stop) = @_;
	my $seq = `twoBitToFa $twoBitFile -seq=$chr -start=$start -end=$stop stdout|grep -v ">"|tr -d "\n"`;
	die "Cannot get the sequence from 2BitToFa" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);	
	return $seq;
}

sub revComp {  ## returns the recoverse complement of a given sequence
	my $seq = shift;
	$seq =~ tr/ATGCatgc-/TACGtacg-/;
	$seq = reverse($seq);
	return $seq;
}

sub fixPhase { 	## This function does the following: If the input sequence is GTGCATCGAAT, phase5Prime is 1, phase3Prime is 1, the function returns gTGCATCGAAt
	my ($seq,$phase5Prime,$phase3Prime) = @_;  
	$seq = uc($seq);	## Make entire sequence upper-case before you do anything to it
	
	if ($phase5Prime != 0) { ## Fix 5' end
		my $s1 = substr($seq,0,$phase5Prime);
		my $s2 = substr($seq,$phase5Prime);
		
		$seq = lc($s1)."|".$s2;
	} else {
		$seq = "|".$seq;
	}
	
	if ($phase3Prime != 0) { ## Fix 3' end
		my $s1 = substr($seq,0,(length($seq) - $phase3Prime));
		my $s2 = substr($seq,(length($seq) - $phase3Prime));
		
		$seq = $s1."|".lc($s2);
	} else {
		$seq = $seq."|";
	}
	return $seq;
}


sub getChromSize {  ## This function returns the chromosome size for a given species/chromosome pair 
	my ($species,$chrInterest) = @_;
	
	my $chromSizeFile = "$twoBitDir/$species/chrom.sizes";
	my $chrSize = "";
	
	open(FICS,"$chromSizeFile") || die "Error opening chromosomes size file '$chromSizeFile' in getChromSize function\n";
	while(my $line = <FICS>) {
		$line =~s/\s+$//;
		my($chr,$size) = (split /\t/,$line);
		if ($chr eq $chrInterest) {
			$chrSize = $size;
			last;
		}
	}
	close FICS;
	return $chrSize;
}

sub chopNAtStartOrEnd {
	my $sequence = shift;
	my $start = getNCharAtEnd($sequence);
	my $seqRev = reverse($sequence);
	my $stop  = getNCharAtEnd($seqRev);

	my $beginString = createString($start,"-");
	my $endString   = createString($stop,"-");
	$sequence = $beginString.substr($sequence,$start,length($sequence)-$start-$stop).$endString;
	return $sequence;
}

sub getNCharAtEnd {
	my $seq = shift;
	my @seqSplit = split(/n|N/,$seq);
	
	my $start = 0;
	foreach my $e(@seqSplit) {
		if ($e ne "") {
			last;
		}
		$start++;
	}
	return $start;
}
