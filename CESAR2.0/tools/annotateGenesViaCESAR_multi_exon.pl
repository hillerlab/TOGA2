#!/usr/bin/env perl

## Virag Sharma, 2018. MPI-CBG and MPI-PKS.
## The CESAR2.0 multi-exon annotation engine

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
my $max_mem = 30;
GetOptions ("v|verbose" => \$verbose, "maxMemory=i" => \$max_mem, "list" =>\$list);

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

my $input         = $ARGV[0];
my $maf_index     = $ARGV[1];
my $gene_pred     = $ARGV[2];
my $reference     = $ARGV[3];
my $species_list  = $ARGV[4];
my $out_dir       = $ARGV[5];
my $two_bit_dir   = $ARGV[6];
my $cesar_dir     = $ARGV[7];

## Check arguments
die "The mafIndex file '$maf_index' does not exist\n" if (! -e $maf_index);
die "The genePred file '$gene_pred' does not exist\n" if (! -e $gene_pred);
die "The genes list file '$input' does not exist, you specified that the input is a list\n" if (! -e $input && $list != 0);

my $clade = "human";

my @species_array = split(",",$species_list);
my %species_interest = map{$_ => 1}@species_array;

## Step 0: check if two bit files are present for all species of interest
foreach my $species(@species_array) {
	my $two_bit_file = "$two_bit_dir/$species/$species.2bit";
	die "The two bit file ('$two_bit_file') for species '$species' is not present. Aborting\n" if (! -e $two_bit_file);
}
my $two_bit_file_ref = "$two_bit_dir/$reference/$reference.2bit";

my %genes_interest = ();
if ($list) {
	open(FI,$input) || die "Error opening genesList file '$input'\n";
	while (my $gene = <FI>) {
		chomp $gene;
		$genes_interest{$gene} = "T";
	}
	close FI;
} else {
	$genes_interest{$input} = "T";
}

## Step1: Get the information about the input transcripts/genes from the genePred file
my %gene_info_hash = ();
open(FI,$gene_pred) || die "Error opening genePred file\n";
while (my $line = <FI>) {
	$line =~s/\s+$//;
	my ($gene,$chr,$strand,$exons) = (split /\t/,$line)[0,1,2,7];
	if (exists $genes_interest{$gene}) {
		$gene_info_hash{$gene} = "$chr#$strand#$exons";
	}
}
close FI;

## check if there is nothing in the geneInfoHash
my @genes_int = keys(%gene_info_hash);
if (scalar(@genes_int) == 0) {
	print "############## ERROR ##############\n";
	print "No information found for the input gene '$input' from the genePred file '$gene_pred'\n";
	print "Could it be that the input is a file and you forgot to specify the '-list' parameter\nCheck and re-run\n";
	print "###################################\n";
	exit;
}

my @temp_files_all = ();  ## create all temp files upfront and store them in @tempFiles array
my $maf_file      = `mktemp /dev/shm/exon.maf.XXXXXXXXX`; chomp $maf_file;
my $cesar_input   = `mktemp /dev/shm/cesarIn.fa.XXXXXX`; chomp $cesar_input;
my $cesar_output  = `mktemp /dev/shm/cesarOut.aln.XXXXXX`; chomp $cesar_output;
push(@temp_files_all,$maf_file,$cesar_input,$cesar_output);
	
foreach my $gene(keys(%gene_info_hash)) {
	print "Processing gene '$gene'\n";
	my($chr_ref,$strand_ref,$all_exons) = (split /#/,$gene_info_hash{$gene});

	print "The gene '$gene' comes on chromosome '$chr_ref', strand '$strand_ref'. Exons are '$all_exons'\n" if ($verbose);
	my @exons_list = split(",",$all_exons);
	@exons_list = reverse(@exons_list) if ($strand_ref eq "-");
	my $noe = scalar(@exons_list);
	
	my $ref_exon_phase_list = get_exon_phase_list(\@exons_list);
	my ($ref_reference_splice_sites,$n_query_seqs) = get_cesar_multi_exon_input($gene,\@exons_list,$ref_exon_phase_list,$chr_ref,$strand_ref,$species_list,$cesar_input,$maf_file);
	
	if ($n_query_seqs == 0) {
		print "Nothing in the maf for the query species (the entire gene could be deleted)|The query locus is not colinear, so CESAR-multi exon is not suitable\n";
		next;
	}

	## Now run CESAR on the input;
	my $cesar_call = "$cesar_dir/cesar $cesar_input --max-memory $max_mem > $cesar_output";
	system($cesar_call) == 0 || die "Error running cesar for the gene '$gene', '$cesar_call' failed\n";
	print "Running cesar call '$cesar_call'\n" if ($verbose);	
	## Split the alignment into pairwise alignments and check if the exon is intact for each species
	open(FIA,"$cesar_output") || die "Error opening cesar_output file $cesar_output\n";
	my @alignment_array = <FIA>;
	close FIA;
		
	for(my $i = 0; $i < scalar(@alignment_array); $i++) {
		if ($alignment_array[$i] =~/referenceExon/) {
			my $ref_sequence   = $alignment_array[$i+1]; chomp $ref_sequence;
			my $query_header   = $alignment_array[$i+2]; chomp $query_header;
			my $query_sequence = $alignment_array[$i+3]; chomp $query_sequence;
				
			$query_header =~s/>//;
			my $species = (split /#/,$query_header)[0];
			# create output subdirectory and remove a potentially existing output file
			`mkdir -p $out_dir/$chr_ref/$species/`;
			`rm -rf $out_dir/$chr_ref/$species/$gene`;

			open(FOS,">>$out_dir/$chr_ref/$species/$gene") || die "Error: cannot write to $out_dir/$chr_ref/$species/$gene file\n";
			## Now split the multi-exon alignment into individual exon alignments
			my $ref_split_exon_alignments = split_multi_exon_to_single_exon($ref_sequence,$query_sequence);
			my @split_exon_alignments = @$ref_split_exon_alignments;
			
			my $k = 1;
			foreach my $exon_aln(@split_exon_alignments) {
				my ($ref_seq_exon,$query_seq_exon) = split(/\n/,$exon_aln);
				my $acc_ref = $ref_reference_splice_sites->{$k}{"acc"};
				my $do_ref  = $ref_reference_splice_sites->{$k}{"do"};
				
				my $exon_start = (split /-/,$exons_list[$k-1])[0];
				my $k_ori = $k;

				if ($ref_seq_exon =~/>>>>/) {  ## this would be the case if there is an intron deletion which is indicated by ">>>" in reference sequence
					my ($ref_seq_exon_new,$query_seq_exon_new,$n_exon_merged) = intron_deletion_correction($ref_seq_exon,$query_seq_exon);
					$do_ref  = $ref_reference_splice_sites->{$k+$n_exon_merged - 1}{"do"};
					$k = $k + $n_exon_merged - 1;
					$ref_seq_exon   = $ref_seq_exon_new;
					$query_seq_exon = $query_seq_exon_new;
				}
				my $exon_end = (split /-/,$exons_list[$k-1])[1];
				
				my $exon_cds = checkExonIntactness($ref_seq_exon,$query_seq_exon,$query_header,$strand_ref,$gene,$k,$ref_exon_phase_list->[$k-1],$acc_ref,$do_ref,$noe,$verbose);
				if ($exon_cds ne "NI") {  ## for exons that are not intact - checkExonIntactness returns "NI";		
					my($chr_query,$exon_start_query,$exon_end_query,$strand_query,$exon_number) = (split /\t/,$exon_cds);
					my $query_key = "$species,$chr_query,$exon_start_query,$exon_end_query,$strand_query";
					my $exon_key = "$gene\_$k";
					$exon_key = "$gene\_$k_ori\_$k" if ($k != $k_ori);
 
					print FOS "$chr_ref\t$exon_start\t$exon_end\t$query_key\t$exon_key\n";
					print "The exon '$k' from gene '$gene' is INTACT in the species '$species'\n" if ($verbose);
				} else {
					print "The exon '$k' from gene '$gene' is NON-INTACT in species '$species'\n" if ($verbose);
				}
				print "DONE\n\n#####\n" if ($verbose);		
				
				$k++ if ($ref_seq_exon !~/>>>>>/);  ## if there was a case of intron deletion, then k has already been incremented.
			}
			close FOS;
			$i = $i+3;
		}
	}
}
safeDie("NA"); ## this call simply deletes all temp files


########################################################
####        End of the main body of the code        ####
####               Subroutines begin                ####
########################################################

sub get_exon_phase_list {		## This function gives the phases (5' and 3') of each exon.
	my $ref_exons_list = shift;
	my @exons_list =  @$ref_exons_list;
	
	my $i = my $length_cum = 0;
	my @exon_length = my @exon_phase_list = ();
	
	foreach my $exon(@exons_list) {
		my ($start,$stop) = (split /-/,$exon)[0,1];
		my $length = abs($start-$stop);
		$length_cum = $length_cum + $length;
		
		$exon_length[$i] = $length_cum;
		$i++;
	}
	
	for (my $i = 0;$i < scalar(@exon_length); $i++) {
		my $exon_length_cum = $exon_length[$i];	
		my $phase_5Prime = my $phase_3Prime = "";
		
		if ($i != 0) {
			$phase_5Prime = 3 - ($exon_length[$i-1]%3);
			$phase_3Prime = ($exon_length[$i]%3);
		} else {
			$phase_5Prime = 0;
			$phase_3Prime = ($exon_length[$i]%3);
		}
		
		$phase_5Prime = 0 if ($phase_5Prime == 3);
		$phase_3Prime = 0 if ($phase_3Prime == 3);
		
		my $phase_string = "$phase_5Prime-$phase_3Prime";
		$exon_phase_list[$i] = $phase_string;
	}
	return \@exon_phase_list;
}

sub get_cesar_multi_exon_input {  ## Returns the input file for CESAR-multi-exon
	my($gene,$ref_all_exons,$ref_exon_phase_list,$chr_ref,$strand_ref,$species_string,$out_file,$tmp_maf) = @_;
	
	my %reference_splice_site = ();
	
	my @all_exons = @$ref_all_exons;
	my $noe = scalar(@all_exons);
	my @species_list = split(",",$species_string);
	my %species_in_maf = ();	

	my $n_query_seqs = 0;
	my $i = 1;
	my %cds_start_hash = my %cds_stop_hash = my %strand_hash = my %chr_hash = my %species_ignore = ();
	open(FO,">$out_file") || die "Error writing to output file '$out_file' in get_cesar_multi_exon_input\n";
	foreach my $ex(@all_exons) {
		my($start,$stop) = (split /-/,$ex)[0,1];
		my($phase_5P,$phase_3P) = (split /-/,$ref_exon_phase_list->[$i-1])[0,1];
		
		my $start_ss = $start - 2;
		my $stop_ss  = $stop + 2;
		
		## get exonic sequence by running twoBitToFa:
		my $exon_sequence_ref = get_sequence_from_2bit($two_bit_file_ref,$chr_ref,$start_ss,$stop_ss);
		$exon_sequence_ref = revComp($exon_sequence_ref) if ($strand_ref eq "-");
		
		my $acc_site = lc(substr($exon_sequence_ref,0,2));
		my $do_site  = lc(substr($exon_sequence_ref,length($exon_sequence_ref)-2,2));
		
		$reference_splice_site{$i}{"acc"} = $acc_site;
		$reference_splice_site{$i}{"do"} = $do_site;

		my $acceptor_profile = "extra/tables/$clade/acc_profile.txt";
		my $donor_profile    = "extra/tables/$clade/do_profile.txt";
		$acceptor_profile    = "extra/tables/$clade/firstCodon_profile.txt" if ($i == 1);
		$donor_profile       = "extra/tables/$clade/lastCodon_profile.txt" if ($i == $noe);			
		$acceptor_profile = "extra/tables/$clade/u12_acceptor_profile.txt" if ($acc_site eq "AC");
		$donor_profile    = "extra/tables/$clade/u12_donor_profile.txt" if ($do_site eq "AT");
	
		## Fix the sequence's split codon bases
		$exon_sequence_ref = substr($exon_sequence_ref,2,length($exon_sequence_ref)-4); 
		$exon_sequence_ref = fix_phase($exon_sequence_ref,$phase_5P,$phase_3P);
		print FO ">$gene\_$i\t$acceptor_profile\t$donor_profile\n$exon_sequence_ref\n";
		
		## Now extract the coordinates
		my $call = "mafExtract $maf_index -region=$chr_ref:$ex stdout|mafSpeciesSubset stdin speciesList=$reference,$species_string species.lst=NULL $tmp_maf";
		system($call) == 0 || die "Error running mafExtract call, '$call' failed\n";
		my ($ref2cdsStartQuery,$ref2cdsStopQuery,$ref2strandQuery,$ref2chrQuery,$ref2moreThanOneStrand,$ref2moreThanOneChr,$ref2speciesPresent) = getCoordinatesFromMaf($tmp_maf,$reference,\@temp_files_all);
		
		my @species_present = @$ref2speciesPresent;
		my %species_present_in_maf = map{$_ => 1}@species_present;

		foreach my $species(@species_list) {
			if (exists $species_present_in_maf{$species}) {
				## Get the attributes for the species that are present in the maf
				my $cds_start      = $ref2cdsStartQuery->{$species};
				my $cds_stop       = $ref2cdsStopQuery->{$species};
				my $strand_query   = $ref2strandQuery->{$species};
				my $chr_query      = $ref2chrQuery->{$species};
				my $strand_plural  = $ref2moreThanOneStrand->{$species};
				my $chr_plural     = $ref2moreThanOneChr->{$species};
	
				print "Exon # $i, the reference coordinates are $chr_ref:$ex\nQuery coordinates for '$species' ==> start is '$cds_start' .. stop is '$cds_stop' .. strand is '$strand_query' .. chromosome is '$chr_query'\n\n" if ($verbose);
				$species_ignore{$species} = "T" if ($strand_plural eq "T" || $chr_plural eq "T");
		
				$cds_start_hash{$species}{$i} = $cds_start;
				$cds_stop_hash{$species}{$i}  = $cds_stop;
				$strand_hash{$species}{$i}    = $strand_query;
				$chr_hash{$species}{$i}       = $chr_query;
				$species_in_maf{$species}     = "T";
			} else {
				$cds_start_hash{$species}{$i} = "NA";
				$cds_stop_hash{$species}{$i}  = "NA";
				$strand_hash{$species}{$i}    = "NA";
				$chr_hash{$species}{$i}       = "NA";
			}
		}
		$i++;	
	}
	print FO "####\n";

	## Now check if all the exons are co-linear. If yes, extract the query sequence and print to the file.
	my @final_species_list = ();
	foreach my $species(@species_list) {
		if (! exists $species_in_maf{$species}) {
			print "Species '$species' ignored because there is nothing in the maf for this species\n" if ($verbose);
			next;
		}

		if (exists $species_ignore{$species}) {
			print "Species '$species' ignored because atleast one of the exons is funky\n" if ($verbose);
			next;
		}
		
		my $query_2bit_file = "$two_bit_dir/$species/$species.2bit"; 
		my @chr_list    = values (%{$chr_hash{$species}});
		my @strand_list = values (%{$strand_hash{$species}});

		my ($chr_query,$chr_count)       = get_unique_count(\@chr_list);
		my ($strand_query,$strand_count) = get_unique_count(\@strand_list);
	
		if ($chr_count == 1 && $strand_count == 1) { ## i.e. the gene is co-linear
			my ($query_sequence,$header) = get_query_sequence($species,\%cds_start_hash,\%cds_stop_hash,$strand_ref,$strand_query,$chr_query,$query_2bit_file,$noe);
			if ($query_sequence eq "No sequence") {
				print "Omitting species '$species' because get_query_sequence did not return a valid query locus.\n" if ($verbose);
			} else {
				print FO ">$header\n$query_sequence\n";
				$n_query_seqs++;	
				push(@final_species_list,$species);
			}
		} else {
			print "Omitting species '$species' because the query locus comes from scaffolds: '@chr_list' and strand: '@strand_list'\n";
		}
	}
	close FO;
	print "CESAR-multi exon mode will be run for the following species '@final_species_list'\n" if ($verbose);
	return (\%reference_splice_site,$n_query_seqs);
}

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

	return(\%cdsStartQuery,\%cdsStopQuery,\%strandQuery,\%chrQuery,\%moreThanOneStrand,\%moreThanOneChr,\@speciesList);
}


sub safeDie {  ## cleans all temp files when the argument is "NA", otherwise cleans all temp files and dies after printing the argument(the error message)

	my $message = shift;
	my $allTempFiles = join(" ",@temp_files_all);
	`rm -rf $allTempFiles`;
	die $message if ($message ne "NA");
}

sub get_sequence_from_2bit {  ## This function gets sequence from2  bit file given the twoBit file, the chrom, the start and the stop coordinate
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

sub fix_phase { 	## This function does the following: If the input sequence is GTGCATCGAAT, phase5Prime is 1, phase3Prime is 1, the function returns gTGCATCGAAt
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

sub get_unique_count {
	## this function checks if a list contains only one unique value
	my $ref_input_list = shift;
	my @input_list = @$ref_input_list;
	
	## Remove "NA" from this list
	@input_list = grep{$_ ne "NA"}@input_list;
	
	my %input_hash = map{$_ => 1}@input_list;
	my @input_list_unique = keys(%input_hash);
	my $count = scalar(@input_list_unique);
	return ($input_list_unique[0],$count);
}

sub get_query_sequence {
	my($species,$ref_start_cds_hash,$ref_stop_cds_hash,$strand_ref,$strand_query,$chr_query,$query_2bit_file,$noe) = @_;
	
	my $cds_start = my $cds_stop = "";

	## if strands do not match, then invert the values:
	if ($strand_query ne $strand_ref) {
		$ref_start_cds_hash = invert_values($ref_start_cds_hash,$species,$noe);
		$ref_stop_cds_hash  = invert_values($ref_stop_cds_hash,$species,$noe);
	}

	if ($strand_query eq "+") {
		$cds_start = $ref_start_cds_hash->{$species}{1};
		$cds_stop  = $ref_stop_cds_hash->{$species}{$noe};
	} else {
		$cds_start = $ref_stop_cds_hash->{$species}{$noe};
		$cds_stop  = $ref_start_cds_hash->{$species}{1};
	}
	
	if ($cds_start eq "NA" || $cds_stop eq "NA") {
		return "No sequence"; ## this implies that either the first exon is deleted or the last exon is deleted or has unaligning sequence
	}                             ## and as result, the query locus cannot be proprely inferred.
	
	## otherwise proceed 
	if ($cds_start > $cds_stop) {
		my $tmp_var = $cds_start;
		$cds_start = $cds_stop;
		$cds_stop  = $tmp_var;
	}

	$cds_start = $cds_start - 500;
	$cds_stop  = $cds_stop + 500;

	print "Query locus for the species '$species' is '$chr_query:$cds_start-$cds_stop', the strand is '$strand_query'\n" if ($verbose);
	my $query_sequence = get_sequence_from_2bit($query_2bit_file,$chr_query,$cds_start,$cds_stop);
	$query_sequence = revComp($query_sequence) if ($strand_query ne $strand_ref);
	my $header = "$species#$cds_start#$cds_stop#$chr_query#$strand_query";
	return ($query_sequence,$header);
}

sub split_multi_exon_to_single_exon {  ## Takes in two sequences which are the aligned-ref and aligned-query (produced by CESAR-multi-exon mode)
	                                   ## and splits the alignment into exonwise mode, so the output resembles the output of CESAR-single-exon mode except a longish query sequence
	my($ref_sequence,$query_sequence) = @_;
	
	$ref_sequence   =~s/\s+$//;
	$query_sequence =~s/\s+$//;
	
	my @seq_split = split(/\s+/,$ref_sequence);
	@seq_split = grep{$_ ne ""}@seq_split;
	
	my $pos_start = 0;
	my @pos_start_list = ();
	foreach my $s(@seq_split) {
		my $pos = index($ref_sequence,$s,$pos_start);
		push(@pos_start_list,$pos);
		$pos_start = $pos + length($s);
	}
	
	my @exon_blocks = ();
	my $k = 1;
	foreach my $x(@pos_start_list) {
		my $seq_ref_exon   = substr($ref_sequence,$x,length($seq_split[$k-1]));
		my $seq_query_exon = substr($query_sequence,$x,length($seq_split[$k-1]));	

		## the following is to ensure that the sequence up and down of the aligning exon in the query genome
		## has no gaps and mimics raw genomic sequence and not some aligned sequence
		my $seq_query_up   = substr($query_sequence,0,$x);
		my $seq_query_down = substr($query_sequence,$x+length($seq_split[$k-1]));

		$seq_query_up   =~s/-//g; 
		$seq_query_down =~s/-//g;
 
		my $seq_ref_up    = create_string(length($seq_query_up)," ");
		my $seq_ref_down  = create_string(length($seq_query_down)," ");
		my $seq_ref_new   = $seq_ref_up.$seq_ref_exon.$seq_ref_down;
		my $seq_query_new = $seq_query_up.$seq_query_exon.$seq_query_down;	
		
		my $new_exon_alignment = "$seq_ref_new\n$seq_query_new";
		push(@exon_blocks,$new_exon_alignment);
		
		$k++;
	}
	return \@exon_blocks;
}

sub intron_deletion_correction {  ### correct for intron-deletions. These are indicated by in ">>>>>>>>>>>>>>>>>>>" reference sequence and
	my($ref_seq,$query_seq) = @_; ### "-------------------" in the query sequence. This function deletes these elements from the sequence.
	
	my @seq_split = split(/>+/,$ref_seq);
	my $exons_merged = scalar(@seq_split);
	
	my $pos_start = 0;
	my @pos_start_list = ();
	foreach my $s(@seq_split) {
		my $pos = index($ref_seq,$s,$pos_start);
		push(@pos_start_list,$pos);
		$pos_start = $pos + length($s);
	}
	
	my $k = 0;
	my $seq_ref_exon_new = my $seq_query_exon_new = "";
	foreach my $x(@pos_start_list) {
		my $seq_ref_exon   = substr($ref_seq,$x,length($seq_split[$k]));
		my $seq_query_exon = substr($query_seq,$x,length($seq_split[$k]));
		$seq_ref_exon_new   = $seq_ref_exon_new.$seq_ref_exon;
		$seq_query_exon_new = $seq_query_exon_new.$seq_query_exon;
		$k++;
	}
	return($seq_ref_exon_new,$seq_query_exon_new,$exons_merged);
}

sub create_string {
	my($length,$char) = @_;
	
	my $string = "";
	my $i = 1;
	while ($i <= $length) {
		$string = $string.$char;
		$i++;
	}
	return $string;
}

sub invert_values {
	my($ref_2_cds_values,$species,$noe) = @_;
	
	my @values_all = ();
	my $j = my $k = 0;
	while ($j < $noe) {
		push(@values_all,$ref_2_cds_values->{$species}{$j+1});
		$j++;
	}

	@values_all = reverse(@values_all);
	while ($k < $noe) {
		$ref_2_cds_values->{$species}{$k+1} = $values_all[$k];
		$k++;
	}
	return $ref_2_cds_values;
}
