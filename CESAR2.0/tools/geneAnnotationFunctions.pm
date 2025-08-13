package geneAnnotationFunctions;
use strict;
use warnings;
use Exporter;
use Scalar::Util::Numeric qw(isint);
use POSIX;

our @ISA = ('Exporter');
our @EXPORT = qw(checkExonIntactness getIntronTable getSLineData); 

my @stopCodonsList = qw(TAA TAG TGA);
my %stopCodons = map{$_ => 1}@stopCodonsList;

sub checkExonIntactness {
	my($refSequence,$querySequence,$header,$refStrand,$gene,$exonN,$phaseString,$accRef,$donorRef,$noe,$verbose) = @_;

	die "Undefined value for strand '$refStrand'\n" if ($refStrand ne "+" && $refStrand ne "-");

	my $querySequenceFull = $querySequence;
	my($alignmentStart,$alignmentStop) = getAlignmentBoundaries($refSequence,$verbose);
	my $refCDS   = substr($refSequence,$alignmentStart,($alignmentStop-$alignmentStart));
	my $queryCDS = substr($querySequence,$alignmentStart,($alignmentStop-$alignmentStart));

	my $queryCDSNoGaps = $queryCDS;
	$queryCDSNoGaps =~s/-//g;
	my $lengthAlignment = length($queryCDSNoGaps);

	my($phase5Prime,$phase3Prime) = (split /-/,$phaseString);

	my $nt5PEnd = my $nt3PEnd = "";
	if ($phase5Prime != 0) {
		if ($phase5Prime == 1) {
			$nt5PEnd = "NN";	
		} elsif ($phase5Prime == 2) {
			$nt5PEnd = "N";
		} else {
			die "Undefined value of phase5Prime '$phase5Prime\n'";
		}
	}

	if ($phase3Prime != 0) {
		if ($phase3Prime == 1) {
			$nt3PEnd = "NN";	
		} elsif ($phase3Prime == 2) {
			$nt3PEnd = "N";
		} else {
			die "Undefined value of phase3Prime '$phase3Prime\n'";
		}
	}

	$refSequence   = uc($nt5PEnd.$refCDS.$nt3PEnd);
	$querySequence = uc($nt5PEnd.$queryCDS.$nt3PEnd);
	
	if ($verbose) {
		print "############\n";
		print "The phase string is '$phaseString'\nThe exon number is '$exonN'\n";
		print "The reference sequence is '$refSequence'\nThe query sequence is     '$querySequence'\n";
		print "############\n";
			
	}

	my $codonPosHash_Ref    = getCodonPosHash($refSequence);
	my $codonPosHash_Query  = getCodonPosHash($querySequence);
	my @refSeqArray = split("",$refSequence);
	my @querySeqArray = split("",$querySequence);
	
	my ($fsCount,$insertWithStop)    = getFrameshifts($refSequence,$querySequence,$codonPosHash_Ref,$codonPosHash_Query,\@querySeqArray,$verbose);
	my $stopCodonCount               = getInframeStopCodons($refSequence,$querySequence,$codonPosHash_Ref,$verbose);
	my $spliceSiteMutations          = getSpliceSiteMutations($querySequenceFull,$alignmentStart,$alignmentStop,$accRef,$donorRef,$exonN,$noe,$verbose);
					
	my $exonCDS = "";
	if ($fsCount == 0 && $insertWithStop == 0 && $stopCodonCount == 0 && $spliceSiteMutations == 0) {
		$exonCDS = getCoordinatesOfIntactExon($header,$alignmentStart,$alignmentStop,$gene,$exonN,$noe,$lengthAlignment,$querySequenceFull,$refStrand);
	} else {
		$exonCDS = "NI";
	}

	return $exonCDS;
}

## Different sub-routines
sub getFrameshifts {
	my($refSequence,$querySequence,$codonPosHash_Ref,$codonPosHash_Query,$ref2QueryArray,$verbose) = @_;
	my $insertions = getIndels($refSequence);
	my $deletions  = getIndels($querySequence);
	
	my %insertionsHash = %$insertions;
	my %deletionsHash  = %$deletions;
	
	my $fsCount = my $insertWithStop = 0;
	my %allFrameshiftsHash = ();  ## create a hash and only look at those indels that cause a frameshift
	
	foreach my $insert(keys(%insertionsHash)) {
		my $insertLength = $insertionsHash{$insert};
		if ($insertLength%3 != 0) {
			$fsCount++; 	
			$allFrameshiftsHash{$insert} = $insertLength;
			print "Insertion --> Happens at position '$insert', the value is '$insertLength'\n" if ($verbose);
		} else {  ## check if a frame-preserving indel has an in-frame stop codon
			my $seqInsert = substr($querySequence,$insert,$insertLength);
			
			my $i = 0;
			while ($i < length($seqInsert)) {
				my $codon = substr($seqInsert,$i,3);
				if (exists $stopCodons{$codon}) {
					print "The framepreserving insertion contains the stop codon '$codon' is present in the sequence '$seqInsert' at position '$i'\n" if ($verbose);
					$insertWithStop++;
				}
				$i = $i + 3;
			}
		}
	}
	
	foreach my $deletion(keys(%deletionsHash)) {
		my $deletionLength = $deletionsHash{$deletion};
		if($deletionLength%3 != 0) {
			$fsCount++;
			$deletionLength = 0 - $deletionLength;
			$allFrameshiftsHash{$deletion} = $deletionLength;
			print "Deletion --> Happens at position '$deletion', the value is '$deletionLength'\n" if ($verbose);
		}
	}
	
	if ($fsCount >= 2) {  ## find if the frameshifts can compensate each other
		print "#####\n\nThe count of frameshifts is '$fsCount', so trying to determine if there are frameshifts that compensated for each other ....\n\n\n" if ($verbose);
		my @allFrameshifts = sort{$a <=> $b}(keys(%allFrameshiftsHash));
		
		my %excludedInsertions = ();
		my %excludedDeletions = ();
		
		my $ref2CompensatedEvents = getCompensatedEvents(\%allFrameshiftsHash,\%excludedInsertions,\%excludedDeletions,$codonPosHash_Ref,$codonPosHash_Query,$ref2QueryArray,$verbose);
		my %compensatedEvents = %$ref2CompensatedEvents;
		my %ignoreCompensatedFS = ();
		
		foreach my $keys(keys(%compensatedEvents)) {
			my @compFSEvents = split(",",$compensatedEvents{$keys});
			push(@compFSEvents,$keys);
			my %compEventsHashLocal = map{$_ => 1}@compFSEvents;
			%ignoreCompensatedFS = (%ignoreCompensatedFS,%compEventsHashLocal);
		}
		
		my @allFrameshifts_NonCompensated = ();
		foreach my $fs(@allFrameshifts) {
			push(@allFrameshifts_NonCompensated,$fs) if (! exists $ignoreCompensatedFS{$fs});
		}
		$fsCount = scalar(@allFrameshifts_NonCompensated);
	}
	
	return ($fsCount,$insertWithStop);
}

## Get insertions and deletions in the sequence
sub getIndels { ## Insertions are gaps in the reference sequence, deletions are gaps in the query sequence
	my $sequence = shift;
	my @seqArray = split(/[A-Z]/,$sequence);
	@seqArray = grep{$_ ne ""}@seqArray;

	my %hashIndel = ();
	my $posStart = 0;
	## Get positions of gaps
	foreach my $element(@seqArray) {	
		my $pos = index($sequence,$element,$posStart);
		my $lengthIndel = length($element);
		$hashIndel{$pos} = $lengthIndel;
		$posStart = $pos + length($element);
	}
	
	return \%hashIndel;
}

sub getCodonPosHash {
	my $sequence = shift;
	my %codonPosHash = ();
	
	my @seqArray = split("",$sequence);
	my $ct = 0;
	my $baseUngapped = 0;
	
	foreach my $base(@seqArray) {
		if ($base ne "-") {
			$baseUngapped++;
			my $codonPos = $baseUngapped%3;
			$codonPos = 3 if ($codonPos == 0);
			$codonPosHash{$ct} = $codonPos;
		} else {
			$codonPosHash{$ct} = "-";
		}
		$ct++;
	}
	return \%codonPosHash;
}

sub getCompensatedEvents {
	my ($ref2allFrameshifts,$ref2ExcludedInsertions,$ref2ExcludedDeletions,$codonPosHash_Ref,$codonPosHash_Query,$ref2QueryArray,$verbose) = @_;		### 
	
	my %fdEvents = %$ref2allFrameshifts;
	my @indelsAll = sort {$a <=> $b} (keys(%fdEvents));  ## This array contains the positions of all frameshifting events
	
	my %excludedInsertions	  = %$ref2ExcludedInsertions;
	my %excludedDeletions	  = %$ref2ExcludedDeletions;

	my %compensationHash = ();		### This is the hash that contains all the events which can be compensated, the key is the first event, the value is the event (or events) which can compensate the first event.
	my %excludeEvents = ();			### Any event which can compensate for another event is dumped into this hash, so that it is not considered again.

	for (my $i = 0;$i < @indelsAll;$i++) {
		my $queryEvent = $indelsAll[$i];	## this is the query, now find events which can compensate for this query event.
		next if (exists $excludeEvents{$queryEvent});
		
		print "Querying '$queryEvent' to find possible compensating events\n" if ($verbose);
		
		my $length = $fdEvents{$queryEvent};	## get the length of the query event	
		my $sumLengthsPCE = $length;
		my @potCompEventList = ();
		
		for (my $j = $i + 1;$j < @indelsAll;$j++) {	## Now start looking for compensatory events. Remember if the index of the query 
													## event is i, the compensatory events should be search from i+1
			my $potCompEvent = $indelsAll[$j];
			my $lengthPCE = $fdEvents{$potCompEvent};
			next if (exists $excludeEvents{$potCompEvent}); 		## the potential compensatory event should not have been "consumed" in compensating for some other event.
				
			$sumLengthsPCE = $sumLengthsPCE + $fdEvents{$potCompEvent};
			push(@potCompEventList,$potCompEvent);
			
			if ($sumLengthsPCE%3 == 0) {	## if the length of query is +2, length of potential compensating event is -5, 
											## $sumLengthsPCE = -3, so the events are compensatory.
				print "This combination looks like a compensatory event '@potCompEventList'\n" if ($verbose);
				my $start = $queryEvent;
				my $stop  = $potCompEvent + abs($lengthPCE) - 1;		## Go until the end of the event. For instance if the compensatory event is an Insertion, 
																		## go until the end of the insertion and check if there is a stop codon in the species reading frame.
				my $startUpstream = $queryEvent - 1;
				my $ancRF = getAncestralReadingFrame($startUpstream,$codonPosHash_Ref);	## What the function returns is the ancestral reading frame
						
				my $speciesRF = "";		
				## Now I convert the ancestral reading frame to the species reading frame
				if ($ancRF eq "-") {
					$speciesRF = 1; ### Ancestral reading frame is zero when we are trying to find the ancestral reading frame for a deletion that occurs right at the beginning of the gene. In this case, startUpstream = -1	
				} else {
					$speciesRF = $ancRF + 1;
					$speciesRF = 1 if ($speciesRF == 4);	
				}
				print "The ancestral reading frame is '$ancRF' and the species reading frame is '$speciesRF'\n" if ($verbose);
					
				my $stopCodonFlag = "";
				print "Now checking for stop codons between '$start' and '$stop'\n" if ($verbose);
				$stopCodonFlag = checkForStopCodons($start,$stop,$ref2ExcludedInsertions,$ref2ExcludedDeletions,$speciesRF,$ref2QueryArray,$verbose);
				
				if ($stopCodonFlag eq "T") { ## if there is no stop codon
					print "The query event is '$queryEvent'. The following events are are compensatory for the query event '@potCompEventList'\n" if ($verbose);
					my %excludeEventsNew = map{$_ => "T"}@potCompEventList; ## then put all the compensatory events in the %excludeEvents hash.
					%excludeEvents = (%excludeEvents,%excludeEventsNew);
					
					my $compEventsString = join(",",@potCompEventList);
					$compensationHash{$queryEvent} = $compEventsString;	## Also populate the CompensationHash accordingly.
				}
				last;		## and leave this loop, since we have already found compensatory event(s) for our query events.
			}
					## otherwise continue to look for compensatory events
		}			## End of the FOR loop that looks for compensatory events
					## IF statement which checks that the queryEvent has not already been "consumed"/ it is a compensatory event for some other upstream event.
	}				## End of the FOR loop that iterates through the list of all Indels.
	
	### Perform the "Return To Ancestral Reading Frame Test"
	foreach my $events( sort {$a <=> $b} keys(%compensationHash)) {	
		## Get eventStop position, or the position of the last of the compensatory events. Remember, compensatory events could be a list of events
		my @tmpCompEvents = split(",",$compensationHash{$events});
		my $stop = $tmpCompEvents[$#tmpCompEvents];
		
		my $startUpstream = $events - 1;
		my $speciesRF = getAncestralReadingFrame($startUpstream,$codonPosHash_Ref);	## What the function returns is the ancestral reading frame. since the species reading frame
																	## has been ancestralized, speciesRF = ancestralRF;
		my $length = abs($fdEvents{$stop});
		
		## Only deletions at the beginning of the gene sequence are excluded from the "Return To Ancestral Reading Frame Test"
		next if ($speciesRF eq "-");
		
		my $start = $events;		## start position is where the first event happens
		my $stopPos = $stop + $length -1;	## stop position is where the last of the compensatory events happen + the length of the event - 1
		my $readingFrameRF = "";
		
		while ($start <= $stopPos) {
			if (! exists $excludedInsertions{$start}) {
				### Reading frame for the reference species
				$readingFrameRF = $codonPosHash_Ref->{$start} if ($codonPosHash_Ref->{$start} ne "-");   ## ## Keep getting the reading frame of the Reference species for every line where there is sequence in the Reference
				### Reading frame for the query species
				$speciesRF++ if ( ($codonPosHash_Query->{$start} ne "-" ) || (exists $excludedDeletions{$start}) ); ## Same as for Reference. Even if there is a ?, we assume that there is sequence here. That is why, the if condition is 
																													## << if ($tmp[13] ne "-") >>  NOT  <<< if ( ($tmp[13] ne "-") AND if ($tmp[13] ne "?") ) >>>												
			}
			
			$start++;
		}
		
		my $speciesRFPrint = $speciesRF%3;
		$speciesRFPrint = 3 if ( $speciesRFPrint == 0);
		
		die "This compensatory event pair ...$events .. $compensationHash{$events} does not return the species reading frame to the ancestral reading frame\n" if ($readingFrameRF ne $speciesRFPrint);
	}
	
	return \%compensationHash;
}


sub getAncestralReadingFrame {
	my ($startUpstream,$codonPosHash_Ref) = @_;
	my $speciesRF = "";
	
	while ($startUpstream >= 0) {		## This is the index in the GeneArray from where I ancestralize the reading frame of the species i.e. a positon just upstream of the indel
		my $codonPos = $codonPosHash_Ref->{$startUpstream};
		
		if  ($codonPos ne "-") {		## Get the line where there is a character in the Reference species. Remember we need 
			$speciesRF = $codonPos;     ## the position where there is a character in the Reference, not in the query species
			last;                       # After all, the function returns the ancestral reading frame
		}
		$startUpstream--;
	}
	
	$speciesRF = 0 if ($speciesRF eq "");
	return $speciesRF;
}

sub checkForStopCodons {									
	my ($start,$stop,$ref2ExcludedInsertions,$ref2ExcludedDeletions,$readingFrame,$ref2QueryArray,$verbose) = @_;
	
	my %excludedInsertions = %$ref2ExcludedInsertions;
	my %excludedDeletions  = %$ref2ExcludedDeletions;
	my @querySeqArray      = @$ref2QueryArray;
	
	my $flagStop = my $geneSeq = "";
	print "Extracting gene sequence between '$start' and '$stop'\n" if ($verbose);
	
	for(my $k = $start;$k < $stop; $k++) {
		if (! exists $excludedInsertions{$k}) {
			$geneSeq = $geneSeq.$querySeqArray[$k];
			$geneSeq = $geneSeq."N" if (exists $excludedDeletions{$k});
		}
	}
	
	$geneSeq =~s/\?/N/g;
	$geneSeq =~s/-//g;
	print "The sequence is '$geneSeq'\n" if ($verbose);
	
	if  ($geneSeq eq "") {  ### This happens when two exon deletions can compensate for each other. For example, an exon deletion of 22nt (equivalent to -1) followed by an exon deletion of 23nt nt(equivalent to -2).
		$flagStop = "T";	### In such cases there will be no readingFrame and no geneSeq. And we have to believe that the events are indeed compensatory
	} else {				
		
		if ($readingFrame == 2) {
			$geneSeq = "N".$geneSeq;
		} elsif ($readingFrame == 3) {
			$geneSeq = "NN".$geneSeq;
		}
	
		my $stopCodon = find_stop_codon($geneSeq,$verbose);
		$flagStop = "T" if ($stopCodon eq "");
	}
	
	return $flagStop;
}

sub find_stop_codon {
	my ($seq,$verbose) = @_;;
	
	my $rem = (length($seq)%3);
	if ($rem == 1) { 
		$seq = $seq."NN";
	} elsif ($rem == 2) {
		$seq = $seq."N";
	}
	
	die "The length of the input sequence '$seq' is not divisible by 3\n" if (length($seq)%3 != 0);
	my $i = 0;
	my $stopCodonPresent = "";
	while ($i < length($seq)) {
		my $codon = substr($seq,$i,3);
		if (exists $stopCodons{$codon}) {
			print "FRAMESHIFT THAT CANNOT BE COMPENSATED!!! The stop codon '$codon' is present in the sequence '$seq' at position '$i'\n" if ($verbose);
			$stopCodonPresent = "T";
			last;
		}
		$i = $i + 3;
	}
	return $stopCodonPresent;
}

sub getInframeStopCodons {
	my($refSequence,$querySequence,$codonPosHash_Ref,$verbose) = @_;
	my $stopCodonCount = 0;
	my $i = 0;
	print "####\nInformation for stop codons:\n\n" if ($verbose);
	while ($i < length($querySequence)) {
		my $codon = substr($querySequence,$i,3);
		my $codonRef = substr($refSequence,$i,3);
		
		if (exists $stopCodons{$codon} && ! exists $stopCodons{$codonRef}) { ## Stop codon in query but not in reference. Otherwise it is not a mutation but most likely a stop codon at the end of the gene
			
			my $string = $codonPosHash_Ref->{$i}.$codonPosHash_Ref->{$i+1}.$codonPosHash_Ref->{$i+2};
			if($string eq "123") {
				print "STOP CODON found at position '$i'. The codon in the reference sequence is '$codonRef', the stop codon in query is '$codon'\n" if ($verbose);
				$stopCodonCount++;
			}
		}
		$i = $i + 3;
	}
	return $stopCodonCount;
}

sub getAlignmentBoundaries {
	my ($refSequence,$verbose) = @_;
	$refSequence =~s/\s+$//;
	my @seqArray = split("",$refSequence);
	
	my $stop  = length($refSequence);
	my $start = 0;
	foreach my $base(@seqArray) {
		if ($base ne " ") {
			last;
		}
		$start++;
	}
	print "The alignment start is '$start' and the stop is '$stop'\n" if ($verbose);
	return ($start,$stop);
}

sub getSpliceSiteMutations {
	my($alignedQuery,$alignmentStart,$alignmentStop,$accRef,$donorRef,$exonN,$noe,$verbose) = @_;
	my $accQuery    = substr($alignedQuery,$alignmentStart-2,2);
	my $donorQuery  = substr($alignedQuery,$alignmentStop,2);
	
	my $spliceSiteMutationCount = 0;
	print "Acceptor site query is '$accQuery', acceptor site reference is '$accRef'\n" if ($verbose);
	print "Donor site query is '$donorQuery', donor site reference is '$donorRef'\n" if ($verbose);

	$spliceSiteMutationCount++ if ($accQuery ne "ag" && $accQuery ne $accRef && $exonN != 1);
	$spliceSiteMutationCount++ if ($donorQuery ne "gt" && $donorQuery ne "gc" && $donorQuery ne $donorRef && $exonN != $noe);
	print "The number of splice site mutations is '$spliceSiteMutationCount'\n" if ($verbose);
	return $spliceSiteMutationCount;
}

sub getCoordinatesOfIntactExon {	
	my($header,$alignmentStart,$alignmentStop,$gene,$exonN,$noe,$lengthAlignment,$querySequenceFull,$strandRef) = @_;
	
	## This is how the header looks like panTro4#20600346#20600543#chr15#-
	my($species,$start,$stop,$chr,$strandQuery) = (split /#/,$header);
	my $alignmentStartRev = length($querySequenceFull) - $alignmentStop;

	$alignmentStart = $alignmentStartRev if ($strandRef ne $strandQuery);
	
	my $exonStart = $start + $alignmentStart;
	my $exonEnd   = $exonStart + $lengthAlignment;
	
	my $strandQueryPrint = "+";
	$strandQueryPrint = "-" if ($strandRef ne $strandQuery);
	
	my $bedLine = "$chr\t$exonStart\t$exonEnd\t$strandQueryPrint\t$exonN";
	return $bedLine;
}

sub getIntronTable {
	my($gene,$refChr,$refStrand,$ref2ExonsList,$ref2SpeciesList,$ref2HashInfo,$verbose) = @_;
	
	my @speciesListArray = @$ref2SpeciesList;
	my @exonsList        = @$ref2ExonsList;
	my %hashInfo         = %$ref2HashInfo;
	my @intronLengthsList = ();

	@exonsList = reverse(@exonsList) if ($refStrand eq "-");
	
	foreach my $species(@speciesListArray) {
		$species =~s/\s+$//;
		
		## Make pairs for exons that should be compared: This is important because we could have an alignment to exon1, e-lines for exon2 and an alignment to exon3. 
		## In this case, we have to check that exon1 and exon3 should not overlap
		my @pairsList = ();
		my $j = "";
		
		for (my $i = 0; $i < scalar(@exonsList)-1; $i++) {
			my $e1 = $i;
			my $e2 = $i + 1;
			
			my $exon1 = $exonsList[$e1];
			my $exon2 = $exonsList[$e2];
			print "Checking the possibility of exonPair --> '$e1' and '$e2' i.e. '$exon1 and '$exon2'\n" if ($verbose);
		
			my($cds1,$chr1,$strand1,$l1) = (split /#/,$hashInfo{$species}{$i});
			if ($chr1 eq "NA" || $strand1 eq "NA") {
				print "Chr value is '$chr1' and strand value is '$strand1', so excluding\n\n" if ($verbose);
				next;
			}
		
			my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$i+1});
			my $pair = "";
		
			if ($chr1 ne "NA" && $strand1 ne "NA" && $chr1 eq $chr2 && $strand1 eq $strand2) {
				$pair = "$e1-$e2"; 
				print "Possible pair $e1 and $e2 will be evaluated. Chr/strands are '$chr1', '$strand1', '$chr2' and '$strand2'\n\n" if ($verbose);
			} else {
				print "Chr value1 is '$chr1' and strand value is '$strand1'; Chr value2 is '$chr2' and strand value is '$strand2' so excluding and checking for downstream pairs\n\n" if ($verbose);
				my $j = $i + 2;
				my $jStart = $j;
			
				print "Scenario2: Evaluating the possibility of '$j' and downstream exons\n\n" if ($verbose);
				if ($j < scalar(@exonsList)) {
					for($j = $jStart; $j < scalar(@exonsList); $j++) {
						my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$j});
						if ($chr1 ne "NA" && $strand1 ne "NA" && $chr1 eq $chr2 && $strand1 eq $strand2) {
							$pair = "$e1-$j";
							print "This is the new pair '$pair'\n\n" if ($verbose);
							$i = $j;
							last;
						}
					}
				}
			}
			print "'$pair' is a possible exon pair\n" if ($pair ne "" && ($verbose) );
			push(@pairsList,$pair) if ($pair ne "");
		}

		if (scalar(@pairsList) > 0) {
			my %pairsListHash = map{$_ => 1} @pairsList;
			## If the exons are deleted, then the above code would give me a list like: 0-5 6-12 18-35 36-37 37-38 38-39
			## Ensure that the list becomes continuous i.e, 0-5,5-6,6-12,12-18,18-35 and so on

			my @pairsListNew = ();
			for(my $k = 0; $k < (scalar(@pairsList)-1); $k++) {
				push(@pairsListNew,$pairsList[$k]);  ## i.e. the original pair
				my ($startOriginal,$stopOriginal) = (split /-/,$pairsList[$k])[0,1];
				
				my($chrOriginal,$strandOriginal) = (split /#/,$hashInfo{$species}{$startOriginal})[1,2];
				
				my $start = (split /-/,$pairsList[$k])[1];
				my $stop  = (split /-/,$pairsList[$k+1])[0];				
				my $pairNew = "$start-$stop";

				## Now check if the "new" exon-pairs that I am adding to the list to ensure continuity, come from the same strand and same chromosome
				my($cds1,$chr1,$strand1,$l1) = (split /#/,$hashInfo{$species}{$start});
				my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$stop});

				push(@pairsListNew,$pairNew) if ( ($chr1 eq $chr2 && $strand1 eq $strand2) && (! exists $pairsListHash{$pairNew} && $start != $stop) );
			}
		
			push(@pairsListNew,$pairsList[$#pairsList]);  ## Do this for the last element
			###
			my ($pairStart,$pairEnd)    = (split /-/,$pairsList[$#pairsList])[0,1];
			my ($chrLast,$strandLast)   = (split /#/,$hashInfo{$species}{$pairStart})[1,2];
			 
			########
			@pairsList = @pairsListNew;
			print "These are the continuous exonPairs '@pairsList'\n" if ($verbose);
			my @list = ();

			## Now compare adjacent exons from the pairs that I have created above
			for (my $i = 0; $i < scalar(@pairsList); $i++) {
				my($e1,$e2) = (split /-/,$pairsList[$i]);
				
				my $exon1 = $exonsList[$e1];
				my $exon2 = $exonsList[$e2];
				###
				my($cds1,$chr1,$strand1,$l1) = (split /#/,$hashInfo{$species}{$e1});
				my($cds2,$chr2,$strand2,$l2) = (split /#/,$hashInfo{$species}{$e2});
		
				my $intronLength = $cds2 - $cds1 - $l1;
				my $keyNew = "$exon1<->$exon2<->$intronLength";
				push(@list,$keyNew);
			}
			
			@list = reverse(@list) if ($refStrand eq "-");

			### Formatting the output
			foreach my $element(@list) {
				my $e1 = my $e2 = my $length = "";
				
				if ($refStrand eq "+") {
					($e1,$e2,$length) = (split /<->/,$element)[0,1,2];
				} else {
					($e2,$e1,$length) = (split /<->/,$element)[0,1,2];
				}
				print "Strand is $refStrand .. exon1 is '$e1', exon2 is '$e2' while the length is '$length'\n\n" if ($verbose);
			
				if ($length > 29) {
					$length = $length - 30;
					my $lengthIntron = $length/2;
					my $l1 = my $l2 = "";
					
					if ($length%2 == 0) {
						$l1 = $lengthIntron + 23;
						$l2 = $lengthIntron + 7;
					} else {
						$l1 = ceil($lengthIntron) + 23;
						$l2 = floor($lengthIntron) + 7;
					}
					
					push(@intronLengthsList,"$gene\tthreePrime\t$refChr\t$e1\t$species\t$l1");
					push(@intronLengthsList,"$gene\tfivePrime\t$refChr\t$e2\t$species\t$l2");
					
				} elsif ($length < 29 && $length >= 0){
					my $lengthIntron = $length/2;
					my $l1 = $lengthIntron;
					my $l2 = $lengthIntron;
				
					if ($length%2 != 0){
						$l1 = ceil($lengthIntron);
						$l2 = floor($lengthIntron);
					}
			
					print "Intron deleted!! The lengths are '$l1' and '$l2'\n\n" if ($verbose);
					push(@intronLengthsList,"$gene\tthreePrime\t$refChr\t$e1\t$species\t$l1");
					push(@intronLengthsList,"$gene\tfivePrime\t$refChr\t$e2\t$species\t$l2");
					
				} elsif ($length < 0) {
						push(@intronLengthsList,"$gene\tthreePrime\t$refChr\t$e1\t$species\t$length");
						push(@intronLengthsList,"$gene\tfivePrime\t$refChr\t$e2\t$species\t$length");
						
						print "Negative intron length !! The length is '$length'\n\n" if ($verbose);
				}
			}
		}
	}
	return \@intronLengthsList;
}

sub getSLineData {
	my $line = shift; 
	
	if ($line =~ /^s (.*)/) {		# sequence part
		# ornAna1.Contig17774                  7736 10 +     18306 CTGGG----GCTGT
		my @f = split(/[ ]+/, $1); 
		if ($#f != 5) {
			die "ERROR in getSLineData: cannot parse 6 elements from $line\n";
		} 
		my ($src, $start, $size, $strand, $srcSize, $seq) = (@f)[0,1,2,3,4,5];
		# few sanity checks
		if (! isint($start) || ! isint($size) || ! isint($srcSize)) {
			die "ERROR in getSLineData: start/size/srcSize are not integers in $line\n";
		}
		if ($strand ne "+" && $strand ne "-") {
			die "ERROR in getSLineData: strand is neither + nor - in $line\n";
		}
		my ($species) = (split(/[\.]/, $src))[0];
		my $chr = substr($src, length($species)+1, length($src));	# given "tupBel1.scaffold_142368.1-170510" $chr is "scaffold_142368.1-170510"

#		print "s line: ($blockNo, $species, $chr, $start, $size, $strand, $srcSize, $seq)  [[$line]]\n" if ($verbose);

		return ($species, $chr, $start, $size, $strand, $srcSize, $seq);
	} else {
		die "ERROR: call getSLineData with no s line: $line\n";
	}
}
