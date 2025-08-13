#!/usr/bin/env python3

"""
A speed-up version of processed_segment.py
Most likely a final solution make
"""

from cesar_wrapper_constants import * ## TODO: Explicit imports!
from cesar_wrapper_executables import * ## TODO: Same!
from collections import defaultdict
from logging import Logger
from math import floor
from shared import nn, parts, safe_div
from typing import (
    Any, Dict, Iterable, List, Optional, Set, Tuple, TypeVar, Union
)
from ucsc_report import (
    exon_aln_entry, exon_aln_header, format_fasta_as_aln, mutation_table, ProjectionPlotter
)

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = ('Michael Hiller', 'Bogdan Kirilenko')

SelenoCysteineTuple = TypeVar('SelenoCysteineTuple', bound='Tuple[int, int, str]')

class Coords:
    """Stores start and stop coordinates"""
    __slots__ = ('start', 'stop')
    def __init__(self, start: int, stop: int) -> None:
        self.start = start
        self.stop = stop

    def tuple(self) -> Tuple[int, int]:
        """Returns a tuple of (start, stop) coordinates"""
        return tuple(sorted((self.start, self.stop)))

    def __repr__(self) -> str:
        return str(self.tuple())


def has_spliceai_data(splice_sites: Dict[int, Dict[str, List[int]]]) -> bool:
    """Checks if any SpliceAI data predictions were provided"""
    return any(
        splice_sites[x][y] for x in splice_sites for y in splice_sites[x]
    )


def is_symbol(base: str) -> bool:
    """Estimates whether the symbol is a valid coding symbol, i.e. letter or gap"""
    return base.isalpha() or base == '-'


def is_complete_codon(codon: str) -> bool:
    """Estimates whether the reference codon is a complete triplet"""
    return sum(x.isalpha() for x in codon) == 3 or codon == GAP_CODON


def strip_noncoding(string: str, uppercase_only: bool = True) -> str:
    """Removes noncoding symbols from the nucleotide string"""
    return ''.join(
        x for x in string if (x.isupper() if uppercase_only else x.isalpha())
    )


def frameshift_codon(codon: str) -> bool:
    """Returns whether a codon contains a frameshift mutation"""
    return codon.count('-') % 3 != 0 and codon != GAP_CODON# or len(codon) % 3 != 0


class ProcessedSegment:
    """
    The main CESAR2.0 results processing class. One instance corresponds to a
    single segment and configuration. Internally, the class does the following:

    """
    __slots__ = (
        'transcript', 'segment', 'configuration', 'matrix', 'u12',
        'mask_terminal_mutations', 'rescue_missing_start', 'rescue_missing_stop',
        'min_gap_size', 'correction_mode', 'min_splice_prob', 'prob_margin',
        'acceptor_flank', 'donor_flank', 
        'correct_ultrashort_introns', 'ignore_alternative_frame',
        'is_paralog', 'is_processed_pseudogene', 'logger',
        'reference', 'query', 'reference_length', 'portion2coords',
        'portion2abs_coords', 'portion2strand', 'exon2portion', 'exon2chrom',
        'exon2chain', 'exon2asmgaps', 'unaligned_exons', 'gap_located_exons',
        'out_of_chain_exons', 'expected_coordinates', 'exon_search_spaces',
        'chains', 'last_aligned_base', 'ref_exon_lengths', 'aln_lengths',
        'reference_length', 'exon2ref_codons', 'exon2codons', 'rel_exon_coords',
        'abs_exon_coords', 'cesar_rel_coords', 'clipped_rel_coords', 'spliceai_sites',
        'splice_site_nucs', 'intron2phase',
        'exon_num', 'has_spliceai_data', 'donor_site_prob',
        'acc_site_prob', 'cesar_donor_support', 'cesar_acc_support',
        'spliceai_donor_support', 'spliceai_acc_support', 'found_in_exp_locus', 'corrected_frame',
        'exon_evidence', 'all_ref_codons', 'all_query_codons', 'ref_aa_seq',
        'query_aa_seq', 'triplet_coordinates', 'codon_coordinates',
        'triplet2ref_codon', 'ref_codon2triplets', 'ref_codons_to_mask',
        'query_codons_to_mask', 'query_atgs', 'split_codon_struct',
        'intron_deleted', 'intersects_gap', 'mutation_list', 'selenocysteine_codons', 
        'exon_presence', 'exon_quality', 'codon2mutations', 'ref_codon2mutations',
        'exon2mutations', 'deletion2mask', 'introns_gained',
        'abs_nuc_ids', 'blosum_self_values', 'blosum_diff_values', 'nuc_ids',
        'blosum_ids', 'nuc_id', 'blosum_id',
        'frameshift_counter', 'spliceai_compensations', 
        'alt_frame_masked', 'comps_to_delete',
        'nonsense_counter', 'indel_counter',
        'acc_mutation_counter', 'donor_mutation_counter', 'compensation_counter',
        'deleted_exon_counter', 'missing_exon_counter', 'mdel_counter',
        'frameshift2compensation', 'alternative_frames', 'stop_updated',
        'loss_status', 'longest_fraction_strict', 'longest_fraction_relaxed',
        'total_intact_fraction', 'missing_portion_of_ref',
        'middle_is_intact', 'middle_is_present',
        'name', 'svg', 'v'
    )

    def __init__(
        self,
        transcript: str,
        segment: int,
        configuration: int,
        cesar_output: List[RawCesarOutput],
        matrix: Dict[str, Dict[str, int]],
        # reference_introns: Dict[int, Exon],
        u12: Dict[int, Set[int]],
        mask_terminal: bool,
        rescue_missing_start: bool,
        rescue_missing_stop: bool,
        min_gap_size: int,
        correction_mode: int,
        min_splice_prob: float,
        prob_margin: float,
        acceptor_flank: int,
        donor_flank: int,
        correct_ultrashort_introns: bool,
        ignore_alternative_frame: bool,
        is_paralog: bool,
        is_processed_pseudogene: bool,
        logger: Logger,
        verbose: bool
    ) -> None:
        self.transcript: str = transcript
        self.segment: int = segment
        self.configuration: int = configuration
        self.matrix: Dict[str, Dict[str, int]] = matrix
        # self.ref_intron_lengths: Dict[int, Exon] = reference_introns
        self.u12: Dict[int, Set[int]] = u12
        self.mask_terminal_mutations: bool = mask_terminal
        self.rescue_missing_start: bool = rescue_missing_start
        self.rescue_missing_stop: bool = rescue_missing_stop
        self.min_gap_size: int = min_gap_size
        self.correction_mode: int = correction_mode
        self.min_splice_prob: float = min_splice_prob
        self.prob_margin: float = prob_margin
        self.acceptor_flank: int = acceptor_flank
        self.donor_flank: int = donor_flank
        self.correct_ultrashort_introns: bool = correct_ultrashort_introns
        self.ignore_alternative_frame: bool = ignore_alternative_frame
        self.is_paralog: bool = is_paralog
        self.is_processed_pseudogene: bool = is_processed_pseudogene
        self.logger: Logger = logger

        self.reference: str = ''
        self.query: str = ''
        self.reference_length: int = 0 ## stores the ungapped reference sequence length
        self.portion2coords: Dict[int, Tuple[int]] = {} ## stores CESAR alignment portion coordinates in the concatenated sequences
        self.portion2abs_coords: Dict[int, Tuple[str, int]] = {} ## stores query coordinate for each CESAR output entity
        self.portion2strand: Dict[int, bool] = {} ## stores strand data for each CESAR output entity
        self.exon2portion: Dict[int, int] = {} ## stores which alignment portion the given exon was aligned in
        self.exon2chrom: Dict[int, str] = {} ## stores the chromosomes the respective exon's alignments refer to
        self.exon2chain: Dict[int, str] = {} ## stores the chains that were fo respective exons' alignment
        self.exon2asmgaps: Dict[int, bool] = {} ## stores whether the exon's search space contained assembly gaps
        self.unaligned_exons: Set[int] = set() ## stores which exons were replace with was_not_aligneds
        self.gap_located_exons: Set[int] = set() ## stores exons which were correspond to a chain gap
        self.out_of_chain_exons: Set[int] = set() ## stores exons which were not aligned
        self.expected_coordinates: Dict[int, Tuple[int]] = {} ## stores expected exon coordinates prior to flank addition
        self.exon_search_spaces: Dict[int, Tuple[int]] = {} ## stores exon coordinates after flank addition
        self.chains: List[str] = [] ## stores all the chains used for all the CESAR alignments
        self.last_aligned_base: int = 0 ## stores the relative coordinate of the last base aligned in CESAR output
        self.ref_exon_lengths: Dict[int, int] = {} ## stores reference exon lengths
        self.aln_lengths: Dict[int, int] = {} ## stores exonwise alignment lengths
        self.reference_length: int = 0 ## overall reference transcript length
        self.exon2ref_codons: Dict[int, Tuple[int]] = {} ## maps query exons to first and last codon of their reference counterpart
        self.rel_exon_coords: Dict[int, Coords] = {} ## maps exon to relative coordinates in the alignment
        self.clipped_rel_coords: Dict[int, Coords] = {} ## maps exons to relative coordinates after clipping termini aligned to assembly gaps
        self.abs_exon_coords: Dict[int, Coords] = {} ## maps reference exons to coordinates in the query
        self.cesar_rel_coords: Dict[int, Coords] = {} ## maps exons to CESAR alignment coordinates
        self.spliceai_sites: Dict[int, Dict[str, Dict[int, float]]] = defaultdict(
            lambda: defaultdict(dict)
        ) ## stores SpliceAI predictions per each exon
        self.splice_site_nucs: Dict[int, str] = defaultdict(dict) ## stores splice site dinucleotides in the (donor, acceptor) format
        self.introns_gained: Dict[int, Tuple[int, int]] = {} ## stores coordinates of gained introns

        ## TODO: Move the code below into a separate method
        phase: int = 0
        for i, aln_portion in enumerate(sorted(cesar_output, key=lambda x: x.exons)):
            # print(f'>portion {i}, ref\n{aln_portion.reference}\n>portion {i}, query\n{aln_portion.query}\n{len(aln_portion.reference)=}, {len(aln_portion.query)=}, {self.exon2portion=}')
            prev_aln_len: int = len(self.reference)
            self.reference += aln_portion.reference
            self.query += aln_portion.query
            curr_aln_len: int = len(self.reference)
            self.chains.append(aln_portion.chain) if aln_portion.chain not in self.chains else None
            self.portion2coords[i] = (prev_aln_len, curr_aln_len)
            self.portion2abs_coords[i] = [
                aln_portion.chrom, aln_portion.start, aln_portion.stop
            ]
            self.portion2strand[i] = aln_portion.strand
            self.expected_coordinates = {
                **self.expected_coordinates, **aln_portion.exon_expected_loci
            }
            self.exon_search_spaces = {
                **self.exon_search_spaces, **aln_portion.exon_search_spaces
            }
            # self.reference_length += len(
            #     [x for x in aln_portion.reference if x.isalpha()]
            # )
            self.gap_located_exons = self.gap_located_exons.union(
                aln_portion.gap_located_exons
            )
            self.out_of_chain_exons = self.out_of_chain_exons.union(
                aln_portion.out_of_chain_exons
            )
            for exon in aln_portion.exons:
                self.exon2portion[exon] = i
                self.exon2chrom[exon] = aln_portion.chrom
                self.exon2chain[exon] = aln_portion.chain
                self.unaligned_exons.add(exon) if aln_portion.was_not_aligned else None
                self.exon2asmgaps[exon] = aln_portion.assembly_gap
            curr_exon: int = min(aln_portion.exons)
            exon_start_encountered: bool = False
            asssembled_start_encountered: bool = False
            asmbl_gap_encountered: bool = False
            coding_nuc_encountered: bool = False
            clipped_start, clipped_stop = None, None
            rel_exon_start, rel_exon_stop = None, None
            for k, base in enumerate(aln_portion.reference):
                j: int = k + prev_aln_len
                query_base: str = aln_portion.query[k]
                if is_symbol(base):
                    if not exon_start_encountered:
                        # continue
                        exon_start_encountered = True
                        rel_exon_start: int = j
                    if query_base.upper() == 'N':
                        asmbl_gap_encountered = True
                    else:
                        if not asssembled_start_encountered:# and (query_base.isalpha() or asmbl_gap_encountered):
                            asssembled_start_encountered = True
                            clipped_start: int = j
                        if asssembled_start_encountered and asmbl_gap_encountered and not coding_nuc_encountered and query_base.isalpha():
                            clipped_start: int = j
                        if asssembled_start_encountered and not asmbl_gap_encountered:
                            clipped_stop: int = j
                        if query_base.isalpha():
                            asmbl_gap_encountered = False
                            coding_nuc_encountered = True
                    # if query_base.upper() != 'N':
                    #     if not asssembled_start_encountered and last_defined_base:
                    #         asssembled_start_encountered = True
                    #         clipped_start: int = j
                    #     if clipped_start and query_base != '-':
                    #         clipped_stop: int = j
                if base in (' ', '>') or j == curr_aln_len - 1:
                    if not exon_start_encountered:
                        continue
                    rel_exon_stop: int = j + int(j == curr_aln_len - 1)
                    if rel_exon_start is None:
                        raise RuntimeError(f'Exon start missing for exon {curr_exon}')
                    self.rel_exon_coords[curr_exon] = Coords(rel_exon_start, rel_exon_stop)
                    self.cesar_rel_coords[curr_exon] = Coords(rel_exon_start, rel_exon_stop)
                    clipped_start = (
                        clipped_start if clipped_start is not None else rel_exon_start
                    )
                    clipped_stop = (
                        clipped_stop if (
                            clipped_stop is not None and clipped_stop != j - 1
                        ) else rel_exon_stop
                    )
                    if clipped_stop - clipped_start < 3:
                        clipped_start, clipped_stop = rel_exon_start, rel_exon_stop
                    self.clipped_rel_coords[curr_exon] = Coords(clipped_start, clipped_stop)
                    abs_exon_start: int = self._abs_coord(clipped_start, portion=i)
                    abs_exon_stop: int = self._abs_coord(clipped_stop, portion=i)
                    self.abs_exon_coords[curr_exon] = Coords(abs_exon_start, abs_exon_stop)

                    exon_start_encountered = False
                    asssembled_start_encountered = False
                    rel_exon_start, rel_exon_stop = None, None
                    clipped_start, clipped_stop = None, None
                    curr_exon += 1
            self._update_spliceai_predictions(aln_portion, num=i)
            if aln_portion.subexon_coordinates:
                self.introns_gained = {**self.introns_gained, **aln_portion.subexon_coordinates}

        self.exon_num: int = max(self.exon2portion) ## contains the exon number in the reference
        self.last_aligned_base: int = self.rel_exon_coords[self.exon_num].stop
        self.has_spliceai_data: bool = has_spliceai_data(self.spliceai_sites)
        self.donor_site_prob: Dict[int, float] = {x: 0.0 for x in range(1, self.exon_num + 1)}
        self.acc_site_prob: Dict[int, float] = {x: 0.0 for x in range(1, self.exon_num + 1)}
        self.cesar_donor_support: Dict[int, bool] = {
            x: False for x in range(1, self.exon_num + 1)
        }
        self.cesar_acc_support: Dict[int, bool] = {
            x: False for x in range(1, self.exon_num + 1)
        }
        self.spliceai_donor_support: Dict[int, bool] = {
            x: False for x in range(1, self.exon_num + 1)
        }
        self.spliceai_acc_support: Dict[int, bool] = {
            x: False for x in range(1, self.exon_num + 1)
        }
        self.found_in_exp_locus: Dict[int, bool] = {
            x: False for x in range(1, self.exon_num + 1)
        }
        self.corrected_frame: Dict[int, bool] = {
            x: (False) for x in range(1, self.exon_num + 1)
        }
        self.intron2phase: Dict[int, int] = {x: 0 for x in range(1, self.exon_num)}
        self.exon_evidence: Dict[int, int] = {x: 2 for x in range(1, self.exon_num + 1)} ## TEMPORARY Replace with the actual exon evidence structure

        self.all_ref_codons: Dict[int, str] = {} ## stores the reference triplets
        self.all_query_codons: Dict[int, str] = {} ## stores the query triplets
        self.ref_aa_seq: Dict[int, str] = {} ## stores translated reference sequence by position
        self.query_aa_seq: Dict[int, str] = {} ## stores translated query sequence by position
        self.triplet_coordinates: Dict[int, List[int]] = {} ## stores the triplet coordinates in the query
        self.codon_coordinates: Dict[int, List[int]] = {} ## stores the reference codon coordinates in the query
        self.triplet2ref_codon: Dict[int, int] = {} ## stores the codon alignment triplet-to-reference codon mapping
        self.ref_codon2triplets: Dict[int, List[int]] = defaultdict(list) ## maps reference codons to corresponding aligned triplets
        self.exon2codons: Dict[int, Tuple[int]] = {} ## stores the exon-to-marginal_codons mapping
        self.ref_codons_to_mask: Dict[int, str] = {} ## stores real values of masked codons in the reference
        self.query_codons_to_mask: Dict[int, str] = {} ## stores real values of masked codons in the query
        self.query_atgs: List[int] = [] ## stores codon numbers for post-start ATG codons in the query
        self.split_codon_struct: Dict[int, Dict[int, int]] = defaultdict(lambda: defaultdict(int)) ## stores the portion corresponding to each affected exon for each split codon

        self.intron_deleted: Dict[int, bool] = {x: False for x in range(1, self.exon_num)} ## intron loss mapper
        self.intersects_gap: Dict[int, bool] = {x: False for x in range(1, self.exon_num + 1)} ## assembly gap intersection status storage
        self.mutation_list: List[Mutation] = [] ## stores the Mutation objects in the order the mutations are encountered
        self.selenocysteine_codons: List[SelenoCysteineTuple] = [] ## stores inframe stop codons corresponding to masked codons in the reference
        self.exon_presence: Dict[int, str] = {} ## stores exon presence statuses
        self.exon_quality: Dict[int, str] = {} ## stores exon quality assignments
        self.codon2mutations: Dict[int, List[int]] = defaultdict(list) ## for fast codon-to-mutation mapping
        self.ref_codon2mutations: Dict[int, List[int]] = defaultdict(list) ## fast mapping for reference codons to mutations
        self.exon2mutations: Dict[int, List[int]] = defaultdict(list) ## for fast exon-to-mutation mapping
        self.deletion2mask: Dict[int, bool] = {x: False for x in range(1, self.exon_num + 1)} ## deletion masking status mapper

        self.abs_nuc_ids: Dict[int, int] = {} ## stores absolute nucleotide identity values per reference exon
        self.blosum_self_values: Dict[int, int] = {} ## stores self BLOSUM score per reference exon
        self.blosum_diff_values: Dict[int, int] = {} ## stores BLOSUM score between reference and query exons per reference exon
        self.nuc_ids: Dict[int, float] = {} ## contains nucleotide identity percentage values per reference exon
        self.blosum_ids: Dict[int, float] = {} ## contains BLOSUM identity percentage values per reference exon
        self.nuc_id: float = 0.0 ## overall nucleotide identity value
        self.blosum_id: float = 0.0 ## overall BLOSUM identity value

        self.frameshift_counter: int = 1
        self.nonsense_counter: int = 1
        self.indel_counter: int = 1
        self.acc_mutation_counter: int = 1
        self.donor_mutation_counter: int = 1
        self.compensation_counter: int = 1
        self.deleted_exon_counter: int = 1
        self.missing_exon_counter: int = 1
        self.mdel_counter: int = 1
        self.frameshift2compensation: Dict[int, int] = {} ## for fast frameshift-to-compensation mapping
        self.alternative_frames: Dict[int, Tuple[int, int, bool]] = {} ## for intervals between compensated frameshifts
        self.spliceai_compensations: List[int] = [] ## for mutations numbers of compensation events involving SpliceAI correction
        self.alt_frame_masked: Dict[int, List[int]] = {} ## to store indices of mutations masked by alternative reading frames
        self.comps_to_delete: List[int] = [] ## a list of compensation entry indices to be removed
        self.stop_updated: bool = False ## switches if rescue_missing_stop is set and alternative stop codon is found

        self.loss_status: str = ''
        self.longest_fraction_strict: float = 0.0
        self.longest_fraction_relaxed: float = 0.0
        self.total_intact_fraction: float = 0.0
        self.missing_portion_of_ref: float = 0.0
        self.middle_is_intact: bool = False
        self.middle_is_present: bool = False

        ## here come the rest of attributes
        self.name: str = (
            f'{self.transcript}#{",".join(self.chains)}'
            # f'_segment{self.segment}_configuration{self.configuration}'
        )
        self.svg: str = ''
        self.v: bool = verbose

    def _to_log(self, msg: str, level: str = 'info') -> None:
        """Report a line to standard output if verbosity is enabled"""
        getattr(self.logger, level)(msg) if self.v else None

    def run(self) -> None:
        """
        The main method of the class; executes the processing pipeline in the
        following order:
        1) Given the exon alignment and SpliceAI-predicted splice site
           coordinates, infers the most likely exon borders;
        2) Processes the raw CESAR alignment results; given the exon borders,
           renders exon codon and amino acid alignments, calculates nucleotide-
           and amino acid-level properties, records all mutations, and classifies
           all exons into intact (I), missing (M), or deleted (D);
        3) Calculates nucleotide and BLOSUM identity values for the whole projection;
        4) Checks the recorded frameshifting mutations for possible
           frame-preserving compensation pairs
        5) Checks whether the projection retains an intact coding sequences after
           exons classified as deleted are removed from it;
        6) Masks the mutations recorded during the sequence processing step
           according to their location in the sequence;
        7) Performs high-level projection feature calculation;
        8) Based on the results of the previous steps, assigns the loss status
           to the projection
        """
        ## define exon borders based on CESAR and SpliceAI results
        self.define_exon_borders()
        ## process the alignment based on defined exon borders:
        ## split the alignment into codons, record all mutations found,
        ## and classify exons; arguably the longest step in the algorithm
        self.process()
        ## classify the exons
        self.classify_exons()
        ## check the frameshifts found for potential compensations
        self.find_compensating_frameshifts()
        ## if there are deleted exons, find those whose deletions are compensated
        self.find_compensating_deletions()
        ## iterate over the mutation list to determine which mutations should be masked
        self.remask_mutations()
        ## compute the intact percentage properties
        self.compute_intact_percentage()
        ## finally, classify the projection by estimating its presence status
        ## also, calculate the overall nucleotide and BLOSUM %id
        self.classify_projection()
        ## if any of compensation events were deprecated, remove them
        self.update_mutation_list()
        ## prepate an SVG plot for the projection
        self.plot_svg()

    def define_exon_borders(self) -> None:
        """
        For each exon, defines its final borders in the alignment.
        By default
        """
        ## if no SpliceAI data were provided, simply report the CESAR predictions
        phase: int = 0
        ## now, iterate over introns, defining donor for 5-flanking exon
        ## and acceptor for 3-flanking exon
        for intron in range(1, self.exon_num):
            upstream_exon: int = intron
            downstream_exon: int = intron + 1
            if upstream_exon in self.unaligned_exons:
                continue
            if downstream_exon in self.unaligned_exons:
                continue
            up_portion: int = self.exon2portion[upstream_exon]
            down_portion: int = self.exon2portion[downstream_exon]
            up_strand: bool = self._exon2strand(upstream_exon)
            down_strand: bool = self._exon2strand(downstream_exon)
            ## get the coordinates for up- and downstream neighbors for the intron
            upstream_rel_start, upstream_rel_stop = self.rel_exon_coords[upstream_exon].tuple()
            # upstream_rel_start, upstream_rel_stop = self.clipped_rel_coords[upstream_exon].tuple()
            upstream_abs_start, upstream_abs_stop = self.abs_exon_coords[upstream_exon].tuple()
            downstream_rel_start, downstream_rel_stop = self.rel_exon_coords[downstream_exon].tuple()
            # downstream_rel_start, downstream_rel_stop = self.clipped_rel_coords[downstream_exon].tuple()
            downstream_abs_start, downstream_abs_stop = self.abs_exon_coords[downstream_exon].tuple()
            ## define the intron phase
            phase: int = self._get_intron_phase(intron)
            ## check both splice sites for their conservation
            donor_intact: bool = self._check_splice_site(
                upstream_exon, upstream_rel_stop, acc=False, check_only=True
            )
            acc_intact: bool = self._check_splice_site(
                downstream_exon, downstream_rel_start, acc=True, check_only=True
            )
            ## record the original intron length
            init_intron_len: int = self._intron_length(
                upstream_rel_stop, downstream_rel_start, upstream_exon, downstream_exon
            )
            ## if precise intron loss has occurred, do not correct the splice sites
            if not donor_intact and not acc_intact:
                intron_deleted, backtrack = self._check_deleted_intron(
                    upstream_rel_stop, downstream_rel_start, upstream_exon, downstream_exon
                )
                ## If ultrashort intron has been introduced instead of recording
                ## precise intron deletion, adjust the downstream exon's start
                if intron_deleted:
                    if backtrack and self.correct_ultrashort_introns:
                        self.rel_exon_coords[downstream_exon].start = upstream_rel_stop
                        self.clipped_rel_coords[downstream_exon].start = upstream_rel_stop
                        self.abs_exon_coords[downstream_exon].start = self._abs_coord(
                            upstream_rel_stop, portion=down_portion#up_portion
                        )
                    self.intron_deleted[intron] = True
                    ## put a placeholder for intron's phase
                    self.intron2phase[intron] = -1
                    continue
            ## otherwise, it's splice site correction time!
            ## first, define the intron's class
            ## start with the upstream neighbor's donor site
            donor_dinuc: str = self.query[
                upstream_rel_stop:self._safe_step(upstream_rel_stop, upstream_exon, 2)
            ].lower()
            ref_donor_canonical_u2: bool = (
                self.u12[upstream_exon]['donor'][0] == 'U2' and
                self.u12[upstream_exon]['donor'][1] in RIGHT_SPLICE_CORR
            )
            ref_donor_u12: bool = self.u12[upstream_exon]['donor'][0] == 'U12'
            ref_donor_canonical_u12: bool = (
                ref_donor_u12 and
                self.u12[upstream_exon]['donor'][1].lower() == RIGHT_SPLICE_CORR_U12
            )
            acc_dinuc: str = self.query[
                self._safe_step(downstream_rel_start, downstream_exon, -2):downstream_rel_start
            ].lower()
            ref_acc_canonical_u2: bool = (
                self.u12[downstream_exon]['acceptor'][0] == 'U2' and
                self.u12[downstream_exon]['acceptor'][1] in LEFT_SPLICE_CORR
            )
            ref_acc_u12: bool = self.u12[downstream_exon]['acceptor'][0] == 'U12'
            ref_acc_canonical_u12: bool = (
                ref_acc_u12 and
                self.u12[downstream_exon]['acceptor'][1].lower() == LEFT_SPLICE_CORR_U12
            )
            is_canonical_u2: bool = ref_acc_canonical_u2 and ref_donor_canonical_u2
            is_u12: bool = ref_acc_u12 and ref_donor_u12
            is_canonical_u12: bool = ref_acc_canonical_u12 and ref_donor_canonical_u12
            is_non_canon_u12: bool = is_u12 and not is_canonical_u12

            donor_sites: Dict[int, float] = self.spliceai_sites[upstream_exon]['donor']
            donor_strand: bool = self._exon2strand(upstream_exon)
            donor_search_stop: Union[int, None] = self.exon_search_spaces[upstream_exon][
                int(donor_strand)
            ]
            if donor_search_stop is not None:
                donor_search_stop = self._rel_coord(donor_search_stop)
            init_donor_prob: float = max(
                donor_sites.get(
                    upstream_rel_stop + int(not up_strand), 0.0
                ),
                donor_sites.get(
                    self._closest_non_gap(upstream_rel_stop + int(not up_strand) - 1, upstream_exon) + 1, 
                    0.0
                )
            )
            donor_supported: bool = init_donor_prob >= self.min_splice_prob
            valid_donor_sites: List[Tuple[int, float]] = [(upstream_rel_stop, init_donor_prob)]
            ## splice sites are corrected if either the CESAR-found site is
            ## mutated or it is not supported by SpliceAI AND correct_all mode is set

            ## define whether the current case complies with the selected
            ## correction mode; mode levels are hierachical,
            ## meaning that for level N all previous levels' prerequisites
            ## are also valid

            ## correction mode is set to 2 or higher:
            ## correct only mutated canonical U2 donors
            donor_mode2_corr: bool = (
                not donor_intact and ref_donor_canonical_u2 and 
                self.correction_mode > 1
            )
            ## correction mode is set to 3 or higher:
            ## correct mutated and unsupported canonical U2 donors
            donor_mode3_corr: bool = (
                not donor_supported  and ref_donor_canonical_u2 and
                self.correction_mode > 2 #== 3
            )
            ## correction mode is set to 4 or higher:
            ## correct mutated and unsupported U2 donors as well as U12 sites
            ## deviating from the GT site
            donor_mode4_corr: bool = (
                ref_donor_canonical_u12 and
                not donor_supported and
                donor_dinuc != RIGHT_SPLICE_CORR_U12 and
                self.correction_mode > 3
            )
            ## correction mode is set to 5:
            ## correct mutated and unsupported canonical U2 and canonical U12 donors
            donor_mode5_corr: bool = (
                ref_donor_canonical_u12 and not donor_supported and self.correction_mode > 4
            )
            ## correction mode is set to 6:
            ## correct mutated and unsupported sites for both U12 and
            ## canonical U2 donors
            donor_mode6_corr: bool = (
                ref_donor_u12 and not donor_supported and self.correction_mode == 6
            )
            ## correction mode is set to 7: correct non-canonical U2 exons
            ## with no SpliceAI support on top of all the previous conditions
            donor_mode7_corr: bool = (
                not ref_donor_u12 and not ref_donor_canonical_u2 and
                not donor_supported and self.correction_mode > 6
            )
            donor_needs_correction: bool = any(
                (
                    donor_mode2_corr, donor_mode3_corr, donor_mode4_corr,
                    donor_mode5_corr, donor_mode6_corr, donor_mode7_corr
                    # donor_mode6_corr
                )
            )
            donor_frame_shift: int = self._get_frameshift(upstream_rel_start, upstream_rel_stop)
            if donor_needs_correction:
                for donor, prob in donor_sites.items():
                    donor -= int(not up_strand)
                    if prob < self.min_splice_prob:
                        continue
                    ## corrections leading to negative exon length are not allowed
                    if donor < upstream_rel_start:
                        self._to_log(
                            f'Donor site at {donor} overlaps the previous exon '
                            f'for exon {intron}'
                        )
                        continue
                    ## neither are corrections leading to exon overlap
                    if (donor + MIN_INTRON_LENGTH) >= downstream_rel_start:
                        self._to_log(
                            f'Donor site at {donor} overlaps the next exon '
                            f'for exon {intron}'
                        )
                        continue
                    ## if a search space (expected locus + exon locus flank)
                    ## is defined for the exon, make sure the correction does not
                    ## exceed the space's downstream boundary
                    if donor_search_stop is not None and donor > donor_search_stop:
                        self._to_log(
                            f'Donor site at {donor} exceeds the expected search space '
                            f'for exon {intron}'
                        )
                        continue
                    ## corrections must preserve the exon's reading frame
                    if not self._site_is_inframe(upstream_rel_stop, donor, donor_frame_shift, donor=True):
                        self._to_log(
                            f'Donor site at {donor} with probability {prob} is out of frame for exon {intron}'
                        )
                        continue
                    ## and should not introduce inframe stop codons
                    if not self._nonsense_free_correction(
                        upstream_rel_stop - phase, donor - phase
                    ):
                        self._to_log(
                            f'Donor site at {donor} introduces a nonsense '
                            f'mutation for exon {intron}'
                        )
                        continue
                    if (
                        not donor_supported and init_donor_prob > 0 and
                        init_donor_prob + self.prob_margin > prob
                    ):
                        self._to_log(
                            f'Donor site at {donor} does not exceed initial '
                            f'probability {init_donor_prob} by the margin of '
                            f'{self.prob_margin}'
                        )
                        continue
                    if upstream_exon in self.introns_gained:
                        if donor <= self.introns_gained[upstream_exon][-1][0]:
                            self._to_log(
                                f'Donor site at {donor} rules out query-specific introns '
                                'defined by SpliceAI'
                            )
                            continue
                    ## once all the conditions are met, add the splice site pair
                    ## to the current pool
                    valid_donor_sites.append((donor, prob))
            ## sort the splice sites by SpliceAI probability in the inverse order
            valid_donor_sites.sort(key=lambda x: -x[1])
            ## now do the same for the downstream exon's acceptor site
            acc_sites: Dict[int, float] = self.spliceai_sites[downstream_exon]['acceptor']
            init_acc_prob: float = max(
                acc_sites.get(
                    self._closest_non_gap(downstream_rel_start - int(down_strand), downstream_exon),
                    0.0
                ),
                acc_sites.get(downstream_rel_start - int(down_strand), 0.0)
            )
            acc_supported: bool = init_acc_prob >= self.min_splice_prob
            acc_strand: bool = self._exon2strand(downstream_exon)
            acc_search_stop: Union[int, None] = self.exon_search_spaces[downstream_exon][
                int(not acc_strand)
            ]

            ## correction mode is set to 2 or higher:
            ## correct only mutated canonical U2 donors
            acc_mode2_corr: bool = (
                not acc_intact and ref_acc_canonical_u2 and
                self.correction_mode > 1
            )
            ## correction mode is set to 3 or higher:
            ## correct mutated and unsupported canonical U2 donors
            acc_mode3_corr: bool = (
                ref_acc_canonical_u2 and not acc_supported and
                self.correction_mode > 2 #== 3
            )
            ## correction mode is set to 4 or higher:
            ## correct mutated and unsupported U2 donors as well as U12 sites
            ## deviating from the GT site
            acc_mode4_corr: bool = (
                ref_acc_canonical_u12 and
                not acc_supported and
                acc_dinuc.lower() != LEFT_SPLICE_CORR_U12 and
                self.correction_mode > 3
            )
            ## correction mode is set to 5 or higher:
            ## correct mutated and unsupported sites for both U12 and
            ## canonical U2 donors
            acc_mode5_corr: bool = (
                ref_acc_canonical_u12 and not acc_supported and self.correction_mode > 4
            )
            ## correction mode is set to 6:
            ## correct mutated and unsupported sites for both U12 and
            ## canonical U2 donors
            acc_mode6_corr: bool = (
                is_u12 and not acc_supported and self.correction_mode == 6
            )
            ## correction mode is set to 7: correct non-canonical U2 exons
            ## with no SpliceAI support on top of all the previous conditions
            acc_mode7_corr: bool = (
                not is_u12 and not ref_acc_canonical_u2 and
                not acc_supported and self.correction_mode > 6
            )
            acc_needs_correction: bool = any(
                (
                    acc_mode2_corr, acc_mode3_corr, acc_mode4_corr,
                    acc_mode5_corr, acc_mode6_corr, acc_mode7_corr
                    # acc_mode6_corr
                )
            )

            if acc_search_stop is not None:
                acc_search_stop = self._rel_coord(acc_search_stop)
            valid_acc_sites: List[Tuple[int, float]] = [(downstream_rel_start, init_acc_prob)]
            acc_frame_shift: int = self._get_frameshift(downstream_rel_start, downstream_rel_stop)
            if acc_needs_correction:
                for acc, prob in acc_sites.items():
                    acc += int(down_strand)
                    if prob < self.min_splice_prob:
                        continue
                    if acc > downstream_rel_stop:
                        self._to_log(
                            (
                                'Acceptor site at %i overlaps the next exon '
                                'for exon %i'
                            ) % (acc, intron + 1)
                        )
                        continue
                    if (acc - MIN_INTRON_LENGTH) <= upstream_rel_stop:
                        self._to_log(
                            (
                                'Acceptor site at %i overlaps the previous exon '
                                'for exon %i'
                            ) % (acc, intron + 1)
                        )
                        continue
                    if acc_search_stop is not None and acc < acc_search_stop:
                        self._to_log(
                            (
                                'Acceptor site at %i exceeds the expected search '
                                'space for exon %i'
                            ) % (acc, intron + 1)
                        )
                        continue
                    if not self._site_is_inframe(downstream_rel_start, acc, acc_frame_shift, donor=False):
                        self._to_log(
                            f'Acceptor site at {acc} with probability {prob} is out of frame '
                            f'for exon {intron+1}'
                        )
                        continue
                    if not self._nonsense_free_correction(
                        acc + (3 - phase) % 3, downstream_rel_start + (3 - phase) % 3
                    ):
                        self._to_log(
                            f'Acceptor site at {acc} introduces a nonsense '
                            f'mutation for exon {intron+1}'
                        )
                        continue
                    if (
                        not acc_supported and init_acc_prob > 0 and
                        init_acc_prob + self.prob_margin > prob
                    ):
                        self._to_log(
                            f'Acceptor site at {acc} does not exceed initial '
                            f'probability {init_acc_prob} by the margin of '
                            f'{self.prob_margin}'
                        )
                        continue
                    if downstream_exon in self.introns_gained:
                        if acc >= self.introns_gained[downstream_exon][0][1]:
                            self._to_log(
                                f'Acceptor site at {acc} rules out query-specific introns '
                                'defined by SpliceAI'
                            )
                            continue
                    valid_acc_sites.append((acc, prob))
            valid_acc_sites.sort(key=lambda x: -x[1])
            ## now, get all the splice site combinations
            intron_combinations: List[Tuple[Tuple[int, float]]] = sorted(
                ((x, y) for x in valid_donor_sites for y in valid_acc_sites),
                key=lambda x: (-x[0][1], -x[1][1])
            )
            ## by default, assign the CESAR-predicted site coordinates to the sites of choice
            best_donor: int = upstream_rel_stop
            best_acc: int = downstream_rel_start
            best_donor_prob: float = init_donor_prob
            best_acc_prob: float = init_acc_prob
            for ((donor, donor_prob), (acc, acc_prob)) in intron_combinations:
                abs_donor: int = self._abs_coord(donor, portion=up_portion)
                abs_acc: int = self._abs_coord(acc, portion=down_portion)
                if donor >= acc:
                    continue
                if getattr(
                    abs_donor, '__le__' if up_strand else '__ge__'
                )(upstream_abs_start if up_strand else upstream_abs_stop):
                    self._to_log(
                        f'Donor {donor} with probability {donor_prob} lies before the start of exon {upstream_exon}'
                    )
                    continue
                if getattr(
                    abs_acc, '__ge__' if down_strand else '__le__'
                )(downstream_abs_stop if down_strand else downstream_abs_start):
                    self._to_log(
                        f'Acceptor {acc} with probability {acc_prob} '
                        f'lies after the end of exon {downstream_exon}'
                    )
                    continue
                if self._alternative_split_contains_stop(donor, acc, phase):
                    self._to_log(
                        f'Alternative intron combination {((donor, donor_prob), (acc, acc_prob))} introduces a stop codon'
                    )
                    continue
                if self._zero_length_exon_introduced(intron, donor, acc):
                    self._to_log(
                        f'Alternative intron combination {((donor, donor_prob), (acc, acc_prob))} results in a zero nucleotide exon'
                    )
                    continue
                best_donor = donor
                best_acc = acc
                best_donor_prob = donor_prob
                best_acc_prob = acc_prob
                break

            ## update the intron support data
            self.donor_site_prob[upstream_exon] = best_donor_prob
            self.acc_site_prob[downstream_exon] = best_acc_prob
            donor_cesar_support: bool = upstream_rel_stop == best_donor
            acc_cesar_support: bool = downstream_rel_start == best_acc
            self.cesar_donor_support[upstream_exon] = donor_cesar_support
            self.cesar_acc_support[downstream_exon] = acc_cesar_support
            donor_spliceai_support: bool = best_donor_prob > self.min_splice_prob
            acc_spliceai_support: bool = best_acc_prob > self.min_splice_prob
            self.spliceai_donor_support[upstream_exon] = donor_spliceai_support
            self.spliceai_acc_support[downstream_exon] = acc_spliceai_support
            self.corrected_frame[upstream_exon] |= (not donor_cesar_support and bool(donor_frame_shift))
            self.corrected_frame[downstream_exon] |= (not acc_cesar_support and bool(acc_frame_shift))

            ## if neither of the sites were corrected and the original intron was
            ## too short, correct the intron and proceed
            upd_intron_len: int = self._intron_length(
                best_donor, best_acc, upstream_exon, downstream_exon
            )
            intron_uncorrected: bool = upd_intron_len == init_intron_len
            intron_too_short: bool = init_intron_len is not None and init_intron_len < MAX_RETAINED_INTRON_LEN
            _istart, _istop = self.rel_exon_coords[upstream_exon].tuple()

            if intron_uncorrected and intron_too_short and self.correct_ultrashort_introns:
                self.rel_exon_coords[downstream_exon].start = upstream_rel_stop
                self.clipped_rel_coords[downstream_exon].start = upstream_rel_stop
                self.abs_exon_coords[downstream_exon].start = self._abs_coord(upstream_rel_stop, portion=up_portion)
                self.intron_deleted[intron] = True
                ## put a placeholder for intron's phase
                self.intron2phase[intron] = -1
                continue

            ## update the intron boundaries
            self._to_log(
                f'Best coordinates for intron {intron} are: donor={best_donor} ({best_donor_prob}), acceptor={best_acc} ({best_acc_prob})'
            )
            if self._abs_coord(best_donor, portion=up_portion) == (upstream_abs_stop if up_strand else upstream_abs_start):
                best_donor = upstream_rel_stop
            self.rel_exon_coords[upstream_exon].stop = best_donor
            if best_donor != upstream_rel_stop:
                self.clipped_rel_coords[upstream_exon].stop = best_donor ## sic!
            self.abs_exon_coords[upstream_exon].stop = self._abs_coord(best_donor, portion=up_portion)
            if self._abs_coord(best_acc, portion=down_portion) == (downstream_abs_start if down_strand else downstream_abs_stop):
                best_acc = downstream_rel_start
            self.rel_exon_coords[downstream_exon].start = best_acc
            if best_acc != downstream_rel_start:
                self.clipped_rel_coords[downstream_exon].start = best_acc
            self.abs_exon_coords[downstream_exon].start = self._abs_coord(best_acc, portion=down_portion)
            ## get the intron's phase
            self.intron2phase[intron] = self._get_intron_phase(intron)

        ## for the first exon, estimate whether start codon is missing
        if self.rescue_missing_start:
            self._rescue_alternative_start()
        ## for the last codon, estimate whether the stop codon is missing
        if self.rescue_missing_stop:
            self._rescue_alternative_stop()

    def process(self) -> None:
        """
        Processes the aligned sequences, inferring the alignment structure
        as well as exon properties.
        """
        ref_codon: str = ''
        query_codon: str = ''
        codon_coords: List[int] = []
        triplet_num: int = 1
        ref_codon_num: int = 1
        aa_num: int = 1
        ins_length: int = 0
        del_length: int = 0
        asmbl_gap: int = 0
        split_portion: int = 0
        initial_ref_stop: str = ''
        codon_spliceai_diff: bool = False
        for exon in range(1, self.exon_num + 1):
            self._to_log(f'Processing exon {exon}')
            portion: int = self.exon2portion[exon]
            aln_len: int = 0
            abs_nuc_id: int = 0
            blosum_self: int = None
            blosum_diff: int = None
            ref_exon_len: int = 0
            cesar_start: int = self.cesar_rel_coords[exon].start
            adj_start: int = self.rel_exon_coords[exon].start
            clipped_start: int = self.clipped_rel_coords[exon].start
            cesar_stop: int = self.cesar_rel_coords[exon].stop
            adj_stop: int = self.rel_exon_coords[exon].stop
            clipped_stop: int = self.clipped_rel_coords[exon].stop
            start: int = min(cesar_start, adj_start)
            stop: int = max(cesar_stop, adj_stop)
            first_codon: int = triplet_num
            self.exon2codons[exon] = (first_codon,)
            first_ref_codon: int = ref_codon_num
            self.exon2ref_codons[exon] = (first_ref_codon,)
            self.triplet2ref_codon[first_codon] = first_ref_codon
            ## logic for clipping unassembled termini: do not account for assembly gaps
            exon_is_unclipped: bool = (
                clipped_start == adj_start and clipped_stop == adj_stop
            )
            contains_term_ns: bool = False
            clipped_start_updated: bool = False
            clipped_stop_updated: bool = False
            if exon in self.introns_gained:
                intron_intervals: List[Tuple[int, int]] = [
                    (self.introns_gained[exon][x-1][1], self.introns_gained[exon][x][0])
                    for x in range(1, len(self.introns_gained[exon]))
                ]
            else:
                intron_intervals: List[Tuple[int, int]] = []
            ## check acceptor site
            _ = self._check_splice_site(exon, adj_start, acc=True, check_only=False)
            self._to_log(f'{exon=}, {cesar_start=}, {cesar_stop=}, {adj_start=}, {adj_stop=}, {clipped_start=}, {clipped_stop=}, {self.abs_exon_coords[exon]=}')
            for i in range(start, stop):
                n1: str = self.reference[i]
                n2: str = self.query[i]
                is_unclipped: bool = exon_is_unclipped or (
                    clipped_start <= i <= clipped_stop and n2.upper() != 'N'
                )
                contains_term_ns |= not is_unclipped
                intronic: bool = any(x[0] <= i < x[1] for x in intron_intervals)
                if intronic:
                    n2 = '-'
                else:
                    aln_len += 1
                spliceai_diff: bool = (
                    (i >= cesar_start and i < adj_start) or
                    (i < cesar_stop and i >= adj_stop)
                )
                codon_spliceai_diff |= spliceai_diff
                if is_symbol(n1) and (n2.isupper() or n2 == '-') and not contains_term_ns and not intronic:
                    abs_nuc_id += int(n1.upper() == n2.upper())
                if n1.isalpha():
                    ref_exon_len += 1
                else:
                    n1 = '-'
                if n1.islower():
                    split_portion += 1
                elif split_portion and n1.isalpha():
                    split_portion = 0
                if spliceai_diff:
                    n2 = '-'
                n1 = n1.upper()
                n2 = n2.upper() ## redundant
                if n2 != '-':
                    codon_coords.append(self._abs_coord(i, portion=portion))
                if n2 == 'N':
                    asmbl_gap += 1
                else:
                    asmbl_gap = 0
                if asmbl_gap >= self.min_gap_size:
                    self.intersects_gap[exon] = True
                ## bases corresponding to gained introns
                if intronic:
                    continue
                ref_codon += n1
                query_codon += n2

                ## if current portion is subjected to terminal undefined region
                ## clipping and contains both ambiguous and defined codons,
                ## updated the absolute coordinates
                if contains_term_ns:
                    if not abs_nuc_id and (n2.upper() != 'N'):
                        adj_clipped_start: int = i + 1
                        self.abs_exon_coords[exon].start = self._abs_coord(
                            adj_clipped_start, portion=portion
                        )
                        clipped_start_updated = True
                    elif abs_nuc_id and i == stop - 1:
                        adj_clipped_stop: int = i - asmbl_gap + 1
                        self.abs_exon_coords[exon].stop = self._abs_coord(
                            adj_clipped_stop, portion=portion
                        )
                        clipped_stop_updated = True

                if not is_complete_codon(ref_codon):
                    continue
                ## if the code reached this point, the full
                ## reference codon has been encountered
                if exon == self.exon_num and self.stop_updated:
                    ## if stop codon has been shifted, record the original stop codon
                    ## but do not compare the resulting reciprocal gap codons
                    if (
                        cesar_stop - 1 in codon_coords and ref_codon in STOPS and
                        query_codon == GAP_CODON
                    ):
                        initial_ref_stop: int = ref_codon
                        continue
                    elif (
                        adj_stop - 1 in codon_coords and ref_codon == GAP_CODON and
                        query_codon in STOPS
                    ):
                        ref_codon = ref_codon[:-3] + initial_ref_stop
                codon_start: int = 0
                subcodon_len_sum: int = 0
                for sub_num, (r, q) in enumerate(process_codon_pair(ref_codon, query_codon)):
                    if '|' in r:
                        r_, r = r.split('|')
                        self.ref_codons_to_mask[triplet_num] = r
                        q_, q = q.split('|')
                        self.query_codons_to_mask[triplet_num] = q
                        codon_extent: int = len(q.replace('-', ''))
                    else:
                        ## check_codon will automatically mask
                        ## any stop codon as XXX_CODON, so the last codon
                        ## should be spared from checking
                        is_last: bool =  i >= self.last_aligned_base
                        ref_is_stop: bool = r in STOPS
                        query_is_stop: bool = q in STOPS
                        regular_stop_found: bool = is_last and ref_is_stop
                        subcodon_len_sum += len(q)
                        if not regular_stop_found:
                            r_, q_ = check_codon(r), check_codon(q)
                            if r != r_:
                                self.ref_codons_to_mask[triplet_num] = r
                            if q != q_:
                                self.query_codons_to_mask[triplet_num] = q
                            # r, q = r_, q_
                            codon_extent: int = len(q.replace('-', ''))
                    subcodon_coords: List[int] = codon_coords[
                        codon_start:codon_start + codon_extent
                    ]
                    codon_start = codon_start + codon_extent
                    if not len(subcodon_coords):
                        subcodon_coords = [
                            self._abs_coord(i - len(query_codon) + 1, portion=portion),
                            self._abs_coord(i + 1, portion=portion)
                        ]
                        codon_coords.extend(subcodon_coords)
                    self.triplet_coordinates[triplet_num] = subcodon_coords
                    self.all_ref_codons[triplet_num] = r_
                    self.all_query_codons[triplet_num] = q_
                    self.triplet2ref_codon[triplet_num] = ref_codon_num
                    self.ref_codon2triplets[ref_codon_num].append(triplet_num)
                    ## record the additional start codons in the query
                    if q_ == START:
                        self.query_atgs.append(triplet_num)
                    ## here go the mutation checks
                    ## first, check whether the start codon is in its place
                    if triplet_num == 1:
                        self._check_for_start_loss()
                    ## if any of the codons were masked, check them for
                    ## frameshift and nonsense mutations
                    if r != r_ or q != q_:
                        self._check_for_frameshift(triplet_num)
                        self._check_for_nonsense_mutation(triplet_num)
                    ## then, check for big indel mutations
                    if r == GAP_CODON:
                        ins_length += bool(not codon_spliceai_diff)#1 if not codon_spliceai_diff else 0
                    else:
                        if ins_length * 3 > BIG_INDEL_SIZE:
                            self._add_indel_mutation(triplet_num, ins_length, insertion=True)
                        ins_length = 0
                    if q == GAP_CODON:
                        del_length +=  bool(not codon_spliceai_diff)#1 if codon_spliceai_diff else 0
                    else:
                        if del_length * 3 > BIG_INDEL_SIZE:
                            self._add_indel_mutation(triplet_num, del_length, insertion=False)
                        del_length = 0
                    if exon == self.exon_num and i >= (self.last_aligned_base - 1) and not self.stop_updated:
                        self._check_for_stop_loss()
                    triplet_num += 1
                for ra, qa in process_and_translate(ref_codon, query_codon):
                    if not (split_portion or contains_term_ns):
                    # if not contains_term_ns:
                        self_score: int = get_blosum_score(ra, ra, self.matrix)
                        if blosum_self is None:
                            blosum_self = self_score
                        else:
                            blosum_self += self_score
                        diff_score: int = get_blosum_score(ra, qa, self.matrix)
                        if blosum_diff is None:
                            blosum_diff = diff_score
                        else:
                            blosum_diff += diff_score
                    self.ref_aa_seq[aa_num] = ra
                    self.query_aa_seq[aa_num] = qa
                    aa_num += 1
                split_portion = 0
                codon_spliceai_diff = False
                self.codon_coordinates[ref_codon_num] = codon_coords
                ref_codon_num += bool(ref_codon != GAP_CODON)
                ref_codon = ''
                query_codon = ''
                codon_coords = []
                contains_term_ns = False
                clipped_start_updated, clipped_stop_updated = False, False
            ## record the split portion of the last codon
            if split_portion:
                latent_triplets: int = floor(len(ref_codon) / 3)
                split_triplet: int = triplet_num + (latent_triplets if latent_triplets > 1 else 0 )
                if ref_exon_len == 1:
                    ## codon can be split between three exons only if mid exon is 1 bp long
                    self.split_codon_struct[split_triplet][exon] = 1
                    self.split_codon_struct[split_triplet][exon+1] = 1
                else:
                    self.split_codon_struct[split_triplet][exon] = split_portion
                    self.split_codon_struct[split_triplet][exon+1] = 3 - split_portion
            self.aln_lengths[exon] = aln_len
            self.abs_nuc_ids[exon] = abs_nuc_id
            self.blosum_self_values[exon] = blosum_self
            self.blosum_diff_values[exon] = blosum_diff
            self.exon2codons[exon] = (
                first_codon, split_triplet if split_portion else triplet_num
            )
            self.exon2ref_codons[exon] = (first_ref_codon, ref_codon_num)
            self.ref_exon_lengths[exon] = ref_exon_len
            ## check donor site befor going to the next exon
            _ = self._check_splice_site(exon, adj_stop, acc=False, check_only=False)
        self.reference_length: int = sum(self.ref_exon_lengths.values())

    def classify_exons(self) -> None:
        """Infers exon presence status (I, M, D) for all exons"""
        for exon in range(1, self.exon_num + 1):
            _nid: int = self.abs_nuc_ids[exon]
            _aln_len: int = self.aln_lengths[exon]
            nid: int = (_nid / _aln_len) * 100
            self.nuc_ids[exon] = nid
            _blosum_self: float = self.blosum_self_values[exon]
            _blosum_diff: float  = self.blosum_diff_values[exon]
            if _blosum_self is None or _blosum_diff is None:
                if exon in self.unaligned_exons:
                    _blosum_self = 0.0
                    _blosum_diff = 0.0
                else:
                    ref_codons: Tuple[int] = sorted(self.exon2ref_codons[exon])
                    first_codon: int = ref_codons[0]
                    first_triplet: int = min(self.ref_codon2triplets[first_codon])
                    last_triplet: int = max(self.ref_codon2triplets[first_codon])
                    ref_only_codon: int = ''.join(
                        self.ref_codons_to_mask.get(
                            x, self.all_ref_codons[x]
                        ) for x in range(first_triplet, last_triplet + 1)
                    )
                    query_only_codon: int = ''.join(
                        self.query_codons_to_mask.get(
                            x, self.all_query_codons[x]
                        ) for x in range(first_triplet, last_triplet + 1)
                    )
                    ref_aa, query_aa = next(process_and_translate(ref_only_codon, query_only_codon)) ## NEEDS TESTING
                    _blosum_self = get_blosum_score(ref_aa, ref_aa, self.matrix)
                    _blosum_diff = get_blosum_score(ref_aa, query_aa, self.matrix)
                self.blosum_self_values[exon] = _blosum_self
                self.blosum_diff_values[exon] = _blosum_diff
            blosum: int = safe_div(nn(_blosum_diff), nn(_blosum_self)) * 100
            self.blosum_ids[exon] = blosum
            exon_quality: str = assess_exon_quality(nid, blosum)
            self.exon_quality[exon] = exon_quality
            self._exon_presence_status(exon)

    def find_compensating_frameshifts(self) -> None:
        """
        Find whether certain frameshifts compensate each other.
        We call a pair/group of frameshifts mutually compensating if:
        1) the sum of frameshift values restores the original frame (e.g., +2 and +1
        or +1 and -1), and
        2) the alternative frame segment contains no inframe stop codons.
        Every compensation event is stored as a separate mutation in the global
        mutation slot containing the range of frameshift mutation affected.
        COMPENSATION mutation are not masked but the affected frameshifts are.
        """
        frameshifts: List[Tuple[int, Mutation]] = [
            (x, y) for x, y in enumerate(self.mutation_list) if
            y.mutation_class in (FS_INS, FS_DEL)
        ]
        if len(frameshifts) < 2: 
            ## a solitary frameshift has nothing to compensate 
            ## or to be compensated by
            return
        ## for each mutation other than the last one, check if it is compensated
        already_compensated: List[int] = [] ## to store frameshifts already compensated
        alt_frame_start: Union[int, None] = None ## first affected codon to measure fixed frame from
        for inum, (_, init_mut) in enumerate(frameshifts[:-1]):
            if inum in already_compensated: 
                ## compensated by the previous mutation, 
                ## no need to compensate anew
                continue
            if isinstance(init_mut.exon, int) and self.exon_presence[init_mut.exon] != 'I': 
                ## do not account for mutations 
                ## in the non-present exons
                continue
            elif any(self.exon_presence[int(x)] != 'I' for x in str(init_mut.exon).split('_')):
                continue
            fs_size: int = int(init_mut.description)
            fs_group: List[int] = [inum]
            init_exon: int = init_mut.exon
            init_exon_corrected: bool = self._exon_was_corrected(init_exon)
            ## for each succeeding mutation, check if it compensates the current group
            for jnum, (_, other_mut) in enumerate(frameshifts[inum+1:], inum+1):
                if jnum in already_compensated:
                    continue
                if isinstance(other_mut.exon, int) and self.exon_presence[other_mut.exon] != 'I':
                    continue
                elif any(self.exon_presence[int(x)] != 'I' for x in str(other_mut.exon).split('_')):
                    continue
                other_exon: int = other_mut.exon
                other_exon_corrected: bool = self._exon_was_corrected(other_exon)
                ## frameshifts already resolved by SpliceAI correction 
                ## cannot compensate frameshifts in other exons
                same_correction_event: bool = self._mutations_in_same_query_codon(init_mut, other_mut)
                spliceai_corr_migrated: bool = self._comp_moved_to_split_codon(
                    init_mut, other_mut
                ) or self._comp_moved_to_next_exon(init_mut, other_mut)
                if (init_exon_corrected or other_exon_corrected) and init_mut.exon != other_mut.exon:
                    if not (same_correction_event or spliceai_corr_migrated):
                        continue
                fs_size += int(other_mut.description)
                fs_group.append(jnum)
                if not (fs_size % 3): 
                    ## the sum of frameshifts divides by 3 with no remainder
                    ## define which codons are affected by the frameshift
                    if alt_frame_start is None:
                        first_affected: int = self._next_complete_codon(frameshifts[fs_group[0]][1].codon, prev=True)
                    else:
                        first_affected = alt_frame_start
                    last_affected: int = self._next_complete_codon(frameshifts[fs_group[-1]][1].codon, prev=False)
                    alt_ref_start: int = self.triplet2ref_codon[first_affected]
                    alt_ref_stop: int = self.triplet2ref_codon[last_affected]
                    alt_start: int = min(self.ref_codon2triplets[alt_ref_start])
                    alt_stop: int = max(self.ref_codon2triplets[alt_ref_stop])
                    ## gather the affected codons, unmask if necessary
                    affected_codons: List[str] = [
                        self.query_codons_to_mask.get(
                            x, self.all_query_codons[x]
                        ) for x in range(alt_start, alt_stop)
                        if any(self.exon_presence[y] == 'I' for y in self._codon2exon(x))
                    ]
                    ## concatenate the codons, remove deletion symbols and split the result into triplets
                    alt_frame: List[str] = parts(''.join(affected_codons).replace('-', ''), 3)
                    ## check the resulting frame for nonsense mutations
                    spliceai_mitigated: bool = any(
                            self._in_modified_seq(
                                frameshifts[x][1].exon, frameshifts[x][1].codon
                            ) for x in fs_group
                        )
                    print(f'{spliceai_mitigated=}, {same_correction_event=}, {spliceai_corr_migrated=}')
                    if not STOPS.intersection(alt_frame) or spliceai_mitigated: 
                        ## create a mutation object and add it to the global pool
                        start_id: str = frameshifts[fs_group[0]][1].mutation_id
                        start_num: str = start_id.split('_')[1]
                        end_id: str = frameshifts[fs_group[-1]][1].mutation_id
                        end_num: str = end_id.split('_')[1]
                        description: str = f'FS_{start_num}-{end_num}'
                        mut_id: str = f'C_{self.compensation_counter}'
                        mutation: Mutation = Mutation(
                            self.transcript,
                            self.exon2chain[init_mut.exon],
                            init_mut.exon,
                            init_mut.codon,
                            init_mut.ref_codon,
                            init_mut.chrom,
                            init_mut.start,
                            init_mut.stop,
                            COMPENSATION,
                            description,
                            False, ## compensated mutation entries are not masked
                            '-', ## since they are not masked, the reason slot is kept empty
                            mut_id
                        )
                        self.mutation_list.append(mutation)
                        self.compensation_counter += 1
                        mut_num: int = len(self.mutation_list) - 1
                        ## update the list of compensated frameshifts
                        already_compensated.extend(fs_group)
                        ## mask the compensated frameshifts in the global pool
                        ## also, save frameshift to compensation mapping
                        for compensated in fs_group:
                            orig_number: int = frameshifts[compensated][0]
                            self.mutation_list[orig_number].is_masked = True
                            self.mutation_list[orig_number].masking_reason = COMPENSATION_REASON
                            self.frameshift2compensation[orig_number] = mut_num
                        muts_to_mask: List[int] = []
                        for k, mut in enumerate(self.mutation_list):
                            if mut.mutation_class != STOP:
                                continue
                            if isinstance(mut.codon, str):
                                first_codon, last_codon = map(int, mut.codon.split('_'))
                                if alt_start < first_codon < alt_stop:
                                    continue
                                if alt_start < last_codon < alt_stop:
                                    continue
                            else:
                                if not alt_start < mut.codon < alt_stop:
                                    continue
                            mut.is_masked = True
                            mut.masking_reason = ALT_FRAME_REASON
                            muts_to_mask.append(k)
                        self.alt_frame_masked[mut_num] = muts_to_mask
                        if alt_frame_start is None:
                            alt_frame_start = first_affected
                        comp_start_codon: int = frameshifts[fs_group[0]][1].ref_codon
                        comp_end_codon: int = frameshifts[fs_group[-1]][1].ref_codon
                        if spliceai_mitigated or same_correction_event or spliceai_corr_migrated:
                            self.spliceai_compensations.append(mut_num)
                        self.alternative_frames[mut_num] = (
                            comp_start_codon, comp_end_codon, spliceai_mitigated or same_correction_event
                        )
                    else:
                        alt_frame_start = None
                    ## once compensating mutation is found, 
                    ## break the inner loop
                    break

    def find_compensating_deletions(self) -> None:
        """
        Find which exon deletion events preserve the original frame. Similarly to
        compensating frameshifts, we call deletions compensating if:
        1) single deleted exon has its length multiple by 3;
        2) several consecutive deletions have their length sum multiple of 3;
        3) deletions affect short terminal exons
        TODO: Where was the FIRST_LAST_DEL_SIZE length threshold inferred from and should it be adjusted?
        """
        if self.exon_num == 1: 
            ## obviously makes no sense for single-exon transcripts
            return
        ## get series of consecutive exon deletions
        deletion_runs: List[int] = get_d_runs(self.exon_presence)
        if not deletion_runs:
            return
        already_compensated: List[int] = []
        for deletion_run in deletion_runs:
            if len(deletion_run) == 1: 
                ## lone deletion can be self-compensated
                lone_del: int = deletion_run[0]
                del_run: int = self.ref_exon_lengths[lone_del]
                if not (del_run % 3):
                    already_compensated.append(lone_del)
                continue
            for inum, init_del in enumerate(deletion_run):
                if init_del in already_compensated: 
                    ## ignore deletion if it has already been compensated
                    continue
                group: List[int] = [init_del]
                group_len: int = self.ref_exon_lengths[init_del]
                for other_del in deletion_run[inum+1:]:
                    group.append(other_del)
                    other_len: int = self.ref_exon_lengths[other_del]
                    group_len += other_len
                    if not (group_len % 3):
                        already_compensated.extend(group)
                        break
        ## for each safely deleted exon, update exon presence status and/or
        ## deletion-class Mutation entries
        for del_num in already_compensated:
            first: bool = del_num == 1
            last: bool = del_num == self.exon_num
            short_enough: bool = self.ref_exon_lengths[del_num] < TERMINAL_EXON_DEL_SIZE
            if (first or last) and short_enough: 
                ## short terminal deleted exons are qualified as missing
                for mut in self.exon2mutations[del_num]:
                    if self.mutation_list[mut].mutation_class == DEL_EXON:
                        self.mutation_list[mut].mutation_class = MISS_EXON
                        self.mutation_list[mut].mutation_id = f'MDEL_{self.mdel_counter}' ## TODO: this disrupt the deletion numeration!!!
                        self.mdel_counter += 1
                        # del deletion2mask[del_num]
                self.exon_presence[del_num] = 'M'
            else: 
                ## simply mask the deletion entry
                for mut in self.exon2mutations[del_num]:
                    if self.mutation_list[mut].mutation_class == DEL_EXON: ## TODO: Swap description and mutation class for deleted exons thorughout the text
                        self.mutation_list[mut].is_masked = True
                        self.deletion2mask[del_num] = True
        self._check_updated_split_codons()

    def remask_mutations(self) -> None:
        """
        For all the mutations in the global mutation pool, assesses whether
        positional masking should be performed or whether any mutations should be unmasked
        """
        muts_to_unmask: List[Tuple[int, str]] = []
        for i, mut in enumerate(self.mutation_list):
            codon: int = mut.codon
            if isinstance(codon, str) and '_' in codon:
                codons: Tuple[int] = tuple(int(x) for x in codon.split('_'))
            else:
                codons: Tuple[int] = (int(codon),)
            codon: int = mut.codon if isinstance(mut.codon, int) else int(mut.codon.split('_')[0])
            exon: Union[int, str] = mut.exon
            if isinstance(exon, str) and '_' in exon:
                first_exon, last_exon = map(int, exon.split('_'))
                exons: List[int] = list(range(first_exon, last_exon + 1))
            else:
                exons: List[int] = [exon]
            is_masked: bool = mut.is_masked
            reason: str = mut.masking_reason
            mut_class: str = mut.mutation_class
            if mut_class == DEL_EXON and self.exon_num == 1:
                continue
            if is_masked or reason != '-' or mut_class in ('COMPENSATION', MISS_EXON):
                continue
            masked_by_deletion: bool = False
            for exon in exons:
                ex_status: str = self.exon_presence[exon]
                if ex_status != 'I' and mut_class not in DEL_MISS:
                    self.mutation_list[i].is_masked = True
                    self.mutation_list[i].masking_reason = (EX_DEL_REASON if ex_status == 'D' else EX_MISS_REASON) #f'Exon is {"deleted" if ex_status == "D" else "missing"}'
                    masked_by_deletion = True
                    continue
            if masked_by_deletion:
                continue
            ## check if mutation lies it C-terminal 10% 
            ## or is masked by alternative start in the N-terminal 10%
            max_codon: int = max(self.all_ref_codons)
            for codon in codons:
                if codon > max_codon:
                    break
                is_deletion: bool = mut_class == DEL_EXON
                is_masked, reason = self._to_mask(codon, is_deletion)
                if not is_masked:
                    break
            ## mask this mutation for now; positional masking will be revised later
            self.mutation_list[i].is_masked, self.mutation_list[i].masking_reason = is_masked, reason
            muts_to_unmask.append(i)
            for exon in exons:
                if self.mutation_list[i].is_masked and mut_class in DEL_MISS:
                    self.deletion2mask[exon] = True
        ## if there are non-masked or non-unmaskable critical mutations left,
        ## unmask compensated frameshifts and positionally masked mutations
        affecting_muts: List[Mutation] = [
            x for x in self.mutation_list if
            x.mutation_class not in SAFE_UNMASKABLE_TYPES and 
            x.masking_reason not in SAFE_UNMASKABLE_REASONS and
            not (x.mutation_class == DEL_EXON and self.deletion2mask.get(x.exon, False))
        ]
        if not affecting_muts:
            return
        unsafe_found: bool = any(
            x.masking_reason != COMPENSATION_REASON for x in affecting_muts
        )
        inactivating_found: bool = any(
            not x.is_masked for x in affecting_muts
        )
        comp_num: int = len(
            [
                x for x in list(set(self.frameshift2compensation.values()))
                if x not in self.spliceai_compensations
            ]
        )
        ## unmask if any of the following cases resolve:
        ## 1) there are two or more alternative frames - DEPRECATED
        ## 2) there is alternative frame and at least one mutation outside of it - POTENTIALLY DEPRECATED
        ## 3) there is at least one unmasked inactivating mutation
        if inactivating_found:# or (unsafe_found and comp_num > 0):# or comp_num > 1:
            for i in muts_to_unmask:
                self.mutation_list[i].is_masked = False
                self.mutation_list[i].masking_reason = '-'
            for fs, comp in self.frameshift2compensation.items():
                if comp in self.spliceai_compensations:
                    continue
                self.mutation_list[fs].is_masked = False
                self.mutation_list[fs].masking_reason = '-'
                if comp not in self.comps_to_delete:
                    self.comps_to_delete.append(comp)
                    for k in self.alt_frame_masked[comp]:
                        self.mutation_list[k].is_masked = False
                        self.mutation_list[k].masking_reason = '-'
            for comp in self.comps_to_delete:
                del self.alternative_frames[comp]

    def compute_intact_percentage(self) -> None:
        """
        Infer the percentage of intact codons and their distribution
        across the sequence.
        First, each codon is classified into one of the following groups:
        * Intact (I): both the codon and respective exon are present in query,
          and the codon is not affected by mutations
        * Missing (M): codon belongs to an exon deemed as missing
        * Deleted (D): respective exon is deleted without affecting the reading frame
        * Lost (L): codon either belongs to deleted exon or is affected by
          detrimental mutation(s)
        If ignore_alternative_frame is set to True, codons confined between 
        compensated frameshifts are ignored for the purpose of these features' computation
        """
        all_codon_status: str = '' ## stores codon presence statuses
        for exon, state in self.exon_presence.items():
            first_codon, last_codon = self.exon2ref_codons[exon]
            # first_codon, last_codon = self.exon2codons[exon]
            for codon_num in range(first_codon, last_codon):
                codon_seq: str = ''.join(
                    self.all_query_codons[x] for x in sorted(self.ref_codon2triplets[codon_num])
                )
                if state == 'D': ## add 'D' to all_codon_status unless deletion is masked or something like that
                    status: str = 'D' if self.deletion2mask[exon] else 'L'
                    all_codon_status += status
                    continue
                if state == 'M': ## codon is missing, simple as
                    all_codon_status += 'M'
                    continue
                ## if alternative frame is not considered for conservation assessment, 
                ## ignore the codons confined between compensared frameshifts
                if self.ignore_alternative_frame:
                    is_in_alt_frame: bool = False
                    for i, (alt_frame_start, alt_frame_end, spliceai) in self.alternative_frames.items():
                        if spliceai:
                            continue
                        if alt_frame_start <= codon_num < alt_frame_end:
                            is_in_alt_frame = True
                            break
                    ## commented is the original intent
                    if is_in_alt_frame:
                        all_codon_status += 'A'
                        continue
                status: str = 'I' ## assign Intact status by default
                ## resulting status depends on the mutation type
                if codon_num in self.ref_codon2mutations: ## resulting status depends on the mutation type
                    for mut_num in self.ref_codon2mutations[codon_num]:
                        mutation: Mutation = self.mutation_list[mut_num]
                        mut_class: str = mutation.mutation_class
                        if mut_class in DEL_MISS: ## already handled above, proceed further
                            continue
                        if mut_class == START_MISSING: ## not considered as inactivating mutation, proceed
                            continue
                        if mut_class == COMPENSATION and mut_num not in self.comps_to_delete: ## compensation events do not count
                            continue
                        if (
                            mut_class in (FS_DEL, FS_INS) and
                            mut_num in self.frameshift2compensation and
                            not self.ignore_alternative_frame
                        ): ## frameshift is compensated, proceed
                            continue
                        if mut_class == BIG_DEL:
                            status = 'D' if all(x == '-' for x in codon_seq) else 'I'
                            continue
                        if mutation.is_masked: ## whatever the reason is, the mutation is masked; proceed
                            continue
                        ## if the algorithm made it to this point, the mutation
                        ## is deleterious and unmasked; this involves nonsenses,
                        ## uncompensated frameshifts, BIG insertions/deletions,
                        ## and splice site mutations
                        status = 'L'
                else:
                    status: str = 'D' if all(x == '-' for x in codon_seq) else 'I'
                all_codon_status += status
        ## for now, sequence properties are calculated as implemented in old CESAR_wrapper.py
        first10: int = len(all_codon_status) // 10 ## TODO: should codon_deletions_omitted used instead of all_codon_status ????
        last10: int = len(all_codon_status) - first10
        ## remove deletions
        codon_deletions_omitted: str = all_codon_status.replace('D', '')
        gene_len: int = len(all_codon_status)
        ## create two copies of status profile: with missing codons removed and
        ## with missing codons regarded as intact
        ## commented is the original intent
        missing_ignored: str = codon_deletions_omitted.replace('A','L').replace('M', '')
        missing_as_intact: str = codon_deletions_omitted.replace('A', 'L').replace('M', 'I')
        ## get the spans of non-lost codons for both profiles
        missing_ignored_spans: List[List[str]] = missing_ignored.split('L')
        missing_as_intact_spans: List[List[str]] = missing_as_intact.split('L')
        ## find the longest span for both profiles
        longest_span_omitted: int = max(len(x) for x in missing_ignored_spans)
        longest_span_as_intact: int = max(
            len(x) for x in missing_as_intact_spans
        )
        ## compute the longest span fraction in both cases
        self.longest_fraction_strict: float = longest_span_omitted / gene_len
        self.longest_fraction_relaxed: float = longest_span_as_intact / gene_len
        self.middle_is_intact: bool = 'L' not in all_codon_status[
            (first10 if self.mask_terminal_mutations else 0):last10
        ]
        self.middle_is_present: bool = 'M' not in all_codon_status[first10:last10]
        ## compute the non-missing fraction of the segment
        non_missing: int = len(all_codon_status) - all_codon_status.count('M')
        self.total_intact_fraction: float = (
            1.0 if not non_missing else 
            (all_codon_status.count('I') + all_codon_status.count('A')) / non_missing
        )

    def classify_projection(self) -> None:
        """
        Classifies the given projection based on the exon and coding sequence
        properties
        """
        self._calc_id_values()
        critical_mutations: List[Mutation] = [
            x for i, x in enumerate(self.mutation_list) if
            not x.is_masked and (
                x.mutation_class != COMPENSATION or i in self.comps_to_delete
            ) or
            x.mutation_class == MISS_EXON
        ]
        ## mutation-less projections with high identity values
        ## are automatically classified as intact
        if (
            not critical_mutations and
            self.longest_fraction_strict > STRICT_FACTION_INTACT_THRESHOLD
        ):
            architecture_non_affecting: List[Mutation] = [
                x for i, x in enumerate(self.mutation_list)
                if not self._architecture_preserving_mutation(x, i)
            ]
            if not architecture_non_affecting:
                self.loss_status = FI
            else:
                self.loss_status = I
            return

        ## projections with all exons marked as deleted are deleted themselves
        if all(x == 'D' for x in self.exon_presence.values()):
            self.loss_status = self._alt_loss_status(L)
            return

        ## projections with low sequence preservation values are classified as lost
        if self.total_intact_fraction < INTACT_CODON_LOSS_THRESHOLD:
            self.loss_status = self._alt_loss_status(L)
            return

        ## projections with low non-deleted codon preservation values are classified as lost
        if self.longest_fraction_relaxed < NON_DEL_LOSS_THRESHOLD:
            self.loss_status = self._alt_loss_status(L)
            return

        ## projections with all exons marked is missing are either missing
        ## or partially missing themselves
        if all(x == M for x in self.exon_presence.values()):
            self.loss_status = self._alt_loss_status(M)
            return
        ## alignment gap sequence clipping sometimes yields highly ambiguous sequences
        ## otherwise classified as (F/P)I; explicitly classify those as missing
        if all(x in ('X', 'x', '*') for x in self.query_aa_seq.values()):
            self.loss_status = self._alt_loss_status(M)
            return

        ## the remaining classification cases are less straightforward
        ## first, get the mutations that affect the following part of classification
        affecting_mutations: List[Mutation] = []
        for mut in critical_mutations:
            if isinstance(mut.exon, str) and '_' in mut.exon:
                first_exon, last_exon = map(int, mut.exon.split('_'))
                exons: List[int] = list(range(first_exon, last_exon + 1))
            else:
                exons: List[int] = [mut.exon]
            if any(self.exon_presence[x] == 'I' for x in exons):
                affecting_mutations.append(mut)
        ## second, get the portion of the reference sequence 
        ## comprised by missing exons
        missing_seq_len: int = sum(
            v for k,v in self.ref_exon_lengths.items() if self.exon_presence[k] == 'M'
        )
        missing_seq_perc: float = missing_seq_len / self.reference_length
        ## third, get the exons which comprise more than 40% of the ref sequence each
        bigger_exons: Set[int] = {
            k for k,v in self.ref_exon_lengths.items() if v / self.reference_length > 0.4
        }

        ## and now, classify the remaining cases
        if self.middle_is_intact:
            ## option 1: middle 80% of the sequence are not affected by mutation/loss
            ## projections satisfying this criteterion can be one of the following:
            ## I, PI, PM, UL
            if self.total_intact_fraction < MIN_INTACT_UL_FRACTION:
                ## highly diverged projection, classified as uncertain loss
                self.loss_status = UL
                return
            if not missing_seq_len:
                ## no missing exons and not mutations in the middle 80%; classified as intact
                self.loss_status = I
                return
            if self.middle_is_present:
                ## missing exons do not comprise the middle 80%; classified as intact
                ## TODO: Isn't it redundant with the previous if-statement?
                self.loss_status = I
                return
            else:
                ## there are missing exons covering the middle 80% of the sequence
                if missing_seq_perc < MAX_MISSING_PM_THRESHOLD:
                    ## less than 50% of the sequence is missing; classified as partially intact
                    self.loss_status = PI
                    return
                ## otherwise, the projection is missing
                self.loss_status = self._alt_loss_status(M)
                return
        else:
            ## option 2: the projection contains at least one affecting mutation
            ## in its middle 80% of the sequence
            ## this means that sequence is either uncertainly lost or lost
            if self.exon_num == 1:
                ## single exon projections require at least two mutations and
                ## longest non-missing intact fraction of less than 60%
                ## to be a confirmed loss
                if len(affecting_mutations) > 1 and self.longest_fraction_strict < STRICT_FACTION_INTACT_THRESHOLD:
                    self.loss_status = self._alt_loss_status(L)
                else:
                    self.loss_status = UL
                return
            else:
                if self.longest_fraction_strict < STRICT_FACTION_INTACT_THRESHOLD:
                    ## the identity value is low enough; suspecting proper loss,
                    ## might be also an uncertain loss
                    ## compute the number of exons affected by mutations
                    mut_affected_exons: Set[int] = set()
                    for aff_mut in affecting_mutations:
                        if isinstance(aff_mut.exon, int):
                            mut_affected_exons.add(aff_mut.exon)
                        else:
                            for aff_ex in aff_mut.exon.split('_'):
                                mut_affected_exons.add(int(aff_ex))
                    del_affected_exons: Set[int] = {
                        k for k, v in self.exon_presence.items() if v == 'D'
                        and not self.deletion2mask[k]
                    }
                    num_affected_exons: int = len(
                        mut_affected_exons.union(del_affected_exons)
                    )
                    max_affected_exon_num: int = get_affected_exon_threshold(
                        sum(x != 'M' for x in self.exon_presence.values())
                    )
                    ## if affected exon number threshold is exceeded,
                    ## report the projection as missing
                    if num_affected_exons >= max_affected_exon_num:
                        self.loss_status = self._alt_loss_status(L)
                        return
                    ## then, check if any of 'big' (>=40% of CDS each) exons are
                    ## deleted or affected by 2+ mutations; this would mark the projection
                    ## as lost
                    for be in bigger_exons:
                        if self.exon_presence[be] == 'D' and not self.deletion2mask[be]:
                            self.loss_status = self._alt_loss_status(L)
                            return
                        if sum(x.exon == be for x in affecting_mutations) > 1:
                            self.loss_status = self._alt_loss_status(L)
                            return
                    ## otherwise, classify the projection as uncertain loss
                    self.loss_status = UL
                    return
                else:
                    ## the intact codon portion is high enough; can be either
                    ## partially missing or uncertain loss
                    ## nothing to say for sure except for that it looks fishy;
                    ## classified as uncertain loss
                    self.loss_status = UL
                    return

    def update_mutation_list(self) -> None:
        """Remove revoked compensation events from the mutation list"""
        new_mut_list: List[int] = []
        for i, mut in enumerate(self.mutation_list):
            ## remove frameshifts participating in certain compensations
            if mut.mutation_class in (FS_INS, FS_DEL) and i in self.frameshift2compensation:
                comp: int = self.frameshift2compensation[i]
                ## do not add deprecated compensations
                if comp in self.spliceai_compensations:
                    continue
            ## as well as compensation entries themselves
            if mut.mutation_class == COMPENSATION:
                if i in self.comps_to_delete:
                    mut.masking_reason = OBSOLETE_COMPENSATION
                if i in self.spliceai_compensations:
                    continue
            ## as well as mutations masked by these compensations
            if mut.mutation_class == STOP:
                remove_from_muts: bool = False
                for comp in self.alt_frame_masked:
                    if i in self.alt_frame_masked[comp]:
                        if comp in self.comps_to_delete:
                            remove_from_muts = True
                        break
                if remove_from_muts:
                    continue
            new_mut_list.append(mut)
        self.mutation_list = new_mut_list

    def plot_svg(self) -> None:
        """Creates an SVG plot for the projection"""
        self.svg = ProjectionPlotter(
            self.name, self.ref_exon_lengths, self.mutation_list, self.exon2chain
        ).plot_svg().replace('\n', '')

    def _abs_coord(self, base: int, exon: int = None, portion: int = None) -> int:
        """
        Given the relative base coordinate in the concatenated alignment,
        return its absolute coordinate in the query genome
        """
        if base is None:
            return None
        if portion is None:
            if exon is None:
                p: int = self._base2portion(base)
            else:
                p: int = self.exon2portion[exon]
        else:
            p: int = portion
        strand: bool = self.portion2strand[p]
        rel_start, rel_stop = self.portion2coords[p]
        gap_num: int = self.query[rel_start:base].count('-')
        abs_start, abs_stop = self.portion2abs_coords[p][1:]
        adj_base: int = base - rel_start
        return (abs_start + adj_base - gap_num) if strand else (abs_stop - adj_base + gap_num)

    def _rel_coord(
        self, base: int, counter: int = 0, exon: int = None, portion: int = None
    ) -> int: ## TODO: Needs testing
        """
        Given the absolute coordinate, reduce it to the relative coordinate
        in the alignment
        """
        if portion is None:
            if exon is None:
                p: int = self._absbase2portion(base)
            else:
                p: int = self.exon2portion[exon]
        else:
            p: int = portion
        strand: bool = self.portion2strand[p]
        ## get the absolute coordinates of the portion
        abs_start, abs_stop = sorted(self.portion2abs_coords[p][1:])
        ## get the respective relative coordinates
        rel_start, rel_stop = sorted(self.portion2coords[p])
        ## locate the raw coordinate
        unadj_coord: int = (base - abs_start if strand else abs_stop - base) + rel_start
        if unadj_coord == 0:
            return unadj_coord
        adj_coord: int = rel_start
        for i, n in enumerate(self.query[rel_start:rel_stop], rel_start):
            if adj_coord == unadj_coord:
                break
            if n.isalpha():
                adj_coord += 1
        gap_num: int = self.query[rel_start:i].count('-')
        return i

    def _closest_non_gap(self, base: int, exon: int) -> int:
        """
        Given the relative coordinate, returns the last previous
        non-gap position in the query
        """
        portion: int = self.exon2portion[exon]
        p_start, _ = self.portion2coords[portion]
        for i in range(base, p_start - 1, -1):
            n: str = self.query[i]
            if n != '-':
                return i
        return base


    def _base2portion(self, base: int) -> int:
        """
        Given the relative base coordinate in the concatenated alignment,
        return the CESAR output piece it referes to
        """
        is_gap: str = self.query[base] == '-'
        for portion, coords in self.portion2coords.items():
            if (base < coords[0] or base >= coords[1]) and is_gap:
                continue
            if (base < coords[0] or base >= coords[1]) and not is_gap:
                continue
            return portion
        raise Exception(
            f'Alignment position {base} does not refere to any CESAR alignment portion'
        )

    def _absbase2portion(self, base: int) -> int:
        """
        Given the absolute base coordinate in the concatenated alignment,
        return the CESAR output piece it referes to
        """
        for portion, coords in self.portion2abs_coords.items():
            if base < coords[1] or base > coords[2]:
                continue
            return portion
        raise Exception(
            f'Alignment position {base} does not refere to any CESAR alignment portion'
        )

    def _base2codon(self, i: int, exon: int, ref_codon: bool = False) -> int:
        """Returns codon number the base belongs to"""
        first, last = self.exon2codons[exon]
        for c in range(first, min(last + 1, max(self.triplet_coordinates) + 1)):
            coords: List[int] = self.triplet_coordinates[c]
            start: int = min(coords)
            end: int = max(coords)
            if start - 1     <= i <= end + 1:
                if ref_codon:
                    return self.triplet2ref_codon[c]
                return c


    def _exon2strand(self, exon: int) -> bool:
        """For a given exon number, returns the strand it was aligned to"""
        return self.portion2strand[self.exon2portion[exon]]

    def _safe_step(self, base: int, exon: int, step: int) -> int:
        """
        For a given relative coordinate and a step value, returns the coordinate
        after step if it lies within the same alignment portion, otherwise
        returns the relative coordinate of the portion's closest boundary
        """
        p: int = self.exon2portion[exon]
        _start, _stop = self.portion2coords[p]
        return max(min(base + step, _stop), _start)

    def _codon2exon(self, codon: int) -> int: ## NEEDS UPDATE
        """For a given codon number, return the (reference) exon it belongs to"""
        if not self.exon2codons:
            ## special condition when procedure is run for a yet unterminated first exon
            yield 1
        for exon, coords in self.exon2codons.items():
            first, last = coords if len(coords) == 2 else coords * 2
            is_split: bool = last in self.split_codon_struct
            if codon >= first and (codon <= last if is_split else codon < last):
                yield exon
            elif exon == max(self.exon2codons) and codon >= first:
                yield exon

    def _codon2chrom(self, codon: int) -> str: ## NEEDS UPDATE
        """
        For a given codon number, returns the chromosome
        it is located on in the query
        """
        exon: int = next(self._codon2exon(codon))
        chrom: str = self.exon2chrom[exon]
        return chrom

    def _codon2coords(self, codon: int) -> Tuple[str, int]:
        """
        Returns reference codon coordinates in query in the (chrom, start, stop) form
        """
        chrom: str = self._codon2chrom(codon)
        coords: List[int] = self.codon_coordinates.get(codon, [])
        if len(coords) > 1:
            return (chrom, min(coords), max(coords))
        else: 
            ## a workaround for exon shorter than one full codon
            ## _codon2exon returns all the exons the codon corresponds to
            exons: List[int] = list(self._codon2exon(codon))
            if codon in self.split_codon_struct:
                for ex in sorted(self.split_codon_struct[codon], reverse=True):
                    if ex in self.abs_exon_coords:
                        exon: int = ex
                        break
            else:
                exon: int = exons.pop()
            coords = self.abs_exon_coords[exon].tuple()
            return (chrom, min(coords), max(coords))
    
    def _exon_was_corrected(self, exon: int) -> bool:
        """Returns whether the exon was SpliceAI-corrected"""
        if self.correction_mode <= 1:
            return False
        start_corrected: int = not self.cesar_acc_support[exon] and exon > 1
        end_corrected: int = not self.cesar_donor_support[exon] and (
            exon < self.exon_num
        )
        return start_corrected or end_corrected

    def _in_modified_seq(self, exon: int, triplet: int) -> bool:
        """Checks if a codon corresponds to a sequence shifted with SpliceAI"""
        first_codon, last_codon = self.exon2codons[exon]
        is_first_codon: bool = triplet == first_codon
        is_last_codon: bool = triplet == last_codon - 1 or (
            triplet in self.split_codon_struct.keys() and exon in self.split_codon_struct[triplet]
        )
        orig_start, orig_end = sorted(map(
            lambda x: self._abs_coord(x, exon=exon), self.cesar_rel_coords[exon].tuple()
        ))
        upd_start, upd_end = self.abs_exon_coords[exon].tuple()
        triplet_bases: List[int] = self.triplet_coordinates[triplet]
        strand: bool = self._exon2strand(exon)
        if orig_start > upd_start:
            if any(upd_start - 3 <= x <= orig_start + 3 for x in triplet_bases):
                return True
            if is_first_codon and strand or is_last_codon and not strand:
                return True
        if orig_start < upd_start:
            if any(orig_start - 3 <= x <= upd_start + 3 for x in triplet_bases):
                return True
            if is_first_codon and strand or is_last_codon and not strand:
                return True
        if orig_end > upd_end:
            if any(upd_end - 3 <= x <= orig_end + 3 for x in triplet_bases):
                return True
            if is_last_codon and strand or is_first_codon and not strand:
                return True
        if orig_end < upd_end:
            if any(orig_end - 3 <= x <= upd_end + 3 for x in triplet_bases):
                return True
            if is_last_codon and strand or is_first_codon and not strand:
                return True
        return False

    def _mutations_in_same_query_codon(self, mut1: Mutation, mut2: Mutation) -> bool:
        """
        For rare cases of convoluted SpliceAI-induced pseudoframeshifts, 
        checks whether the two triplets are in fact part of one codon
        """
        if isinstance(mut1.exon, str):
            first_exon: int = max(map(int, mut1.exon.split('_')))
        else:
            first_exon: int = mut1.exon
        if isinstance(mut2.exon, str):
            last_exon: int = min(map(int, mut2.exon.split('_')))
        else:
            last_exon: int = mut2.exon
        if last_exon != first_exon + 1:
            return False
        full_query_seq: str = ''
        first_codon: int = mut1.codon
        last_codon: int = mut2.codon
        for codon in range(first_codon, last_codon + 1):
            if codon in self.query_codons_to_mask:
                query_codon_seq: str = self.query_codons_to_mask[codon]
            else:
                query_codon_seq: str = self.all_query_codons[codon]
            full_query_seq += query_codon_seq.replace('-', '')
        return len(full_query_seq) == 3

    def _comp_moved_to_split_codon(self, mut1: Mutation, mut2: Mutation) -> bool:
        """
        Checks for an even rarer event of SpliceAI-corrected frameshift
        migrating to the next exon
        """
        if isinstance(mut1.exon, str):
            first_exon: int = max(map(int, mut1.exon.split('_')))
        else:
            first_exon: int = mut1.exon
        if isinstance(mut2.exon, str):
            last_exon: int = min(map(int, mut2.exon.split('_')))
        else:
            last_exon: int = mut2.exon
        if last_exon != first_exon + 1:
            return False
        last_codon: int = mut2.ref_codon
        if self.exon2ref_codons[last_exon][0] != last_codon:
            return False
        for triplet in self.ref_codon2triplets[last_codon]:
            if triplet in self.split_codon_struct:
                triplet_split_struct: Dict[int, int] = self.split_codon_struct[triplet]
                if first_exon in triplet_split_struct and last_exon in triplet_split_struct:
                    return True
        return False

    def _comp_moved_to_next_exon(self, mut1: Mutation, mut2: Mutation) -> bool:
        """
        The final line of defense: check whether the SpliceAI-induced frame-restoring
        FS mutation 'leaked' further into the next exon
        """
        if isinstance(mut1.exon, str):
            first_exon: int = max(map(int, mut1.exon.split('_')))
        else:
            first_exon: int = mut1.exon
        if isinstance(mut2.exon, str):
            last_exon: int = min(map(int, mut2.exon.split('_')))
        else:
            last_exon: int = mut2.exon
        if last_exon != first_exon + 1:
            return False
        exon_start, exon_end = sorted(self.abs_exon_coords[mut2.exon].tuple())
        mut2_start, mut2_end = sorted((mut2.start, mut2.stop))
        prev_was_corrected: bool = self._exon_was_corrected(mut1.exon)
        if not prev_was_corrected:
            return False
        return (mut2_start <= exon_start <= mut2_end) or (mut2_start <= exon_start <= mut2_end)

    def _min_annot_codon(self) -> int:
        """Returns the first non-gap codon in the query"""
        for i, coords in self.triplet_coordinates.items():
            if len(coords):
                return i

    def _next_complete_codon(self, num: int, prev: bool = False) -> int:
        """
        Returns the number of the first frame-preserving codon 
        previous/following to a given one
        """
        if num == 1 and prev:
            return num
        if num == max(self.triplet2ref_codon) and not prev:
            return num
        search_range: Iterable[int] = (
            range(num - 1, 0, -1) if prev else
            range(num + 1, max(self.triplet2ref_codon) + 1)
        )
        max_triplet: int = max(self.triplet2ref_codon)
        for i in search_range:
            if i not in self.query_codons_to_mask and i not in self.ref_codons_to_mask:
                return max(i, 1) if prev else min(i, max_triplet)
        return 1 if prev else max_triplet

    def _merged_exon_streak(self, exon: int) -> Tuple[int]: ## TODO: Tested only on minimal examples; test further
        """Returns the iterable of exons merged in query for the given exons"""
        streak: List[int] = [1]
        for i, is_lost in self.intron_deleted.items():
            if not is_lost:
                if exon in streak:
                    return tuple(streak)
                elif i == exon:
                    return (i,)
                streak = []
            # streak.append(i)
            streak.append(i+1)
        if exon in streak:
            return tuple(streak)

    def _get_intron_phase(
        self, intron: int, ref: bool = True
    ) -> int:
        """Returns intron phase for the given intron number in the query"""
        prev_intron_phase: int = 0
        total_seq: str = ''
        merged_exons: List[int] = [intron]
        for x in range(intron - 1, 0, -1):
            _prev_phase: int = self.intron2phase.get(x, 0)
            if _prev_phase < 0:
                merged_exons.append(x)
                continue
            prev_intron_phase = _prev_phase
            break
        for ex in merged_exons[::-1]:
            prev_exon_start, prev_exon_end = self.rel_exon_coords[ex].tuple()
            total_seq += strip_noncoding(
                getattr(self, 'reference' if ref else 'query')[prev_exon_start:prev_exon_end],
                uppercase_only=not ref
            )
        last_exon_length: int = len(total_seq)
        return (last_exon_length - 3 + prev_intron_phase) % 3

    def _intron_gap_in_seq(self, start: int, end: int) -> bool:
        """Checks if a given sequence contains at least 30 consecutive gaps"""
        return GAP_CODON * 10 in self.reference[start:end]

    def _intron_length(
        self, donor: int, acceptor: int, donor_exon: int, acc_exon: int
    ) -> Union[int, None]:
        """
        Calculates the intron length; 
        returns None if donor and acceptor are located on different contigs
        """
        donor_chrom: str = self.portion2abs_coords[self.exon2portion[donor_exon]][0]
        acc_chrom: str = self.portion2abs_coords[self.exon2portion[acc_exon]][0]
        if donor_chrom != acc_chrom:
            return None
        abs_donor: int = self._abs_coord(donor, exon=donor_exon)
        abs_acc: int = self._abs_coord(acceptor, exon=acc_exon)
        start, end, = sorted((abs_donor, abs_acc))
        return end - start

    def _mut_number_in_seq(
        self,
        subexons: List[Tuple[int, int, float, float]],
        exon: int
    ) -> Tuple[int, int]:
        """
        Computes the number of nonsense
        """
        exon_streak: Tuple[int] = self._merged_exon_streak(exon)
        first_ex: int = exon_streak[0]
        last_ex: int = exon_streak[-1]
        seq_start: int = self.rel_exon_coords[first_ex].start
        seq_end: int = self.rel_exon_coords[last_ex].stop
        prev_phase: int = self.intron2phase.get(first_ex - 1, 0)
        next_phase: int = self.intron2phase.get(last_ex, 0)
        # print(f'{exon_streak=}, {prev_phase=}, {next_phase=}, {self.intron2phase=}')
        if prev_phase < 0 or next_phase < 0:
            raise ValueError(
                'Trying to access the phase value for a deleted intron'
            )
        seq_start += (3 - prev_phase) % 3
        seq_end -= next_phase

        orig_ref_codon: str = ''
        corr_ref_codon: str = ''
        orig_query_codon: str = ''
        corr_query_codon: str = ''
        gained_mut_num: int = 0
        lost_mut_num: int = 0
        for i in range(seq_start, seq_end):
            orig_ref_codon += self.reference[i]
            corr_ref_codon += self.reference[i] if any(x[0] <= i < x[1] for x in subexons) else ''
            orig_query_codon += self.query[i]
            corr_query_codon += self.query[i] if any(x[0] <= i < x[1] for x in subexons) else ''
            orig_ref_complete: bool = is_complete_codon(orig_ref_codon)
            corr_ref_complete: bool = is_complete_codon(corr_ref_codon)
            orig_query_complete: bool = is_complete_codon(orig_query_codon)
            corr_query_complete: bool = is_complete_codon(corr_query_codon)
            orig_ref_fs: bool = frameshift_codon(orig_ref_codon)
            corr_ref_fs: bool = frameshift_codon(corr_ref_codon)
            orig_query_fs: bool = frameshift_codon(orig_query_codon)
            corr_query_fs: bool = frameshift_codon(corr_query_codon)
            if orig_ref_complete or i == seq_end - 1:
                if orig_ref_fs:
                    if orig_ref_codon != corr_ref_codon:
                        lost_mut_num += 1
            if corr_ref_complete or i == seq_end - 1:
                if corr_ref_fs:
                    if orig_ref_codon != corr_ref_codon:
                        gained_mut_num += 1
                corr_ref_codon = ''
            if orig_ref_complete:
                orig_ref_codon = ''
            if orig_query_complete or i == seq_end - 1:
                if orig_query_fs or strip_noncoding(orig_query_codon).upper() in STOPS:
                    if orig_query_codon != corr_query_codon:
                        lost_mut_num += 1
            if corr_query_complete or i == seq_end - 1:
                if corr_query_fs or strip_noncoding(corr_query_codon).upper() in STOPS:
                    if orig_query_codon != corr_query_codon:
                        gained_mut_num += 1
                corr_query_codon = ''
            if orig_query_complete:
                orig_query_codon = ''
        return gained_mut_num, lost_mut_num

    def _get_frameshift(self, start: int, end: int):
        """
        Returns the difference in frame phase between reference and query sequences 
        confined within two coordinates
        """
        ## strip gained intron sequenes
        ref_seq: str = ''.join(
            [
                self.reference[i] for i in range(start,end) if 
                self.query[i].isupper() or self.query[i] == '-'
            ]
        )
        query_seq: str = ''.join(
            [
                self.query[i] for i in range(start,end) if 
                self.query[i].isupper() or self.query[i] == '-'
            ]
        )
        ## count overall deletion number in both sequences
        ref_del: str = ref_seq.count('-') % 3
        query_del: str = query_seq.count('-') % 3
        return ref_del - query_del

    def _site_is_inframe(
        self, cesar: int, site: int, offset: Optional[int] = 0, donor: Optional[bool]=True
    ) -> bool:
        """
        For a given CESAR-predicted exon start/stop and a SpliceAI-predicted
        acceptor/donor site, return whether they lie within the same reading frame
        Extension/contraction are defined as follows:
        ## donor and corrected site is to the right from alignment -> extension
        ## donor and corrected site is to the left from alignment boundary -> contraction
        ## acceptor and corrected site is to the right from alignment boundary -> contraction
        ## acceptor and corrected site is to the left from alignment boundary -> extension
        """
        ## get the distance between two coordinates in coding symbols
        distance: str = len(
            strip_noncoding(self.query[slice(*sorted((cesar, site)))], uppercase_only=False)
        )
        ## define which side the frame is corrected from
        side_sign: int = 1 if donor else -1
        ## define which direction the sequence is shifted to
        shift_sign: int = 1 if site > cesar else -1
        ## define difference modifier (i.e., whether the sequence was added or subtracted)
        dist_sign: int = side_sign * shift_sign
        return (dist_sign * distance + offset) % 3 == 0

    def _update_spliceai_predictions(self, p: RawCesarOutput, num: int) -> None:
        """
        Recalculates the raw SpliceAI predictions for easier use during the
        self.process() run
        """
        for exon in p.exons:
            space_start, space_stop = p.exon_search_spaces[exon]
            if space_start is None:
                self.spliceai_sites[exon]['donor'] = {
                    self._rel_coord(x, portion=num):y for x, y in p.spliceai_sites['donor'].items()
                }
                self.spliceai_sites[exon]['acceptor'] = {
                    self._rel_coord(x, portion=num):y for x, y in p.spliceai_sites['acceptor'].items()
                }
                continue
            aln_start, aln_stop = self.abs_exon_coords[exon].tuple()
            start: int = min(space_start, aln_start)
            stop: int = max(space_stop, aln_stop)
            for key in p.spliceai_sites:
                for coord, prob in p.spliceai_sites[key].items():
                    if not start - 1 <= coord <= stop + 1:
                        continue
                    rel_coord: int = self._rel_coord(coord, portion=num)
                    self.spliceai_sites[exon][key][rel_coord] = prob

    def _check_splice_site(
        self,
        exon: int,
        base: int,
        acc: bool = True,
        check_only: bool = False
    ) -> bool: ## TODO: Backtracking at erroneous introns results in overlaps -> fix it
        """
        For a given splice site at a given position, check if the site is intact;
        returns a boolean value indicating whether splice site is intact and does
        not need correction
        Arguments are:
        ::exon:: is the number of exon the splice site corresponds to;
        ::base:: corresponds to coordinate of the respective marginal
            nucleotide in the query sequence in the alignment;
        ::acc:: designates that the splice site in question is acceptor site; the
            opposite means the donor site;
        ::check_only:: prevents adding a splice site mutation to the global pool
        """
        ## return True for non-spliced sites
        if exon == 1 and acc:
            if not check_only:
                self.splice_site_nucs[exon]['A'] = 'NA'
            return True
        if exon == self.exon_num and not acc:
            if not check_only:
                self.splice_site_nucs[exon]['D'] = 'NA'
            return True
        if exon in self.unaligned_exons:
            if not check_only:
                self.splice_site_nucs[exon]['A' if acc else 'D'] = '??'
            return True
        ## assess whether the dinucleotide at respective site is intact
        dinucleotide: str = (
            self.query[self._safe_step(base, exon, -2):base] if acc else
            self.query[base:self._safe_step(base, exon, 2)]
        ).lower()
        is_intact: bool = dinucleotide in (
            LEFT_SPLICE_CORR if acc else RIGHT_SPLICE_CORR
        )
        if is_intact:
            if not check_only:
                self.splice_site_nucs[exon][('A' if acc else 'D')] = dinucleotide
            return True
        ## infer special conditions under which affected splice site is treated as valid
        is_intact = False
        to_mask: bool = False
        reason: str = '-'
        donor_exon: int = exon if not acc else exon - 1
        acc_exon: int = exon if acc else exon + 1
        side: str = 'acceptor' if acc else 'donor'
        site_is_u12: bool = self.u12[exon][side][0] == 'U12'
        site_is_non_cano: bool = (
            not site_is_u12 and self.u12[exon][side][1].lower() not in (
                LEFT_SPLICE_CORR if acc else RIGHT_SPLICE_CORR
            )
        )
        ref_site: str = self.u12[exon][side][1]
        donor_coord: int = self.rel_exon_coords[exon-1].stop if acc else base
        acc_coord: int = base if acc else self.rel_exon_coords[exon+1].start
        if acc:
            prev_site_mutated: Union[int, None] = next(
                (
                    x for x, mut in enumerate(self.mutation_list)if
                    mut.mutation_class == SSM_D and mut.exon == exon - 1
                ), None
            )
        else:
            prev_site_mutated = None
        if 'n' in dinucleotide:
            to_mask, reason = True, 'Ambiguous symbol'
        elif len(dinucleotide) < 2:
            is_intact = True
            dinucleotide = (
                (dinucleotide + '?' * (2 - len(dinucleotide))) if not acc else
                ('?' * (2 - len(dinucleotide)) + dinucleotide)
            )
            to_mask, reason = True, 'Insufficient data'
        elif (
            self._check_deleted_intron(donor_coord, acc_coord, donor_exon, acc_exon)[0] and 
            prev_site_mutated is not None
        ):
            to_mask, reason = True, 'Intron deletion'
            self.mutation_list[prev_site_mutated].is_masked = to_mask
            self.mutation_list[prev_site_mutated].masking_reason = reason
            dinucleotide: str = '--'
            old_prev_name: str = self.mutation_list[prev_site_mutated].description.split('->')[0]
            new_prev_name: str = f'{old_prev_name}->{dinucleotide}'
            self.mutation_list[prev_site_mutated].description = new_prev_name
        elif site_is_u12:
            to_mask, reason = True, U12_REASON
        elif site_is_non_cano:
            to_mask, reason = True, NON_CANON_U2_REASON
        if check_only:
            return is_intact
        ## after this point, create a mutation instance
        self.splice_site_nucs[exon][('A' if acc else 'D')] = dinucleotide
        mut_name: str = f'{ref_site}->{dinucleotide}'
        counter: int = self.donor_mutation_counter if not acc else self.acc_mutation_counter
        mut_id: str = f'SSM{"A" if acc else "D"}_{counter}'
        chain: int = self.exon2chain[exon]
        codon: int = min(
                max(self.exon2codons[exon][int(not acc)] - int(not acc), 1),
                max(self.triplet2ref_codon)
        )
        ref_codon: int = self.triplet2ref_codon[codon]
        splice_site_mut: Mutation = Mutation(
            self.transcript,
            chain,
            exon,
            codon,
            ref_codon,
            *self._codon2coords(codon),
            SSM_A if acc else SSM_D,
            mut_name,
            to_mask,
            reason,
            mut_id
        )
        mut_num: int = len(self.mutation_list)
        self.mutation_list.append(splice_site_mut)
        self.codon2mutations[codon].append(mut_num)
        self.ref_codon2mutations[ref_codon].append(mut_num)
        if acc:
            self.acc_mutation_counter += 1
        else:
            self.donor_mutation_counter += 1


    def _check_deleted_intron(
        self, donor: int, acc: int, donor_exon: int, acc_exon: int
    ) -> Tuple[bool, bool]:
        """
        Checks whether the intron was deleted in the query, Returns a tuple of
        two boolean values:
        * output[0] defines whether the intron was lost
        * output[1] defines whether a short intron was erroneously introduced by
          CESAR, and acceptor site needs backtracking
        """
        if self.reference[donor+1] == '>' and self.reference[acc-1] == '>':
            return (True, False)
        donor_portion: int = self.exon2portion[donor_exon]
        acc_portion: int = self.exon2portion[acc_exon]
        if donor_portion != acc_portion:
            return (False, False)
        elif acc - donor < MAX_RETAINED_INTRON_LEN and self.correct_ultrashort_introns:
            if self._nonsense_free_correction(donor, acc):
                return (True, True)
        return (False, False)

    def _nonsense_free_correction(self, upstream: int, downstream: int) -> bool:
        """Returns whether the spice site extension contains no inframe stop codons"""
        ## correction resulted in exon sequence contraction -> no codon gain
        if downstream <= upstream:
            return True
        for codon in parts(self.query[upstream:downstream], 3):
            if codon.upper() in STOPS:
                return False
        return True

    def _alternative_split_contains_stop(
        self, donor: int, acc: int, phase: int, central: int = None
    ) -> bool:
        """
        Check whether the newly selected splice sites create a nonsense mutation
        """
        if not phase:
            return False
        pre_split_counter: int = 0
        for i in range(donor - 1, -1, -1):
            pre_split_counter += int(is_symbol(self.reference[i]))
            if pre_split_counter == phase:
                break
        post_split_counter: str = 0
        for j in range(acc, len(self.query)):
            post_split_counter += int(is_symbol(self.reference[j]))
            if post_split_counter == 3 - phase - int(central is not None):
                break
        new_split_codon: str = (self.query[i:donor] + self.query[acc:j]).replace('-', '')
        for codon in parts(new_split_codon, 3):
            if codon.upper() in STOPS:
                return True
        return False

    def _zero_length_exon_introduced(self, intron: int, don_: int, acc_: int) -> bool:
        """
        Checks if adjusted intron boundaries render any of the flanking exons to
        length zero
        """
        up_exon: int = intron
        down_exon: int = intron + 1
        # up_start: int = self.rel_exon_coords[up_exon].start
        down_stop: int = self.rel_exon_coords[down_exon].stop
        return (don_ - up_exon) <= 0 or (down_stop - acc_) <= 0

    def _rescue_alternative_start(self) -> int:
        """Searches for an alternative upstream start codon for the first exon"""
        first_codon_start: int = self.rel_exon_coords[1]
        if self.query[first_codon_start:first_codon_start+3].upper() == START:
            return
        if self.reference[first_codon_start:first_codon_start+3].upper() != START:
            return
        for i in range(first_codon_start - 3, 0, -3):
            codon: str = self.query[i:i+3]
            if codon.upper() == START:
                self.rel_exon_coords[1].start = i
                self.abs_exon_coords[1].start = self._abs_coord(i, exon=1)
                return

    def _rescue_alternative_stop(self) -> int:
        """Searches for an alternative downstream stop codon for the last exon"""
        last_codon_stop: int = self.rel_exon_coords[self.exon_num].stop
        last_query_triplet: str = self.query[last_codon_stop-3:last_codon_stop]
        if last_query_triplet in STOPS:
            return
        if last_query_triplet.upper() == NNN_CODON:
            return
        if self.reference[last_codon_stop-3:last_codon_stop] not in STOPS:
            return
        for i in range(last_codon_stop, len(self.query), 3):
            codon: str = self.query[i:i+3]
            if codon.upper() in STOPS:
                self.rel_exon_coords[self.exon_num].stop = i + 3
                self.clipped_rel_coords[self.exon_num].stop = i + 3
                self.abs_exon_coords[self.exon_num].stop = self._abs_coord(i + 3, exon=self.exon_num)
                self.stop_updated = True
                return

    def _check_for_start_loss(self) -> None:
        """
        Checks if the codon differs from the start codon; if so, adds a Mutation entry
        """
        if 1 in self.unaligned_exons:
            return
        query_codon: str = self.all_query_codons[1]
        if query_codon == START:
            return
        ref_codon: str = self.all_ref_codons[1]
        if ref_codon != START:
            return
        start_loss: Mutation = Mutation(
            self.transcript,
            self.exon2chain[1],
            1, ## exon number
            1, ## codon number
            1, ## reference codon number
            *self._codon2coords(1),
            START_MISSING,
            f'{ref_codon}->{query_codon}',
            True, ## start loss is not considered to be inactivating and is masked by default
            'Missing start masked',
            'START_1'
        )
        mut_num: int = len(self.mutation_list)
        self.mutation_list.append(start_loss)
        self.codon2mutations[1].append(mut_num)
        self.ref_codon2mutations[1].append(mut_num)

    def _check_for_frameshift(
        self, codon: int, ref_codon: str = None, query_codon: str = None, exon: int = None
    ) -> None:
        """
        Checks whether the provided codon contains a frameshifting mutation
        """
        if ref_codon is None:
            ref_codon: str = self.ref_codons_to_mask.get(
                codon, self.all_ref_codons[codon]
            )
        if query_codon is None:
            query_codon: str = self.query_codons_to_mask.get(
                codon, self.all_query_codons[codon]
            )
        ref_gap_num: int = ref_codon.count('-')
        query_gap_num: int = query_codon.count('-')
        delta: int = ref_gap_num - query_gap_num
        if not delta % 3:
            return
        if exon is None:
            if codon in self.split_codon_struct:
                split_structure: Dict[int, int] = self.split_codon_struct[codon]
                sub_start: int = 0
                for ex, portion in split_structure.items():
                    ref_sub: str = ref_codon[sub_start:sub_start + portion]
                    query_sub: str = query_codon[sub_start:sub_start + portion]
                    sub_ref_gap_num: int = ref_sub.count('-')
                    sub_query_gap_num: int = query_sub.count('-')
                    if sub_ref_gap_num != sub_query_gap_num:
                        exon = ex
                        break
                    sub_start += portion
                if exon is None:
                    ## if mutations affects both exons sharing the split codon,
                    ## arbitrarily assign it to the downstream exon
                    exon = ex ## TODO: Should it be assigned to both/all three exons instead?
            else:
                exon: int = next(self._codon2exon(codon))
        if exon in self.unaligned_exons:
            return
        ref_codon_num: int = self.triplet2ref_codon[codon]
        description: str = str(delta)
        mutation_class: str = FS_INS if delta > 0 else FS_DEL
        mutation_id: str = f'FS_{self.frameshift_counter}'
        frameshift: Mutation = Mutation(
            self.transcript,
            self.exon2chain[exon],
            exon,
            codon,
            ref_codon_num,
            *self._codon2coords(codon),
            mutation_class,
            description,
            False, ## Point mutations are added unmasked by default
            '-', ## with no masking reasons provided
            mutation_id
        )
        mut_num: int = len(self.mutation_list)
        self.mutation_list.append(frameshift)
        self.frameshift_counter += 1
        self.codon2mutations[codon].append(mut_num)
        self.ref_codon2mutations[ref_codon_num].append(mut_num)

    def _check_for_nonsense_mutation(
        self,
        codon: int,
        ref_codon: str = None,
        query_codon: str = None,
        exon: int = None
    ) -> None:
        """Checks if the provided codon is a premature stop codon"""
        if ref_codon is None:
            ref_codon: str = self.ref_codons_to_mask.get(
                codon, self.all_ref_codons[codon]
            )
        if query_codon is None:
            query_codon: str = self.query_codons_to_mask.get(
                codon, self.all_query_codons[codon]
            )
        if query_codon not in STOPS:
            return
        if ref_codon in STOPS:
            return
        if ref_codon == NNN_CODON:
            self.selenocysteine_codons.append((exon, codon, query_codon))
            return
        if self.stop_updated:
            return
        if exon is None:
            if codon in self.split_codon_struct:
                split_structure: Dict[int, int] = self.split_codon_struct[codon]
                sub_start: int = 0
                for ex, portion in split_structure.items():
                    ref_sub: str = ref_codon[sub_start:sub_start + portion]
                    query_sub: str = query_codon[sub_start:sub_start + portion]
                    ref_sub: str = ref_codon[0:sub_start] + query_sub + ref_codon[sub_start+portion:]
                    if ref_sub in STOPS:
                        exon = ex
                        break
                    sub_start += portion
                if exon is None:
                    ## if mutations affects both exons sharing the split codon,
                    ## arbitrarily assign it to the downstream exon
                    exon = ex ## TODO: Should it be assigned to both/all three exons instead?
            else:
                exon: int = next(self._codon2exon(codon))
        if exon in self.unaligned_exons:
            return
        ref_num_codon: int = self.triplet2ref_codon[codon]
        if ref_codon == NNN_CODON:
            is_masked, reason = True, 'Reference codon masked'
        else:
            is_masked, reason = False, '-'
        description: str = f'{ref_codon}->{query_codon}'
        mut_id: str = f'STOP_{self.nonsense_counter}'
        nonsense: Mutation = Mutation(
            self.transcript,
            self.exon2chain[exon],
            exon,
            codon,
            ref_num_codon,
            *self._codon2coords(codon),
            STOP,
            description,
            is_masked,
            reason,
            mut_id
        )
        mut_num: int = len(self.mutation_list)
        self.mutation_list.append(nonsense)
        self.nonsense_counter += 1
        self.codon2mutations[codon].append(mut_num)
        self.ref_codon2mutations[ref_num_codon].append(mut_num)

    def _add_indel_mutation(self, codon: int, length: int, insertion: bool) -> None:
        """
        Adds a BIG_INDEL (large insertion/deletion) entry to the mutation pool
        """
        first_codon: int = codon - length
        first_exon: int = next(self._codon2exon(first_codon))
        last_exon: int = next(self._codon2exon(codon))
        for i in range(first_exon, last_exon + 1):
            if i in self.unaligned_exons:
                return
        first_ref_codon: int = self.triplet2ref_codon[first_codon] ## TODO: Implement
        last_ref_codon: int = self.triplet2ref_codon[codon]
        first_chain: str = self.exon2chain[first_exon]
        last_chain: str = self.exon2chain[last_exon]
        chain_label: str = f'{first_chain}_{last_chain}'
        first_chrom: str = self.exon2chrom[first_exon]
        last_chrom: str = self.exon2chrom[last_exon]
        chrom_label: str = f'{first_chrom}_{last_chrom}'
        first_codon_strand: bool = self._exon2strand(first_exon)
        last_codon_strand: bool = self._exon2strand(last_exon)
        min_annot_codon: int = self._min_annot_codon()
        if first_codon < min_annot_codon:
            start = (min if first_codon_strand else max)(self.triplet_coordinates[min_annot_codon])
        else:
            start = min(self.triplet_coordinates[first_codon])
        stop: int = (max if last_codon_strand else min)(self.triplet_coordinates[codon])
        exon_label: str = f'{first_exon}_{last_exon}'
        codon_label: str = f'{first_codon}_{codon}'
        ref_codon_label: str = f'{first_ref_codon}_{last_ref_codon}'
        description: str = f'+{length*3}' if insertion else f'-{length*3}'
        mut_id: str = f'BI_{self.indel_counter}'
        indel: Mutation = Mutation(
            self.transcript,
            # chain_label,
            str(first_chain),
            exon_label,
            codon_label,
            ref_codon_label,
            chrom_label,
            start,
            stop,
            BIG_INS if insertion else BIG_DEL,
            description,
            True, ## for compliance with the default no_fpi mode of TOGA 1.0 #False,
            '-',
            mut_id
        )
        mut_num: int = len(self.mutation_list)
        self.mutation_list.append(indel)
        self.indel_counter += 1
        for c in range(first_codon, codon + 1):
            self.codon2mutations[c].append(mut_num)
            ref_c: int = self.triplet2ref_codon[c]
            if mut_num not in self.ref_codon2mutations[ref_c]:
                self.ref_codon2mutations[ref_c].append(mut_num)

    def _check_for_stop_loss(self) -> None:
        """For the last codon, determine if it is a valid stop codon """
        last_codon: int = max(self.all_query_codons)
        query_codon: str = self.query_codons_to_mask.get(
            last_codon, self.all_query_codons[last_codon]
        )
        if query_codon in STOPS:
            return
        ref_codon: str = self.ref_codons_to_mask.get(
            last_codon, self.all_ref_codons[last_codon]
        )
        reason: str = (
            'Non-canonical stop in reference' if ref_codon not in STOPS else 'Missing stop masked'
        )
        exon: int = self.exon_num
        if exon in self.unaligned_exons:
            return
        last_ref_codon: int = self.triplet2ref_codon[last_codon]
        stop_loss: Mutation = Mutation(
            self.transcript,
            self.exon2chain[self.exon_num],
            self.exon_num,
            last_codon,
            last_ref_codon,
            *self._codon2coords(last_codon),
            STOP_MISSING,
            f'{ref_codon}->{query_codon}',
            True, ## for now, stop loss is not considered to be an inactivating mutation,
            reason,
            'STOP_LOSS_1'
        )
        mut_num: int = len(self.mutation_list)
        self.mutation_list.append(stop_loss)
        self.codon2mutations[last_codon].append(mut_num)
        self.ref_codon2mutations[last_ref_codon].append(mut_num)

    def _exon_presence_status(self, exon: int) -> None:
        """
        Classify processed exon into missing (M), deleted (D), or intact (I).
        Exons are classified by the following logic:
        I. Exons for which the predefined locus is available and does intersect
        the CESAR alignment locus are automatically classified as intact;
        II. Exons which do not have a predefined locus in the chain/LASTZ data or
        are aligned outside of this locus by CESAR are classified as intact if
        they pass the identity thresholds (45% of nucleotide identity,
        45% of  BLOSUM score);
        III. Exons which do not meet these criteria are still considered intact
        if their boundaries are supported by SpliceAI predictions;
        IV. Exons which do not meet any of these requirements are considered
        missing if their alignment locus intersects the assembly gap or deleted
        otherwise

        See ref. paper, figure S<> for visual representation and additional data.
        """
        abs_nid: int = 0
        ref_exon_len: int = 0
        blosum_self: int = 0
        blosum_diff: int = 0
        merged_exons: Tuple[int] = self._merged_exon_streak(exon)
        for ex in merged_exons:
            _nid: int = self.abs_nuc_ids[ex]
            _exon_len: int = self.ref_exon_lengths[ex]
            abs_nid += _nid
            ref_exon_len += _exon_len
            _blosum_self: float = self.blosum_self_values[ex]
            _blosum_diff: float  = self.blosum_diff_values[ex]
            if _blosum_self is None or _blosum_diff is None:
                if ex in self.unaligned_exons:
                    _blosum_self: float = 0.0
                    _blosum_diff: float = 0.0
                else:
                    ref_codons: Tuple[int] = sorted(self.exon2ref_codons[exon])
                    first_codon: int = ref_codons[0]
                    first_triplet: int = min(self.ref_codon2triplets[first_codon])
                    last_triplet: int = max(self.ref_codon2triplets[first_codon])
                    ref_only_codon: int = ''.join(
                        self.ref_codons_to_mask.get(
                            x, self.all_ref_codons[x]
                        ) for x in range(first_triplet, last_triplet + 1)
                    )
                    query_only_codon: int = ''.join(
                        self.query_codons_to_mask.get(
                            x, self.all_query_codons[x]
                        ) for x in range(first_triplet, last_triplet + 1)
                    )
                    ref_aa, query_aa = next(process_and_translate(ref_only_codon, query_only_codon)) ## NEEDS TESTING
                    _blosum_self = get_blosum_score(ref_aa, ref_aa, self.matrix)
                    _blosum_diff = get_blosum_score(ref_aa, query_aa, self.matrix)
            blosum_self += _blosum_self
            blosum_diff += _blosum_diff
        nid: float = (abs_nid / ref_exon_len) * 100
        blosum: float = safe_div(nn(blosum_diff), nn(blosum_self)) * 100
        first_codon, last_codon = self.exon2codons[exon]
        first_ref_codon, last_ref_codon = self.exon2ref_codons[exon]
        ## Checkpoint 0: Classify the gap-aligning exons as Missing before going any further
        if len(merged_exons) == 1:
            query_seq: str = self._exon_seq(exon, ref=False)
            total_len: int = len(query_seq)
            defined_len: int = total_len - query_seq.upper().count('N')
            if (defined_len / total_len) < 0.1:
                self._to_log(
                    'Exon %s has only %.3f %% of its sequence defined; classifying it as Missing' % (
                        exon, (defined_len / total_len) * 100
                    )
                )
                self.exon_presence[exon] = 'M'
                mutation: Mutation = Mutation(
                    self.transcript,
                    self.exon2chain[exon],
                    exon,
                    f'{first_codon}_{last_codon}',
                    f'{first_ref_codon}_{last_ref_codon}',
                    self.exon2chrom[exon],
                    *self.abs_exon_coords[exon].tuple(),
                    MISS_EXON,
                    '-',
                    False, ## missing exon entries are not masked
                    '-', ## therefore, the reason slot is kept empty
                    f'MIS_{self.missing_exon_counter}'
                )
                self.missing_exon_counter += 1
                return
        ## Checkpoint 1: Exon resides in its chain-defined locus
        is_present: bool = self.expected_coordinates[exon][0] is not None
        is_real: bool = exon not in self.unaligned_exons
        is_defined: bool = exon not in self.gap_located_exons
        has_aln_data: bool = is_present and is_real and is_defined
        if has_aln_data:
            overlaps_locus: bool = intersection(
                *self.abs_exon_coords[exon].tuple(), *self.expected_coordinates[exon]
            ) > 0
        else:
            overlaps_locus: bool = False
        if overlaps_locus:# or not is_defined:
            self.found_in_exp_locus[exon] = True
        if (has_aln_data and overlaps_locus):
            self._to_log(f'Exon {exon} is found in its expected locus')
            self.exon_presence[exon] = 'I'
            return
        ## Checkpoint 2: For exons with no chained-defined locus,
        ## check if they exceed quality thresholds
        ## TODO: The comparison can be moved to a separate function 
        ## to be executed within process()
        if nid >= MIN_ID_THRESHOLD and blosum >= MIN_BLOSUM_THRESHOLD:
            self._to_log(f'Exon {exon} has sufficient identity values')
            self.exon_presence[exon] = 'I'
            return
        ## Checkpoint 3: If both splice site locations are supported by
        ## SpliceAI predictions, the exon is still worth saving
        sai_acc_support: bool = (
            self.spliceai_acc_support[min(merged_exons)] or
            min(merged_exons) == max(merged_exons) == 1 != self.exon_num
        )
        sai_donor_support: bool = (
            self.spliceai_donor_support[max(merged_exons)] or
            max(merged_exons) == min(merged_exons) == self.exon_num != 1
        )
        if sai_acc_support and sai_donor_support:
            self._to_log(f'Exon boundaries are supported by SpliceAI for exon {exon}')
            self.exon_presence[exon] = 'I'
            return
        ## Checkpoint 4: Exons reaching this point are absent from the alignment;
        ## those which lie in the vicinity of alignment gaps are classified as
        ## missing, otherwise they are marked as deleted
        first_codon, last_codon = self.exon2codons[exon]
        first_ref_codon, last_ref_codon = self.exon2ref_codons[exon]
        if self.intersects_gap[exon] or self.exon2asmgaps[exon] or exon in self.out_of_chain_exons:
            if self.intersects_gap[exon] or self.exon2asmgaps[exon]:
                self._to_log(
                    (
                        'Search space for exon %i contains assembly gap; '
                        'marking exon %i as missing'
                    ) % (exon, exon)
                )
            else:
                self._to_log(
                    (
                        'Exon %i falls outside of the subchain; '
                        'marking exon %i as missing'
                    ) % (exon, exon)
                )
            self.exon_presence[exon] = 'M'
            mutation: Mutation = Mutation(
                self.transcript,
                self.exon2chain[exon],
                exon,
                f'{first_codon}_{last_codon}',
                f'{first_ref_codon}_{last_ref_codon}',
                self.exon2chrom[exon],
                *self.abs_exon_coords[exon].tuple(),
                MISS_EXON,
                '-',
                False, ## missing exon entries are not masked
                '-', ## therefore, the reason slot is kept empty
                f'MIS_{self.missing_exon_counter}'
            )
            self.missing_exon_counter += 1
        else:
            self._to_log(f'Exon {exon} is deleted from the {self.name} projection')
            is_terminal: bool = exon == 1 or exon == self.exon_num
            is_short: bool = self.ref_exon_lengths[exon] < TERMINAL_EXON_DEL_SIZE
            deletion_masked: bool = False
            reason_for_masking: str = '-'
            presence_class: str = 'D'
            mut_class: str = DEL_EXON
            mut_id: str = f'DEL_{self.deleted_exon_counter}'
            if is_terminal and is_short:
                presence_class = 'M'
                mut_class = MISS_EXON
                mut_id = f'MDEL_{self.mdel_counter}'
            elif not self.ref_exon_lengths[exon] % 3 and self.exon_num > 1:
                deletion_masked: bool = True
                reason_for_masking: str = 'Frame-preserving deletion'
                self.deletion2mask[exon] = True
            else:
                self.deletion2mask[exon] = False
            self.exon_presence[exon] = presence_class
            mutation: Mutation = Mutation(
                self.transcript,
                self.exon2chain[exon],
                exon,
                f'{first_codon}_{last_codon}',
                f'{first_ref_codon}_{last_ref_codon}',
                self.exon2chrom[exon],
                *self.abs_exon_coords[exon].tuple(),
                mut_class,
                '-',
                deletion_masked,
                reason_for_masking,
                mut_id
            )
            if 'MDEL' in mut_id:
                self.mdel_counter += 1
            else:
                self.deleted_exon_counter +=1
        self.mutation_list.append(mutation)
        mut_num: int = len(self.mutation_list) - 1
        for codon in range(first_codon, last_codon):
            self.codon2mutations[codon].append(mut_num)
            ref_codon: int = self.triplet2ref_codon[codon]
            if mut_num not in self.ref_codon2mutations[ref_codon]:
                self.ref_codon2mutations[ref_codon].append(mut_num)
        self.exon2mutations[exon].append(mut_num)

    def _get_splice_pairs(self) -> List[Tuple[int, int]]:
        """
        Infer pairs of exons to be spliced given the known exon presence status
        Exons are considered to be a spliced pair if they follow each other or
        are separated by one or more deleted exons.
        Note that missing exons disrupt splice pairs since it their effect on
        actual splicing cannot be estimated.
        For example, in the following transcript:
        1 2 3 4 5 6 7 8 9 10
        I I D D I D M I D I
        , splice pairs are (1,2), (2,5), (5,8), and (8,10)
        """
        output: List[Tuple[int, int]] = []
        from_: int =  None
        for exon, status in self.exon_presence.items():
            if status == 'M': 
                ## unset the pair start since M exons disrupt pairs
                from_ = None
                continue
            if status == 'I':
                if from_: 
                    ## add a pair of spliced exons
                    output.append((from_, exon))
                ## set a new pair start to the current exon
                from_ = exon
        return output

    def _check_updated_split_codons(self) -> None:
        """
        After exon classification and deleted exon exclusion,
        checks novel spliced exons for potential mutations in the split codons
        """
        if self.exon_num == 1: ## irrelevant for single-exon transcripts
            return
        splice_pairs: List[Tuple[int, int]] = self._get_splice_pairs()
        for from_, to_ in splice_pairs:
            if from_ + 1 == to_: ## the site has been already analysed during .process()
                continue
            ## reconstruct the respective split codons
            donor_codon_num: int = self.exon2codons[from_][1]
            donor_codon: str = self.query_codons_to_mask.get(
                donor_codon_num, self.all_query_codons[donor_codon_num]
            )
            donor_portion: int = self.split_codon_struct[donor_codon_num][from_]
            donor_split: str = ''
            for i in donor_codon:
                donor_split += i
                if sum([x.isalpha() for x in donor_split]) != donor_portion:
                    continue
                break
            acc_codon_num: int = self.exon2codons[to_][0]
            acc_codon: str = self.query_codons_to_mask.get(
                acc_codon_num, self.all_query_codons[acc_codon_num]
            )
            acc_portion: int = self.split_codon_struct[acc_codon_num][to_]
            acc_split: str = ''
            for j in acc_codon[::-1]:
                acc_split += j
                if sum([x.isalpha() for x in acc_split]) != acc_portion:
                    continue
                break
            acc_split = acc_split[::-1]
            ## compare the resulting codon to the latter split codon in the reference
            ref_codon: str = self.all_ref_codons[acc_codon_num]
            ## check the resulting codon
            restored_codon: str = (donor_split + acc_split).replace('-', '')
            if restored_codon == '':
                continue
            for sub in parts(restored_codon, 3):
                curr_mut_num: int = len(self.mutation_list)
                self._check_for_frameshift(
                    acc_codon_num, ref_codon, restored_codon, to_
                )
                self._check_for_nonsense_mutation(
                    acc_codon_num, ref_codon, restored_codon, to_
                )
                if len(self.mutation_list) > curr_mut_num:
                    for m in range(curr_mut_num, len(self.mutation_list)):
                        self.mutation_list[m].is_masked = False
                        self.mutation_list[m].masking_reason = ALT_MASKING_REASON

    def _ref_len_till_codon(self, codon: int) -> int:
        """
        Returns the relative position of a given codon in the ungapped reference
        sequence
        """
        len_till_codon: int = 0
        for c in range(1, codon + 1): 
            ## TODO: +1 added for the sake of compliance with 1.0; 
            ## I believe this is likely a bug on the 1.0's side
            codon_seq: str = self.ref_codons_to_mask.get(c, self.all_ref_codons[c])
            len_till_codon += sum(x.isalpha() for x in codon_seq)
        return len_till_codon

    def _to_mask(self, codon: int, is_deletion: bool) -> Tuple[bool, str]:
        """
        Defines whether mutation should be masked for inactivated status check
        """
        ## calculate the first and last 10% of the reference sequence
        first_ten: int = round(self.reference_length / 10)
        last_ten: int = self.reference_length - first_ten
        ## find the position of the affected codon in the reference sequence
        len_till_codon: int = self._ref_len_till_codon(codon)
        ## define whether the codon falls into first or last 10% of the reference
        codon_in_first_ten: bool = len_till_codon <= first_ten
        codon_in_last_ten: bool = len_till_codon >= last_ten
        portion: float = len_till_codon / self.reference_length * 100
        portion = round(portion if codon_in_first_ten else 100 - portion, 1)

        ## RECALC FOR THE SAKE OF CONSISTENCY WITH 1.0; TO BE REVERTED UPON RELEASE
        first_ten: int = max(self.ref_codon2triplets) // 10
        last_ten: int = max(self.ref_codon2triplets) - first_ten
        codon_in_first_ten: bool = self.triplet2ref_codon[codon] <= first_ten
        codon_in_last_ten: bool = self.triplet2ref_codon[codon] >= last_ten
        portion: float = len_till_codon / self.reference_length * 100
        portion = round(portion if codon_in_first_ten else 100 - portion, 1)

        ## if mask_terminal_mutation is set, mask all mutations
        ## in the first and 10% of the sequence; mutations in last 10% are masked regardless
        first_ten_condition: bool = (
            self.mask_terminal_mutations and codon_in_first_ten and not is_deletion
        )
        if first_ten_condition or codon_in_last_ten:
            ## mask the mutation if it falls within the first/last 10%,
            ## calculate the actual portion it falls in
            terminus: str = 'N' if codon_in_first_ten else 'C'
            return True, f'{terminus}-terminal {portion}%'
        to_mask: bool = any(
            (
                x  for x in self.query_atgs if x > codon and
                self.triplet2ref_codon[x] <= first_ten
            )
        )
        return (
            to_mask,
            'Alternative start found' if to_mask else '-'
        )

    def _restored_frame_mut(self, mut: Mutation) -> bool:
        comp_mut: bool = mut.mutation_class == COMPENSATION
        exon_nums: Iterable[int] = [mut.exon] if isinstance(mut.exon, int) else map(int, mut.exon.split('_'))
        exon_frame_corr: bool = any(self.corrected_frame[x] for x in exon_nums)
        is_in_alt_frame: bool = exon_frame_corr and mut.masking_reason == ALT_FRAME_REASON
        is_comp: bool = exon_frame_corr and mut.masking_reason == COMPENSATION_REASON
        return comp_mut or is_comp or is_in_alt_frame

    def _architecture_preserving_mutation(self, mut: Mutation, num: int) -> bool:
        """
        Returns whether mutation preserves the gene architecture, i.e.:
        1) is a known U12/non-canonical U2 site;
        2) was caused by splice site shift;
        3) corresponds to intron gain or precise intron deletion;
        """
        if mut.masking_reason == 'Intron deletion':
            return True
        if mut.mutation_class == INTRON_GAIN:
            return True 
        if self._restored_frame_mut(mut):
            return True
        if mut.masking_reason in SAFE_SPLICE_SITE_REASONS:
            return True
        if mut.mutation_class in BIG_INDEL:
            return True
        if mut.masking_reason == COMPENSATION_REASON: ## TODO: Would this be correct?
            comp_num: int = self.frameshift2compensation[num]
            ## check if it's a deprecated (unmasked) compensation   
            if comp_num in self.comps_to_delete:
                return False
            ## check if it's a SpliceAI-corrected compensation
            if comp_num in self.spliceai_compensations:
                return True
            frameshifts: List[Mutation] = [
                x for i,x in enumerate(self.mutation_list)
                if self.frameshift2compensation.get(i, -1) == comp_num
            ]
            return any(
                self._in_modified_seq(x.exon, x.codon) for x in frameshifts
            )
        return False

    def _alt_loss_status(self, status: str) -> str:
        if self.is_paralog:
            return PG
        if self.is_processed_pseudogene:
            return PP
        return status

    def _calc_id_values(self) -> None:
        """Calculates the nucleotide and BLOSUM %identity for the projection"""
        abs_id: int = sum(self.abs_nuc_ids.values())
        reference_length: int = sum(self.ref_exon_lengths.values())
        self.nuc_id: float = (abs_id / reference_length) * 100
        blosum_self: int = sum(
            x if x is not None else 0.0 for x in self.blosum_self_values.values()
        )
        blosum_diff: int = sum(
            x if x is not None else 0.0 for x in self.blosum_diff_values.values()
        )
        self.blosum_id: float = safe_div(nn(blosum_diff), nn(blosum_self)) * 100

    def _exon_seq(self, exon: int, ref: bool) -> str:
        """For a given exon number, return the aligned nucleotide sequence"""
        cesar_start, cesar_stop = self.cesar_rel_coords[exon].tuple()
        adj_start, adj_stop = self.rel_exon_coords[exon].tuple()
        if ref:
            seq: str = self.reference[cesar_start:cesar_stop]
            if adj_start < cesar_start:
                seq = '-' * (cesar_start - adj_start) + seq
            if adj_stop > cesar_stop:
                if exon == self.exon_num and self.stop_updated:
                    init_stop: str = seq[-3:]
                    seq = seq[:-3] + '-' * (adj_stop - cesar_stop - 3) + init_stop
                else:
                    seq += '-' * (adj_stop - cesar_stop)
        else:
            start: int = min(cesar_start, adj_start)
            stop: int = max(cesar_stop, adj_stop)
            seq: str = self.query[start:stop]
            if adj_start != cesar_start:
                if adj_start > cesar_start:
                    seq = '-' * (adj_start - cesar_start) + seq[(adj_start - cesar_start):]
                else:
                    seq = seq[:(cesar_start-start)].upper() + seq[(cesar_start-start):]
            if adj_stop != cesar_stop:
                if adj_stop < cesar_stop:
                    seq = seq[:-(cesar_stop - adj_stop)] + '-' * (cesar_stop - adj_stop)
                else:
                    if exon == self.exon_num and self.stop_updated:
                        seq = seq[:-(adj_stop - cesar_stop + 3)] + seq[-(adj_stop - cesar_stop):].upper()
                    else:
                        seq = seq[:-(adj_stop - cesar_stop)] + seq[-(adj_stop - cesar_stop):].upper()
        return seq

    def _cds_seq(self) -> str:
        """
        Retrieves the complete nucleotide coding sequence for the query 
        """
        cds: str = ''
        for exon in range(1, self.exon_num + 1):
            exon_seq: str = self._exon_seq(exon, ref=False)
            exon_seq = strip_noncoding(exon_seq, uppercase_only=True)
            cds += exon_seq
        return cds

    def _query_protein_seq(self) -> str:
        """Returns the reference protein sequence corrected for all compensated frameshifts"""
        frame_changed: bool = False
        ## check whether there are any uncompensated frameshifts;
        ## if there are any, report the original protein sequence
        for mut in self.mutation_list:
            if mut.mutation_class not in FS_INDELS:
                continue
            if mut.masking_reason != COMPENSATION_REASON:
                frame_changed = False
                break
            frame_changed = True
        if frame_changed:
            codon_seq: str = ''.join(
                self.query_codons_to_mask.get(x, self.all_query_codons[x]).upper() 
                for x in self.all_query_codons
            ).replace('-', '')
            if not (len(codon_seq) % 3):
                return ''.join(AA_CODE.get(x, 'X') for x in parts(codon_seq, 3))
        return ''.join(x for x in self.query_aa_seq.values() if x != '-')

    def _intact_exon_portion(self) -> Tuple[int, float]:
        """
        Reports the percent of reference exons not affected by any mutation in the query
        """
        affected_exons: Set[int] = set()
        for mut in self.mutation_list:
            if mut.mutation_class in SAFE_UNMASKABLE_TYPES:
                continue
            if mut.masking_reason in SAFE_UNMASKABLE_REASONS:
                continue
            if isinstance(mut.exon, str):
                first, last = map(int, mut.exon.split('_'))
                for x in range(first, last + 1):
                    affected_exons.add(x)
            else:
                affected_exons.add(mut.exon)
        num: int = len(affected_exons)
        percent: float = round((num / self.exon_num) * 100, 2)
        return (num, percent)

    def bed12(self, raw: bool = True, browser: bool = False) -> str:
        """
        Returns projection represenntation in BED format.
        By default, formats basic projection§ properties as a BED12 entry.
        Fragmented projections are split into separate entries for each input chain.

        If ::raw:: is set to True (default value), reports all the exons with
        anyhow defined coordinates; otherwise, missing and deleted exons are
        removed from the final BED12 block layout.

        If ::chain_id:: is set to True (default feature), populates the 'score'
        column with an actual chain ID used for reconstructing the fragment.
        Note that valid BED scores range between 0 and 1000 while chain IDs
        range from 1 to potentially infinity, therefore scores should be replaced
        with a proper placeholder if the BED file is further used for UCSC browser
        upload or any other BED-handling software manipulations.
        """
        color: str = CLASS_TO_COL[self.loss_status]
        if browser:
            exon_alns: List[str] = [
                self._exon_aln_for_browser(exon) for exon in range(1, self.exon_num + 1)
            ]
        for chain in self.chains:
            exons: List[int] = [
                exon for exon in range(1, self.exon_num + 1) if self.exon2chain[exon] == chain
            ]
            if not exons and len(self.chains) > 1:
                continue
            strand: bool = self.portion2strand[self.exon2portion[exons[0]]]
            chrom: str = self.exon2chrom[exons[0]]
            exons = exons if strand else exons[::-1]
            lengths: List[int] = []
            starts: List[int] = []
            ex_num: int = len(exons)
            _stop: Union[int, None] = None
            last_exon: int = 0
            for exon in exons:
                portion: int = self.exon2portion[exon]
                if not raw and self.exon_presence[exon] != 'I':
                    ex_num -= 1
                    continue
                exon_seq: str = self._exon_seq(exon, ref=False).replace('-', '')
                if exon in self.unaligned_exons or not len(exon_seq):
                    ex_num -= 1
                    continue
                start, stop = self.abs_exon_coords[exon].tuple()
                if stop - start == 0:
                    ex_num -= 1
                    continue
                if exon in self.introns_gained:
                    sub_exons: List[int] = [
                        sorted(
                            map(lambda x: self._abs_coord(x, portion=portion), x[:2])
                        ) for x in self.introns_gained[exon]
                    ]
                    rel_start, rel_stop = self.rel_exon_coords[exon].tuple()
                    if not strand:
                        sub_exons = sub_exons[::-1]
                    ex_num += len(sub_exons) - 1
                else:
                    sub_exons: List[int] = [(start, stop)]
                for i, (sub_start, sub_stop) in enumerate(sub_exons, start=1):
                    if i == len(sub_exons):
                        sub_stop = stop
                    if _stop is None:
                        cds_start: int = sub_start
                    if sub_start == _stop:
                        ex_num -= 1
                        lengths = lengths[:-1]
                        prev_start: int = starts[-1] + cds_start
                        if stop - prev_start == 0:
                            ex_num -= 1
                            continue
                        length: int = sub_stop - prev_start
                        if length <= 0:
                            raise ValueError(
                                f'Non-positive length observed in exon {exon}'
                            )
                        lengths.append(length)
                    else:
                        length: int = sub_stop - sub_start
                        if length <= 0:
                            raise ValueError(
                                f'Non-positive length observed in exon {exon}'
                            )
                        lengths.append(length)
                        starts.append(sub_start - cds_start)
                    _stop = sub_stop
                last_exon = exon
            if not ex_num: ## possible only if all the exons were not aligned or are completely deleted
                self._to_log('All exons for current fragment have been deleted')
                continue
            cds_stop: int = self.abs_exon_coords[last_exon].tuple()[1]
            lengths = ','.join(map(str, lengths)) + ','
            starts = ','.join(map(str, starts)) + ','
            score: str = str(chain) if browser else '0'
            output: List[Any] = [
                chrom, cds_start, cds_stop,
                f'{self.name}', score,
                '+' if strand else '-', cds_start, cds_stop, color,
                ex_num, lengths, starts
            ]
            if browser:
                protein_aln: str = self._aa_aln_for_browser()
                plot: str = self.svg.replace('\n', '')
                mut_table: str = mutation_table(self.mutation_list)
                exon_aln_field: str = ''.join(exon_alns)
                cds_field: str = self._cds_seq()
                alt_prot_field: str = self._query_protein_seq()
                affected_num, affected_percent = self._intact_exon_portion()
                output.extend(
                    [
                        CLASS_TO_NAME[self.loss_status],
                        self.longest_fraction_strict ,
                        self.longest_fraction_relaxed, 
                        self.total_intact_fraction,
                        '0.0',
                        int(self.middle_is_intact),
                        int(self.middle_is_present),
                        protein_aln, plot, mut_table, 
                        exon_aln_field, cds_field, alt_prot_field,
                        affected_num, affected_percent
                    ]
                )
                pass
            yield '\t'.join(map(str, output))

    def _aa_aln_for_browser(self) -> str:
        """
        Generates a protein alignment line for the UCSC report
        """
        ref_seq: str = ''.join(self.ref_aa_seq.values())
        query_seq: str = ''.join(self.query_aa_seq.values())
        return format_fasta_as_aln(ref_seq, query_seq, protein=True)

    def _orth_label(self) -> str:
        """Returns orthology label for exon alignment"""
        if self.is_paralog:
            return PARALOG
        if self.is_processed_pseudogene:
            return PROC_PSEUDOGENE
        return ORTHOLOG

    def _exon_aln_for_browser(self, exon: int) -> str:
        ref_seq: str = self._exon_seq(exon, ref=True).replace('>', '-').replace(' ', '-')
        query_seq: str = self._exon_seq(exon, ref=False).replace('>', '-')
        ## TODO: add splice sites
        chrom: str = self.exon2chrom[exon]
        start, end = self.abs_exon_coords[exon].tuple()
        exp_start, exp_end = self.expected_coordinates[exon]
        found_in_exp: bool = self.found_in_exp_locus[exon]
        nuc_id: float = self.nuc_ids[exon]
        blosum: float = self.blosum_ids[exon]
        has_gap: bool = self.intersects_gap[exon] or self.exon2asmgaps[exon]
        aln_class: str = self.exon_quality[exon]
        header: str = exon_aln_header(
            exon, chrom, start, end, exp_start, exp_end, 
            found_in_exp, nuc_id, blosum, has_gap, aln_class
        )
        return exon_aln_entry(ref_seq, query_seq, header)

    def bdb(self) -> str:
        """
        Returns the alignment processing results in the legacy format
        """
        r_post: str = ' | REFERENCE'
        q_post: str = ' | QUERY'
        chains: str = ','.join(self.chains)
        prot_name: str = f'>{self.name} | {chains} | PROT'
        codon_name: str = f'>{self.name}| {chains} | CODON'
        ref_exon_template: str = f'>{self.name} | {{}} | {{}} | reference_exon'
        query_exon_template: str = (
            f'>{self.name} | {{}} | {{}} | {{}}:{{}}-{{}} | {{}} | {{}} | N/A | '
            'N/A | exp:N/A-N/A | N/A | False | query_exon'
        )
        ref_prot: str = prot_name + r_post
        query_prot: str = prot_name + q_post
        ref_codon: str = codon_name + r_post
        query_codon: str = codon_name + q_post
        output: str = [
            self.name, ref_prot, ''.join(self.ref_aa_seq.values()),
            query_prot, ''.join(self.query_aa_seq.values()),
            ref_codon, ' '.join(self.all_ref_codons.values()),
            query_codon, ' '.join(self.all_query_codons.values())
        ]
        for e in range(1, self.exon_num + 1):
            nuc_id: float = round(self.nuc_ids[e], 2)
            blosum_id: float = round(self.blosum_ids[e], 2)
            ref_exon_name: str = ref_exon_template.format(e, self.exon2chain[e])
            query_exon_name: str = query_exon_template.format(
                e, self.exon2chain[e], self.exon2chrom[e],
                *self.abs_exon_coords[e].tuple(), nuc_id, blosum_id
            )
            ref_seq: str = self._exon_seq(e, ref=True).replace('>', '-').replace(' ', '-')
            query_seq: str = self._exon_seq(e, ref=False).replace('>', '-').upper()
            output.extend([ref_exon_name, ref_seq, query_exon_name, query_seq])
        return '\n'.join(output) + '\n'

    def codon_fasta(self) -> str:
        """
        Returns codon alignment in FASTA format,
        with codons separated by whitespaces
        """
        chains: str = ','.join(self.chains)
        r_name: str = f'>{self.name}| {chains} | CODON | REFERENCE'
        q_name: str = f'>{self.name}| {chains} | CODON | QUERY'
        return '\n'.join(
            [
                r_name,
                ' '.join(self.all_ref_codons.values()),
                q_name,
                ' '.join(self.all_query_codons.values())
            ]
        )

    def aa_fasta(self) -> str:
        """
        Returns amino acid alignment in FASTA format
        """
        chains: str = ','.join(self.chains)
        r_name: str = f'>{self.name} | PROT | REFERENCE'
        q_name: str = f'>{self.name} | PROT | QUERY'
        return '\n'.join(
            [
                r_name, ''.join(self.ref_aa_seq.values()),
                q_name, ''.join(self.query_aa_seq.values())
            ]
        )

    def exon_fasta(self) -> str:
        """
        Returns exon nucleotide alignment in FASTA format
        """
        output: List[str] = []
        ref_exon_template: str = f'>{self.name} | {{}} | {{}} | reference_exon'
        query_exon_template: str = (
            f'>{self.name} | {{}} | {{}} | {{}}:{{}}-{{}} | {{}} | {{}} | '
            f'{{}}:{{}} | {{}} | {{}} | query_exon'
        )
        orth: str = self._orth_label()
        for e in range(1, self.exon_num + 1):
            nuc_id: float = round(self.nuc_ids[e], 2)
            blosum_id: float = round(self.blosum_ids[e], 2)
            ref_exon_name: str = ref_exon_template.format(e, self.exon2chain[e])
            found_in_exp: str = 'INC' if self.found_in_exp_locus[e] else 'EXCL'
            exp_start, exp_stop = self.expected_coordinates[e]
            if exp_start is not None:
                exp_coords: str = f'{exp_start}-{exp_stop}'
            else:
                exp_coords: str = 'NA-NA'
            query_exon_name: str = query_exon_template.format(
                e, self.exon2chain[e], self.exon2chrom[e],
                *self.abs_exon_coords[e].tuple(), nuc_id, blosum_id,
                self.exon2chrom[e], exp_coords, found_in_exp, orth
            )
            ref_seq: str = self._exon_seq(e, ref=True).replace('>', '-').replace(' ', '-')
            query_seq: str = self._exon_seq(e, ref=False).replace('>', '-')
            output.extend([ref_exon_name, ref_seq, query_exon_name, query_seq])
        return '\n'.join(output)

    def exon_meta(self) -> Iterable[str]:
        """
        Returns a generator object of tab-separated exon metadata lines
        """
        for ex in range(1, self.exon_num + 1):
            chain: str = self.exon2chain[ex]
            chrom: str = self.exon2chrom[ex]
            start, stop = self.abs_exon_coords[ex].tuple()
            strand: str = '+' if self._exon2strand(ex) else '-'
            was_aligned: bool = 'UNALIGNED' if ex in self.unaligned_exons else 'ALIGNED'
            has_gap: str = 'INTERSECTS_GAP' if self.intersects_gap[ex] else 'GAP_FREE'
            exp_start, exp_stop = self.expected_coordinates[ex]
            if exp_start is None:
                exp_coords: str = f'{chrom}:NA-NA'
            else:
                exp_coords: str = f'{chrom}:{exp_start}-{exp_stop}'
            found_in_exp: str = 'INCL' if self.found_in_exp_locus[ex] else 'EXCL'
            # intersects_gap: str = ''
            start_from_cesar: bool = self.cesar_acc_support[ex]
            stop_from_cesar: bool = self.cesar_donor_support[ex]
            start_from_cesar = (
                'START_ALIGNED' if start_from_cesar or ex == 1 else 'START_UNALIGNED'
            )
            stop_from_cesar = (
                'STOP_ALIGNED' if stop_from_cesar or ex == self.exon_num else 'STOP_UNALIGNED'
            )
            start_from_sai: bool = self.spliceai_acc_support[ex]
            stop_from_sai: bool = self.spliceai_donor_support[ex]
            start_from_sai = (
                'FIRST_EXON' if ex == 1 else 'ACC_SUPPORTED' if start_from_sai
                else 'ACC_UNSUPPORTED'
            )
            stop_from_sai = (
                'DONOR_SUPPORTED' if stop_from_sai else 'LAST_EXON'
                if ex == self.exon_num else 'DONOR_UNSUPPORTED'
            )
            acc_prob: bool = self.acc_site_prob[ex]
            donor_prob: bool = self.donor_site_prob[ex]
            chain_support: str = (
                'CHAIN_UNSUPPORTED' if ex in self.out_of_chain_exons else 'CHAIN_SUPPORTED'
            )
            out_line: List[str] = [
                self.name, ex, chain, chrom, start, stop,
                strand, self.exon_presence[ex], was_aligned,
                self.exon_quality[ex], start_from_cesar, stop_from_cesar,
                start_from_sai, acc_prob, stop_from_sai, donor_prob,
                exp_coords, found_in_exp, has_gap, chain_support,
                self.nuc_ids[ex], self.blosum_ids[ex]
            ]
            yield '\t'.join(map(str, out_line))

    def mutation_file(self) -> str:
        """Returns a TSV-formatted list of mutations found"""
        output: str = ''
        for exon, subexons in self.introns_gained.items():
            portion: int = self.exon2portion[exon]
            introns: List[Tuple[int, int]] = [
                (
                    self._abs_coord(subexons[i-1][1], portion=portion), 
                    self._abs_coord(subexons[i][0], portion=portion)
                )
                for i in range(1, len(subexons))
            ]
            chrom: str = self.exon2chrom[exon]
            strand: bool = self._exon2strand(exon)
            intron_gain_counter: int = 1
            for intron in introns:
                first_codon: int = self._base2codon(intron[0], exon)
                last_codon: int = self._base2codon(intron[1], exon)
                codon_label: str = f'{first_codon}_{last_codon}'
                first_ref_codon: int = self._base2codon(intron[0], exon, True)
                last_ref_codon: int = self._base2codon(intron[1], exon, True)
                ref_codon_label: str = f'{first_ref_codon}_{last_ref_codon}'
                mutation_id: str = f'GAIN_{intron_gain_counter}'
                mut_record: Mutation = Mutation(
                    self.transcript,
                    self.exon2chain[exon],
                    exon,
                    codon_label,
                    ref_codon_label,
                    chrom,
                    *sorted(intron),
                    INTRON_GAIN,
                    '-',
                    True, ## intron gains are not considered as inactivating mutations
                    INTRON_GAIN_MASK_REASON,
                    mutation_id
                )
                self.mutation_list.append(mut_record)
                intron_gain_counter += 1
        mutations: List[Mutation] = sorted(
            self.mutation_list, key=lambda x: (
                x.exon if isinstance(x.exon, int) else int(x.exon.split('_')[0]),
                x.codon if isinstance(x.codon, int) else int(x.codon.split('_')[0])
            )
        )
        for mut in mutations:
            output += str(mut) + '\n'
        return output

    def transcript_meta(self) -> str:
        """
        Returns the intact/present portion stats as defined in the older
        CESAR wrapper output formatted in a TSV manner
        """
        return '\t'.join(
            [
                self.name,
                self.loss_status,
                str(self.nuc_id),
                str(self.blosum_id),
                str(self.longest_fraction_strict),
                str(self.longest_fraction_relaxed),
                str(self.total_intact_fraction),
                'MIDDLE_IS_INTACT' if self.middle_is_intact else 'MIDDLE_IS_AFFECTED',
                'MIDDLE_IS_PRESENT' if self.middle_is_present else 'MIDDLE_IS_MISSING'
            ]
        )

    def cds_nuc(self) -> str:
        """Returns stripped query nucleotide sequence in FASTA format"""
        header: str = f'>{self.name}'
        seq: str = self._cds_seq()
        return f'{header}\n{seq}'

    def splice_site_table(self) -> str:
        """
        Returns the tab-separated four column table containing splice site
        dinucleotides in the "projection-exon-acceptor-donor" format
        """
        out_str: str = ''
        for exon in self.splice_site_nucs:
            acc: str = self.splice_site_nucs[exon]['A']
            # acc_class: str = self.u12[exon]['acceptor'][0]
            donor: str = self.splice_site_nucs[exon]['D']
            # donor_class: str = self.u12[exon]['donor'][0]
            out_str += f'{self.name}\t{exon}\t{acc}\t{donor}'
            if exon < self.exon_num:
                out_str += '\n'
        return out_str

    def selenocysteine_codon_table(self) -> str:
        out_str: str = ''
        if not self.selenocysteine_codons:
            return out_str
        sc_num: int = len(self.selenocysteine_codons)
        for i, (exon, codon_num, codon) in enumerate(self.selenocysteine_codons, start=1):
            if exon is None:
                exon = next(self._codon2exon(codon_num))
            strand: bool = self._exon2strand(exon)
            chrom: str = self.exon2chrom[exon]
            coords: List[int] = self.triplet_coordinates[codon_num]
            start: int = min(coords) - (0 if strand else 1)
            end: int = max(coords) + (1 if strand else 0)
            out_line: str = '\t'.join(
                map(str, (self.name, exon, codon_num, chrom, start, end, codon))
            )
            out_str += out_line
            if i < sc_num:
                out_str += '\n'
        return out_str


    def splice_site_shifts(self) -> str:
        """
        Returns a tab-separated table containing data on SpliceAI-guided splice site shifts
        """
        out_str: List[str] = []
        if self.correction_mode <= 1:
            return ''
        for ex in range(1, self.exon_num + 1):
            final_start, final_end = self.abs_exon_coords[ex].tuple()
            coord_str: str = f'{self.exon2chrom[ex]}:{final_start}-{final_end}'
            strand: bool = self._exon2strand(ex)
            strand_str: str = '+' if strand else '-'
            p: int = self.exon2portion[ex]
            abs_cesar_start, abs_cesar_end = sorted(
                map(lambda x: self._abs_coord(x, exon=ex, portion=p), self.cesar_rel_coords[ex].tuple())
            )
            acc_from_sai: bool = self.spliceai_acc_support[ex]
            donor_from_sai: bool = self.spliceai_donor_support[ex]
            start_shifted: bool = final_start != abs_cesar_start
            end_shifted: bool = final_end != abs_cesar_end
            acc_shifted: bool = acc_from_sai and (start_shifted if strand else end_shifted)
            donor_shifted: bool = donor_from_sai and (end_shifted if strand else start_shifted)
            if acc_shifted:
                side: str = 'acceptor'
                intron_class: str = self.u12[ex][side][0]
                diff: int = (final_start - abs_cesar_start) if strand else (final_end - abs_cesar_end)
                diff_str: str = ('+' if diff > 0 else '') + str(diff)
                prob: bool = self.acc_site_prob[ex]
                dinuc: str = self.splice_site_nucs[ex]['A']
                out_line: Tuple[str] = (
                    self.name, str(ex), coord_str, strand_str, 
                    side, diff_str, intron_class, str(prob), dinuc
                )
                out_str.append('\t'.join(out_line))
            if donor_shifted:
                side: str = 'donor'
                intron_class: str = self.u12[ex][side][0]
                diff: int = (final_end - abs_cesar_end) if strand else (final_start  - abs_cesar_start)
                diff_str: str = ('+' if diff > 0 else '') + str(diff)
                prob: bool = self.donor_site_prob[ex]
                dinuc: str = self.splice_site_nucs[ex]['D']
                out_line: Tuple[str] = (
                    self.name, str(ex), coord_str, strand_str, 
                    side, diff_str, str(intron_class), str(prob), dinuc
                )
                out_str.append('\t'.join(out_line))
        out_str: str = '\n'.join(out_str)
        return out_str