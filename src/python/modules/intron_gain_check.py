"""
Uses CESAR2 alignment results and SpliceAI predictions
to search for additional introns in the query data
"""

from cesar_wrapper_constants import *
from cesar_wrapper_executables import CesarInput, RawCesarOutput, fast_seq_id
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
# from functools import reduce
from logging import Logger
from math import floor
# from Operator import mul
from shared import nn, parts,   safe_div
from typing import Any, Dict, List, Optional, Tuple, TypeVar, Union

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = ('Michael Hiller')

ExCand = TypeVar('ExCand', bound='ExonCandidate')
PotExGroup = TypeVar('PotExGroup', bound='PotentialExonGroup')

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


class PotentialExonGroup:
    """Stores coordinates for potential exons resulting from intron gain search"""
    __slots__ = (
        'exon_tuples', 'acceptor_phases', 'donor_phases',
        'redeemed_muts', 'intron_meta'
    )

    def __init__(self) -> None:
        self.exon_tuples: List[Tuple[int, int, float, float]] = []
        self.acceptor_phases: List[int] = []
        self.donor_phases: List[int] = []
        self.redeemed_muts: int = 0
        self.intron_meta: List[IntronMeta] = []

    def add(self, exon: Tuple[int, int, float, float]) -> None:
        self.exon_tuples.append(exon)

    def len(self) -> int:
        return len(self.exon_tuples)

    def group_acc(self) -> int:
        return min(x[0] for x in self.exon_tuples)

    def group_donor(self) -> int:
        return max(x[1] for x in self.exon_tuples)

    def group_donor_prob(self) -> float:
        donor: int = self.group_donor()
        return [x for x in self.exon_tuples if x[1] == donor][0][3]

    def group_acc_prob(self) -> float:
        acc: int = self.group_acc()
        return [x for x in self.exon_tuples if x[0] == acc][0][2]

    def all_sites_exceed_threshold(self, prob: float) -> bool:
        acc: int = self.group_acc()
        donor: int = self.group_donor()
        return all(
            (x[2] >= prob if x[2] != acc else True) and 
            (x[3] >= prob if x[3] != donor else True) 
            for x in self.exon_tuples
        )

    def acceptor_phase(self) -> int:
        if not self.acceptor_phases:
            return 0
        return self.acceptor_phases[-1]

    def donor_phase(self) -> int:
        if not self.donor_phases:
            return 0
        return self.donor_phases[-1]

    def copy(self) -> PotExGroup:
        return deepcopy(self)

    def __repr__(self) -> str:
        return (
            f'PotentialExonGroup[exons={self.exon_tuples}, '
            f'acceptor phase={self.acceptor_phases}, donor phase={self.donor_phases}]'
        )


@dataclass
class IntronMeta:
    __slots__ = ('gap_support', 'mutation_support', 'acc_prob', 'donor_prob')
    gap_support: bool
    mutation_support: bool
    donor_prob: float
    acc_prob: float

    def __repr__(self) -> str:
        gap_support: str = 'GAP_SUPPORTED' if self.gap_support else 'GAP_UNSUPPORTED'
        mut_support: str = 'MUTATION_SUPPORTED' if self.mutation_support else 'MUTATION_UNSUPPORTED'
        return '\t'.join(map(str, (gap_support, mut_support, self.donor_prob, self.acc_prob)))


@dataclass
class ExtraRoundCesarInput:
    """Stores input data for extra round of CESAR alignment"""
    __slots__ = (
        'ref_seqs', 'query_seqs', 'exon_seqs',
        'introns', 'upstream_seqs', 'downstream_seqs', 'ref_phases',
        'needs_realigment', 'el2exons', 'orig_ref', 'orig_query', 'orig_id',
        'intron_evidence'
    )
    ref_seqs: Dict[int, str]
    query_seqs: Dict[int, str]
    exon_seqs: Dict[int, List[str]]
    introns: Dict[int, List[str]]
    upstream_seqs: Dict[int, Tuple[str, str]]
    downstream_seqs: Dict[int, Tuple[str, str]]
    ref_phases: Dict[int, Tuple[int, int]]
    needs_realigment: Dict[int, bool]
    el2exons: Dict[int, List[int]]
    orig_ref: Dict[int, str]
    orig_query: Dict[int, str]
    orig_id: Dict[int, float]
    intron_evidence: Dict[int, List[IntronMeta]]
    # exon2extra_seq: Tuple[int, int]
    # subexon_coords: Dict[int, List[Tuple[int, int, float, float]]]


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


def strip_noncoding(string: str) -> str:
    """Removes noncoding symbols from the nucleotide string"""
    return ''.join(x for x in string if x.isalpha())


def frameshift_codon(codon: str) -> bool:
    """Returns whether a codon contains a frameshift mutation"""
    return codon.count('-') % 3 != 0 and codon != GAP_CODON# or len(codon) % 3 != 0


def binary_streaks(all_values: List[Any], heads: List[Any]) -> List[List[Any]]:
    """
    """
    all_streaks: List[List[Any]] = []
    curr_streak: List[Any] = []
    heads_streak: bool = False
    for x in all_values:
        is_heads: bool = x in heads
        # if is_heads != heads_streak:
        if is_heads != heads_streak or is_heads and heads_streak:
            if curr_streak:
                all_streaks.append(curr_streak)
            curr_streak = []
            heads_streak = is_heads
        curr_streak.append(x)
    if curr_streak:
        all_streaks.append(curr_streak)
    return all_streaks


def calc_intron_score(
    coverage_score: float,
    intron_num: int,
    # gained_muts: int,
    # lost_muts: int
    # frameshift_num_diff: int,
    # frame_diff: int
    lost_muts: int,
    corr_frame: int
) -> float:
    """Calculates the score for intron layout prioritization"""
    # return coverage_score ** (intron_num * ((2.5 ** gained_muts) / (2 ** lost_muts)))
    # return coverage_score ** (intron_num / (2 ** lost_muts))
    frame_modif: int = int(corr_frame != 0)
    # return coverage_score ** (intron_num / ( frameshift_num_diff * (2 ** frame_diff)))
    power_numerator: int = (5 ** frame_modif) #intron_num * (5 ** frame_modif)
    power: float = power_numerator / (2 ** lost_muts)
    return coverage_score ** power


def introduce_intron_sequences(
    ref_seq: str, 
    query_seq: str, 
    introns: List[str] 
    # leading_seq: Optional[int] = 0, 
    # trailing_seq: Optional[int] = 0
) -> Tuple[str, str, List[Tuple[int, int]], bool]:
    """
    Replaces all precise intron deletion instances with lowercase intron sequences
    from the list
    """
    intron_coords: List[Tuple[int, int]] = []
    intron_start: int = 0
    intron_end: int = 0
    prev_coding: bool = False
    prev_del: bool = False
    intron_started: bool = False
    for i, (r, q) in enumerate(zip(ref_seq, query_seq)):
        curr_lower: bool = r.islower()
        curr_coding: bool = r.isupper() or r == '-'
        intron_del: bool = q == '>'
        if (intron_del or curr_lower and prev_coding) and not intron_started:
            intron_start = i
            intron_started = True
        elif (not intron_del and prev_del or curr_coding and not prev_coding) and intron_started:
            intron_end = i
            intron_coords.append((intron_start, intron_end))
            intron_started = False
        prev_coding: bool = curr_coding
        prev_del: bool = intron_del
    ## TODO:
    ## 1) Coordinates obviously do not correspond to actual sequence after ungabunga
    ## 2) Coordinates o
    upd_ref_seq, upd_query_seq = '', ''
    prev_intron_stop: int = 0
    if len(intron_coords) != len(introns):
        raise Exception(
            'ERROR: Number of intron intervals does not correspond to number '
            'of itnron sequences provided'
        )
    last_intron: int = len(introns) - 1
    subexon_coords: List[Tuple[int, int]] = []
    subexon_start: int = 0
    subexon_end: int = len(ref_seq)
    for j, (start, end) in enumerate(intron_coords):
        upd_ref_seq += ref_seq[prev_intron_stop:start]
        upd_query_seq += query_seq[prev_intron_stop:start].upper()
        subexon_end = len(upd_ref_seq)
        subexon_coords.append((subexon_start, subexon_end))
        intron_seq: str = introns[j]
        gap_seq: str = '-' * len(intron_seq)
        upd_ref_seq += gap_seq
        upd_query_seq += intron_seq.lower()
        subexon_start = len(upd_ref_seq)
        if j == last_intron:
            upd_ref_seq += ref_seq[end:]
            upd_query_seq += query_seq[end:].upper()
        prev_intron_stop = end
    subexon_coords.append((subexon_start, len(upd_ref_seq)))
    upd_ref_seq = upd_ref_seq.upper()
    frameshifted_aln: bool = upd_query_seq.count(' ') % 3 != 0
    upd_query_seq = upd_query_seq.replace(' ', '-')
    
    return upd_ref_seq, upd_query_seq, subexon_coords, frameshifted_aln


def get_all_indels(seq: str) -> List[int]:
    """
    """
    all_indels: List[int] = []
    prev_indel: bool = False
    indel_len: int = 0
    for i in seq:
        if i == '-':
            if not prev_indel:
                prev_indel = True
            indel_len += 1
        elif i != '-' and prev_indel:
            all_indels.append(indel_len)
            prev_indel = False
            indel_len = 0
    return all_indels


def strip_intron_deletion_runs(ref: str, query: str) -> Tuple[str, str]:
    """
    Removes symbols corresponding to precise intron deletions from both reference 
    and query lines of CESAR2.0 output
    """
    upd_ref, upd_query = '', ''
    for r,q in zip(ref, query):
        if r == '>':
            continue
        upd_ref += r
        upd_query += q
    return upd_ref, upd_query


def frameshift_in_cesar_data(ref: str, query: str) -> bool:
    """
    Given reference and query fields from CESAR2 output, 
    calculates the total frame phase difference between the sequences 
    """
    if '>' in ref:
        ref, query = strip_intron_deletion_runs(ref, query)
    indels_in_ref: str = ref.replace(GAP_CODON, '').count('-')
    indels_in_query: str = query.replace(GAP_CODON, '').count('-')
    # print(f'{indels_in_ref=}, {indels_in_query=}')
    return abs(indels_in_ref - indels_in_query) % 3 != 0 


class IntronGainChecker:
    __slots__ = (
        'portions', 'min_intron_gain_score', 
        # 'min_intron_prob_gapped', 'min_intron_prob_ungapped', 
        'min_intron_prob_trusted', 'min_intron_prob_supported', 'min_intron_prob_unsupported',
        'max_intron_num', 'lower_splice_threshold', 'logger',
        'rel_exon_coords', 'abs_exon_coords', 'exon2portion', 'intron2phase',
        'introns_gained', 'intron_meta'
    )

    def __init__(
        self,
        cesar_output: List[RawCesarOutput],
        min_intron_gain_score: float,
        # min_intron_prob_gapped: float,
        # min_intron_prob_ungapped: float,
        min_intron_prob_trusted: float,
        min_intron_prob_supported: float,
        min_intron_prob_unsupported: float,
        max_intron_num: int,
        logger: Logger
    ) -> None:
        self.portions: List[RawCesarOutput] = sorted(cesar_output, key=lambda x: x.exons)
        self.rel_exon_coords: Dict[int, Coords] = {}
        self.abs_exon_coords: Dict[int, Coords] = {}
        self.min_intron_gain_score: float = min_intron_gain_score
        # self.min_intron_prob_gapped: float = min_intron_prob_gapped
        # self.min_intron_prob_ungapped: float = min_intron_prob_ungapped
        # self.lower_splice_threshold: float = min(min_intron_prob_gapped, min_intron_prob_ungapped)
        self.min_intron_prob_trusted: float = min_intron_prob_trusted
        self.min_intron_prob_supported: float = min_intron_prob_supported
        self.min_intron_prob_unsupported: float = min_intron_prob_unsupported
        self.lower_splice_threshold: float = min(
            min_intron_prob_trusted, min_intron_prob_supported, min_intron_prob_unsupported
        )
        self.max_intron_num: int = max_intron_num
        self.exon2portion: Dict[int, int] = {}
        self.intron2phase: Dict[int, int] = {}
        self.introns_gained: Dict[int, List[PotExGroup]] = {}
        self.intron_meta: Dict[int, List[IntronMeta]] = {}
        self.logger: Logger = logger

        for i, aln_portion in enumerate(self.portions):
            curr_aln_len: int = len(aln_portion.reference)
            curr_exon: int = min(aln_portion.exons)
            exon_start_encountered: bool = False
            rel_exon_start, rel_exon_stop = None, None
            for j, base in enumerate(aln_portion.reference):
                # j: int = k + prev_aln_len
                query_base: str = aln_portion.query[j]
                if is_symbol(base):
                    if not exon_start_encountered:
                        # continue
                        exon_start_encountered = True
                        rel_exon_start: int = j
                if base in (' ', '>') or j == curr_aln_len - 1:
                    if not exon_start_encountered:
                        continue
                    rel_exon_stop: int = j + int(j == curr_aln_len - 1 and base.isalpha())
                    if rel_exon_start is None:
                        raise RuntimeError(f'Exon start missing for exon {curr_exon}')
                    self.rel_exon_coords[curr_exon] = Coords(rel_exon_start, rel_exon_stop)
                    abs_exon_start: int = self._abs_coord(aln_portion, rel_exon_start)
                    abs_exon_stop: int = self._abs_coord(aln_portion, rel_exon_stop)
                    self.abs_exon_coords[curr_exon] = Coords(abs_exon_start, abs_exon_stop)
                    intron_phase: int = self._get_query_intron_phase(aln_portion, curr_exon)
                    ## record intron phase
                    if base == '>':
                        ## mark deleted introns' phases as negative numbers for certain checks
                        self.intron2phase[curr_exon] = -intron_phase if intron_phase % 3 else -3
                    else:
                        self.intron2phase[curr_exon] = intron_phase
                    exon_start_encountered = False
                    rel_exon_start, rel_exon_stop = None, None
                    self.exon2portion[curr_exon] = i
                    curr_exon += 1
        ## TODO: Simple copy-paste from the processed_segment.py for now
        ## Adjust!

    def _to_log(self, msg: str, level: str = 'info'):
        """
        Adds the message to the log
        """
        if not hasattr(self, 'logger'):
            return
        getattr(self.logger, level)(msg)

    def record_intron_gains(self) -> None:
        """
        """
        for portion in self.portions:
            for exon in sorted(portion.exons):
                _ = self._check_for_intron_gain(portion, exon)

    def return_updated_cesar_output(self) -> List[Union[RawCesarOutput, CesarInput]]:
        """
        Returns the updated input data:
        * Exon groups with no potential exon gains found are returned as they
          were provided;
        * Exon groups suspected for intron gains are returned in the following way:
        1) Exons not affected by intron gains are split into intron gain are returned
           as RawCesarOutput projects with their data updated and alignment sequences
           clipped accordingly;
        2) Exons affected by intron gain are returned as CesarInput object to
           be subjected for additional round of alignment
        """
        out_list: List[Union[RawCesarOutput, CesarInput]] = []
        total_aln_sum: int = 0
        for l, portion in enumerate(self.portions):
            exons: List[int] = portion.exons
            ## case 1: none of the exons in this portion were affected by intron gain;
            ## return the original alignment chunk as-is
            if not any(x in self.introns_gained for x in exons):
                out_list.append(portion)
                total_aln_sum += len(portion.reference)
                continue
            ## case 2: alignment contains a single reference exon that was
            ## affected by intron gain; extract query subexons, flip reference
            ## and query sequences and produce a CesarInput object
            ## TODO: Devise a container class
            split_exons: List[int] = [x for x in exons if x in self.introns_gained]
            full_sequence_len: int = len(portion.reference)
            ref_seqs: Dict[int, str] = {}
            query_seqs: Dict[int, str] = {}
            exon_seqs: Dict[int, List[str]] = {}
            introns: Dict[int, List[str]] = {}
            upstream_seqs: Dict[int, Tuple[str, str]] = {}
            downstream_seqs: Dict[int, Tuple[str, str]] = {}
            needs_realignment: Dict[int, bool] = {}
            el2exons: Dict[int, List[int]] = {}
            exon2orig_ref: Dict[int, str] = {}
            exon2orig_query: Dict[int, str] = {}
            exon2orig_id: Dict[int, float] = {}
            # exon2extra_seq: Dict[int, Tuple[int, int]] = {}
            el2ref_phases: Dict[int, Tuple[int, int]] = {}
            streaks: List[List[int]] = binary_streaks(exons, split_exons)
            last_streak: int = len(streaks) - 1
            el_counter: int = 0
            intron_meta: Dict[int, IntronMeta] = {}
            for j, streak in enumerate(streaks):
                if any(x in split_exons for x in streak): ## TODO: split exons streak are now expected to consist of one exon only; simplify the code
                    last_exon: int = len(streak) - 1
                    for k, exon in enumerate(streak):
                        intron_meta[exon] = self.intron_meta[exon]
                        exon_start, exon_end = self.rel_exon_coords[exon].tuple()
                        orig_ref: str = portion.reference[exon_start:exon_end]
                        orig_query: str = portion.query[exon_start:exon_end]
                        orig_id: float = fast_seq_id(orig_ref, orig_query)
                        exon2orig_ref[exon] = orig_ref
                        exon2orig_query[exon] = orig_query
                        exon2orig_id[exon] = orig_id
                        _ref_seqs: List[str] = []
                        _exon_seqs: List[str] = []
                        _introns: List[str] = []
                        for subexon_coords in self.introns_gained[exon]:
                            ## TODO:
                            ## 1) think how to organize output (list of ExtraRound objects?)
                            ## 2) How to set the element counter for each candidate group
                            first_phase: int = (3 - self._get_frame(exon-1)) % 3
                            last_phase: int = self._get_frame(exon)
                            prev_phase: int = first_phase
                            next_phase: int = last_phase
                            subexon_seqs: List[str] = [
                                strip_noncoding(portion.query[slice(*x[:2])])
                                for x in subexon_coords
                            ]
                            # subexon_seqs: List[str] = [
                            #     portion.query[slice(*x[:2])] for x in subexon_coords
                            # ]
                            last_seq: int = len(subexon_seqs) - 1
                            formatted_subexons: List[str] = []
                            for i, seq in enumerate(subexon_seqs):
                                if i == 0:
                                    prev_streak_affected: bool = j > 0 and any(x in split_exons for x in streaks[j-1])
                                    if j == 0 or prev_streak_affected:
                                        up_start_coord: int = 0 if j == 0 else self.rel_exon_coords[exon - 1].stop
                                        up_end_coord: int = subexon_coords[i][0]
                                        upstream_ref: str = portion.reference[up_start_coord:up_end_coord]
                                        upstream_query: str = portion.query[up_start_coord:up_end_coord]
                                    else:
                                        upstream_ref: str = ''
                                        upstream_query: str = ''
                                upstream_seqs[el_counter] = (upstream_ref, upstream_query)
                                if i == last_seq and j == last_streak:
                                    down_start_coord: int = subexon_coords[i][1]
                                    down_end_coord: int = full_sequence_len
                                    downstream_ref: str = portion.reference[down_start_coord:down_end_coord]
                                    downstream_query: str = portion.query[down_start_coord:down_end_coord]
                                else:
                                    downstream_ref: str = ''
                                    downstream_query: str = ''
                                downstream_seqs[el_counter] = (downstream_ref, downstream_query)
                                seq_len: int = len(seq)
                                stripped_seq_len: int  = len(strip_noncoding(seq))
                                # if i == last_seq:
                                #     next_phase = last_phase
                                # else:
                                # phase_shift_prev: int = seq[:prev_phase].count('-')
                                # next_phase = (stripped_seq_len - prev_phase + phase_shift_prev) % 3
                                next_phase = (stripped_seq_len - prev_phase) % 3
                                formatted_seq: str = (
                                    seq[:prev_phase].lower() +
                                    seq[prev_phase:seq_len-next_phase] +
                                    seq[seq_len-next_phase:].lower()
                                )
                                formatted_subexons.append(strip_noncoding(formatted_seq))
                                prev_phase = (3 - next_phase) % 3
                            intron_sequences: List[str] = [
                                strip_noncoding(
                                    portion.query[subexon_coords[x-1][1]:subexon_coords[x][0]]
                                )
                                for x in range(1, len(subexon_coords))
                            ]
                            target_seq: str = strip_noncoding(
                                portion.reference[subexon_coords[0][0]:subexon_coords[-1][1]]
                            )
                            target_left_phase: int = (3 - self.intron2phase.get(exon-1,0)) % 3
                            target_right_phase: int = self.intron2phase.get(exon, 0)
                            _ref_seqs.append(target_seq)
                            _exon_seqs.append(formatted_subexons)
                            _introns.append(intron_sequences)
                        ref_seqs[el_counter] = _ref_seqs
                        exon_seqs[el_counter] = _exon_seqs
                        introns[el_counter] = _introns
                        needs_realignment[el_counter] = True
                        el2exons[el_counter] = [exon]
                        el2ref_phases[el_counter] = (target_left_phase, target_right_phase)
                        el_counter += 1
                else:
                    if j == 0:
                        untouched_seq_start: int = 0
                    else:
                        untouched_seq_start: int = self.rel_exon_coords[min(streak) - 1].stop
                    if j == len(streaks) - 1:
                        untouched_seq_end: int = full_sequence_len
                    else:
                        untouched_seq_end: int = self.rel_exon_coords[max(streak) + 1].start
                    untouched_ref_seq: str = portion.reference[untouched_seq_start:untouched_seq_end]
                    untouched_query_seq: str = portion.query[untouched_seq_start:untouched_seq_end]
                    ref_seqs[el_counter] = untouched_ref_seq
                    query_seqs[el_counter] = untouched_query_seq
                    needs_realignment[el_counter] = False
                    el2exons[el_counter] = streak
                    el_counter += 1
            output_instance: ExtraRoundCesarInput = ExtraRoundCesarInput(
                ref_seqs, query_seqs, exon_seqs, introns,
                upstream_seqs, downstream_seqs, el2ref_phases,
                needs_realignment, el2exons, 
                exon2orig_ref, exon2orig_query, exon2orig_id,
                intron_meta
                #exon2extra_seq#, all_subexons
            )
            out_list.append(output_instance)
            total_aln_sum += len(portion.reference)
        return out_list


    def _abs_coord(self, portion: RawCesarOutput, base: int) -> int:
        """
        Given the relative base coordinate in the alignment chunk,
        return its absolute coordinate in the query genome
        """
        if base is None:
            return None
        strand: bool = portion.strand
        gap_num: int = portion.query[:base].count('-')
        abs_start, abs_stop = portion.start, portion.stop
        return (abs_start + base - gap_num) if strand else (abs_stop - base + gap_num)

    def _rel_coord(
        self,
        portion: RawCesarOutput,
        base: int
    ) -> int:
        """
        """
        abs_start, abs_stop = portion.start, portion.stop
        strand: bool = portion.strand
        unadj_coord: int = (base - abs_start if strand else abs_stop - base)
        if unadj_coord == 0:
            return unadj_coord
        adj_coord: int = 0
        for i, n in enumerate(portion.query):
            if adj_coord == unadj_coord:
                break
            if n.isalpha():
                adj_coord += 1
        return i

    def _get_query_intron_phase(
        self, portion: RawCesarOutput, intron: int
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
                portion.reference[prev_exon_start:prev_exon_end]
                # portion.query[prev_exon_start:prev_exon_end]
            )
        last_exon_length: int = len(total_seq)
        # print(f'{intron=}, {total_seq=}')
        # print(f'{intron=}, {merged_exons=}, {last_exon_length=}, {prev_intron_phase=}')
        # print(f'{portion.reference[prev_exon_start:prev_exon_end]=}')
        # print(f'{strip_noncoding(portion.reference[prev_exon_start:prev_exon_end])=}')
        # print(f'{portion.query[prev_exon_start:prev_exon_end]=}')
        # print(f'{strip_noncoding(portion.query[prev_exon_start:prev_exon_end])=}')
        return (last_exon_length - 3 + prev_intron_phase) % 3

    def _check_for_intron_gain(self, portion: RawCesarOutput, exon: int) -> Any:
        """
        """
        if portion.was_not_aligned:
            return
        start, stop = self.rel_exon_coords[exon].tuple()
        abs_start, abs_stop = self.abs_exon_coords[exon].tuple()
        strand: bool = portion.strand
        ## define frame phase for both donor and acceptor sides
        acc_phase: int = self.intron2phase.get(exon - 1, 0)
        donor_phase: int = self.intron2phase.get(exon, 0)
        adj_start: int = start - int(strand)
        # adj_stop: int = stop + int(strand)
        adj_stop: int = stop + int(not strand)
        alt_accs: Dict[int, float] = {}

        for alt_acc, alt_acc_prob in portion.spliceai_sites['acceptor'].items():
            alt_acc_: int = self._rel_coord(portion, alt_acc)# + int(strand))
            internal_site: bool = adj_start <= alt_acc_ < adj_stop
            exceeds_threshold: bool = alt_acc_prob >= self.lower_splice_threshold
            if internal_site and (exceeds_threshold or alt_acc_ == adj_start):
                alt_accs[alt_acc_] = alt_acc_prob
        if start not in alt_accs:
            alt_accs[start] = 0.0
        alt_donors: Dict[int, float] = {}
        for alt_donor, alt_donor_prob in portion.spliceai_sites['donor'].items():
            alt_donor_: int = self._rel_coord(portion, alt_donor)# - int(not strand))
            internal_site: bool = adj_start < alt_donor_ <= adj_stop
            exceeds_threshold: bool = alt_donor_prob >= self.lower_splice_threshold
            if internal_site and (exceeds_threshold or alt_donor == adj_stop):
                alt_donors[alt_donor_] = alt_donor_prob
        if stop not in alt_donors:
            alt_donors[stop] = 0.0
        potential_exons: List[Tuple[int, int, float, float]] = []
        # print(f'{alt_accs=}, {alt_donors=}')
        ## now, get all the plausible acceptor-donor pairs
        for acc, acc_prob in sorted(alt_accs.items()):
            if acc != start:
                acc = self._adjust_splice_site(portion, acc, True)
            for donor, donor_prob in sorted(alt_donors.items()):
                if donor != stop:
                    donor = self._adjust_splice_site(portion, donor, False)
                if acc >= donor:
                    continue
                orig_seq: str = portion.query[slice(*self.rel_exon_coords[exon].tuple())].replace('-', '')
                mod_seq: str = portion.query[acc:donor].replace('-', '')
                # print(f'{portion.query[slice(*self.rel_exon_coords[exon].tuple())]=}, {mod_seq=}, {exon=}')
                if orig_seq and not mod_seq:
                    continue
                # if len(portion.query) >= 3 and len(mod_seq) < 3:
                #     continue
                new_group: Tuple[int, int, float, float] = (
                    acc, donor, acc_prob, donor_prob
                )
                potential_exons.append(new_group)
        # print(f'{potential_exons=}')

        ## all exon groups should start with the original acceptor site;
        ## add respective candidates as primers to
        potential_exon_groups: List[PotentialExonGroup] = []
        for pot_exon in potential_exons:
            if pot_exon[0] != start:
                continue
            pot_acc, pot_donor = pot_exon[:2]
            new_group: PotentialExonGroup = PotentialExonGroup()
            new_group.add(pot_exon)
            new_group.acceptor_phases.append((3 - acc_phase) % 3)
            new_group.donor_phases.append(
                (
                    len(strip_noncoding(portion.query[pot_acc:pot_donor])) - new_group.acceptor_phase()
                ) % 3
            )
            potential_exon_groups.append(new_group)
        potential_exon_groups.sort(key=lambda x: x.group_donor())

        closest_donor: int = potential_exon_groups[0].group_donor()
        # print(f'{potential_exons=}, {exon=}')
        for i, pot_exon in enumerate(potential_exons):
            acc, donor, acc_prob, donor_prob = pot_exon
            ## exons starting with an original acceptor site have been
            ## already added as group starters; dismiss them
            if acc == start:
                continue
            ## for each active group, check if the new exon can be added to it
            if acc <= closest_donor:
                self._to_log(f'SKIPPING: Acceptor {acc} (exon {exon}) lies upstream to the most upstream known donor ({closest_donor})')
                continue
            new_groups: List[PotentialExonGroup] = []
            groups_to_remove: List[PotentialExonGroup] = []
            for group in potential_exon_groups:
                group_acc: int = group.group_acc()
                group_donor: int = group.group_donor()
                last_subexon: bool = group_donor == stop
                ## exons cannot overlap each other
                if acc <= group_donor:
                    self._to_log(f'SKIPPING: Acceptor {acc} (exon {exon}) lies upstream to the group donor ({group_donor})')
                    break
                if len(portion.query[group_donor:acc].replace('-', '')) <= MIN_INTRON_LENGTH:
                    self._to_log(f'SKIPPING: SKIPPING: Acceptor {acc} and group donor ({group_donor}) result in zero base long intron for exon {exon}')
                    continue
                new_intron_len: int = acc - group_donor
                ## introduced introns cannot be shorter than the set minimal lenght
                if new_intron_len <= MIN_INTRON_LENGTH:
                    self._to_log(f'SKIPPING: Intron at {group_donor}-{acc} (exon {exon}) is too short (mininal length: 30 bp)')
                    continue
                ## introduced introns must have sufficient SpliceAI support;
                ## intron location in the gapped alignment region might serve
                ## as a stipulation in terms of minimal probability
                has_gap: bool = self._intron_gap_in_seq(portion, group_donor, acc)
                # _, nonsense_lost, init_frame, corr_frame = self._candidate_intron_summary(
                #     portion, [(group_donor, acc)], exon
                # )
                # has_muts: bool = bool(nonsense_lost) or (init_frame and corr_frame < )
                has_muts: bool = self._has_mutations(
                    portion, exon, group_donor, acc
                )
                both_evidence_sources: bool = has_gap and has_muts
                any_evidence_source: bool = has_gap or has_muts
                self._to_log(f'Evidence for intron at {group_donor}-{acc} in exon {exon}: gap support - {has_gap}, mutation support - {has_muts}')
                ## most trusted intron candidates: have both mutations and alignment gaps
                if both_evidence_sources and (
                    (acc_prob < self.min_intron_prob_trusted) or 
                    (group.group_donor_prob() < self.min_intron_prob_trusted and not last_subexon)
                ):
                    self._to_log(
                        f'SKIPPING: Gap- and  mutation-supported intron at {group_donor}-{acc} (exon {exon}) has insufficient SpliceAI support'
                    )
                    continue
                ## second line introns: have either mutations and assembly gaps
                if any_evidence_source and (
                    (acc_prob < self.min_intron_prob_supported) or 
                    (group.group_donor_prob() < self.min_intron_prob_supported and not last_subexon)
                ):
                    self._to_log(
                        f'SKIPPING: Gap-/mutation-supported intron at {group_donor}-{acc} (exon {exon}) has insufficient SpliceAI support'
                    )
                    continue
                ## third line introns: no mutation or gap support but the boundaries are still SpliceAI-supported
                if not any_evidence_source and (
                    (acc_prob < self.min_intron_prob_unsupported) or 
                    (group.group_donor_prob() < self.min_intron_prob_unsupported and not last_subexon)
                ):
                    self._to_log(
                        f'SKIPPING: Gap- and mutation-independent intron at {group_donor}-{acc} (exon {exon}) has insufficient SpliceAI support'
                    )
                    continue
                # if not has_gap and (
                #     (acc_prob < self.min_intron_prob_ungapped) or 
                #     (group.group_donor_prob() < self.min_intron_prob_ungapped and not last_subexon   )
                # ):
                #     print(
                #         f'SKIPPING: Gap-independent intron at {group_donor}-{acc} (exon {exon}) has insufficient SpliceAI support'
                #     )
                #     continue
                # if has_gap and (
                #     (acc_prob < self.min_intron_prob_gapped) or 
                #     (group.group_donor_prob() < self.min_intron_prob_gapped and not last_subexon)
                # ):
                #     print(
                #         f'SKIPPING: Gap-guided intron at {group_donor}-{acc} (exon {exon}) has insufficient SpliceAI support'
                #     )
                #     continue
                ## introduced introns cannot split inframe stop codons
                if self._alternative_split_contains_stop(
                    portion, group_donor, acc, group.donor_phase()
                ):
                    self._to_log(f'SKIPPING: Split stop codon created by intron at {group_donor}-{acc} (exon {exon})')
                    continue
                ## abide to the maximum allowed intron number 
                if group.len() == (self.max_intron_num - 1) and donor != stop:
                    self._to_log('SKIPPING: Maximum intron number exceeded in the layout (exon {exon})')
                    groups_to_remove.append(group)
                    continue
                ## if the code has reached this point, the exon can extend the current group
                ## it it ends at original donor site, terminate the group
                self._to_log(f'Adding new intron: donor - {group_donor} ({group.group_donor_prob()}), acceptor - {acc} ({acc_prob})')
                modif_group: Tuple[Tuple[int, int, float, float]] = group.copy()
                intron_meta: IntronMeta = IntronMeta(
                    has_gap, has_muts, group.group_donor_prob(), acc_prob
                )
                modif_group.add(pot_exon)
                modif_group.intron_meta.append(intron_meta)
                modif_group.acceptor_phases.append(
                    (3 - new_group.donor_phase()) % 3
                )
                modif_group.donor_phases.append(
                    (
                        len(strip_noncoding(portion.query[acc:donor])) -
                        modif_group.acceptor_phase()
                    ) % 3
                )
                new_groups.append(modif_group)
            potential_exon_groups= [x for x in potential_exon_groups if x not in groups_to_remove]
            potential_exon_groups.extend(new_groups)
            potential_exon_groups.sort(key=lambda x: x.group_donor())
            closest_donor = potential_exon_groups[0].group_donor()

        ## for each group, find the number of mutations removed with the
        sorted_groups: List[PotExGroup] = []
        init_exon_len: int = len(strip_noncoding(portion.query[start:stop]))
        if len(potential_exon_groups) == 1:
            return
        group2score: Dict[int, float] = {}
        for j, group in enumerate(potential_exon_groups):
            if group.group_acc() != start or group.group_donor() != stop:
                continue
            if len(group.exon_tuples) == 1:
                continue
            ## get the number of mutations in the intronic sequence
            # mut_num: int = 0
            # splice_prob_sum: int = 0
            intron_len_sum: int = 0
            intron_num: int = 0
            gained_muts: int = 0
            lost_muts: int = 0
            for i in range(1, len(group.exon_tuples)):
                prev_exon: Tuple[int, int, float, float] = group.exon_tuples[i-1]
                next_exon: Tuple[int, int, float, float] = group.exon_tuples[i]
                intron_start: int = prev_exon[1]
                intron_end: int = next_exon[0]
                intron_len_sum += len(
                    strip_noncoding(
                        portion.query[intron_start:intron_end]
                    )
                )
                intron_num += 1
            gained_nonsense, init_frame, corr_frame = self._candidate_intron_summary(
                portion, group.exon_tuples, exon
            )
            if gained_nonsense:
                continue
            if corr_frame:# and corr_frame != init_frame:
                continue
            coverage_score: float = (init_exon_len - intron_len_sum) / init_exon_len
            if coverage_score < self.min_intron_gain_score: ## TODO: Placeholder!
                continue
            # print(f'{group=}')
            # for _e in group.exon_tuples:
            #     print(f'{portion.reference[slice(*_e[:2])]=}')
            #     print(f'{portion.query[slice(*_e[:2])]=}')
            group2score[j] = (corr_frame, coverage_score, intron_num)
        if not group2score:
            return
        self.introns_gained[exon] = [
            x.exon_tuples for i, x in enumerate(potential_exon_groups) if i in group2score
        ]
        self.intron_meta[exon]= [
            x.intron_meta for i, x in enumerate(potential_exon_groups) if i in group2score
        ]

    def _alternative_split_contains_stop(
        self, portion: RawCesarOutput, donor: int, acc: int, phase: int, central: int = None
    ) -> bool:
        """
        Check whether the newly selected splice sites create a nonsense mutation
        """
        if not phase:
            return False
        pre_split_counter: int = 0
        if donor > 1:
            for i in range(donor - 1, 0, -1):
                pre_split_counter += int(is_symbol(portion.reference[i]))
                if pre_split_counter == phase:
                    break
        else:
            i = 1
        post_split_counter: str = 0
        if acc < len(portion.query) - 1:
            for j in range(acc, len(portion.query)):
                post_split_counter += int(is_symbol(portion.reference[j]))
                if post_split_counter == 3 - phase - int(central is not None):
                    break
        else:
            j = len(portion.query) - 1
        # central_portion: int = ## TODO: Accommodate for ultrashort internal codons
        new_split_codon: str = strip_noncoding(portion.query[i:donor] + portion.query[acc:j])
        for codon in parts(new_split_codon, 3):
            if codon.upper() in STOPS:
                return True
        return False

    def _intron_gap_in_seq(
        self, portion: RawCesarOutput, start: int, end: int
    ) -> bool:
        """
        """
        return portion.reference[start:end].count(GAP_CODON) >= 10

    def _adjust_splice_site(
        self, portion: RawCesarOutput, base: int, acceptor: bool, reverse: bool = False
    ) -> int:
        """
        Adjusts splice site position in 1-based into zero-based coordinates,
        ignoring the indel symbols
        """
        strand: bool = portion.strand
        if acceptor and not strand:
            return base
        if not acceptor and strand:
            return base
        last_base: int = min(base, len(portion.query)-1)
        search_range: Iterable[int] = (
            range(last_base, len(portion.query)) if acceptor else range(last_base, -1, -1)
        )
        for i in search_range:
            n: str = portion.query[i]
            if n != '-':
                return i + 1 if acceptor else i - 1
        return base

    def _get_terminal_gaps(self, portion: RawCesarOutput, base: int, forward: True) -> int:
        """Gets the length of terminal gap symbols for an exon"""
        search_range: Iterable[int] = (
            range(base, len(portion.query)) if forward else range(base, -1, -1)
        )
        for i in search_range:
            n: str = portion.query[i]
            if n != '-':
                return abs(i - base)
        return len(portion.query) - base if forward else 0 


    def _merged_exon_streak(self, exon: int) -> Tuple[int]: ## TODO: Tested only on minimal examples; test further
        """Returns the iterable of exons merged in query for the given exons"""
        streak: List[int] = [1]
        for i, phase in self.intron2phase.items():
            if phase >= 0:
                if exon in streak:
                    return tuple(streak)
                elif i == exon:
                    return (i,)
                streak = []
            # streak.append(i)
            streak.append(i+1)
        if exon in streak:
            return tuple(streak)

    def _get_frame(self, intron: int) -> int:
        """
        Returns the correct frame for the given intron, taken intron deletions into account
        """
        if intron not in self.intron2phase:
            return 0
        phase: int = self.intron2phase[intron]
        if phase > 0:
            return phase
        else:
            # merged_exons: List[int] = self._merged_exon_streak(intron)
            return -phase % 3

    def _has_mutations(self, portion: RawCesarOutput, exon: int, start: int, end: int) -> bool:
        """
        Indicates whether exon alignment contains frameshifting mutations
        and whether potential intron sequence contains any nonsense mutations
        """
        exon_start, exon_end = self.rel_exon_coords[exon].tuple()
        ## get frameshifts, if any are found
        ref_seq: str = ''.join(
            [
                portion.reference[i] for i in range(exon_start, exon_end) if 
                portion.query[i].isupper() or portion.query[i] == '-'
            ]
        )
        query_seq: str = ''.join(
            [
                portion.query[i] for i in range(exon_start, exon_end) if 
                portion.query[i].isupper() or portion.query[i] == '-'
            ]
        )
        ref_del: str = ref_seq.count('-') % 3
        query_del: str = query_seq.count('-') % 3
        has_frameshift: bool = ref_del - query_del

        ## get inframe stops within the intron boundaries;
        ## get the absolute value of the leading intron if it is deleted
        acc_phase: int = (3 - abs(self.intron2phase.get(exon - 1, 0))) % 3
        # init_ref_seq: str = portion.reference[exon_start:exon_end]
        # init_query_seq: str = portion.query[exon_start:exon_end]
        pre_intron_seq: str = portion.query[exon_start:start].replace('-', '')
        # print(f'{pre_intron_seq=}')
        pre_intron_phase: int = (len(pre_intron_seq) - acc_phase) % 3
        split_nuc: str = pre_intron_seq[-pre_intron_phase:] if pre_intron_phase else ''
        intron_seq: str = split_nuc + portion.query[start:end].replace('-', '')
        # print(f'{intron_seq=}')
        has_nonsense_mut: bool = any(x.upper() in STOPS for x in parts(intron_seq, 3))
        # print(f'{has_frameshift=}, {has_nonsense_mut=}, {acc_phase=}, {pre_intron_phase=}')
        return bool(has_frameshift) or bool(has_nonsense_mut)

    def _candidate_intron_summary(
        self,
        portion: RawCesarOutput,
        subexons: List[Tuple[int, int, Optional[float], Optional[float]]],
        exon: int
    ) -> Tuple[int, int, int, int, bool]:
        """
        Summarizes the mutation and frame data for a resulting layout
        """
        first_subexon_start: int = subexons[0][0]
        last_subexon_end: int = subexons[-1][1]
        last_exon: bool = exon == max(self.rel_exon_coords.keys())

        ## get the sequence of the focatl exon
        exon_start, exon_end = self.rel_exon_coords[exon].tuple()
        ## get the query introns' phases
        exon_prev_phase: int = self.intron2phase.get(exon - 1, 0)
        exon_next_phase: int = self.intron2phase.get(exon, 0)
        ## at this point, the code does not care about flanking introns' deletion;
        ## use absolute values of their intron phases to keep the calculations consistent
        exon_left_offset: int = (3 - abs(exon_prev_phase)) % 3
        exon_right_offset: int = abs(exon_next_phase)
        ## get initial sequences; strip them of full codon deletions for convenience
        # init_ref_seq: str = portion.reference[exon_start:exon_end].replace(GAP_CODON, '')
        init_query_seq: str = portion.query[exon_start:exon_end].replace(GAP_CODON, '')
        corr_query_seq: str = ''.join([
            portion.query[i] for i in range(exon_start, exon_end) 
            if any(x[0] <= i < x[1] for x in subexons) or 
            i < first_subexon_start or i >= last_subexon_end
        ]).replace(GAP_CODON, '')
        ## get the initial and the corrected frame phases
        init_frame: int = len(init_query_seq.replace('-', '')) % 3
        corr_frame: int = (len(corr_query_seq.replace('-', '')) - exon_right_offset - exon_left_offset) % 3
        ## check the corrected exon frame for inframe stop codons
        query_cds: str = corr_query_seq.replace('-', '')
        query_cds_in_frame: str = query_cds[exon_left_offset:len(query_cds) - exon_right_offset]
        # query_split_stripped_codons: List[str] = parts(corr_query_seq.replace('-', ''), 3)
        query_split_stripped_codons: List[str] = parts(query_cds_in_frame, 3)
        stripped_split_codon_num: int = len(query_split_stripped_codons)
        gained_nonsense: int = sum(
            x.upper() in STOPS for x in query_split_stripped_codons[:stripped_split_codon_num-(int(last_exon))]
        )
        
        ## check whether the exon is a part of a merged exon run
        exon_streak: Tuple[int] = self._merged_exon_streak(exon)
        if len(exon_streak) > 1:
            ## if it is
            last_exon: bool = exon == max(self.rel_exon_coords.keys())
            first_ex: int = exon_streak[0]
            last_ex: int = exon_streak[-1]
            seq_start: int = self.rel_exon_coords[first_ex].start
            seq_end: int = self.rel_exon_coords[last_ex].stop
            prev_phase: int = self.intron2phase.get(first_ex - 1, 0)
            next_phase: int = self.intron2phase.get(last_ex, 0)
            ## introns flanking the streak must not be deleted, 
            ## otherwise the respective exons must be included in the merged exon streak
            if prev_phase < 0 or next_phase < 0:
                raise ValueError(
                    'Trying to access the phase value for a deleted intron'
                )
            left_offset: int = (3 - prev_phase) % 3
            right_offset: int = next_phase
            last_exon: bool = max(exon_streak) == max(self.rel_exon_coords.keys())
            corr_query_seq: str = ''.join([
                portion.query[i] for i in range(seq_start, seq_end) 
                if any(x[0] <= i < x[1] for x in subexons) or 
                i < first_subexon_start or i >= last_subexon_end
            ]).replace(GAP_CODON, '')
            query_cds: str = corr_query_seq.replace('-', '')
            query_cds_in_frame: str = query_cds[left_offset:len(query_cds)-right_offset]
            # query_split_stripped_codons: List[str] = parts(corr_query_seq.replace('-', ''), 3)
            query_split_stripped_codons: List[str] = parts(query_cds_in_frame, 3)
            stripped_split_codon_num: int = len(query_split_stripped_codons)
            run_gained_nonsense: int = sum(
                x.upper() in STOPS for x in query_split_stripped_codons[:stripped_split_codon_num-(int(last_exon))]
            )
            gained_nonsense = max(gained_nonsense, run_gained_nonsense)
        return gained_nonsense, init_frame, corr_frame

