#!/usr/bin/env python3

"""Contains executables used in the CESAR wrapper"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from _io import TextIOWrapper
# from CESAR_wrapper import *
from collections import defaultdict
from common import chain_extract_id, parts
from dataclasses import dataclass
from datetime import datetime
# from GLP_values import *
from modules.cesar_wrapper_constants import *
# from SeqOp import extract, intersection, reverse_complement
from modules.shared import (
    CONTEXT_SETTINGS, intersection, nn, reverse_complement
)
# from string_splitter import split_at_
from sys import stderr
from typing import Any, Dict, Iterable, List, Set, TextIO, Tuple, Union

import click
import ctypes

__author__ = 'Yury V. Malovichko'
__version__ = '0.5'
__year__ = '2023'

@dataclass
class Exon:
    """
    Minimal data class for exon block storage
    """
    __slots__ = ['num', 'start', 'stop']
    num: int
    start: int
    stop: int

    def length(self) -> int:
        return self.stop -  self.start

class ExonDict(dict):
    """
    An auxiliary class to store {exon_number : Exon} mapping
    """
    def max(self):
        if not super().keys():
            return None
        maxval: int = max(super().keys())
        return super().__getitem__(maxval)

    def min(self):
        if not super().keys():
            return None
        minval: int = min(super().keys())
        return super().__getitem__(minval)


@dataclass
class AnnotationEntry:
    """
    Minimal data class for progenitor transcript/gene properties' storage
    """
    __slots__ = (
        'name', 'chrom', 'start', 'stop', 'strand', 'exon_number', 'exons'
    )
    name: str
    chrom: str
    start: int
    stop: int
    strand: bool
    exon_number: int
    exons: ExonDict[int, Exon]


@dataclass
class Segment:
    """
    Stores data on LASTZ/HMM segments both before and after sequence extension
    """
    __slots__ = [
        'name', 'chrom',
        'start', 'stop', 'strand',
        'min_exon', 'max_exon', 'exons'
    ]
    name: str
    chrom: str
    start: int
    stop: int
    strand: bool
    min_exon: int
    max_exon: int
    exons: ExonDict[int, Exon]

    def __repr__(self):
        return '\t'.join(map(str, (self.__getattribute__(x) for x in self.__slots__)))


@dataclass
class Exon2BlockMapper:
    """
    Maps reference exons to alignment chain blocks; used for segment search space
    extension and exon loci elucidation
    """
    __slots__ = [
        'chainid', 'prob',
        'tstart', 'tstop', 'qstart', 'qstop',
        'tchrom', 'qchrom', 'tstrand', 'qstrand',
        'blocks', #'e2b', 'fe2b',
        'e2c', 'out_of_chain', 'missing', 'gap_located',
        'dangling', 'spanning_chain'#, 'aa_qual'
    ]
    chainid: str
    # prob: float
    tstart: int
    tstop: int
    qstart: int
    qstop: int
    tchrom: str
    qchrom: str
    tstrand: bool
    qstrand: bool
    blocks: Dict[str, Tuple[int, int]]
    e2c: Dict[int, Tuple[int]]
    out_of_chain: Set[int]
    missing: Set[int]
    gap_located: Set[int]
    dangling: Set[int]
    spanning_chain: bool
    # aa_qual: Dict[int, bool]

    def get_exon_coords(self, e: int) -> Tuple[int, int]:
        # starts: List[int] = [self.blocks[x[0]][2] for x in self.e2b[e]]
        # stops: List[int] = [self.blocks[x[0]][3] for x in self.e2b[e]]
        # start, stop = sorted([min(starts), max(stops)])
        # return start, stop
        # return self.e2c[e]
        return self.chainid, self.qchrom, *(self.e2c[e])

    def codirected(self) -> bool:
        return self.tstrand == self.qstrand

    def get_suitable_blocks(self, start: int, stop: int) -> List[Tuple[int]]:
        """
        Given a start and a stop coordinate in the query,
        return all the aligned blocks confined between them
        """
        out_list: List[Tuple[int]] = []
        for name, block in self.blocks.items():
            if '_' in name:
                continue
            q_start, q_stop = block[2:]
            if q_start < start or q_stop > stop:
                continue
            out_list.append(block)
        return out_list


def extract_introns(annot: AnnotationEntry) -> ExonDict:
    """
    Returns a dictionary of introns for a given AnnotationEntry object
    """
    output: ExonDict = ExonDict()
    for num in range(1, annot.exon_number):
        intron_start: int = annot.exons[num if annot.strand else num+1].stop + 1
        intron_stop: int = annot.exons[num+1 if annot.strand else num].start
        intron: Exon = Exon(num, intron_start, intron_stop)
        output[num] = intron
    return output


def define_max_space(
    blocks: List[Tuple[int]], min_size: int, max_size: int
) -> Tuple[int, int]:
    """
    For a given block set, select the largest interblock range so that
    min_size <= x <= max_size
    """
    start, stop = (None, None)
    if not blocks:
        return (start, stop)
    curr_max_size: int = 0
    for i, f_block in enumerate(blocks):
        for l_block in blocks[i+1:]:
            space_size: int = l_block[2] - f_block[3] - 2
            if space_size > curr_max_size and min_size <= space_size <= max_size:
                curr_max_size = space_size
                start = f_block[3] + 1
                stop = l_block[2] - 1
    return (start, stop)


@dataclass
class ProjectionCoverageData:
    """
    An auxiliary data class containing data on marginal exons and fragment boundaries
    Attributes are:
    * :min_exon: and :max_exon: are minimal and maximal exon numbers for which
      the respective exons were found to be covered by the chain
    * :first_exon: and :last_exon: are minimal and maximal exon numbers for
      which the respective exons are expected to be found in this fragment after
      its search space extrapolation;
    * :start: and :stop: are absolute boundary coordinates after search space
      extrapolation
    """
    __slots__ = [
        'min_exon', 'max_exon', 'first_exon', 'last_exon', 'start', 'stop'
    ]
    min_exon: int
    max_exon: int
    first_exon: int
    last_exon: int
    start: int
    stop: int


def consistent_coords(
    coords: Tuple[str, Union[int, None]], chrom_len: int
) -> Tuple[str, Union[int, None]]:
    """
    Modifies coordinates tuple so that start and end coordinates lie within
    the [0; chromosome length] interval
    """
    start_coord: int = nn(min(coords[2], chrom_len)) if coords[2] is not None else None
    stop_coord: int = nn(min(coords[3], chrom_len)) if coords[3] is not None else None
    return (*coords[:2], start_coord, stop_coord)

class ProjectionGroup:
    """
    A projection-defining object comprising of a chain mapper and/or a LASTZ segment
    """
    __slots__ = [
        'mapper', 'segment', 'annot', 'start', 'stop',
        'chrom_len', 'first_exon', 'last_exon',
        'chrom', 'chain', 'acc_flank', 'donor_flank', 'out_of_chain_exons',
        'max_space_size', 'spliceai_sites', 'strand', 'shifted_exons',
        'evidence_source', 'exon_expected_loci', 'exon_search_spaces',
        'exon_coords', 'missing', 'gap_located', 'discordant_cases',
        'alternative_loci', 'cesar_exon_grouping', 'assembly_gaps', 'v'
    ]
    def __init__(
        self,
        mapper: Exon2BlockMapper,
        segment: Segment,
        annot: AnnotationEntry,
        chain: str,
        start: int,
        stop: int,
        chrom_len: int,
        first_exon: int,
        last_exon: int,
        acc_flank: int,
        donor_flank: int,
        out_of_chain_exons: Set[int],
        max_space_size: int,
        spliceai_sites: Dict[str, List[int]]
    ) -> None:
        """
        An auxiliary class for defining exon coordinates and CESAR alignment groups
        """
        self.mapper: Exon2BlockMapper = mapper
        self.segment: Segment = segment
        self.annot: AnnotationEntry = annot
        self.chain: str = chain
        self.start: int = start
        self.stop: int = stop
        self.chrom_len: int = chrom_len
        self.first_exon: int = first_exon
        self.last_exon: int = last_exon
        self.acc_flank: int = acc_flank
        self.donor_flank: int = donor_flank
        self.out_of_chain_exons: Set[int] = out_of_chain_exons
        self.max_space_size: int = max_space_size
        self.spliceai_sites: Dict[str, List[int]] = spliceai_sites
        self.strand: bool = (
            (self.mapper.qstrand == self.annot.strand) if self.mapper
            else self.segment.strand
        )
        self.chrom: str = self.mapper.qchrom if self.mapper else self.segment.chrom
        self.shifted_exons: Dict[str, List[int]] = {
            x: [] for x in ['segment_up', 'segment_down', 'chain_up', 'chain_down']
        }
        self.evidence_source: Dict[int, int] = {}
        self.exon_expected_loci: Dict[int, Tuple[int]] = {}
        self.exon_search_spaces: Dict[int, Tuple[int]] = {}
        self.exon_coords: Dict[int, Tuple[str, int]] = {}
        self.missing: Set[int] = set()
        self.gap_located: Set[int] = {x for x in self.mapper.gap_located} if self.mapper else set()
        self.discordant_cases: List[int] = []
        self.alternative_loci: List[int] = []
        self.cesar_exon_grouping: List[List[int]] = []
        self.assembly_gaps: List[Tuple[str, int]] = []
        # self.v: bool = verbose
        self.run()

    def run(self) -> None:
        """
        Primary working method; defines the trusted evidence source and expected
        coordinates for each exon, then groups the exons into CESAR alignment
        groups
        """
        self.get_shifted_exons()
        self.reconcile_exon_sources()
        self.resolve_discordant_cases()
        self.group_exons()

    def get_shifted_exons(self) ->  None:
        """
        Detects exons shifted in their position between chain and segment data
        """
        if not self.mapper or not self.segment:
            return
        self.shifted_exons = dict(
            zip(
                self.shifted_exons.keys(),
                detect_shifted_exons(self.segment, self.mapper)
            )
        )

    def reconcile_exon_sources(self) -> None:
        """
        For each exon, defines its expected loci; multiple sources of evidences
        are reconciled, discordant cases are set aside for separate examination
        """
        exon_range: Iterable[int] = (
            range(self.first_exon, self.last_exon + 1) if self.strand else
            range(self.last_exon, self.first_exon - 1, -1)
        )
        for num in exon_range:
            click.echo(f'Analysing exon {num}')
            chain_data_reliable, seg_data_reliable = False, False
            exon_in_chain: bool = False
            if self.mapper:
                exon_in_chain = num not in self.mapper.missing
                exon_has_block: bool = exon_in_chain and num not in self.mapper.gap_located
                chain_not_shifted: bool = (
                    num not in self.shifted_exons['chain_up'] and
                    num not in self.shifted_exons['chain_down']
                )
                chain_data_reliable = exon_has_block and chain_not_shifted
                exon_in_chain = exon_in_chain and chain_not_shifted
            if self.segment:
                seg_locus: Exon = self.segment.exons.get(num, None)
                segment_not_shifted: bool = (
                    num not in self.shifted_exons['segment_up'] and
                    num not in self.shifted_exons['segment_down']
                )
                seg_data_reliable: bool = bool(seg_locus) and segment_not_shifted
                if seg_data_reliable:
                    seg_data_reliable = num in self.segment.exons
            if seg_data_reliable: ## segment already has a suitable exon locus
                click.echo(f'Exon {num} is a part of the segment')
                locus_length: int = seg_locus.length()
                e_length: Exon = self.annot.exons[num].length()
                if chain_data_reliable:
                    block_start, block_stop = self.mapper.get_exon_coords(num)[2:]
                    click.echo(
                        f'Adjusting exon {num} coordinates '
                        'by the respective chain block'
                    )
                    inter_block_locus: int = intersection(
                        seg_locus.start,
                        seg_locus.stop,
                        block_start,
                        block_stop
                    )
                    if inter_block_locus > 0:
                        click.echo(
                            f'Adjusting the search locus for exon {num}'
                        )
                        self.exon_coords[num] = (
                            self.chrom,
                            self.strand,
                            min(seg_locus.start, block_start),
                            max(seg_locus.stop, block_stop)
                        )
                        self.evidence_source[num] = 3 ## concordant chain and segment predictions
                    else:
                        click.echo(
                            'WARNING: LASTZ locus and chain block '
                            f'for exon {num} in segment do not intersect!'
                        )
                        click.echo(
                            'Downgrading to the segment locus'
                        )
                        self.discordant_cases.add(num)
                        # self.exon_coords[num] = (locus.start, locus.stop)
                else:
                    click.echo(
                        'No chain block available for '
                        f'exon {num}, using segment locus'
                    )
                    self.exon_coords[num] = (
                        self.chrom,
                        self.strand,
                        seg_locus.start,
                        seg_locus.stop
                    )
                    self.evidence_source[num] = 1 ## segment-based location prediction
            elif chain_data_reliable or exon_in_chain:
                click.echo(
                    f'Boundaries for exon {num} are '
                    'inferred from the corresponding chain block'
                )
                # print(f'{self.mapper.get_exon_coords(num)=}')
                self.exon_coords[num] = (
                    self.chrom,
                    self.strand,
                    *self.mapper.get_exon_coords(num)[2:]
                )
                self.evidence_source[num] = 2 ## chain-based location prediction
            else: ## exon has neither LASTZ nor chain evidence; mark it as missing
                click.echo(
                    f'Exon {num} has no defined locus for alignment'
                )
                ## if SpliceAi data are present, try finding an optimal locus
                ## confined between a donor and an acceptor site within the
                ## already defined search_locus_start and search_locus_stop
                # if self.spliceai_sites:
                #     spliceai_coords: List[Tuple[int]] = self._get_spliceai_coords(
                #         *self.exon_coords[num][2:]
                #     )
                #     for (sai_start, sai_stop) in spliceai_coords:
                #         if (sai_stop - sai_start) > self.annot.exons[num].length():
                #             click.echo(
                #                 f'Inferring coordinates for exon {num} '
                #                 'from the SpliceAI predictions'
                #             )
                #             self.exon_coords[num][2] = sai_start
                #             self.exon_coords[num][3] = sai_stop
                #             self.evidence_source[num] = 4 ## SPLICEAI evidence source
                #             continue
                self.missing.add(num)
                self.evidence_source[num] = 0 ## exon location is undefined prior to CESAR run

    def resolve_discordant_cases(self) -> None:
        """
        For exons with data from both LASTZ chained alignments and LASTZ segments
        present and contradicting each other, selects the more reliable evidence
        source
        """
        if not self.discordant_cases:
            return
        click.echo('Analysing the discordant cases')
        min_exon: int = min(
            list(filter(
                lambda x: x not in discordant_cases, self.exon_coords.keys())
            )
        )
        max_exon: int = max(
            list(filter(
                lambda x: x not in discordant_cases, self.exon_coords.keys())
            )
        )
        draft_start: int = (
            self.exon_coords[min_exon][0] if self.strand \
            else self.exon_coords[max_exon][0]
        )
        draft_stop: int = (
            self.exon_coords[max_exon][1] if self.strand \
            else self.exon_coords[min_exon][1]
        )
        for d in sorted(discordant_cases):
            ref_exon: Exon = self.annot.exons[d]
            segment_exon: Exon = self.segment.exons[d]
            chain_block: Tuple[int, int] = self.mapper.get_exon_coords(d)[2:]
            segment_inter: int = intersection(
                segment_exon.start,
                segment_exon.stop,
                draft_start,
                draft_stop
            )
            block_inter: int = intersection(
                chain_block[0],
                chain_block[1],
                draft_start,
                draft_stop
            )
            ## case 1: one evidence is clearly out of the draft transcript locus
            if segment_inter >= 0 and block_inter < 0:
                click.echo(
                    f'Preferring LASTZ evidence for exon {d} since chain '
                    'block lies further from the rest of the segment'
                )
                self.exon_coords[d] = (
                    self.chrom,
                    self.strand,
                    segment_exon.start,
                    segment_exon.stop
                )
                self.evidence_source[d] = 1
                # alternative_loci[d].append(chain_block)
            elif segment_inter < 0 and block_inter >= 0:
                click.echo(
                    f'Preferring chain evidence for exon {d} since segment '
                    'locus lies further from the current '
                )
                self.exon_coords[d] = (
                    self.chrom,
                    self.strand,
                    *chain_block
                )
                self.evidence_source[d] = 2
                # alternative_loci[d].append((segment_exon.start, segment_exon.stop))
            ## case 2: exon loci lie outside of the concordant locus
            ## go for the locus lying closer to the concordant  (draft) locus
            elif segment_inter < 0 and block_inter < 0:
                click.echo(
                    'Adding the proximal locus to the primary segment '
                    'and putting the other to alternative stash '
                    f'for exon {d}'
                )
                if d < min_exon:
                    dist_from_segment: int = (
                        draft_start - segment_exon.stop
                        if self.strand else
                        segment_exon.start - draft_stop
                    )
                    dist_from_block: int = (
                        draft_start - chain_block[1]
                        if self.strand else
                        chain_block[0] - draft_stop
                    )
                    self.exon_coords[d] = (
                        self.chrom,
                        self.strand,
                        *(
                            (segment_exon.start, segment_exon.stop) if
                            dist_from_segment < dist_from_block else
                            chain_block
                        )
                    )
                    self.alternative_loci[d].append(
                        self.chrom,
                        self.strand,
                        *(
                            chain_block if
                            dist_from_segment < dist_from_block else
                            (segment_exon.start, segment_exon.stop)
                        )
                    )
                else:
                    dist_from_segment: int = (
                        segment_exon.start - draft_stop
                        if self.strand else
                        draft_start - segment_exon.stop
                    )
                    dist_from_block: int = (
                        chain_block[0] - draft_stop
                        if self.strand else
                        draft_start - chain_block[1]
                    )
                    pick_segment: bool = dist_from_segment < dist_from_block
                    self.exon_coords[d] = (
                        self.chrom,
                        self.strand,
                        *(
                            (segment_exon.start, segment_exon.stop) if
                            pick_segment else chain_block
                        )
                    )
                    self.alternative_loci[d].append(
                        self.chrom,
                        self.strand,
                        *(
                            chain_block if pick_segment else
                            (segment_exon.start, segment_exon.stop)
                        )
                    )
                self.evidence_source[d] = 4 ## both loci
            else:
                ## first, intersect both LASTZ locus and block with exons having defined coordinates
                ## if any of these two intersect at least one exon, prioritize the other
                intersected_by_segment: Set[int] = set()
                intersected_by_block: Set[int] = set()
                for e in self.exon_coords:
                    if e == d or e in self.missing:
                        continue
                    if intersection(
                        segment_exon.start, segment_exon.stop,
                        self.exon_coords[e][0], self.exon_coords[e][1]
                    ) > 0:
                        intersected_by_segment.add(e)
                    if intersection(
                        chain_block[0], chain_block[1],
                        self.exon_coords[e][0], self.exon_coords[e][1]
                    ) > 0:
                        intersected_by_block.add(e)
                if intersected_by_block and intersected_by_segment:
                    raise Exception(
                        f'For exon {d}, both the chain block and the LASTZ hit '
                        'correspond to the already defined exon(s)'
                    )
                elif intersected_by_block:
                    click.echo(
                        f'Preferring LASTZ evidence for exon {d} since '
                        'chain block corresponds to another exon'
                    )
                    self.exon_coords[d] = (
                        self.chrom,
                        self.strand,
                        segment_exon.start,
                        segment_exon.stop
                    )
                    self.evidence_source[d] = 1
                elif intersected_by_segment:
                    click.echo(
                        f'Preferring chain evidence for exon {d} '
                        'since LASTZ locus corresponds to another exon'
                    )
                    self.exon_coords[d] = (
                        self.chrom,
                        self.strand,
                        *chain_block
                    )
                    self.evidence_source[d] = 2
                else: ## non-overlapping stray hits within the draft frame!
                    click.echo(
                        f'A complex case encountered for exon {d}; '
                        'setting chain block as primary coordinates, adding '
                        'the LASTZ locus to the alternative stash'
                    )
                    ## by default, we trust the chain data more
                    self.exon_coords[d] = (
                        self.chrom,
                        self.strand,
                        *chain_block
                    )
                    self.alternative_loci[d].append(
                        (
                            self.chrom,
                            self.strand,
                            segment_exon.start,
                            segment_exon.stop
                        )
                    )
                    self.evidence_source[d] = 4

    def group_exons(self) -> None:
        """
        Groups exons for further exonwise CESAR alignment.
        The following logic applies:
        1) Exons with no defined alignment loci form a separate group with
            alignment coordinates defined by boundaries of neighboring exons with
            defined loci;
        1a) Consecutive exons with missing loci are grouped together;
        1b) Missing exons shorter than three bases are grouped with any
            neighboring exons, be it with or without coordinates;
        2) Exons with defined loci form separate groups by default;
        2a) Exons separated by less than 2*flank size + minimum intron size
            are grouped together;
        2b) Exons with 'dangling' ends (terminal portions longer than quarter
            exon length not covered by any chain block) are grouped with their
            neighbors from the respective side;
        The described procedure applies to both primary and
        alternative exon loci.
        """
        ## TODO: reusing self.exon_coords was a bad idea; think how to
        exon_range: Iterable[int] = (
            range(self.first_exon, self.last_exon + 1) if self.strand else
            range(self.last_exon, self.first_exon - 1, -1)
        )
        if not self.segment and self.mapper.spanning_chain:
            self.cesar_exon_grouping.append(list(exon_range))
            group_start: int = min(
                self.exon_coords[x][2] for x in self.exon_coords if self.exon_coords[x][2] != self.exon_coords[x][3]
            )
            group_stop: int = max(
                self.exon_coords[x][3] for x in self.exon_coords if self.exon_coords[x][2] != self.exon_coords[x][3]
            )
            for num in exon_range:
                self.exon_coords[num] = (self.chrom, self.strand, group_start, group_stop)
            return
        # prev_dangling: bool = False
        ## expected loci are 'raw' chain-derived coordinates;
        ## these will not be anyhow updated further
        self.exon_expected_loci = {
            k:consistent_coords(v, self.chrom_len) for k,v in self.exon_coords.items()
            if all(z is not None for z in v)
        }
        ## exon search spaces account for provided donor/acceptor flanks;
        ## these loci are defined only for exons with defined alignment coordinates
        ## and will be further slightly updated to account for potential flank
        ## intersection
        self.exon_search_spaces = {
            k:consistent_coords(v, self.chrom_len) for k,v in self.exon_coords.items()
            if all(z is not None for z in v)
        }
        # self.exon_search_spaces = self.add_flanks_to_exon_coords()
        ## missing exons unite those which lie out of the defined projection portion
        ## and those for which the original mapper failed to defined plausible coordinates
        all_missing: Set[int] = self.missing.union(self.out_of_chain_exons)

        ## TODO: The following procedure must be (nearly) overlap-proof:
        ## 1) Find overlaps between search spaces of all defined exons 
        ##     a) Do not consider gap-aligned exons at this point
        ##     b) Define shared coordinates for each such group; 
        ##     c) Define 
        for num in exon_range:
            vacant_group: bool = bool(self.cesar_exon_grouping) ## a shortcut to define if any group has been already added to the grouping list
            not_yet_grouped: bool = vacant_group and num not in self.cesar_exon_grouping[-1] ## a shortcut to define if an exon has already been added to the last group
            already_grouped: bool = vacant_group and not not_yet_grouped
            if already_grouped:
                if num in all_missing:
                    click.echo(f'Exon {num} has been already grouped; skipping')
                    continue
                else:
                    click.echo(
                        f'Exon {num} has been already grouped; '
                        f'checking for the trailing intron length'
                    )
            last_group_contains_missing: bool = vacant_group and any(
                map(lambda x: x in all_missing, self.cesar_exon_grouping[-1])
            )
            is_short: bool = self.annot.exons[num].length() < 3
            if num in all_missing: ## estimate expected loci for missing exons or group them with neighbors
                click.echo(f'Exon {num} has no defined coordinates')
                prevs: List[int] = [
                    x for x in range(self.first_exon, num) if x not in all_missing
                ]
                nexts: List[int] = [
                    x for x in range(num + 1, self.last_exon + 1) if
                    x not in all_missing
                ]
                prev: int = None if not prevs else max(prevs)
                next_: int = None if not nexts else min(nexts)
                next_neighbor: int = next_ if self.strand else prev
                if not prev and not next_: ## case 0: a single-exon segment with its only exon undefined; must not ever happen
                    self.exon_coords[num] = (self.start, self.stop)
                elif num == self.first_exon: ## case 1: first exon
                    click.echo(f'Exon {num} is the first exon')
                    # print(f'{self.strand=}, {next_=}, {self.exon_coords[next_]=}, {self.donor_flank=}, {self.acc_flank=}')
                    if self.strand:
                        self.exon_coords[num] = (
                            self.chrom,
                            self.strand,
                            nn(self.start),
                            nn(self.exon_coords[next_][2] - 1 - self.acc_flank)
                        )
                    else:
                        self.exon_coords[num] = (
                            self.chrom,
                            self.strand,
                            nn(
                                min(self.exon_coords[next_][3] + 1 + self.acc_flank, self.chrom_len)
                            ),
                            nn(min(self.stop, self.chrom_len))
                        )
                    # print(f'UNDEFINED: {self.exon_coords[num]=}')
                elif num == self.last_exon: # case 2: last exon
                    click.echo(f'Exon {num} is the last exon')
                    # print(f'{self.exon_coords[prev]=}, {self.strand=}, {self.start=}, {self.stop=}')
                    if self.strand:
                        self.exon_coords[num] = (
                            self.chrom,
                            self.strand,
                            nn(
                                min(self.exon_coords[prev][3] + 1 + self.donor_flank, self.chrom_len)
                            ),
                            nn(min(self.stop, self.chrom_len))
                        )
                    else:
                        self.exon_coords[num] = (
                            self.chrom,
                            self.strand,
                            nn(min(self.chrom_len, self.start)),
                            nn(
                                min(self.exon_coords[prev][2] - 1 - self.donor_flank, self.chrom_len)
                            )
                        )
                    # print(f'{self.exon_coords[num]=}')
                else: # case 3: middle exon
                    ## all the space between the neighboring defined loci will
                    ## be used as alignment locus
                    click.echo(f'Exon {num} is a middle exon')
                    if self.strand:
                        search_locus_start: int = nn(
                            self.start if not prev else
                            self.exon_coords[prev][3] + 1 + self.acc_flank
                        )
                        search_locus_stop: int = nn(
                            self.stop if not next_ else
                            self.exon_coords[next_][2] - 1 - self.donor_flank
                        )
                    else:
                        search_locus_start: int = nn(
                            self.start if not next_ else
                            self.exon_coords[next_][3] + 1 + self.donor_flank
                        )
                        search_locus_stop: int = nn(
                            self.stop if not prev else
                            self.exon_coords[prev][2] - 1 - self.acc_flank
                        )
                    search_locus_start = min(search_locus_start, self.chrom_len)
                    search_locus_stop = min(search_locus_stop, self.chrom_len)
                    self.exon_coords[num] = (
                        self.chrom, self.strand,
                        search_locus_start, search_locus_stop
                    )
                # introns[num] = (None, None) ## intron is undefined in this case

                self.exon_search_spaces[num] = (self.chrom, self.strand, None, None)
                if vacant_group:
                    prev_group_start: int = min(
                        (
                            self.exon_coords[x][2]
                            for x in self.cesar_exon_grouping[-1]
                            if x in self.exon_coords.keys()
                        )
                    ) #- (self.acc_flank if self.strand else self.donor_flank)
                    prev_group_stop: int = max(
                        (
                            self.exon_coords[x][3]
                            for x in self.cesar_exon_grouping[-1]
                            if x in self.exon_coords.keys()
                        )
                    ) #+ (self.donor_flank if self.strand else self.acc_flank)
                else:
                    prev_group_start, prev_group_stop = 0, 0

                this_group_start: int = self.exon_coords[num][2] #- (
                #     self.acc_flank if self.strand else self.donor_flank
                # )
                this_group_stop: int = self.exon_coords[num][3] #+ (
                #     self.donor_flank if self.strand else self.acc_flank
                # )
                # if vacant_group and prev_dangling: ## previous exon has defined borders but its 3'-end is not covered by a chain block
                #     click.echo(f'Previous exon has an unaligned 3\'-terminus; adding exon {num} to the last defined group')
                #     self.cesar_exon_grouping[-1].append(num)
                if last_group_contains_missing: ## last added group comprises of exons with undefined loci
                    click.echo(f'Previous group contains missing exons; adding exon {num} to the last defined group')
                    self.cesar_exon_grouping[-1].append(num)
                elif vacant_group and intersection(
                    *sorted((this_group_start, this_group_stop)),
                    *sorted((prev_group_start, prev_group_stop))
                ) > 0:
                    ## predicted borders intersect with the previous group;
                    ## in this case, the following exon is most likely to overlap
                    ## with the last group as well, so it is also added
                    click.echo(
                        f'Search space of exon {num} intersects with that of the '
                        f'previous groups; adding exon {num} to the last defined group'
                    )
                    self.cesar_exon_grouping[-1].append(num)
                    if next_neighbor and next_neighbor not in self.cesar_exon_grouping[-1]:
                        self.cesar_exon_grouping[-1].append(next_neighbor)
                # elif is_short:
                #     ## short exons are grouped with the next neighboring exon
                #     click.echo(
                #         f'Exon {num} is shorter than the defined length threshold; '
                #         'adding it to the last defined group'
                #     )
                #     if next_neighbor:
                #         if vacant_group and next_neighbor in self.cesar_exon_grouping[-1]:
                #             self.cesar_exon_grouping[-1].append(num)
                #         else:
                #             self.cesar_exon_grouping.append([num, next_neighbor])
                #     else:
                #         self.cesar_exon_grouping[-1].append(num)
                else:
                    self.cesar_exon_grouping.append([num])
                # print([x for e in self.cesar_exon_grouping[-1] for x in self.exon_coords[e]])
                # prev_dangling = False

                # next_exon: Union[int, None] = None
                # # print(f'{num=}, {self.first_exon=}, {self.last_exon=}, {self.strand=}')
                # if num != self.last_exon and self.strand:
                #     next_exon = num + 1
                # if num != self.first_exon and not self.strand:
                #     next_exon = num - 1
                # if next_exon in all_missing:
                #     next_exon = None
                # # next_exon: int = num + 1 if self.strand else num -1 
                # if next_exon is not None and next_exon not in self.cesar_exon_grouping[-1]:
                #     next_group_start, next_group_stop = sorted(self.exon_coords[next_exon][2:])
                #     if next_exon not in all_missing:
                #         next_group_start -= (self.acc_flank if self.strand else self.donor_flank)
                #         next_group_stop += (self.donor_flank if self.strand else self.acc_flank)
                #     print(f'{prev_group_start=}, {prev_group_stop=}, {this_group_start=}, {this_group_stop=}, {next_group_start=}, {next_group_stop=}')
                #     intersection_to_next: int = intersection(
                #         *sorted((this_group_start, this_group_stop)),
                #         next_group_start, next_group_stop
                #     )
                #     if intersection_to_next > 0:
                #         click.echo(
                #             f'Next exon {next_exon} is added to the current group due to '
                #             f'search space intersection {intersection_to_next}'
                #         )
                #         self.cesar_exon_grouping[-1].append(next_exon)

                group: List[int] = self.cesar_exon_grouping[-1]
                if len(group) == 1:
                    continue
                group_start: int = nn(
                    min(
                        min(self.exon_coords[x][2] for x in group if x in self.exon_coords),
                        self.chrom_len
                    )
                )
                group_stop: int = nn(
                    min(
                        max(self.exon_coords[x][3] for x in group if x in self.exon_coords),
                        self.chrom_len
                    )
                )
                for e in group:
                    self.exon_coords[e] = (
                        self.chrom, self.strand, group_start, group_stop
                    )
                # print([x for e in self.cesar_exon_grouping[-1] for x in self.exon_coords[e]])
                continue

            ## exon has defined coordinates; check for conditions under which
            ## exon should be grouped with its neighbors, place it a separate
            ## group of one exon otherwise
            click.echo(f'Exon {num} has a defined expected locus')
            # print(f'{self.exon_coords[num]=}')
            # print(f'{self.exon_search_spaces[num]=}')
            # self.exon_search_spaces[num] = self.exon_coords[num]
            if self.strand and num < self.last_exon:
                next_: int = num + 1
            elif not self.strand and num > self.first_exon:
                next_: int = num - 1
            else: ## a terminal exon reached; there is no next exon, so set a placeholder
                next_: int = 0
                # if not vacant_group or not_yet_grouped:
                #     self.cesar_exon_grouping.append([num])
                # continue ## nothing to do otherwise, and there's no intron
            prev_: int = num - 1 if self.strand else num + 1

            ## correct the raw coordinates:
            if num == self.first_exon and num not in self.gap_located:
                click.echo('Updating search locus coordinate for the first exon')
                orig_start: int = self.exon_search_spaces[num][2 if self.strand else 3]
                flanked_start: int = nn(
                    min(
                        orig_start - (self.acc_flank if self.strand else -self.acc_flank),
                        self.chrom_len
                    )
                )
                self.exon_search_spaces[num] = (
                    *self.exon_search_spaces[num][:2], flanked_start, self.exon_search_spaces[num][3]
                ) if self.strand else (
                    *self.exon_search_spaces[num][:3], flanked_start
                )
            if num == self.last_exon and num not in self.gap_located:
                click.echo(f'Updating search locus coordinate for the last exon; {self.chrom_len=}')
                # click.echo(f'Search space before adjustment: {self.exon_search_spaces[num]}')
                orig_stop: int = self.exon_search_spaces[num][3 if self.strand else 2]
                flanked_stop: int = nn(
                    min(
                        orig_stop + (self.donor_flank if self.strand else -self.donor_flank),
                        self.chrom_len
                    )
                )
                self.exon_search_spaces[num] = (
                    *self.exon_search_spaces[num][:3], flanked_stop
                ) if self.strand else (
                    *self.exon_search_spaces[num][:2], flanked_stop, self.exon_search_spaces[num][3]
                )
                # click.echo(f'Search space after adjustment: {self.exon_search_spaces[num]}')
            else:
                self.exon_search_spaces[num] = (
                    *self.exon_search_spaces[num][:2],
                    nn(min(self.exon_search_spaces[num][2], self.chrom_len)),
                    nn(min(self.exon_search_spaces[num][3], self.chrom_len))
                )

            ## reset grouping flags
            vacant_group: bool = bool(self.cesar_exon_grouping)
            not_yet_grouped: bool = vacant_group and num not in self.cesar_exon_grouping[-1]

            if next_ and next_ not in all_missing: # next exon has defined borders
                intron_start: int = self.exon_coords[num][3] + 1 + (
                    self.donor_flank if self.strand else self.acc_flank
                )
                intron_stop: int = self.exon_coords[next_][2] - 1 - (
                    self.acc_flank if self.strand else self.donor_flank
                )
                next_is_short: bool = self.annot.exons[next_].length() < 3
                # introns[num] = (intron_start, intron_stop)
                intron_length: int = intron_stop - intron_start + 1
                short_intron: bool = intron_length < MIN_INTRON_LENGTH
                if short_intron or is_short or next_is_short:# or dangling_down:
                    ## intron elimination is suspected, or current exon is
                    ## too short, or downstream terminus is not covered
                    if is_short:
                        click.echo(
                            f'Exon {num} is merged with {next_} due to '
                            'being shorter than 3 bp'
                        )
                    elif short_intron:
                        click.echo(
                            f'Exons {num} and {next_} are to be merged due to '
                            f'insuffiecient intron length ({intron_length})'
                        )

                    if not vacant_group or not_yet_grouped: ## append a whole new group
                        self.cesar_exon_grouping.append([num, next_])
                    else: ## current exon is already in the last group
                        self.cesar_exon_grouping[-1].append(next_)

                elif not vacant_group or not_yet_grouped:  ## no need for grouping the exons together
                    self.cesar_exon_grouping.append([num])

            else: ## next exon's coordinates are missing
                ## if 3'-end is dangling, grouping will be better resolved
                ## when processing the next exon
                # introns[num] = (None, None)
                if not vacant_group or not_yet_grouped:
                    self.cesar_exon_grouping.append([num])

            ## update raw coordinates; these are guaranteed to be defined for
            ## exons with defined loci
            # print(f'{next_=}, {next_ not in all_missing=}, {next_ not in self.gap_located=}, {num not in self.gap_located =}')
            if next_ and next_ not in all_missing:# and next_ not in self.gap_located and num not in self.gap_located:
                # print('Updating search locus coordinates for the two neighboring exons')
                raw_curr_exon_stop: int = self.exon_coords[num][3]
                raw_next_exon_start: int = self.exon_coords[next_][2]
                curr_exon_stop: int = self.exon_search_spaces[num][3]
                next_exon_start: int = self.exon_search_spaces[next_][2]
                # if curr_exon_stop >= next_exon_start:
                #     self.exon_search_spaces[num] = (
                #         *self.exon_search_spaces[num][:3], upd_flanked_stop
                #     )
                #     self.exon_search_spaces[next_] = (
                #         *self.exon_search_spaces[next_][:2], upd_next_flanked_start,
                #         self.exon_search_spaces[next_][3]
                #     )
                # else:
                # if curr_exon_stop < next_exon_start:
                # if curr_exon_stop >= next_exon_start:
                flanked_stop: int  = curr_exon_stop + (
                    self.donor_flank if self.strand else self.acc_flank
                )
                next_flanked_start: int = next_exon_start - (
                    self.acc_flank if self.strand else self.donor_flank
                )
                flanked_stop = nn(min(flanked_stop, self.chrom_len))
                next_flanked_start = nn(min(next_flanked_start, self.chrom_len))
                # print(f'{curr_exon_stop=}, {next_exon_start=}, {raw_curr_exon_stop=}, {raw_next_exon_start=}, {flanked_stop=}, {next_flanked_start=}')
                ## TODO: REVISE THE LOGIC; MAKE SURE THAT IT COVERS ALL THE CASES
                # upd_flanked_stop: int = min(
                #     max(
                #         min(flanked_stop, nn(next_flanked_start - 1)),
                #         raw_curr_exon_stop
                #     ),
                #     raw_next_exon_start - 1
                # )
                upd_flanked_stop: int = max(
                    min(flanked_stop, nn(next_flanked_start - 1)),
                    raw_curr_exon_stop
                )
                # upd_next_flanked_start: int = max(
                #     min(
                #         max(min(flanked_stop + 1, self.chrom_len), next_flanked_start),
                #         raw_next_exon_start
                #     ),
                #     raw_curr_exon_stop + 1
                # )
                upd_next_flanked_start: int = min(
                    max(min(flanked_stop + 1, self.chrom_len), next_flanked_start),
                    raw_next_exon_start
                )
                # print(f'{upd_flanked_stop=}, {upd_next_flanked_start=}')
                if num not in self.gap_located:
                    self.exon_search_spaces[num] = (
                        *self.exon_search_spaces[num][:3],
                        upd_flanked_stop
                    )
                if next_ not in self.gap_located:
                    self.exon_search_spaces[next_] = (
                        *self.exon_search_spaces[next_][:2],
                        upd_next_flanked_start,
                        self.exon_search_spaces[next_][3]
                    )
                    # print(f'{num=}, {self.exon_search_spaces[num]=}, {next_=}, {self.exon_search_spaces[next_]=}')
                else:
                    pass ## TODO: What to do here? We still need to catch all the marginal cases

            ## exons within each group should share the same coordinates
            # print(self.cesar_exon_grouping)
            group: List[int] = self.cesar_exon_grouping[-1]
            # if len(group) == 1:
                # print(f'{self.exon_coords[num]=}')
                # continue
            group_start: int = min(
                self.exon_coords[x][2] for x in group if x in self.exon_coords
            )
            group_stop: int = max(
                self.exon_coords[x][3] for x in group if x in self.exon_coords
            )
            ## in few cases, group borders can overlap with those of the previous
            ## group upon adding a new exon; thus, check for the potential
            ## intron deletion
            if len(self.cesar_exon_grouping) > 1:
                prev_group: List[int] = self.cesar_exon_grouping[-2]
                prev_start, prev_stop = self.exon_coords[prev_group[-1]][2:]
                prev_intron_length: int = group_start - prev_stop
                if prev_intron_length < MIN_INTRON_LENGTH:
                    click.echo(
                        'Previous intron is too short for groups comprising '
                        f'of exons {",".join(map(str, group))}'
                    )
                    if all(x in all_missing for x in prev_group) and prev_intron_length >= 0:
                        click.echo(
                            'Previous group comprises of exon with undefined loci; '
                            'adjusting the coordinates'
                        )
                        new_prev_stop: int = nn(min(group_start - 1, self.chrom_len))
                        for e in prev_group:
                            self.exon_coords[e] = (
                                self.chrom, self.strand, prev_start, new_prev_stop
                            )
                    else:
                        click.echo(
                            'Previous group comprises of exons with defined loci; '
                            'combining the two groups'
                        )
                        self.cesar_exon_grouping = self.cesar_exon_grouping[:-2]
                        prev_group.extend(group)
                        self.cesar_exon_grouping.append(prev_group)
                        group = prev_group.copy()
                        group_start = prev_start

            group: List[int] = self.cesar_exon_grouping[-1]
            group_start: int = min(
                self.exon_coords[x][2] for x in group if x in self.exon_coords
            )
            group_stop: int = max(
                self.exon_coords[x][3] for x in group if x in self.exon_coords
            )
            for e in group:
                self.exon_coords[e] = (
                    self.chrom,
                    self.strand,
                    nn(min(group_start, self.chrom_len)),
                    nn(min(group_stop, self.chrom_len))
                )

        ## final adjustments performed after all the exons are grouped
        if self.mapper:
            for i, group in enumerate(self.cesar_exon_grouping):
                if not all(x in all_missing for x in group):
                    # continue
                    ## for groups containing exons with undefined loci, check if their
                    ## search spaces could be further reduced

                    ## say, two exons with undefined loci are confined between two of those
                    ## which correspond to the aligned blocks

                    ## exons ---------exon1----------ex2--e3---------------------exon4------
                    ## chain =======[block1]====[bl2]=========[bloc3]===========[block4]====

                    ## here, exon 1 and exon 4 have defined loci, while exons 2 and 3 are
                    ## deemed missing and thus form a multi-exon CESAR group
                    ## in this case, one can suggest that exons 2 and 3 are likely confined
                    ## between blocks 2 and 3
                    click.echo(
                        f'Shrinking search space for exon group {i}'
                    )
                    # print(f'{i=}, {group=}, {self.exon_coords[group[0]]=}')
                    group_start, group_stop = self.exon_coords[group[0]][2:]
                    # print(f'{group_start=}, {group_stop=}')
                    suitable_blocks: List[Tuple[int]] = self.mapper.get_suitable_blocks(
                        group_start, group_stop
                    )
                    group_len: int = sum(
                        [self.annot.exons[x].stop - self.annot.exons[x].start for x in group]
                    ) + MIN_INTRON_LENGTH * (len(group) - 1)
                    max_inter_coords: Tuple[Union[int, None]] = define_max_space(
                        suitable_blocks, group_len, self.max_space_size
                    )
                    # print(f'{group_start=}, {group_stop=}, {max_inter_coords=}, {group_len=}, {self.max_space_size=}')
                    if max_inter_coords[0] is None:
                        continue
                    if (max_inter_coords[1] - max_inter_coords[0]) < (group_stop - group_start):
                        click.echo(
                            f'Previous coords for group {i}: '
                            f'{self.chrom}:{group_start}-{group_stop} '
                            f'(size = {group_stop-group_start}), '
                            f'shrinked coords: '
                            f'{self.chrom}:{max_inter_coords[0]}-{max_inter_coords[1]} '
                            f'(size = {max_inter_coords[1]-max_inter_coords[0]})'
                        )
                        for e in group:
                            self.exon_coords[e] = (self.chrom, self.strand, *max_inter_coords)
                else:
                    pass ## TODO: Move all the expected coordinate adjustments here
                    #     pass
                    # for j, exon in enumerate(group):
                    #     pass
                        # if exon == self.first_exon:
                        #     pass
                        # if exon == self.last_exon:
                        #     pass
                        # next_: int = exon + (1 if self.strand else - 1)
                        # if next_ in group and next not in all_missing:
                        #     pass
                        # else:
                        #     pass


        ## alternative exon locations are inserted instead of the primary
        ## location in the respective groups provided they fit by their
        ## coordinates, otherwise they are put into a separate group
        for a_num in self.alternative_loci: ## TODO: add alternative groups to cesar_exon_grouping
            groups_containing_exon: List[List[int]] = [
                x for x in self.cesar_exon_grouping if a_num in x
            ]
            for a, (alt_start, alt_stop) in enumerate(self.alternative_loci[a_num]):
                subst: Tuple[int, int] = (a_num, a)
                alt_len: int = alt_stop - alt_start
                group_found: bool = False
                for group in groups_containing_exon:
                    if len(group) == 1: ## a redundant replacement
                        continue
                    group_start: int = self.exon_coords[group[0]][0]
                    group_stop: int = self.exon_coords[group[-1]][1]
                    intersection_size: int = intersection(
                        alt_start, alt_stop, group_start, group_stop
                    )
                    ## alternative exon lies further than MIN_INTRON_LENGTH from
                    ## the rest of the group -> skipping
                    if intersection_size < -MIN_INTRON_LENGTH:
                        continue
                    ## primary exon lies lies not on the group boundary,
                    ## alternative locus is not fully enclosed -> skipping
                    is_group_edge: bool = a_num == group[0] or a_num == group[-1]
                    if not is_group_edge and alt_len > intersection_size:
                        continue
                    new_group: List[Union[int, Tuple[int]]] = group.copy()
                    i: int = new_group.index(a_num)
                    new_group[i] = subst
                    self.cesar_exon_grouping.append(new_group)
                    group_found = True
                if not group_found:
                    self.cesar_exon_grouping.append([subst])

    def add_flanks_to_exon_coords(self) -> Dict[str, Tuple[str, bool, int, int]]:
        """Adds flanks to exon coordinates"""
        all_missing: Set[int] = self.missing.union(self.out_of_chain_exons)
        output: Dict[str, Tuple[str, bool, int, int]] = {}
        for k,v in self.exon_coords:
            if k in all_missing:
                output[k] = v
                continue
            if k in self.gap_located:
                output[k] = v
                continue
            chrom, strand, start, end = v
            if self.strand:
                start -= (self.acc_flank if k != self.first_exon else 0)
                end += (self.donor_flank if k != self.last_exon else 0)
            else:
                start -= (self.donor_flank if k != self.last_exon else 0)
                end += (self.acc_flank if k != self.first_exon else 0)
            search_coords: Tuple[str, bool, int, int] = consistent_coords(
                (chrom, strand, start, end), self.chrom_len
            )
            output[k] = search_coords
        return output


def all_projection_exons(
    projection: Tuple[Union[Exon2BlockMapper, None], Union[Segment, None]]
) -> List[int]:
    """
    Returns the list of exon numbers present in the chain-segment combination
    """
    output: List[int] = []
    chain, segment = projection
    if chain is not None:
        output.extend(list(chain.e2c.keys()))
    if segment is not None:
        output.extend(list(segment.exons.keys()))
    return output

## OBSOLETE
# @dataclass
# class CesarInput:
#     __slots__ = [
#         'input_str', 'projection', 'group', 'chain', 'exons', 'memory',
#         'chrom', 'start', 'stop', 'strand', 'asmbl_gaps', 'was_not_aligned',
#         'evidence', 'contains_max'
#     ]
#
#     input_str: str
#     projection: str
#     group: int
#     exons: List[int]
#     memory: float
#     chain: str
#     chrom: str
#     start: int
#     stop: int
#     strand: bool
#     asmbl_gaps: List[Tuple[int]]
#     was_not_aligned: bool
#     evidence: Dict[int, int]
#     contains_max: bool

## OBSOLETE
# def table_to_cesar_group(table: Iterable[Iterable[str]], min_gap_size: int = 10) -> Dict[int, Dict[int, str]]:
#     """
#     Given the NumPy record of the transcript preprocessing,
#     return the restored CESAR input
#     """
#     out_dict: Dict[int, Dict[int, str]] = defaultdict(dict)
#     num_group: int = 0
#     for i in range(table.shape[0]):
#         line: List[str] = [x.decode('utf-8') for x in table[i]]
#         print(line)
#         input_str: str = ''
#         tr, chain, group, exons, mem, chrom, start, stop, strand = line[:9]
#         group: int = int(group)
#         exons: List[int] = [int(x) for x in exons.split(',')]
#         mem: float = float(mem)
#         start, stop = int(start), int(stop)
#         strand: bool = bool(int(strand))
#         exon_headers: List[str] = line[9].split('|')
#         exon_seqs: List[str] = line[10].split('|')
#         if len(exon_headers) != len(exon_seqs):
#             raise Exception(
#                 'CESAR input data is likely corrupted; unequal number of exon '
#                 f'headers and sequences for group {group} at line {i}'
#             )
#         for j in range(len(exon_headers)):
#             input_str += f'{exon_headers[j]}\n{exon_seqs[j]}\n'
#         space_header, space_seq = line[11:13]
#         input_str += f'####\n{space_header}\n{space_seq}'
#         was_not_aligned: bool = bool(int(line[13]))
#         evidence: Dict[int, int] = dict([map(int, x.split('|')) for x in line[14].split(',')])
#         contains_max: bool = bool(int(line[15]))
#         gaps: List[Tuple[int]] = find_gaps(
#             space_seq, min_gap_size, start if strand else stop, strand, chrom
#         )
#         out_dict[tr][num_group] = CesarInput(
#             input_str, tr, num_group, exons, mem, chain,
#             chrom, start, stop, strand, gaps, was_not_aligned, evidence, contains_max
#         )
#         num_group += 1
#
#     return out_dict

@dataclass
class CesarInput:
    __slots__ = [
        'segment', 'projection', 'group', 'chain', 'exons', 'memory',
        'exon_headers', 'exon_seqs', 'query_seq',
        'chrom', 'start', 'stop', 'strand', 'u12_data', 'spliceai_data',
        'expected_coords', 'search_space_coords',
        'gap_located_exons', 'out_of_chain_exons', 'intersects_asmbl_gaps',
        'evidence', 'contains_last_exon', 'acceptor_flank', 'donor_flank',
        'will_be_aligned'
    ]

    segment: str
    projection: str
    group: int
    exons: List[int]
    memory: float
    exon_headers: List[str]
    exon_seqs: List[str]
    query_seq: str
    chain: str
    chrom: str
    start: int
    stop: int
    strand: bool
    u12_data: Dict[int, Set[int]]
    spliceai_data: Dict[str, Dict[int, float]]
    expected_coords: Dict[int, Tuple[Union[int, None]]]
    search_space_coords: Dict[int, Tuple[Union[int, None]]]
    gap_located_exons: List[int]
    out_of_chain_exons: List[int]
    intersects_asmbl_gaps: bool
    contains_last_exon: bool
    acceptor_flank: bool
    donor_flank: bool
    # evidence: Dict[int, int]
    will_be_aligned: bool

def table_to_cesar_group(table: Iterable[Iterable[str]]) -> Dict[int, Dict[int, CesarInput]]:
    """
    Given the NumPy record of the transcript preprocessing,
    return the restored CESAR input
    """
    out_dict: Dict[int, Dict[int, CesarInput]] = defaultdict(dict)
    num_group: int = 0
    for i in range(table.shape[0]):
        line: List[str] = [x.decode('utf-8') for x in table[i]]
        input_str: str = ''
        # segment: int = int(line[0])
        segment: int = line[0] ## PLACEHOLDER; CORRECT FURTHER
        projection: int = int(line[1])
        exon_group: int = int(line[2])
        chain: str = line[3]
        chrom: str = line[4]
        start, stop = map(int, line[5:7])
        strand: bool = line[7] == '+'
        ex_nums: List[int] = sorted(map(int, line[8].split(',')))

        exon_names: List[str] = line[9].split('|')
        exon_seqs: List[str] = line[10].split('|')
        if len(exon_names) != len(exon_seqs):
            raise Exception(
                'CESAR input data is likely corrupted; unequal number of exon '
                f'headers and sequences for group {group} at line {i}'
            )
        for x in zip(exon_names, exon_seqs):
            input_str += '\n'.join(x) + '\n'
        query_header: str = f'{chrom}:{start}-{stop}'
        query_seq: str = line[11]
        input_str += (query_header + '\n' + query_seq) + '\n'

        u12_data: Dict[int, Set[int]] = defaultdict(set)
        acc_u12_status: str = line[12]
        for acc_stat in acc_u12_status.split(','):
            ex, stat = acc_stat.split(':')
            if stat == '1':
                ex = int(ex)
                u12_data[ex].add(0)
        donor_u12_status: str = line[13]
        for donor_stat in donor_u12_status.split(','):
            ex, stat = donor_stat.split(':')
            if stat == '1':
                ex = int(ex)
                u12_data[ex].add(1)

        ## TODO: For an upcoming U12 handling update
        # intron_data: Dict[int, Dict[str, Tuple[str]]] = {
        #     {x: {'donor': (,), 'acceptor': (,)}} for x in ex_nums
        # }
        # intron_line: str = line[???]
        # if intron_data:
        #     for intron in intron_data.split(';'):
        #         intron_num, intron_meta = intron.split(':')
        #         intron_num = int(intron_num)
        #         intron_class, donor_dinuc, acc_dinuc = intron_meta.split('#')
        #         donor_exon: int = intron_num
        #         acc_exon: int = intron_num + 1
        #         intron_data[donor_exon]['donor'] = (intron_class, donor_dinuc)
        #         intron_data[acc_exon]['acceptor'] = (intron_class, acc_dinuc)

        spliceai_data: Dict[str, Dict[int, float]] = {'acceptor': {}, 'donor': {}}
        spliceai_acc_data: str = line[14]
        if spliceai_acc_data:
            for pred in spliceai_acc_data.split(','):
                site, prob = pred.split(':')
                site = int(site)
                prob = float(prob)
                spliceai_data['acceptor'][site] = prob
        spliceai_donor_data: str = line[15]
        if spliceai_donor_data:
            for pred in spliceai_donor_data.split(','):
                site, prob = pred.split(':')
                site = int(site)
                prob = float(prob)
                spliceai_data['donor'][site] = prob

        exp_coords: Dict[int, Tuple[int]] = {}
        exp_coords_data: str = line[16]
        for locus in exp_coords_data.split(','):
            ex, coords = locus.split(':')
            ex = int(ex)
            exp_start, exp_stop = coords.split('-')
            exp_start = None if exp_start == 'None' else int(exp_start)
            exp_stop = None if exp_stop == 'None' else int(exp_stop)
            exp_coords[ex] = (exp_start, exp_stop)

        search_space_coords: Dict[int, Tuple[int]] = {}
        search_space_data: str = line[17]
        for locus in search_space_data.split(','):
            ex, coords = locus.split(':')
            ex = int(ex)
            search_start, search_stop = coords.split('-')
            search_start = None if search_start == 'None' else int(search_start)
            search_stop = None if search_stop == 'None' else int(search_stop)
            search_space_coords[ex] = (search_start, search_stop)

        gap_located_exons: List[int] = sorted(
            map(int, line[18].split(','))
        ) if line[18] else []
        out_of_chain_exons: List[int] = sorted(
            map(int, line[19].split(','))
        ) if line[19] else []
        intersects_gaps: bool = line[20] == '1'
        contains_last_exon: bool = line[21] == '1'
        acc_flank: int = int(line[22])
        donor_flank: int = int(line[23])
        will_be_aligned: bool = line[24] == '1'
        mem: float = float(line[25])

        out_dict[segment][num_group] = CesarInput(
            segment, projection, exon_group, ex_nums, mem,
            exon_names, exon_seqs, query_seq,
            chain, chrom, start, stop, strand, u12_data, spliceai_data,
            exp_coords, search_space_coords, gap_located_exons, out_of_chain_exons,
            intersects_gaps, contains_last_exon, acc_flank, donor_flank,
            will_be_aligned
        )
        num_group += 1

    return out_dict

def a_table_to_cesar_group(table: Iterable[Iterable[str]]) -> Dict[int, Dict[int, CesarInput]]:
    """
    Given the NumPy record of the transcript preprocessing,
    return the restored CESAR input
    """
    out_dict: Dict[int, Dict[int, CesarInput]] = defaultdict(dict)
    num_group: int = 0
    for i in range(table.shape[0]):
        line: List[str] = [x.decode('utf-8') for x in table[i]]
        input_str: str = ''
        # segment: int = int(line[0])
        segment: int = line[0] ## PLACEHOLDER; CORRECT FURTHER
        projection: int = int(line[1])
        exon_group: int = int(line[2])
        chain: str = line[3]
        chrom: str = line[4]
        start, stop = map(int, line[5:7])
        strand: bool = line[7] == '+'
        ex_nums: List[int] = sorted(map(int, line[8].split(',')))

        exon_names: List[str] = line[9].split('|')
        exon_seqs: List[str] = line[10].split('|')
        if len(exon_names) != len(exon_seqs):
            raise Exception(
                'CESAR input data is likely corrupted; unequal number of exon '
                f'headers and sequences for group {group} at line {i}'
            )
        for x in zip(exon_names, exon_seqs):
            input_str += '\n'.join(x) + '\n'
        query_header: str = f'{chrom}:{start}-{stop}'
        query_seq: str = line[11]
        input_str += (query_header + '\n' + query_seq) + '\n'

        # u12_data: Dict[int, Set[int]] = defaultdict(set)
        # acc_u12_status: str = line[12]
        # for acc_stat in acc_u12_status.split(','):
        #     ex, stat = acc_stat.split(':')
        #     if stat == '1':
        #         ex = int(ex)
        #         u12_data[ex].add(0)
        # donor_u12_status: str = line[13]
        # for donor_stat in donor_u12_status.split(','):
        #     ex, stat = donor_stat.split(':')
        #     if stat == '1':
        #         ex = int(ex)
        #         u12_data[ex].add(1)

        ## TODO: For an upcoming U12 handling update
        intron_data: Dict[int, Dict[str, Tuple[str]]] = {
            x: {'donor': tuple(), 'acceptor': tuple()} for x in ex_nums
        }
        splice_site_line: str = line[12]
        if splice_site_line:
            for ex in splice_site_line.split(';'):
                ex_num, ex_meta = ex.split(':')
                ex_num = int(ex_num)
                ex_meta = (None if x == 'None' else x for x in ex_meta.split(','))
                acc_class, acc_site, donor_class, donor_site = ex_meta
                # donor_exon: int = intron_num
                # acc_exon: int = intron_num + 1
                intron_data[ex_num]['donor'] = (donor_class, donor_site)
                intron_data[ex_num]['acceptor'] = (acc_class, acc_site)

        spliceai_data: Dict[str, Dict[int, float]] = {'acceptor': {}, 'donor': {}}
        spliceai_acc_data: str = line[13]
        if spliceai_acc_data:
            for pred in spliceai_acc_data.split(','):
                site, prob = pred.split(':')
                site = int(site)
                prob = float(prob)
                spliceai_data['acceptor'][site] = prob
        spliceai_donor_data: str = line[14]
        if spliceai_donor_data:
            for pred in spliceai_donor_data.split(','):
                site, prob = pred.split(':')
                site = int(site)
                prob = float(prob)
                spliceai_data['donor'][site] = prob

        exp_coords: Dict[int, Tuple[int]] = {}
        exp_coords_data: str = line[15]
        for locus in exp_coords_data.split(','):
            ex, coords = locus.split(':')
            ex = int(ex)
            exp_start, exp_stop = coords.split('-')
            exp_start = None if exp_start == 'None' else int(exp_start)
            exp_stop = None if exp_stop == 'None' else int(exp_stop)
            exp_coords[ex] = (exp_start, exp_stop)

        search_space_coords: Dict[int, Tuple[int]] = {}
        search_space_data: str = line[16]
        for locus in search_space_data.split(','):
            ex, coords = locus.split(':')
            ex = int(ex)
            search_start, search_stop = coords.split('-')
            search_start = None if search_start == 'None' else int(search_start)
            search_stop = None if search_stop == 'None' else int(search_stop)
            search_space_coords[ex] = (search_start, search_stop)

        gap_located_exons: List[int] = sorted(
            map(int, line[17].split(','))
        ) if line[17] else []
        out_of_chain_exons: List[int] = sorted(
            map(int, line[18].split(','))
        ) if line[18] else []
        intersects_gaps: bool = line[19] == '1'
        contains_last_exon: bool = line[20] == '1'
        acc_flank: int = int(line[21])
        donor_flank: int = int(line[22])
        will_be_aligned: bool = line[23] == '1'
        mem: float = float(line[24])

        out_dict[segment][num_group] = CesarInput(
            segment, projection, exon_group, ex_nums, mem,
            exon_names, exon_seqs, query_seq,
            chain, chrom, start, stop, strand, intron_data, spliceai_data,
            exp_coords, search_space_coords, gap_located_exons, out_of_chain_exons,
            intersects_gaps, contains_last_exon, acc_flank, donor_flank,
            will_be_aligned
        )
        num_group += 1

    return out_dict

def table_to_intron_data(table: Iterable[Iterable[str]]) -> Dict[int, Set[int]]:
    """Infers U12 exons from the HDF5 table record. Likely a provisional function"""
    out_dict: Dict[int, Set[int]] = defaultdict(set)
    for i in range(table.shape[0]):
        line: Iterable[str] = table[i][:]
        exon_num: int = line[1]
        out_dict[exon_num].add(0) if line[2] == 1 else None
        out_dict[exon_num].add(1) if line[3] == 1 else None
    return out_dict


@dataclass
class CesarExonEntry:
    """
    Stores minimally processed CESAR2 alignment results; the alignment is bound
    to specific location in query, both reference and query sequences are
    stripped of non-coding symbols but not split into codons
    """
    __slots__ = [
        'chain', 'chrom', 'start', 'stop', 'name', 'num', 'strand',
        'seq', 'hit_seq', 'splice_sites', 'acc_site', 'donor_site', 'trailing_intron',
        'evidence', 'intersects_gap', 'was_not_aligned'
    ]
    chain: str
    chrom: str
    start: int
    stop: int
    name: str
    num: int
    strand: str
    seq: str
    hit_seq: str
    # splice_sites: Dict[int, Dict[str, Tuple[str;]]]
    acc_site: str
    donor_site: str
    trailing_intron: str
    evidence: int
    intersects_gap: bool
    was_not_aligned: bool

    def __str__(self):
        return '\t'.join(map(str, (self.__getattribute__(x) for x in self.__slots__[1:6])))


@dataclass
class RawCesarOutput:
    """
    A data class to store the raw CESAR 2.0 output. Contains minimally processed
    results as well as auxiliary data on alignment location in the query
    """
    __slots__ = [
        'chain', 'chrom', 'start', 'stop', 'strand',
        'reference', 'query', 'exons', 'exon_expected_loci', 'exon_search_spaces',
        'spliceai_sites', 'gap_located_exons', 'out_of_chain_exons',
        'was_not_aligned', 'assembly_gap', 'subexon_coordinates'
    ]
    chain: str
    chrom: str
    start: int
    stop: int
    strand: bool
    reference: str
    query: str
    exons: List[int]
    exon_expected_loci: Dict[int, Tuple[int]]
    exon_search_spaces: Dict[int, Tuple[int]]
    spliceai_sites: Dict[str, List[int]]
    gap_located_exons: Set[int]
    out_of_chain_exons: Set[int]
    was_not_aligned: bool
    assembly_gap: bool
    subexon_coordinates: Dict[int, List[Tuple[int, int]]]


@dataclass
class Mutation:
    """
    A slightly extended version of named tuples from older inact_mut_check.py .
    Stores basic mutation properties.
    """
    __slots__ = [
        'transcript', 'chain', 'exon', 'codon', 'ref_codon',
        'chrom', 'start', 'stop',
        'mutation_class', 'description', 'is_masked', 'masking_reason',
        'mutation_id'
    ]
    transcript: str
    chain: Union[int, str]
    exon: int
    codon: int
    ref_codon: int
    chrom: str
    start: int
    stop: int
    mutation_class: str
    description: str
    is_masked: bool
    masking_reason: str
    mutation_id: str

    def __str__(self) -> str:
        output: List[str] = [
            f'{self.transcript}#{self.chain}'
        ] + [
            str(self.__getattribute__(x)) for x in self.__slots__[2:]
        ]
        output[9] = 'MASKED' if self.is_masked else 'NOT_MASKED'
        return '\t'.join(output)


def fa_name(name: str) -> str:
    """
    Checks if FASTA format name contains a '>' symbol at the beginning, adds it otherwise
    """
    name = name.lstrip().lstrip('\t')
    if not name:
        raise ValueError('Empty sequence name encountered')
    if name[0] != '>':
        name = '>' + name
    return name


def find_gaps(seq: str, min_length: int) -> bool:#List[Tuple[int, int]]:
    """
    Reports if assembly gap of specified length have been recorded
    in the provided sequence
    """
    return 'N' * min_length in seq.upper()


def parse_extr_exon_fasta(raw_fasta: str) -> Dict[str, Dict[int, str]]:
    """
    An ad hoc FASTA parser for CESAR preprocessing step. Given a raw FASTA string
    containing exon sequences for one or more transcripts, returns a
    {transcript: {exon_number: exon_sequence}} nested dictionary. Exon naming
    convention is {transcript_name}|{exon_number}.
    """
    output_dict: Dict[str, Dict[int, str]] = defaultdict(dict)
    has_entry: bool = False
    curr_seq: str = ''
    for line in raw_fasta.split('\n'):
        if not line:
            continue
        if line[0] == '>':
            if has_entry:
                output_dict[transcript][exon] = curr_seq
                curr_seq = ''
            has_entry = True
            transcript, exon = line.split('|')
            transcript = transcript.lstrip('>')
            exon = int(exon)
        else:
            curr_seq += line
    if has_entry:
        output_dict[transcript][exon] = curr_seq
    return output_dict


def prepare_exons(
    raw_exons: Dict[int, str],
    exon_flanks: Dict[int, Tuple[int, int]],
    mask_stops: bool = False
) -> Dict[int, str]:
    """
    """
    s_sites: List[Any] = []
    exon_seqs: Dict[int, str] = {}
    max_exon: int = max(raw_exons.keys())
    phase: int = 0
    for num, raw_seq in sorted(raw_exons.items(), key=lambda x: x[0]):
        raw_seq = raw_seq.upper()
        is_first: bool = num == 1
        is_last: bool = num == max_exon
        # stderr.write(f'>{num}\n{raw_seq}\n')
        ## process splice sites if any were captured during sequence extraction
        acc_flank, donor_flank = exon_flanks[num]
        acc_: str = raw_seq[:acc_flank]
        donor_: str = raw_seq[-donor_flank:]
        if is_first:
            acc_ = 'NN'
            # acc_flank = 0
        elif len(acc_) < 2:
            acc_ = 'N' * (2 - len(acc_)) + acc_
        if is_last:
            donor_ = 'NN'
            # donor_flank = 0
        elif len(donor_) < 2:
            donor_ = donor_ + 'N' * (2 - len(donor_))
        s_sites.append(acc_.lower())
        s_sites.append(donor_.lower())

        exon_seq: str = raw_seq[acc_flank:-donor_flank]
        ## assess exon's phase
        prev_phase: int = (3 - phase) if phase > 0 else 0
        next_phase: int = (len(exon_seq) - prev_phase) % 3
        inframe_end: int = len(exon_seq) - next_phase
        ## check sequence for inframe stop codons
        seq_to_check: str = exon_seq[:inframe_end]
        if num > 1 and phase:
            seq_to_check = exon_seqs[num-1][-phase:] + seq_to_check
        exon_len: int = len(exon_seq) + phase
        for i in range(0, len(seq_to_check), 3):
            triplet: str = seq_to_check[i:i+3]
            if triplet.upper() in STOPS and not (i + 3 >= exon_len and is_last):
                if not mask_stops:
                    raise Exception(
                        'Inframe stop codon found in the reference sequence. '
                        'If this is a known readthrough/selenocysteine/repurposed '
                        'codon, consider setting --mask_inframe_stops flag; '
                        'otherwise, remove transcript {tr} from TOGA input.'
                    ) ## TODO: Add transcript name slot to AnnotationEntry
                if phase and i < phase:
                    exon_seqs[num-1] = (
                        exon_seqs[num-1][:-phase] + 'n' * phase
                    )
                if phase and i < phase and i + 3 > phase:
                    exon_seq = 'n' * prev_phase + exon_seq[prev_phase:]
                else:
                    exon_seq = exon_seq[:i-phase] + 'NNN' + exon_seq[i-phase+3:]
        ## phase exon sequence for CESAR input
        exon_seq = (
            exon_seq[:prev_phase].lower() + exon_seq[prev_phase:inframe_end] +
            exon_seq[inframe_end:].lower()
        )
        phase = next_phase
        exon_seqs[num] = exon_seq
    return exon_seqs, s_sites


def get_exons(
    annot: AnnotationEntry,
    chrom_seq: str,
    ss_size: int = SS_SIZE,
    exon_flank: int = FLANK_SPACE,
    mask_stops: bool = False
) -> Tuple[Any, Any]:
    """
    Extract exons sequences for reference.
    Taken from original CESAR_wrapper
    """
    # exons_raw: List[Tuple[int, int]] = [(x.start, x.stop) for x in annot.exons.values()]
    # if not annot.strand:
    #     exons_raw = exons_raw[::-1]
    # exons_raw.reverse() if annot.strand == '-' else None
    # exons_pos: Dict[int, Tuple[int, int]] = {}  # contain exon_num : positions
    s_sites: List[Any] = []
    exon_flanks: Dict[int, Dict[str, int]] = {}
    # for num, exon in enumerate(exons_raw, 1):
    #     # start, end if strand == + end, start otherwise
    #     # need to know start and end to extract from 2bit file
    #     exons_pos[num] = (
    #         (int(exon[0]), int(exon[1]))
    #         if annot.strand
    #         else (int(exon[1]), int(exon[0]))
    #     )
    exons_pos: Dict[int, Tuple[int, int]] = {
        k:(v.start, v.stop) for k,v in sorted(
            annot.exons.items(), key=lambda x: x[0]
        )
    }

    max_exon_num: int = max(exons_pos.keys())
    _all_positions: List[int] = sorted(flatten(exons_pos.values()))
    gene_borders = {_all_positions[0], _all_positions[-1]}
    # extract sequences
    exon_seqs: Dict[int, str] = {}  # exon number: sequence dict
    chrom: str = annot.chrom
    # if chrom not in chrom_seq:
    #     raise KeyError(f"Error! Cannot find chrom {chrom} in 2bit file {twobit}")
    # click.echo("\nExons sequences ####\n")
    phase: int = 0
    for num, pos in exons_pos.items():
        is_first_exon = num == 1
        is_last_exon = num == max_exon_num

        # for twoBitToFa start must be < end
        # determine search start and end
        # do not subtract/add SS_SIZE if gene border: no splice sites then
        min_pos, max_pos = sorted(pos)
        start: int = min_pos - ss_size if min_pos not in gene_borders else min_pos
        end: int = max_pos + ss_size if max_pos not in gene_borders else max_pos
        # print(num, start, end)
        # get exon 10 bp flanks:
        left_brd_: int = max(min_pos - exon_flank, 0)
        right_brd_: int = min(max_pos + exon_flank, len(chrom_seq))
        left_flank_coord: Tuple[int, int] = (left_brd_, min_pos)
        right_flank_coord: Tuple[int, int] = (max_pos, right_brd_)
        left_flank: str = chrom_seq[left_brd_ : min_pos].upper()
        right_flank: str = chrom_seq[max_pos : right_brd_].upper()
        # correct for strand:
        left_flank = left_flank if annot.strand else reverse_complement(right_flank)
        right_flank = right_flank if annot.strand else reverse_complement(left_flank)
        # placeholder in case we could not extract flanks
        left_flank = left_flank if len(left_flank) > 0 else 'X'
        right_flank = right_flank if len(right_flank) > 0 else 'X'

        exon_seq_w_ss: str = chrom_seq[start:end].upper()
        if not annot.strand: # revert if negative strand
            exon_seq_w_ss = reverse_complement(exon_seq_w_ss)
        # trim splice sites
        ## TODO: Surely can be simplified
        if not is_first_exon and not is_last_exon: # both splice sites are present
            exon_seq = exon_seq_w_ss[ss_size:-ss_size]
            acc_: str = exon_seq_w_ss[:ss_size]
            don_: str = exon_seq_w_ss[-ss_size:]
        elif is_first_exon and is_last_exon: # no splice sites
            exon_seq = exon_seq_w_ss
            acc_: str = 'NN'
            don_: str = 'NN'
        elif is_first_exon: ## only donor site available
            exon_seq = exon_seq_w_ss[:-ss_size]
            acc_: str = 'NN'
            don_: str = exon_seq_w_ss[-ss_size:]
        elif is_last_exon: ## onlu acceptor site available
            exon_seq = exon_seq_w_ss[ss_size:]
            acc_: str = exon_seq_w_ss[:ss_size]
            don_: str = 'NN'
        else: ## legacy branch, kept for the consistency reasons
            raise RuntimeError('Unreachable branch reached')
        s_sites.append(acc_)
        s_sites.append(don_)
        ## define unsplit inframe portion of the exon sequence
        prev_phase: int = (3 - phase) if phase > 0 else 0
        next_phase: int = (len(exon_seq) - prev_phase) % 3
        inframe_end: int = len(exon_seq) - next_phase
        ## check sequence for inframe stop codons
        seq_to_check: str = (
            (exon_seqs[num-1][-phase:] if num > 1 else '') + exon_seq[:inframe_end]
        )
        exon_len: int = len(exon_seq) + phase
        for i in range(0, len(seq_to_check), 3):
            triplet: str = seq_to_check[i:i+3]
            if triplet.upper() in STOPS and not (i + 3 >= exon_len and is_last_exon):
                if not mask_stops:
                    raise Exception(
                        'Inframe stop codon found in the reference sequence. '
                        'If this is a known readthrough/selenocysteine/repurposed '
                        'codon, consider setting --mask_inframe_stops flag; '
                        'otherwise, remove transcript {tr} from TOGA input.'
                    ) ## TODO: Add transcript name slot to AnnotationEntry
                if phase and i < phase:
                    exon_seqs[num-1] = (
                        exon_seqs[num-1][:-phase] + 'n' * phase
                    )
                if phase and i < phase and i + 3 > phase:
                    exon_seq = 'n' * prev_phase + exon_seq[prev_phase:]
                else:
                    exon_seq = exon_seq[:i-phase] + 'NNN' + exon_seq[i-phase+3:]
        ## phase exon sequence for CESAR input
        exon_seq = (
            exon_seq[:prev_phase].lower() + exon_seq[prev_phase:inframe_end] +
            exon_seq[inframe_end:].lower()
        )
        phase = next_phase
        ## check for start and stop codon presence and inframe stop codons
        ## TODO: Implement
        exon_seqs[num] = exon_seq
        exon_flanks[num] = {'L': left_flank, 'R': right_flank}
        # click.echo(f'Exon {num} in the range {pos}; sequence:\n{exon_seq}' + '\n')
    return exons_pos, exon_seqs, s_sites, exon_flanks


def make_matrix(matrix: str = Union[str, TextIO]) -> Dict[str, Dict[str, int]]:
    """
    Given an alignment matrix file, parses it as
    """
    matrix_dict: Dict[str, Dict[str, int]] = defaultdict(dict)
    if isinstance(matrix, TextIOWrapper):
        lines: List[str] = matrix.readlines()
    else:
        with open(matrix, 'r') as h:
            lines: List[str] = h.readlines()
    header: List[str] = lines[0].split()
    for line in lines[1:]:
        data: List[str] = line.split()
        if len(data) - 1 != len(header):
            raise Exception(
                'Score matrix is incomplete. '
                'Check that the number of rows corresponds to the number of columns.'
            )
        key: str = data[0]
        for i in range(1, len(data)):
            subkey: str = header[i-1]
            value: int = int(data[i])
            matrix_dict[key][subkey] = value

    return matrix_dict


def parseU12(file_content: TextIO, progenitor: str = None) -> Dict[int, Set[int]]:
    """
    Parses the three-column file containing external U12 exon evidence
    """
    if not progenitor:
        print(
            'WARNING: Source transcript was not set for the U12 parser; '
            'returning empty object.'
        )
    progenitor = progenitor.split('#')[0]
    out_dict: Dict[int, Set[int]] = defaultdict(set)
    for line in file_content.readlines():
        data: List[str] = line.strip().split('\t')
        transcript: str = data[0].split('#')[0]
        if transcript == progenitor:
            exon_num: int = int(data[1])
            flank_num: int = int(data[2] == 'D')
            out_dict[exon_num].add(flank_num)

    return out_dict


def infer_u12_sites(
    sites: List[str],
    u12: Dict[int, Set[int]]
) -> Dict[int, Set[int]]:
    """
    Updates the U12 splice site dict by checking the terminal nucleotides exonwise.
    TODO: The current exon processing procedure is rather cumbersome. Consider
    adding the functionality encapsulated in this function into get_exons or a similar
    function
    """
    sites_split: List = parts(sites[1:-1], 2)
    for exon, intron in enumerate(sites_split):
        don_, acc_ = map(lambda x: x.lower(), intron)
        if (don_ not in DONOR_SITE and don_ != 'nn') or (acc_ not in ACCEPTOR_SITE and acc_ != 'nn'):
            u12[exon+1].add(1)
            u12[exon+2].add(0)
        ## TODO: I am not sure if the approach above works; sites should be treated separately
        ## better ask Bogdan
        # u12[exon+1].add(1) if don_ not in DONOR_SITE else None
        # u12[exon+2].add(0) if acc_ not in ACCEPTOR_SITE else None

    return u12


def check_ref_exons(
    exon_seqs: Dict[int, str],
    mask_stops: bool = True,
    stops: Set[str] = STOPS
) -> Tuple[List[str], Set[int]]:
    """Check if the reference sequence is correct.

    Should start with ATG and end with a stop.
    Mask_stops controls handling of inframe stops.
    """
    sec_codons: Set[int] = set()  # start coordinates of inframe TGA codons
    gene_seq: str = "".join([exon_seqs[i] for i in exon_seqs.keys()])
    codons: List[str] = parts(gene_seq, n=3)  # split a seq of letters in chunks of len == 3
    if codons[0] != "ATG":
        sys.stderr.write("Input is corrupted! Reference sequence should start with ATG!\n")
    elif codons[-1] not in STOPS:
        sys.stderr.write("Input is corrupted! Reference sequence should end with a stop codon!\n")
    stop_codons = [(n, c) for n, c in enumerate(codons[:-1]) if c in STOPS]
    if len(stop_codons) == 0:  # no stop codons -> nothing else to do
        return exon_seqs, set()
    # there are stop codons in reference sequence:
    sys.stderr.write("Warning! There are inframe stop codons!\n")
    for stop in stop_codons:
        click.echo(f"Codon num {stop[0] + 1} - {stop[1]}\n")
        codons[stop[0]] = NNN_CODON if mask_stops else codons[stop[0]]
        if stop[1] == "TGA":
            # maybe a sec codon
            sec_codons.add(stop[0])

    ## if inframe stop codons were found but masking was set to False,
    ## show error message and die
    if not mask_stops:
        sys.stderr.write(">>>STOP_CODON>>>\n")
        die("Inframe stop codons were found but masking was disabled; exiting", 0)

    safe_seq: str = "".join(codons)
    stop_masked: Dict[int, str] = {}
    prev_index: int = 0
    ## TODO: terminal nucleotide formatting can be performed right at this point
    for num, exon_seq in exon_seqs.items():
        exon_len: int = len(exon_seq)
        stop_masked[num] = safe_seq[prev_index : prev_index + exon_len]
        prev_index += exon_len
    return stop_masked, sec_codons


def prepare_exons_for_cesar(exon_seqs: Dict[int, str]) -> Dict[int, str]:
    """
    Formats exon sequences for compatibility with CESAR 2
    Off-frame terminal nucleotides are converted to lowercase:
    ATG|TTT|A CT|GTA|AAG|TGC|C TT|AGT|TGA
    ATG|TTT|a ct|GTA|AAG|TGC|c tt|AGT|TGA
    """
    left_pointer: int = 0  # init value
    cesar_input_exons: Dict[int, str] = {}  ## stores results

    for k, exon_seq in exon_seqs.items():
        # apply left pointer
        if left_pointer != 3:
            exon_seq = (
                exon_seq[:left_pointer].lower() + exon_seq[left_pointer:].upper()
            )
        # define number of letters to lowercase at the right side
        right_pointer: int = (len(exon_seq) - left_pointer) % 3
        # re-define left pointer
        left_pointer = 3 - right_pointer
        # apply right-side pointer
        if right_pointer != 0:
            exon_seq = (
                exon_seq[:-right_pointer] + exon_seq[-right_pointer:].lower()
            )
        # save prepared exon into special dict
        cesar_input_exons[k] = exon_seq
    return cesar_input_exons

def add_exon_headers(
    exon_numbers: List[int],
    source: str,
    u12_sites: Dict[int, Set[int]],
    common_acceptor: str,
    common_donor: str,
    first_acceptor: str,
    last_donor: str,
    u12_acceptor: str,
    u12_donor: str
) -> Dict[int, str]:
    """
    Add the exon headers based on exon numbers and splice site configurations
    """
    out_name: str = f'>{source}_exon'
    if len(exon_numbers) == 1:
        return {1: f'{out_name}1'}
    names: Dict[int, str] = {}
    for exon in exon_numbers:
        name: str = f'{out_name}{exon}'
        if exon == 1:
            name += '\t' + first_acceptor
        else:
            name += (
                '\t' + u12_acceptor if \
                exon in u12_sites and 0 in u12_sites[exon] else \
                '\t' + common_acceptor
            )
        if exon == len(exon_numbers):
            name += '\t' + last_donor
        else:
            name += (
                '\t' + u12_donor if \
                exon in u12_sites and 1 in u12_sites[exon] else \
                '\t' + common_donor
            )
        names[exon] = name
    return names


def transcriptwise_subchains(
    chain: str, tr2coords: Dict[str, Tuple[int, int]]
) -> List[Tuple[int]]:
    """
    Converts a string-like alignment representation into a list
    of coordinate tuples in the (ref_start, ref_end, query_start, query_end).
    All blocks are sorted and numbered according to the referenve coordinates.
    Chain gaps are introduced as similar tuples named after flanking blocks
    """
    ####### TODO: NEEDS TESTING FOR FUNCTIONALITY, COMPATIBILITY, AND PERFORMANCE SPEED
    sorted_blocks: List[Tuple[int]] = []
    tr2blocks: Dict[str, Dict[str, Tuplep[int]]] = defaultdict(dict)
    ## split chain string into separate lines
    # split_chain: List[str] = chain.split('\n')
    # init1 = datetime.now()
    # split_chain = chain.split('\n')
    # stderr.write(f'Pythonic chain splitting time: {str(datetime.now() - init1)}\n')
    init1 = datetime.now()
    split_chain: Iterable[str] = split_at_(chain, '\n')
    stderr.write(f'Cythonic chain splitting time: {str(datetime.now() - init1)}\n')
    ## extract all the necessary data from the chain header
    # chain_header: List[str] = split_chain[0].split()
    chain_header: List[str] = next(split_chain).split()
    ref_chr: str = chain_header[2]
    ref_chr_size: int = int(chain_header[3])
    ref_strand: bool = chain_header[4] == '+'
    ref_start: int = int(chain_header[5])
    ref_end: int = int(chain_header[6])
    if not ref_strand:
        t_chain_start: int = ref_chr_size - ref_end
        t_chain_stop: int = ref_chr_size - ref_start
        ref_start = ref_chr_size - ref_start
        ref_end = ref_chr_size - ref_end
    else:
        t_chain_start, t_chain_stop = ref_start, ref_end
    query_chr: str = chain_header[7]
    query_chr_size: int = int(chain_header[8])
    query_strand: bool = chain_header[9] == '+'
    query_start: int = int(chain_header[10])
    query_end: int = int(chain_header[11])
    if not query_strand:
        query_start = query_chr_size - query_start
        query_end = query_chr_size - query_end
    chain_meta: Tuple[int] = (
        ref_chr, query_chr, t_chain_start, t_chain_stop, ref_strand, query_strand
    )
    ## sort the transcripts
    sorted_transcripts: Tuple[str] = tuple(
        sorted(
            (k for k in tr2coords), key=lambda x: tr2coords[x][:2]
        ) ## sort by coordinates
    )
    abs_start: int = sorted(tr2coords[sorted_transcripts[0]])[0]
    abs_stop: int = sorted(tr2coords[sorted_transcripts[-1]])[1]
    curr_tr: int = 0
    # b_num: int = 1 if ref_strand else len(split_chain) - 1
    b_num: int = 1
    gap_num: str = ''
    gap_coord_tuple: Tuple[int] = (0,0,0,0)
    r_start: int = ref_start
    q_start: int = query_start
    last_tr: str = sorted_transcripts[-1]

    for i, block in enumerate(split_chain):
        # print(f'{block=}')
        if not block:
            break
        if '\t' in block:
            size, dt, dq = map(int, block.split('\t'))
        else:
            size, dt, dq = int(block), 0, 0
        r_end: int = r_start + size * (1 if ref_strand else -1)
        q_end: int = q_start + size * (1 if query_strand else -1)
        coord_tuple: Tuple[int] = (
            *sorted((r_start, r_end)), *sorted((q_start, q_end))
        )
        ## reset the current start values
        r_start = r_end
        q_start = q_end
        ## first, check if the block intersects the transcript array
        up_block: bool = coord_tuple[1] < abs_start
        down_block: bool = coord_tuple[0] > abs_stop

        ## if iterator has not reached the first transcript's search space,
        ## there is no need to check for further intersections;
        ## however, in order to ensure that every transcript has a space to align,
        ## create a provisional gap object
        if (ref_strand and up_block) or (not ref_strand and down_block):
            r_end = r_start + dt * (1 if ref_strand else -1)
            q_end = q_start + dq * (1 if query_strand else -1)
            # if not gap_coord_tuple:
            gap_num = '0_0'
            gap_coord_tuple = (r_start, r_end, q_start, q_end)
            r_start = r_end
            q_start = q_end
            continue
        ## if iterator goes past the last transcript's locus, the process halts;
        ## however, for trans
        if ((ref_strand and down_block) or (not ref_strand and up_block)):
            print('The loop has passed the last transcript; skipping')
            discard_next_block: bool = True
            for tr in sorted_transcripts[curr_tr:]:
                if gap_num in tr2blocks[tr] and b_num not in tr2blocks[tr]:
                    discard_next_block = False
                    break
            if discard_next_block:
                for tr in sorted_transcripts[curr_tr:]:
                    if not tr2blocks[tr]:
                        tr2blocks[tr][gap_num] = gap_coord_tuple
                break
        ## check if the block intersects any of the transcripts
        for tr_num, tr in enumerate(sorted_transcripts[curr_tr:], start=curr_tr):
            tr_start, tr_end = sorted(tr2coords[tr])
            # print(f'{tr=}, {tr_start=}, {tr_end=}, {tr_end-tr_start=}')
            if ref_strand:
                if tr_end < coord_tuple[0]:
                    ## the transcript lies upstream to the aligned block;
                    ## update the transcript pointer and proceed further
                    curr_tr += 1
                    # print(f'{gap_num=}, {tr2blocks[tr]=}, {gap_num in tr2coords[tr]=}')
                    if gap_num and gap_num in tr2blocks[tr]:
                        tr2blocks[tr][str(b_num)] = coord_tuple
                    if not tr2blocks[tr]:
                        tr2blocks[tr][gap_num] = gap_coord_tuple
                    # print(f'Incrementing current transcript number after {tr}')
                    continue
                if tr_start >= coord_tuple[1]:
                    ## somehow the transcript iterator moved downstream
                    ## to the current block; as long as transcripts are
                    ## sorted properly, there is no need to check the remaining
                    ## transcripts
                    # print(f'Breaking at block {b_num} for transcript {tr}, positive strand')
                    break
            else:
                ## logic is essentially the same except for now the iterator
                ## goes upstream along the chromosome
                if tr_start >= coord_tuple[1]:
                    curr_tr += 1
                    continue
                if tr_end < coord_tuple[0]:
                    # print(f'Breaking at block {b_num}, negative strand')
                    break
            tr2blocks[tr][str(b_num)] = coord_tuple
        if curr_tr >= len(sorted_transcripts):
            break
        ## if this is the last block, there is no chain gap to follow
        # if i == len(split_chain) - 1:
        #     break
        ## otherwise, create a gap block
        gap_num: str = f'{b_num}_{b_num+1}' if ref_strand else f'{b_num-1}_{b_num}'
        r_end: int = r_start + dt * (1 if ref_strand else -1)
        q_end: int = q_start + dq * (1 if query_strand else -1)
        gap_coord_tuple: Tuple[int] = (r_start, r_end, q_start, q_end)
        ## and do the same trick
        for tr_num, tr in enumerate(sorted_transcripts[curr_tr:], start=curr_tr):
            tr_start, tr_end = tr2coords[tr]
            if ref_strand:
                if tr_end < gap_coord_tuple[0]:
                    curr_tr += 1
                    continue
                # if tr_start >= gap_coord_tuple[1]:
                #     print(f'Breaking at gap {gap_num}, positive strand')
                #     break
            else:
                if tr_start >= coord_tuple[1]:
                    curr_tr += 1
                    continue
                # if tr_end < coord_tuple[0]:
                #     print(f'Breaking at gap {gap_num}, negative strand')
                #     break
            tr2blocks[tr][gap_num] = gap_coord_tuple
            if b_num not in tr2blocks[tr]:
                tr2blocks[tr][str(b_num)] = coord_tuple
        r_start = r_end
        q_start = q_end
        if curr_tr >= len(sorted_transcripts):
            # print(f'Breaking the outer loop at block {b_num}')
            if gap_num in tr2blocks[tr]:
                tr2blocks[tr][b_num] = coord_tuple
            break
        b_num += 1# if ref_strand else -1
    return tr2blocks, chain_meta

def intersect_exons_to_blocks( ## TODO: clearly must be moved to Exon2BlockMapper methods
    exons: Dict[int, Exon],
    blocks: Dict[str, List[int]],
    chainid: str,
    # chain_prob: float,
    t_chain_start: int,
    t_chain_stop: int,
    t_chrom: str,
    q_chrom: str,
    q_chrom_size: int,
    t_strand: str,
    q_strand: bool,
    codirected: bool,
    acc_flank: int,
    donor_flank: int,
    extra_flank_modifier: float,
    verbose: bool
) -> Exon2BlockMapper:
    """
    An adapted version of the original intersect_exons_blocks_gaps function
    from CESAR_wrapper.py

    Given a list of the Exon entities and a (sub)chain structure,
    returns exon-to-block/gap mapping
    """
    # allt_starts, allt_stops, allq_starts, allq_stops = [], [], [], []
    # print(f'{blocks=}')
    ## TODO:
    ## 1) Reduce the number of for-loops; sort the blocks and add the break statement (smart counter);
    ## 1a) Benchmark against the last wrapper run; make sure no results have changed;
    ## 2) Consider increasing the number of flanking blocks up to 2 (MT3 case);
    ## 2a) Benchmark against the last wrapper run; check the number of missing/deleted exons, also benchmark vs the mm10 annot
    # for block in blocks.values():
    #     allt_starts.append(block[0])
    #     allt_stops.append(block[1])
    #     allq_starts.append(block[2])
    #     allq_stops.append(block[3])
    # global_block: Tuple[int] = (
    #     min(allt_starts),
    #     max(allt_stops),
    #     min(allq_starts),
    #     max(allq_stops)
    # )

    ## for exons which are not entirely covered by coding blocks, store which
    ## flank is uncovered; notation is (upstream_uncovered, downstream_uncovered)
    dangling_ends: Dict[int, Tuple[bool]] = {e: (False, False) for e in exons}
    # missing: Set[int] = {x for x in exons if not ex_num2blocks.get(x, [])}
    missing: Set[int] = set()
    out_of_chain: Set[int] = set()
    # out_of_chain: Set[int] = set()
    gap_located: Set[int] = set()
    exon2raw_loci: Dict[int, Tuple[int]] = {} ## for raw, untrimmed projection coordinates
    exons2coordinates: Dict[int, Tuple[int]] = {} ## for clipped/extrapolated projection coordinates
    ## first, estimate if all the exons are covered entirely by the blocks
    # for e, exon in exons.items():
    #     exon_start, exon_stop = min(exon.start, exon.stop), max(exon.start, exon.stop)
    #     exon_size: int = exon_start - exon_stop
    #     global_block_intersection: int = intersection(
    #         exon_start,
    #         exon_stop,
    #         global_block[0],
    #         global_block[1]
    #     )
    #     if global_block_intersection < exon_size:
    #         up_dangling: bool = exon_start < global_block[0]
    #         down_dangling: bool = exon_stop > global_block[0]
    #         dangling_ends[e] = (up_dangling, down_dangling)

    ## then, find the exon-to-block intersections blockwise
    ex_num2blocks: Dict[int, List[int]] = defaultdict(list)
    # fl_ex_num2blocks: Dict[int, List[int]] = defaultdict(list)
    # aa_ex_num2blocks: Dict[int, List[int]] = defaultdict(list)
    init0 = datetime.now()
    sorted_block_keys: List[str] = sorted(
        blocks.keys(), key=lambda x: (blocks[x][0], blocks[x][1])
    )
    block_sort_point = datetime.now()
    stderr.write(f'Block sorting time: {str(block_sort_point - init0)}\n')
    tstart: int = blocks[sorted_block_keys[0]][0]
    tstop: int = blocks[sorted_block_keys[-1]][1]
    qstart: int = min(x[2] for x in blocks.values())
    qstop: int = max(x[3] for x in blocks.values())
    curr_block: int = 0 ## a 'smart pointer' to track the current block
    sorted_exon_keys: List[int] = sorted(
        exons.keys(), key=lambda x: exons[x].start
    )

    for e_num in sorted_exon_keys:
        init = datetime.now()
        exon: Exon = exons[e_num]
        extra_flank: int = int(exon.length() * extra_flank_modifier)
        min_coord, max_coord = 0, 0
        for b_pointer, b_num in enumerate(sorted_block_keys[curr_block:], start=curr_block):
            this_block: List[str] = blocks[b_num]
            if this_block[0] > exon.stop:
                inter_point = datetime.now()
                stderr.write(f'Intersection time: {str(inter_point - init)}\n')
                break
            if not (this_block[1] - this_block[0]):
                continue
            if not (this_block[3] - this_block[2]):
                continue
            inter_size: int = intersection(exon.start, exon.stop, *this_block[:2])
            if inter_size < 0:
                continue
            curr_block = b_pointer
            ex_num2blocks[e_num].append(b_num)
        loop_end_point = datetime.now()
        stderr.write(f'Regular loop end time: {str(loop_end_point - init)}\n')

        ## infer the exons coordinates from the chain blocks

        ## case 0: if there are no blocks of non-zero length corresponding
        ## to the exons, mark it as missing and proceed further
        init2 = datetime.now()
        if not ex_num2blocks[e_num]:
            click.echo(
                f'Exon {exon.num} is missing from the chain'
            ) if verbose else None
            if exon.stop < t_chain_start or exon.start > t_chain_stop:
                # print(f'{t_chain_start=}, {t_chain_stop=}, {qstart=}, {qstop=}, {exon.start=}, {exon.stop=}')
                out_of_chain.add(e_num)
            # else:
            missing.add(e_num)
            continue
        up_block_name: str = ex_num2blocks[e_num][0]
        up_block: List[int] = blocks[up_block_name] if up_block_name else None
        down_block_name: str = ex_num2blocks[e_num][-1]
        down_block: List[int] = blocks[down_block_name] if down_block_name else None

        ## case 1: exon is enclosed within only one block
        if up_block_name == down_block_name:
            ## case 1.1: enclosing chain part is an interblock gap
            if '_' in up_block_name:
                click.echo(
                    f'Exon {exon.num} corresponds to an unaligned chain gap'
                ) if verbose else None
                # missing.add(e)
                if MAX_CHAIN_GAP_SIZE >= abs(up_block[3] - up_block[2]) >= exon.length(): ## TODO: Needs testing
                    min_coord, max_coord = up_block[2:]
                    gap_located.add(e_num)
                else:
                    missing.add(e_num)
            ## case 1.2: enclosing chain part is an aligned block
            else:
                min_coord = up_block[2]
                max_coord = up_block[3]
                trim_up: int  = exon.start - up_block[0]
                trim_down: int = up_block[1] - exon.stop
                click.echo(
                    f'Exon {exon.num} is fully enclosed within '
                    f'block {up_block_name}'
                    if trim_up >= 0 and trim_down >= 0  else
                    f'Exon {exon.num} is partially covered by a single '
                    f'block {up_block_name}'
                ) if verbose else None
                min_coord += trim_up if codirected else trim_down
                max_coord -= trim_down if codirected else trim_up
        ## case 2: there are more than one block corresponding to the exon
        else:
            # print(f'FLAG {blocks=}')
            # print(f'FLAG {up_block_name=}, {down_block_name=}')
            new_up_block_name, next_up_block_name = (
                up_block_name.split('_') if '_' in up_block_name else ('', '') ## TODO: replace with a concise function split_gap_name()
            )
            new_up_block: List[int] = blocks[new_up_block_name] if new_up_block_name else None
            next_up_block: List[int] = blocks[next_up_block_name] if next_up_block_name else None
            next_down_block_name, new_down_block_name = (
                down_block_name.split('_') if '_' in down_block_name else ('', '')
            )
            new_down_block: List[int] = blocks[new_down_block_name] if new_down_block_name else None
            next_down_block: List[int] = blocks[next_down_block_name] if next_down_block_name else None
            ## case 2.1: any of two marginal blocks are interblock gaps
            ## in this case, extend the exon coordinate to the next block
            if new_up_block_name:
                side: str = 'Left' if codirected else 'Right'
                click.echo(
                    f'{side} flank for exon {exon.num} is a chain gap; '
                    f'extending coordinates to block {new_up_block_name}'
                ) if verbose else None
                trim_up: int = next_up_block[0] - exon.start
                dangling_up: bool = trim_up > (
                    acc_flank if codirected else donor_flank
                )
                extra_flank_up: int = extra_flank * dangling_up
                if codirected:
                    min_coord = max(0, next_up_block[2] - trim_up - extra_flank_up)
                else:
                    max_coord = next_up_block[3] + trim_up + extra_flank_up
            if new_down_block_name:
                side: str = 'Right' if codirected else 'Left'
                click.echo(
                    f'{side} flank for exon {exon.num} is a chain gap; '
                    f'extending coordinates to block {new_down_block_name}'
                )
                click.echo(','.join(map(str, blocks[new_down_block_name])))
                trim_down: int = next_down_block[1] - exon.stop
                dangling_down: bool = trim_down < -1 * (
                    donor_flank if codirected else acc_flank
                )
                extra_flank_down: int = extra_flank * dangling_down
                if codirected:
                    max_coord = max(0, next_down_block[3] - trim_down - extra_flank_down)
                else:
                    min_coord = next_down_block[2] + trim_down + extra_flank_down
            ## case 2.2: at least one of the marginal blocks is a true
            ## aligned block; trim the coordinates by the exon coordinates
            ## from the reference annotation
            click.echo(
                f'Exon {exon.num} intersects the following blocks: '
                f'{up_block_name}, {down_block_name}'
            ) if verbose else None
            trim_up: int = up_block[0] - exon.start
            trim_down: int = down_block[1] - exon.stop
            dangling_up: bool = trim_up > (
                acc_flank if codirected else donor_flank
            )
            dangling_down: bool =  trim_down < -1 * (
                donor_flank if codirected else acc_flank
            )
            # if codirected:
                # print(f'{blocks[up_block_name][2]=}, {blocks[down_block_name][3]=}')
            # else:

                # print(f'{blocks[down_block_name][2]=}, {blocks[up_block_name][3]=}')
            # print(f'{exon.start = }, {exon.stop = }, {trim_up = }, {trim_down = }')
            if not min_coord:
                side: str = 'left' if codirected else 'right'
                click.echo(
                    f'Trimming exon {exon.num} coordinates '
                    f'from the {side} side'
                ) if verbose else None
                min_coord = max(
                    up_block[2] - trim_up if codirected else
                    down_block[2] + trim_down, 0
                )
            if not max_coord:
                side: str = 'right' if codirected else 'left'
                click.echo(
                    f'Trimming exon {exon.num} coordinates '
                    f'from the {side} side'
                ) if verbose else None
                max_coord = max(
                    down_block[3] - trim_down if codirected else
                    up_block[3] + trim_up, 0
                )
        ## a short sanity check: order the coordinates
        min_coord, max_coord = sorted([min_coord, max_coord])
        ## exon loci longer than five times the exon length are
        ## suspected for chain distruption
        if (max_coord - min_coord) / exon.length() >= 5 and e_num not in gap_located:
            ## calculate the alignment sum on each side of the expected locus
            ## select the side which accumulates more aligned nucleotides,
            ## use the respective border coordinate as an anchor
            upstream_blocks: List[Tuple[int]] = [
                y for x, y in blocks.items() if y[3] < min_coord and '_' not in x
            ]
            downstream_blocks: List[Tuple[int]] = [
                y for x, y in blocks.items() if y[2] > max_coord and '_' not in x
            ]
            upstream_aln_sum: int = (
                0 if not upstream_blocks else sum([x[3] - x[2] for x in upstream_blocks])
            )
            downstream_aln_sum: int = (
                0 if not downstream_blocks else sum([x[3] - x[2] for x in upstream_blocks])
            )
            if upstream_aln_sum > downstream_aln_sum:
                max_coord = min_coord + int(1.1 * exon.length())
            else:
                min_coord = max(0, max_coord - int(1.1 * exon.length()))
        click.echo(
            f'Resulting coordinates for exon {exon.num} are: '
            f'{min_coord}-{max_coord}'
        ) if verbose else None
        min_coord = nn(min(min_coord, q_chrom_size))
        max_coord = nn(min(max_coord, q_chrom_size))
        exons2coordinates[e_num] = (min_coord, max_coord)
        # dangling_ends[e_num] = cov_status
        exon_proc_point = datetime.now()
        stderr.write(f'Exon processing after intersection time: {str(exon_proc_point - init2)}\n')
        stderr.write(f'Full exon processing time: {str(exon_proc_point - init)}\n')

    spanning_chain: bool = all(e in missing or e in gap_located for e in exons)
    # print(f'{out_of_chain=}')
    # print(f'GAP LOCATED: {gap_located=}')

    init3 = datetime.now()
    mapper: Exon2BlockMapper = Exon2BlockMapper(
        chainid,
        # chain_prob,
        # *global_block,
        tstart,
        tstop,
        qstart,
        qstop,
        t_chrom,
        q_chrom,
        t_strand,
        q_strand,
        blocks,
        # ex_num2blocks,
        # fl_ex_num2blocks,
        exons2coordinates,
        out_of_chain,
        missing,
        gap_located,
        dangling_ends,
        spanning_chain
        # aa_sat
    )
    instance_prep_time = datetime.now()
    stderr.write(f'Instance preparation time: {str(instance_prep_time - init3)}\n')

    return mapper#, curr_block


def detect_shifted_exons(segment: Segment, chain: Exon2BlockMapper) -> Any:
    """
    Detect if segment and chain are shifted against each other
    Return indices of exons for which chain and segment locations coincide
    Note that:
    1) chain/segment can be up- and downstream shifted at the same time
       ==chain==chain============chain===chain================chain===chain=====
       ===============seg==seg===segme===segme===seg==segm======================
       ===============exon==exo==exon====exon====ex===exon======================
    2) chain and segment can be mutually shifted from different sides
       =========exon-1==exon-1========exon==exon==========exon+1==exon+1========
       =========chain======================chain================================
       ===============================segm========================segment=======
    3) chain and segment CANNOT be mutually upshifted or mutually downshifted
    """
    min_seg: int = segment.exons.min() if segment.strand else segment.exons.max()
    max_seg: int = segment.exons.max() if segment.strand else segment.exons.min()
    all_blocks: List[int] = sorted([x for x in chain.e2c if x not in chain.missing])
    ## inferring which strand the chain-projected exons lie on might be complicated, use this proxy
    strand: bool = chain.e2c[max(all_blocks)][0] > chain.e2c[min(all_blocks)][1]
    min_block: int = min(all_blocks) if strand else max(all_blocks)
    max_block: int = max(all_blocks) if strand else min(all_blocks)
    ## define max for the progenitor transcript
    max_exon: int = max(min_seg.num, max_seg.num, min_block, max_block) + 1

    segment_upstream_shift: List[int] = [
        x for x in segment.exons if
        segment.exons[x].stop <= chain.get_exon_coords(min_block)[2] and
        (x > min_block if strand else x < min_block)
    ]
    if segment_upstream_shift:
        segment_upstream_shift.extend(
            list(range(1, min(segment_upstream_shift))) if chain.codirected()
            else list(range(max(segment_upstream_shift) + 1, max_exon))
        )

    segment_downstream_shift: List[int] = [
        x for x in segment.exons if
        segment.exons[x].start >= chain.get_exon_coords(max_block)[3] and
        (x < max_block if strand else x > max_block)
    ]
    if segment_downstream_shift:
        segment_downstream_shift.extend(
            list(range(max(segment_downstream_shift) + 1, max_exon))
            if chain.codirected() else list(range(1, min(segment_downstream_shift)))
        )

    chain_upstream_shift: List[int] = [
        x for x in all_blocks if chain.get_exon_coords(x)[3] < min_seg.start and
        (x > min_seg.num if segment.strand else x < min_seg.num)
    ]
    if chain_upstream_shift:
        chain_upstream_shift.extend(
            list(range(1, min(chain_upstream_shift))) if segment.strand
            else list(range(max(chain_upstream_shift) + 1, max_exon))
        )

    chain_downstream_shift: List[int] = [
        x for x in all_blocks if chain.get_exon_coords(x)[2] >= max_seg.stop and
        (x < max_seg.num if segment.strand else x > max_seg.num)
    ]
    if chain_downstream_shift:
        chain_downstream_shift.extend(
            list(range(max(chain_downstream_shift) + 1, max_exon))
            if segment.strand else list(range(1, min(chain_downstream_shift)))
        )
    ## think of proper return value
    return (
        sorted(segment_upstream_shift),
        sorted(segment_downstream_shift),
        sorted(chain_upstream_shift),
        sorted(chain_downstream_shift)
    )


def dump_for_cesar(
    que_names: List[str],
    que_seqs: List[str],
    ref_names: List[str],
    ref_seqs: List[str],
    path: str
) -> None:
    """
    Writes input for CESAR2.0 to <path>
    Arguments are:
    :que_names: is a list of query sequence names (with donor-acceptor profiles)
    :que_seq: is a list of query sequences
    :ref_names: is a list of reference sequence names (with donor-acceptor profiles)
    :ref_seqs: is a list of reference sequences
    :path: is a file path to write the output to
    NOTE that sequences and their names should be presented in the corresponding order
    Notation used here assumes that queries are written to the file before the references
    """
    if len(que_names) != len(que_seqs):
        raise AttributeError('Inconsistent number of query names and sequences')
    if len(ref_names) != len(ref_seqs):
        raise AttributeError('Inconsistent number of reference names and sequences')
    out_line: str = ''
    for i in range(len(que_names)):
        name: str = fa_name(que_names[i])
        out_line += name + '\n' + que_seqs[i] + '\n'
    out_line += '####' + '\n'
    for j in range(len(ref_names)):
        name: str = fa_name(ref_names[j])
        out_line += name + '\n' + ref_seqs[j] + '\n'

    if path is None:
        return out_line

    with open(path, 'w') as h:
        h.write(out_line)


def parse_cesar(cesar_out: str) -> List[Union[str, int]]:
    """
    Parses raw CESAR2.0 output; returns clipped aligned sequences and
    alignment coordinates in the reference sequence.

    :cesar_out: is a string representation of CESAR results
    (i.e., output file .read() content or CESAR output to stdout)

    NOTE that TOGA uses CESAR2.0 for single-reference alignments, therefore
    the function assumes that the output contains only four lines.
    """
    data: List[str] = cesar_out.split('\n')
    ref_line: str = data[1]
    que_line: str = data[3]
    coords: List[Tuple[int, int]] = []
    aligned: bool = False
    deleted_intron: bool = False
    q_gap: int = 0
    ref_seq: str = ''
    que_seq: str = ''
    trailing_intron_seq: str = ''
    exon_num = 1
    for i in range(len(ref_line)):
        r_pos: str = ref_line[i]
        q_pos: str = que_line[i]
        if r_pos == '>':# and deleted_intron: ## intron deletion marked by '>' signs in reference; assign it to the next exon
            trailing_intron_seq += q_pos
            deleted_intron = True
        elif (q_pos.isupper() or q_pos == '-') and not aligned: ## aligned segment start
            start: int = i - q_gap
            aligned = True
            deleted_intron = False
            # ref_seq += r_pos
            # que_seq += q_pos
        if aligned:
            if q_pos.islower() or i == len(ref_line) - 1 or r_pos == '>': #(r_pos == '>' and not deleted_intron): ## exon sequence stop
                stop: int = i - 1 - q_gap
                if i == len(ref_line) - 1 and (q_pos.isupper() or q_pos == '-'):
                    ref_seq += r_pos
                    que_seq += q_pos
                ## save the splice site sequences from query
                acceptor_site: str = que_line[max(0, start - 2):start]
                donor_site: str = que_line[min(i, len(que_line)):i+2]
                coords.append(
                    (start, stop, ref_seq, que_seq, acceptor_site, donor_site, trailing_intron_seq)
                )
                exon_num += 1
                aligned = False
                ref_seq = ''
                que_seq = ''
                trailing_intron_seq = ''
                if r_pos == '>':
                    deleted_intron = True
                    trailing_intron_seq += q_pos
            else:
                if q_pos == '-': ## alignment gap in query -> to be subtracted from the sequence stop
                    q_gap += 1
                ref_seq += r_pos
                que_seq += q_pos

    return coords


def get_trailing_from_cesar_output(cesar_out: str) -> str:
    """
    Returns trailing sequence from raw CESAR output string
    """
    data: List[str] = cesar_out.split('\n')
    ref_line: str = data[1]
    que_line: str = data[3]
    trailing: str = ''
    for i, r_pos in reversed(list(enumerate(ref_line))):
        if r_pos != ' ':
            break
        q_pos: str = que_line[i]
        trailing += q_pos

    return trailing[::-1]


def scan_seq_till_stop(seq: str) -> str:
    """
    For a given string, report the sequence preceding the first stop codon found;
    report an empty sequence if no stop codons were found
    """
    for i in range(len(seq)):
        if seq[i-2:i+1].upper() in STOPS:
            return seq[:i+1].upper()
    return ''


def extend_aln_till_stop(entry: CesarExonEntry, seq: str) -> CesarExonEntry:
    """
    Given a CESAR2.0 output storage entity and a downstream sequence,
    inserts the extension with the new found stop codon into the alignment
    """
    extension: str = scan_seq_till_stop(seq)
    adj_ext: str = '-' * (3 - len(extension) % 3) + extension
    entry.seq = entry.seq[:-3] + '-' * (len(adj_ext) - 3) + entry.seq[-3:]
    entry.hit_seq = entry.hit_seq[:-3] + adj_ext
    entry.stop += len(extension)
    return entry


def check_codon(codon: str) -> str:
    """Masks stop and frameshift-containing codons, returns input codon otherwise"""
    if codon in STOPS:
        # mask stop
        return XXX_CODON
    elif codon == GAP_CODON:
        # whole codon is missing
        return codon
    elif codon.count("-") > 0:
        # FS
        return XXX_CODON
    elif len(codon) % 3:
        # FS; rarely occurs when SpliceAI correction restores the shifted frame
        return XXX_CODON
    else:  # normal codon
        return codon


def get_blosum_score(
    aa1: str,
    aa2: str,
    matrix: Dict[str, Dict[str, int]],
    ins_pen: int = INS_PEN,
    del_pen: int = DEL_PEN
) -> int:
    """
    Returns the protein alignment score for a pair of amino acids.
    Arguments are:
    :aa1: and :aa2: are amino acid symbols to compare
    :matrix: is a dict-of-dicts representation of the amino acid scoring matrix
    :ins_pen: and :del_pen: are insertion and deletion penalties, respectively
    """
    if aa1 == '-' and aa2 == '-':
        return 0
    if aa1 != '-' and aa2 != '-':
        return matrix[aa1][aa2]
    elif aa1 == '-' and aa2 != '-':
        return del_pen
    elif aa1 != '-' and aa2 == '-':
        return ins_pen
    else:
        raise ValueError(
            'CESAR2.0 output is likely corrupted'
        )


def process_codon_pair(ref_codon: str, query_codon: str) ->  Tuple[Tuple[str, str]]:
    """
    Codon-wise version of extract_codon_data from CESAR_wrapper.py .
    Splits a pair of aligned nucleotide sequnces into codons for future codon alignment,
    assuming that the reference sequence contains exactly one full codon
    """
    if len(ref_codon) == len(query_codon) == 3:
        ## the simplest case; return codons as they are
        return ((ref_codon, query_codon),)
    elif len(ref_codon) % 3 == 0:
        ## frame-preserving indels; split the sequences into trimers
        ## and return reference-query codons as a zip object
        ref_subcodons: List[str] = parts(ref_codon, 3)
        query_subcodons: List[str] = parts(query_codon, 3)
        return zip(ref_subcodons, query_subcodons)
    elif len(ref_codon) < 3:
        ## reference sequence provided does not contain any full codons;
        ## an unlikely situation in the current scheme but let's keep
        ## this functionality for later convenience
        return ((f'{XXX_CODON}|{ref_codon}', f'{XXX_CODON}|{query_codon}'),)
    else:
        ## codon partitioning contains uncompensated indels,
        ## most likely insertions in the query sequence
        ## return the latter three nucleotides in each sequence as they come
        ## and the leading nucleotide(s) as was_not_aligned codons
        frame_out: int = len(ref_codon) % 3
        ref_int: str = ref_codon[frame_out:]
        query_int: str = query_codon[frame_out:]
        fs_ref: str = f'{GAP_CODON}|{ref_codon[:frame_out]}'
        fs_query: str = f'{XXX_CODON}|{query_codon[:frame_out]}'
        if len(ref_codon) < 6:
            return ((fs_ref, fs_query), (ref_int, query_int))
        else:
            return (
                (fs_ref, fs_query), *zip(parts(ref_int, 3), parts(query_int, 3))
            )


def process_and_translate(ref_codon: str, query_codon: str) -> Tuple[Tuple[str, str]]:
    """
    Codon-wise version of translate_codons from CESAR_wrapper.py
    Infers information on encoded amino acids from a pair of aligned nucleotide
    sequences, assuming that the reference sequence contains exactly one full codon
    """
    ref_gapless: str = ref_codon.replace('-', '')
    query_gapless: str = query_codon.replace('-', '')
    ## for reference, the function scraps amino acid data from all the nucleotides present
    ## unlike reference, query sequence may contain an arbitrary number of
    ## codons fo varying completeness; any difference in codon number should be compensated
    ## by leading deletions
    query_codon_number: int = len(query_codon) // 3
    ref_aa: List[str] = list(
            '-' * (query_codon_number - 1) + (
            AA_CODE.get(ref_gapless, 'X') if ref_gapless else '-'
        )
    )
    ## for query, codons are translated ONLY if the frame is preserved
    ## TODO: I disagree with this notion, and this should be discussed with both Bogdan and Michael
    ## imagine the following alignment:
    ## ref A A - A
    ## que A - A A
    ## we surely rescue the reference codon
    ## but the query encodes a complementary amino acid since the frameshifts are properly compensated
    if len(query_codon) % 3:
        query_aa: List[str] = list('X' * query_codon_number)
    else:
        ## Once again, I do not agree; consider the following
        ## A - A - - A
        ## - B - B B -
        ## why is it -B- + BB- and not BBB?
        query_aa: List[str] = [AA_CODE.get(x, 'X') for x in parts(query_codon, 3)]

    return zip(ref_aa, query_aa)


def assess_exon_quality(nuc: float, blosum: float) -> str:
    """
    Assesses exon quality based on nucleotide and BLOSUM identity.
    Classes are:
    * High quality (HQ): PID >=
    * Average quality (AQ):
    * Low quality (LQ):
    * NA: exon has identity below all thresholds, most likely missing or deleted
    """
    if nuc >= HQ_PID and blosum >= HQ_BLOSUM:
        return 'HQ'
    elif nuc >= A_T_PID and blosum >= A_T_BLOSUM:
        return 'AQ'
    elif nuc >= LO_T_PID and blosum >= LO_T_BLOSUM:
        return 'LQ'
    else:
        return 'NA'


def get_d_runs(exon_states: Dict[int, str]) -> List[List[int]]:
    """
    Of the ordered exon presence state list, select consecutively deleted exons'
    streaks
    """
    out_list: List[List[int]] = []
    curr_group: List[int] = []
    for i, state in exon_states.items():
        if state == 'D':
            curr_group.append(i)
        elif curr_group:
            out_list.append(curr_group)
            curr_group = []
    if curr_group:
        out_list.append(curr_group)
    return out_list


def get_affected_exon_threshold(exon_num: int) -> int:
    """
    Returns the maximum number of affected exons to be allowed for a projection
    to be classified as uncertain rather than confirmed loss:
    1) For single-exon transcripts, this is (obviously) 1;
    2) For transcripts consisting of ten or less exons, this makes 2;
    3) For transcripts longer than 10 exons, the maximum portion is
       20% of exon number;
    """
    if exon_num == 1:
        return 1
    if exon_num <= 10:
        return 2
    return exon_num / 5


def dfs(
    groups: Dict[int, List[CesarExonEntry]],
    group: int,
    path: List[CesarExonEntry]#,
    # strand: bool
) -> Iterable[List[CesarExonEntry]]:
    """
    A depth-first search-like function for reconstructing all the possible
    combinations (=alternative segments) from the CesarExonEntry exon objects.
    Input it:
    :groups: is a dictionary of group_number:[CesarExonEntry] pairs
        where value contains all the exon entries belonging to the group
    :group: is a current exon group (pass 1 if running in your own code)
    :path: is an existing combination of exon entries from the previous groups
    :strand: is a strand the segment is located on
    """
    last_group: bool = group == max(groups)
    for alt in groups[group]:
        new_path: List[CesarExonEntry] = path.copy()
        # for entry in sorted(alt, key=lambda x: x.num if strand else -x.num):
        for entry in sorted(alt, key=lambda x: x.num):
            new_path.append(entry)
        if last_group:
            yield new_path
        else:
            new_group: int = group + 1
            yield from dfs(groups, new_group, new_path)#, strand)


def ddfs(
    groups: Dict[int, List[RawCesarOutput]],
    group: int,
    path: List[RawCesarOutput]#,
    # strand: bool
) -> Iterable[List[RawCesarOutput]]:
    """
    A PROVISIONAL FUNCTION VERSION adapted for the updated ProcessedSegment from
    processed_segment.py; to be merged with the primary function

    A depth-first search-like function for reconstructing all the possible
    combinations (=alternative segments) from the CesarExonEntry exon objects.
    Input it:
    :groups: is a dictionary of group_number:[CesarExonEntry] pairs
        where value contains all the exon entries belonging to the group
    :group: is a current exon group (pass 1 if running in your own code)
    :path: is an existing combination of exon entries from the previous groups
    :strand: is a strand the segment is located on
    """
    last_group: bool = group == max(groups)
    for alt in groups[group]:
        new_path: List[CesarExonEntry] = path.copy()
        # for entry in sorted(alt, key=lambda x: x.num if strand else -x.num):
        # for entry in sorted(alt, key=lambda x: x.num):
        #     new_path.append(entry)
        new_path.append(alt)
        if last_group:
            yield new_path
        else:
            new_group: int = group + 1
            yield from ddfs(groups, new_group, new_path)#, strand)


def cesar_memory_check(block_lengths: List[int], query_length: int) -> float:
    """
    Estimates memory consumption for CESAR2.0 . Arguments are:
    :block_lengths: is a list of reference (exon) sequence sizes
    :query_length: is a length of the query sequence
    """
    num_states, rlength, extra = 0, 0, 100000
    for block_length in block_lengths:
        num_codons = block_length // 3
        num_states += 6 + 6 * num_codons + 1 + 2 + 2 + 22 + 6
        rlength += block_length
    MEM: int = (
        (num_states * 4 * 8)
        + (num_states * query_length * 4)
        + (num_states * 304)
        + (2 * query_length + rlength) * 8
        + (query_length + rlength) * 2 * 1
        + extra
    )
    # convert bytes to GB
    GB: float = MEM / 1000000000
    return GB



##############################
### Below are the functions directly borrowed from the TOGA 1.0 CESAR_wrapper.py
### These are likely to be replaced with Cython/alternative C-related solutions
##############################
chain_coords_conv_lib_path: str = os.path.join(LOCATION, 'chain_coords_converter_slib.so')
ch_lib = ctypes.CDLL(chain_coords_conv_lib_path)
ch_lib.chain_coords_converter.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_char_p),
]
ch_lib.chain_coords_converter.restype = ctypes.POINTER(ctypes.c_char_p)

def get_chain(chain_file, chain_id):
    """Return chain string according the parameters passed."""
    chain = None  # to calm IDE down
    if chain_file.endswith(".bst"):
        # we have bdb file; extract with BDB extractor
        chain = chain_extract_id(chain_file, chain_id)
        return chain
    elif chain_file.endswith(".gz"):  # a gzipped chain file was given
        # gzip and redirect scteam to chain_filter_by_id binary
        extract_by_id_cmd = (
            f"gzip -dc {chain_file} | ./modules/chain_filter_by_id stdin {chain_id}"
        )
        try:  # check that output is OK
            chain = subprocess.check_output(extract_by_id_cmd, shell=True).decode(
                "utf-8"
            )
        except subprocess.CalledProcessError:
            # die if the command died
            die(
                "Error! Process {extract_by_id_cmd} died! Please check if input data is correct",
                1,
            )
        return chain

    else:  # just a chain file, extract the chain we need
        # the same as above, but without gzip stage
        extract_by_id_cmd = f"./modules/chain_filter_by_id {chain_file} {chain_id}"
        try:  # also check if output is OK, die otherwise
            chain = subprocess.check_output(extract_by_id_cmd, shell=True).decode(
                "utf-8"
            )
        except subprocess.CalledProcessError:
            die(
                "Error! Process {extract_by_id_cmd} died! Please check if input data is correct",
                1,
            )
        return chain


def range_corrector(g_range):
    """Swap start and end if start > end."""
    chrom, start_end = g_range.split(":")
    start_end_split = start_end.split("-")
    start, end = int(start_end_split[0]), int(start_end_split[1])
    if start < end:
        return g_range
    else:
        return f"{chrom}:{end}-{start}"


def chain_cut(chain_str, gene_range, gene_flank, extra_flank=0):
    """Call chain_cut binary.

    Project reference gene coordinates to query through a chain.
    Also add flanks if shift is > 0.
    """
    # need to get genomic region for the gene
    # also need to translate python data types to C
    # to call the shared library; I do it 2 times here
    # for shift = 0 and shifts = 2 (add flanks around gene)
    c_chain = ctypes.c_char_p(chain_str.encode())
    c_shift_2 = ctypes.c_int(2)
    c_shift_0 = ctypes.c_int(0)
    granges_num = 1
    c_granges_num = ctypes.c_int(granges_num)  # we need only one grange to analyze
    granges_arr = (ctypes.c_char_p * (granges_num + 1))()  # granges_num + 1
    granges_bytes = [gene_range.encode("utf-8")]
    # need to do this tricks to pass strings array to C
    granges_arr[:-1] = granges_bytes
    granges_arr[granges_num] = None
    raw_ch_conv_s2 = ch_lib.chain_coords_converter(
        c_chain, c_shift_2, c_granges_num, granges_arr
    )
    chain_coords_conv_out_s2 = []  # keep lines here
    # convert C output to python-readable type
    for i in range(granges_num + 1):
        chain_coords_conv_out_s2.append(raw_ch_conv_s2[i].decode("utf-8"))
    # chain', 'chr5', '+', '137889395', '148245211', 'chr18', '+', '34409342', '44120958
    chain_data = chain_coords_conv_out_s2[0].split(" ")
    t_strand = True if chain_data[2] == "+" else False
    q_strand = True if chain_data[7] == "+" else False
    t_size = int(chain_data[3])
    q_size = int(chain_data[8])

    # re-define arrays to avoid segfault
    c_chain = ctypes.c_char_p(chain_str.encode())
    granges_arr = (ctypes.c_char_p * (granges_num + 1))()  # granges_num + 1
    granges_bytes = [gene_range.encode("utf-8")]
    granges_arr[:-1] = granges_bytes
    granges_arr[granges_num] = None

    raw_ch_conv_s0 = ch_lib.chain_coords_converter(
        c_chain, c_shift_0, c_granges_num, granges_arr
    )
    chain_coords_conv_out_s0 = []  # keep lines here

    # convert C output to python-readable type
    for i in range(granges_num + 1):
        chain_coords_conv_out_s0.append(raw_ch_conv_s0[i].decode("utf-8"))
    # another approach to detect range
    # sometimes blocks go so far
    # ------------------genegene-------------------
    # block-------------blockblock------------block

    # to avoid very huge query sequences program controls it's size
    search_region_shift_str = range_corrector(
        chain_coords_conv_out_s2[1].split("\t")[1]
    )
    search_region_abs_str = range_corrector(chain_coords_conv_out_s0[1].split("\t")[1])

    chrom = search_region_shift_str.split(":")[0]
    search_reg_shift = [
        int(x) for x in search_region_shift_str.split(":")[1].split("-")
    ]
    search_reg_abs = [int(x) for x in search_region_abs_str.split(":")[1].split("-")]
    search_reg_flanked = [
        search_reg_abs[0] - gene_flank,
        search_reg_abs[1] + gene_flank,
    ]

    # define actual starts and ends
    act_start = (
        search_reg_shift[0]
        if search_reg_shift[0] > search_reg_flanked[0]
        else search_reg_flanked[0]
    )
    act_end = (
        search_reg_shift[1]
        if search_reg_shift[1] < search_reg_flanked[1]
        else search_reg_flanked[1]
    )

    if extra_flank > 0:  # add extra flanks if required
        act_start = act_start - extra_flank if act_start - extra_flank > 0 else 0
        act_end = (
            act_end + extra_flank if act_end + extra_flank < q_size else q_size - 1
        )

    act_search_range = f"{chrom}:{act_start}-{act_end}"
    # ext_search_range = f"{chrom}:{}-{}"
    del raw_ch_conv_s0  # not sure if this is necessary
    del raw_ch_conv_s2  # but just in case
    return (
        act_search_range,
        search_region_shift_str,
        (t_strand, t_size, q_strand, q_size),
    )


def fast_seq_id(ref: str, query: str) -> float:
    """
    An ad-hoc identity calculator for intron gain-bound CESAR run results
    """
    len_sum: int = 0
    match_sum: int = 0
    for n1, n2 in zip(ref, query):
        ref_is_coding: bool = n1.isalpha() or n1 == '-'
        query_is_upper: bool = n2.isupper() or n2 == '-'
        if not (ref_is_coding and query_is_upper):
            # print(f'{n1=}, {n2=}, IGNORED')
            continue
        len_sum += 1
        match_sum += int(n1.upper() == n2)
        # print(f'{n1=}, {n2=}')
    return  match_sum / len_sum
