#!/usr/bin/env python3

"""
A collection of executables related to data preprocessing for CESAR
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from modules.cesar_wrapper_constants import (
    MAX_CHAIN_GAP_SIZE, MAX_CHAIN_INTRON_LEN, MIN_INTRON_LENGTH, STOPS
)
from collections import defaultdict
from common import chain_extract_id, parts
from dataclasses import dataclass
from shared import intersection, nn
from typing import Dict, Iterable, List, Set, Tuple, Union

import ctypes
import logging
import subprocess

__author__ = 'Yury V. Malovichko'
__version__ = '1.0'
__year__ = '2023'

@dataclass
class Exon:
    """
    Minimal data class for exon block storage
    """
    __slots__ = ('num', 'start', 'stop')
    num: int
    start: int
    stop: int

    def length(self) -> int:
        return self.stop -  self.start
    
    def coords(self) -> Tuple[int, int]:
        return tuple(sorted((self.start, self.stop)))

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
    __slots__ = (
        'chainid', 'prob',
        'tstart', 'tstop', 'qstart', 'qstop',
        'tchrom', 'qchrom', 'tstrand', 'qstrand',
        'blocks', #'e2b', 'fe2b',
        'e2c', 'out_of_chain', 'missing', 'gap_located',
        'spanning_chain', 'init_cov'
    )
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
    init_cov: Dict[int, int]
    # dangling: Set[int]
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
    logger: logging.Logger,
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
    # missing: Set[int] = {x for x in exons if not ex_num2blocks.get(x, [])}
    missing: Set[int] = set()
    out_of_chain: Set[int] = set()
    # out_of_chain: Set[int] = set()
    gap_located: Set[int] = set()
    exon2init_cov: Dict[int, Union[int, None]] = {} ## for coding block coverage before trimming/extension
    exons2coordinates: Dict[int, Tuple[int]] = {} ## for clipped/extrapolated projection coordinates

    ## then, find the exon-to-block intersections blockwise
    ex_num2blocks: Dict[int, List[int]] = defaultdict(list)
    # fl_ex_num2blocks: Dict[int, List[int]] = defaultdict(list)
    # aa_ex_num2blocks: Dict[int, List[int]] = defaultdict(list)
    sorted_block_keys: List[str] = sorted(
        blocks.keys(), key=lambda x: (blocks[x][0], blocks[x][1])
    )
    if len(sorted_block_keys) == 1 and '_' in sorted_block_keys[0]:
        upstream_chain_gap: Union[str, None] = sorted_block_keys[0]
        downstream_chain_gap: Union[str, None] = sorted_block_keys[0]
    else:
        if '_' in sorted_block_keys[0]:
            upstream_chain_gap: Union[str, None] = sorted_block_keys.pop(0)
        else:
            upstream_chain_gap: Union[str, None] = None
        if '_' in sorted_block_keys[-1]:
            downstream_chain_gap: Union[str, None] = sorted_block_keys.pop(-1)
        else:
            downstream_chain_gap: Union[str, None] = None
    
    # print(f'{blocks=}')
    tstart: int = blocks[sorted_block_keys[0]][0]
    tstop: int = blocks[sorted_block_keys[-1]][1]
    qstart: int = min(x[2] for x in blocks.values())
    qstop: int = max(x[3] for x in blocks.values())
    curr_block: int = 0 ## a 'smart pointer' to track the current block
    sorted_exon_keys: List[int] = sorted(
        exons.keys(), key=lambda x: exons[x].start
    )
    ex2del: Dict[int, bool] = defaultdict(bool)

    for e_num in sorted_exon_keys:
        exon: Exon = exons[e_num]
        extra_flank: int = int(exon.length() * extra_flank_modifier)
        coords: Tuple[int, int] = exon.coords()
        min_coord, max_coord = 0, 0
        init_cov: int = 0
        for b_pointer, b_num in enumerate(sorted_block_keys[curr_block:], start=curr_block):
            this_block: List[str] = blocks[b_num]
            if this_block[0] > exon.stop:
                break
            if not (this_block[1] - this_block[0]):
                continue
            inter_size: int = intersection(exon.start, exon.stop, *this_block[:2])
            if inter_size < 0:
                continue
            if not (this_block[3] - this_block[2]):
                ex2del[e_num] = True
                continue
            curr_block = b_pointer
            ex_num2blocks[e_num].append(b_num)

        ## infer the exons coordinates from the chain blocks

        ## case 0: if there are no blocks of non-zero length corresponding
        ## to the exons, mark it as missing and proceed further
        if not ex_num2blocks[e_num]:
            logger.info(
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
                logger.info(
                    f'Exon {e_num} corresponds to an unaligned chain gap'
                ) if verbose else None
                # missing.add(e)
                if up_block[1] - up_block[0] == 0 or up_block[3] - up_block[2] == 0:
                    logger.info(
                        f'Exon {e_num} corresponds to a chain insertion; marking it as missing'
                    )
                    missing.add(e_num)
                elif MAX_CHAIN_GAP_SIZE >= abs(up_block[3] - up_block[2]) >= exon.length(): ## TODO: Needs testing
                    logger.info(
                        f'Exon {e_num} corresponds to a chain gap of plausible size'
                    )
                    min_coord, max_coord = up_block[2:]
                    init_cov = 0
                    gap_located.add(e_num)
                else:
                    logger.info(
                        f'Exon {e_num} corresponds to a chain gap either too small or too large '
                        'to define its coordinates; marking it as missing'
                    )
                    missing.add(e_num)
            ## case 1.2: enclosing chain part is an aligned block
            else:
                min_coord = up_block[2]
                max_coord = up_block[3]
                raw_min_coord, raw_max_coord = up_block[2:]
                trim_up: int  = exon.start - up_block[0]
                trim_down: int = up_block[1] - exon.stop
                min_coord += trim_up if codirected else trim_down
                max_coord -= trim_down if codirected else trim_up
                if trim_up >= 0 and trim_down >= 0:
                    logger.info(
                        f'Exon {exon.num} is fully enclosed within '
                        f'block {up_block_name}'
                    ) if verbose else None
                    init_cov = exon.length()
                else:
                    logger.info(
                        f'Exon {exon.num} is partially covered by a single '
                        f'block {up_block_name}'
                    ) if verbose else None
                    init_cov = intersection(
                        *sorted(coords),
                        *sorted(up_block[:2])
                    )
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
            for _block in ex_num2blocks[e_num]:
                if '_' in _block:
                    continue
                inter_size: int =  intersection(*coords, *sorted(blocks[_block][:2]))
                init_cov += inter_size if inter_size > 0 else 0
            ## case 2.1: any of two marginal blocks are unaligned gaps gaps
            ## in this case, extend the exon coordinate to the next block
            if new_up_block_name:
                side: str = 'Left' if codirected else 'Right'
                logger.info(
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
                logger.info(
                    f'{side} flank for exon {exon.num} is a chain gap; '
                    f'extending coordinates to block {new_down_block_name}'
                )
                logger.info(','.join(map(str, blocks[new_down_block_name])))
                trim_down: int = next_down_block[1] - exon.stop
                dangling_down: bool = trim_down < -1 * (
                    donor_flank if codirected else acc_flank
                )
                extra_flank_down: int = extra_flank * dangling_down
                if codirected:
                    max_coord = max(0, next_down_block[3] - trim_down + extra_flank_down)
                else:
                    min_coord = next_down_block[2] + trim_down - extra_flank_down
            ## case 2.2: at least one of the marginal blocks is a true
            ## aligned block; trim the coordinates by the exon coordinates
            ## from the reference annotation
            logger.info(
                f'Exon {exon.num} intersects the following blocks: '
                f'{up_block_name}, {down_block_name}'
            ) if verbose else None
            trim_up: int = up_block[0] - exon.start
            trim_down: int = down_block[1] - exon.stop
            if trim_up > 2000 or trim_down < -2000:
                missing.add(e_num)
            else:
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
                    logger.info(
                        f'Trimming exon {exon.num} coordinates '
                        f'from the {side} side'
                    ) if verbose else None
                    min_coord = max(
                        up_block[2] - trim_up if codirected else
                        down_block[2] + trim_down, 0
                    )
                if not max_coord:
                    side: str = 'right' if codirected else 'left'
                    logger.info(
                        f'Trimming exon {exon.num} coordinates '
                        f'from the {side} side'
                    ) if verbose else None
                    max_coord = max(
                        down_block[3] - trim_down if codirected else
                        up_block[3] + trim_up, 0
                    )
        ## a short sanity check: if an exon ended up with reverted coordinates,
        ## record it as missing and proceed
        if min_coord >= max_coord:
            logger.warning(
                'Negative size locus for exon %i after coordinate trimming; marking it as missing' % e_num
            )
            missing.add(e_num)
            min_coord, max_coord = 0, 0
        # min_coord, max_coord = sorted([min_coord, max_coord])      
        ## exon loci longer than five times the exon length are
        ## suspected for chain distruption
        elif (max_coord - min_coord) / exon.length() >= 5 and e_num not in gap_located:
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
        min_coord = nn(min(min_coord, q_chrom_size))
        max_coord = nn(min(max_coord, q_chrom_size))
        logger.info(
            f'Resulting coordinates for exon {exon.num} are: '
            f'{min_coord}-{max_coord}'
        ) if verbose else None
        exons2coordinates[e_num] = (min_coord, max_coord)
        exon2init_cov[e_num] = init_cov
        # dangling_ends[e_num] = cov_status

    # spanning_chain: bool = all(e in missing or e in gap_located for e in exons)
    spanning_chain: bool = all(
        e in gap_located or (ex2del[e] and e in missing) for e in exons
    )
    if spanning_chain:
        if upstream_chain_gap:
            tstart = blocks[upstream_chain_gap][0]
            qstart = blocks[upstream_chain_gap][2]
        if downstream_chain_gap is not None:
            tstop = blocks[downstream_chain_gap][1]
            qstop = blocks[downstream_chain_gap][3]
        qstart, qstop = sorted((qstart, qstop))

    # print(f'{out_of_chain=}')
    # print(f'GAP LOCATED: {gap_located=}')
    if not spanning_chain and len(exons) > 1:
        for i, curr_exon in enumerate(sorted_exon_keys[1:], start=1):
            if curr_exon in out_of_chain or curr_exon in missing or curr_exon in gap_located:
                continue
            prev_exon: int = sorted_exon_keys[i-1]
            if prev_exon in out_of_chain or prev_exon in missing or prev_exon in gap_located:
                continue
            # print(f'{prev_exon=}, {curr_exon=}')
            ref_intron_start: int = exons[prev_exon].coords()[1]
            ref_intron_end: int = exons[curr_exon].coords()[0]
            ref_intron_len: int = ref_intron_end - ref_intron_start
            if codirected:
                query_intron_start: int = exons2coordinates[prev_exon][1]
                query_intron_end: int = exons2coordinates[curr_exon][0]
            else:
                query_intron_start: int = exons2coordinates[curr_exon][1]
                query_intron_end: int = exons2coordinates[prev_exon][0]
            query_intron_len: int = query_intron_end - query_intron_start
            if query_intron_len >= ref_intron_len * 5 and query_intron_len >= MAX_CHAIN_INTRON_LEN:
                logger.warning(
                    'Intron %i is too long (%i) and likely contains insertion' % (
                        min(prev_exon, curr_exon), query_intron_len
                    )
                )
                ref_strand: bool = curr_exon > prev_exon
                query_strand: bool = ref_strand == codirected
                prev_range: Iterable[int] = (
                    range(1, prev_exon + 1) if ref_strand else range(max(exons), prev_exon - 1, -1)
                )
                sum_before: int = sum(
                    exons[e].length() for e in prev_range if e not in out_of_chain and e not in missing
                )
                next_range: Iterable[int] = (
                    range(curr_exon, max(exons) + 1) if ref_strand else range(1, curr_exon + 1)
                )
                sum_after: int = sum(
                    exons[e].length() for e in next_range if e not in out_of_chain and e not in missing
                )
                if sum_before != sum_after:
                    loosing_range: Iterable[int] = prev_range if sum_before < sum_after else next_range
                    logger.warning(
                        'Marking exons [%s] as missing due to intronic insertion location' % ','.join(map(str, loosing_range))
                    )
                    for e in loosing_range:
                        if e in exons2coordinates:
                            del exons2coordinates[e]
                        missing.add(e)

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
        exon2init_cov,
        spanning_chain
        # aa_sat
    )

    return mapper#, curr_block


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
        'min_exon', 'max_exon', 'first_exon', 'last_exon', 
        'start', 'stop', 'init_start', 'init_stop', 'qstart', 'qstop'
    ]
    min_exon: int
    max_exon: int
    first_exon: int
    last_exon: int
    start: int
    stop: int
    init_start: int
    init_stop: int
    qstart: int
    qstop: int


class ProjectionGroup:
    """
    A projection-defining object comprising of a chain mapper and/or a LASTZ segment
    """
    __slots__ = [
        'mapper', 'annot', 'start', 'stop', 'init_start', 'init_stop',
        'qstart', 'qstop', 'chrom_len', 'first_exon', 'last_exon',
        'chrom', 'chain', 'acc_flank', 'donor_flank', 'out_of_chain_exons',
        'max_space_size', 'strand', 'shifted_exons',
        'evidence_source', 'exon_expected_loci', 'exon_search_spaces',
        'exon_coords', 'missing', 'gap_located', 'all_missing', 'discordant_cases',
        'cesar_exon_grouping', 'group_coords', 'assembly_gaps', 'spliceai_sites',
        'logger', 'v'
    ]
    def __init__(
        self,
        mapper: Exon2BlockMapper,
        annot: AnnotationEntry,
        chain: str,
        start: int,
        stop: int,
        init_start: int,
        init_stop: int,
        qstart: int,
        qstop: int,
        chrom_len: int,
        first_exon: int,
        last_exon: int,
        acc_flank: int,
        donor_flank: int,
        out_of_chain_exons: Set[int],
        max_space_size: int, 
        logger: logging.Logger
    ) -> None:
        """
        An auxiliary class for defining exon coordinates and CESAR alignment groups
        """
        self.logger: logging.Logger = logger
        self.mapper: Exon2BlockMapper = mapper
        self.annot: AnnotationEntry = annot
        self.chain: str = chain
        self.start: int = start
        self.stop: int = stop
        self.init_start: int = init_start
        self.init_stop: int = init_stop
        self.qstart: int = qstart
        self.qstop: int = qstop
        self.chrom_len: int = chrom_len
        self.first_exon: int = first_exon
        self.last_exon: int = last_exon
        self.acc_flank: int = acc_flank
        self.donor_flank: int = donor_flank
        self.max_space_size: int = max_space_size
        self.strand: bool = mapper.qstrand == self.annot.strand
        self.chrom: str = self.mapper.qchrom if self.mapper else self.segment.chrom
        self.shifted_exons: Dict[str, List[int]] = {
            x: [] for x in ['segment_up', 'segment_down', 'chain_up', 'chain_down']
        }

        self.out_of_chain_exons: Set[int] = out_of_chain_exons
        self.missing: Set[int] = {*mapper.missing}
        # self.gap_located: Set[int] = {*mapper.gap_located}
        # self.all_missing: Set[int] = self.missing.union(self.out_of_chain_exons)
        # self.missing: Set[int]
        self.all_missing: Set[int] = {*self.missing, *self.out_of_chain_exons}
        self.gap_located: Set[int] = {*mapper.gap_located}

        self.spliceai_sites: Dict[str, Dict[int, float]] = {'donor': {}, 'acceptor': {}}

        self.evidence_source: Dict[int, int] = {}
        self.exon_expected_loci: Dict[int, Tuple[int]] = {}
        self.exon_search_spaces: Dict[int, Tuple[int]] = {}
        self.exon_coords: Dict[int, Tuple[str, int]] = {}
        self.populate_coord_slots(mapper)
        
        self.cesar_exon_grouping: List[List[int]] = []
        self.group_coords: Dict[int, Tuple[int, int]] = {}
        self.assembly_gaps: List[Tuple[str, int]] = []
        # self.v: bool = verbose
        if mapper.spanning_chain:
            self.spanning_chain_grouping()
        else:
            defined_groups: List[List[int]] = self.group_defined_exons()
            self.final_grouping(defined_groups)
        # self.run()

    def _to_log(self, msg: str = '', lvl: str = 'info') -> None:
        getattr(self.logger, lvl)(msg)

    def _die(self, msg: str) -> None:
        self._to_log(msg, 'critical')
        sys.exit(1)

    def _consistent_coords(self, coords: Tuple[str, bool, int, int]) -> Tuple[str, bool, int, int]:
        """
        Modifies a coordinate tuple to ensure that no numeric value 
        exceeds chromosome length or falls below zero
        """
        start_coord: int = nn(min(coords[2], self.chrom_len)) if coords[2] is not None else None
        stop_coord: int = nn(min(coords[3], self.chrom_len)) if coords[3] is not None else None
        return (*coords[:2], start_coord, stop_coord)

    def _exon2group(self, exon: int, groups: List[List[int]]) -> List[int]:
        """Returns the group the exon belongs to"""
        for group in groups:
            if exon in group:
                return group
        return []

    def _coords_for_defined_group(self, group: List[int]) -> Tuple[int, int]:
        """
        Returns group search space coordinates
        """
        group_start: int = min(
            self.exon_search_spaces[x][2] for x in group
            if x not in self.all_missing
        )
        group_end: int = max(
            self.exon_search_spaces[x][3] for x in group
            if x not in self.all_missing
        )
        return group_start, group_end

    def populate_coord_slots(self, mapper: Exon2BlockMapper) -> None:
        """Transfers coordinates from the Mapper object to the respective owned slots"""
        chrom: str = self.chrom
        strand: bool = self.strand
        for exon in range(self.first_exon, self.last_exon + 1):
            if exon in self.all_missing:
                # self.exon_coords[exon]
                # self.exon_expected_loci[exon]
                # self.exon_search_spaces[exon]
                continue
            exon_start, exon_end = mapper.get_exon_coords(exon)[2:]
            coords: Tuple[str, bool, int, int] = self._consistent_coords(
                (chrom, strand, exon_start, exon_end)
            )
            self.exon_coords[exon] = coords
            self.exon_expected_loci[exon] = coords
            # if exon in self.gap_located:
            #     self.exon_search_spaces[exon] = coords
            #     continue
            if strand:
                # search_start: int = exon_start - (self.acc_flank if exon != self.first_exon else 0)
                # search_end: int = exon_end + (self.donor_flank if exon != self.last_exon else 0)
                search_start: int = exon_start - self.acc_flank
                search_end: int = exon_end + self.donor_flank
            else:
                # search_start: int = exon_start - (self.donor_flank if exon != self.last_exon else 0)
                # search_end: int = exon_end + (self.acc_flank if exon != self.first_exon else 0)
                search_start: int = exon_start - self.donor_flank
                search_end: int = exon_end + self.acc_flank
            search_coords: Tuple[str, bool, int, int] = self._consistent_coords(
                (chrom, strand, search_start, search_end)
            )
            self.exon_search_spaces[exon] = search_coords

    def spanning_chain_grouping(self) -> None:
        """
        If transcript is projected via a spanning chain, 
        gather all exons as a single CESAR group and exit
        """
        self.cesar_exon_grouping.append(
            list(range(self.first_exon, self.last_exon + 1))
        )
        group_start: int = min(
            self.exon_search_spaces[x][2] for x in self.exon_coords if self.exon_coords[x][2] != self.exon_coords[x][3]
        )
        group_end: int = max(
            self.exon_search_spaces[x][3] for x in self.exon_coords if self.exon_coords[x][2] != self.exon_coords[x][3]
        )
        self.group_coords[0] = (group_start, group_end)

    def group_defined_exons(self) -> List[List[int]]:
        """
        Groups exons with chain-defined loci in the query 
        according to their coordinate/search space intersection
        """
        init_exon_range: List[int] = [
            x for x in range(self.first_exon, self.last_exon + 1)
            if x not in self.all_missing
        ]
        init_exon_range.sort(key=lambda x: self.exon_search_spaces[x][2:])
        self._to_log(
            'Defining search spaces and grouping for well-defined exons: %s' % ','.join(
                map(str, init_exon_range)
            )
        )
        if not self.strand:
            init_exon_range = init_exon_range[::-1]
        groups: List[List[int]] = []
        pointer: int = 0
        present_exon_num: int = len(init_exon_range)
        while pointer <= present_exon_num - 1:
            curr_exon: int = init_exon_range[pointer]
            # print(f'Current exon is {curr_exon}, pointer is {pointer}')
            curr_group: List[int] = [pointer]
            curr_start, curr_end = self.exon_search_spaces[curr_exon][2:]
            for next_pointer in range(pointer+1, present_exon_num):
                next_exon: int = init_exon_range[next_pointer]
                next_start, next_end = self.exon_search_spaces[next_exon][2:]
                inter: int = intersection(
                    curr_start, curr_end, next_start, next_end
                )
                if inter >=  -MIN_INTRON_LENGTH:
                    self._to_log(
                        'Search space for exon %i intersects that of exon %i' % (curr_exon, next_exon)
                    )
                    curr_group.append(next_pointer)
                    # print(f'{curr_start=}, {next_start=}, {min(curr_start, next_start)=}')
                    curr_start = min(curr_start, next_start)
                    # print(f'{curr_end=}, {next_end=}, {max(curr_end, next_end)=}')
                    curr_end = max(curr_end, next_end)
            last_intersected: int = max(curr_group)
            group_to_add: List[int] = init_exon_range[pointer:last_intersected+1]
            self._to_log('Exons %s are added as a single group' % ','.join(map(str, group_to_add)))
            groups.append(group_to_add)
            pointer = last_intersected + 1
            # print(f'Updated pointer is {pointer}')
        ## phase two: for exons within each group, update the expected coordinates
        for group in groups:
            sorted_by_coord: List[int] = sorted(
                group, key=lambda x: self.exon_search_spaces[x][2:]
            )
            # print(f'{group=}, {sorted_by_coord=}')
            for i, exon in enumerate(sorted_by_coord):
                raw_curr_exon_start, raw_curr_exon_end = self.exon_expected_loci[exon][2:]
                flanked_curr_exon_start, flanked_curr_exon_end = self.exon_search_spaces[exon][2:]
                # print(f'{raw_curr_exon_start=}, {raw_curr_exon_end=}, {flanked_curr_exon_start=}, {flanked_curr_exon_end=}')
                if i < len(group) - 1:
                    next_exon: int = sorted_by_coord[i+1]
                    raw_next_exon_start, raw_next_exon_end = self.exon_expected_loci[next_exon][2:]
                    flanked_next_exon_start, flanked_next_exon_end = self.exon_search_spaces[next_exon][2:]
                    # print(f'{raw_next_exon_start=}, {raw_next_exon_end=}, {flanked_next_exon_start=}, {flanked_next_exon_end=}')
                    upd_curr_flanked_end: int = max(
                        raw_curr_exon_end,
                        min(flanked_curr_exon_end, flanked_next_exon_start - 1)
                    )
                    upd_next_flanked_start: int = min(
                        raw_next_exon_start,
                        max(flanked_next_exon_start, flanked_curr_exon_end + 1)
                    )
                    # if exon not in self.gap_located:
                    self.exon_search_spaces[exon] = (
                        *self.exon_search_spaces[exon][:2],
                        flanked_curr_exon_start,
                        upd_curr_flanked_end
                    )
                    # if next_exon not in self.gap_located:
                    self.exon_search_spaces[next_exon] = (
                        *self.exon_search_spaces[next_exon][:2],
                        upd_next_flanked_start,
                        flanked_next_exon_end
                    )
        return groups

    def final_grouping(self, defined_groups: List[List[int]]) -> None:
        """Resolves grouping for undefined exons and defines group coordinates"""
        final_groups: List[List[int]] = [
            [x for x in y] for y in defined_groups
        ]
        # print(f'{final_groups=}, {self.exon_search_spaces=}')
        all_missing_exons: List[int] = sorted(
            x for x in range(self.first_exon, self.last_exon + 1) if x in self.all_missing
        )
        if not self.strand:
            all_missing_exons = all_missing_exons[::-1]
        curr_group: List[int] = []
        for exon in all_missing_exons:
            for final_group in final_groups:
                if exon_belongs_to_group(exon, final_group):
                    final_group.append(exon)
                    break
            else:
                prev_exon: int = exon - (1 if self.strand else -1)
                if curr_group and prev_exon not in curr_group:
                    final_groups.append(curr_group)
                    curr_group = []
                curr_group.append(exon)
        if curr_group:
            final_groups.append(curr_group)
        ## step 2: define coordinates
        final_groups = sorted(final_groups, key=lambda x: min(x) if self.strand else -min(x))
        for i, group in enumerate(final_groups):
            group = sorted(group, reverse=not self.strand)
            if any(x not in self.all_missing for x in group):
                # first_exon: int = (min if self.strand else max)(group)
                # group_start: int = self.exon_search_spaces[first_exon][2]
                # last_exon: int = (max if self.strand else min)(group)
                # group_end: int = self.exon_search_spaces[last_exon][3]
                # self.cesar_exon_grouping.append(group)
                # self.group_coords[i] = (group_start, group_end)
                group_start, group_end = self._coords_for_defined_group(group)
            else:
                group_start, group_end = None, None
                if self.first_exon in group:
                    if self.strand:
                        group_start = nn(min(self.start, self.chrom_len))
                    else:
                        group_end = nn(min(self.stop, self.chrom_len))
                if self.last_exon in group:
                    if self.strand:
                        group_end = nn(min(self.stop, self.chrom_len))
                    else:
                        group_start = nn(min(self.start, self.chrom_len))
                prevs: List[int] = [
                        x for x in range(self.first_exon, min(group)) 
                        if x not in self.all_missing
                    ]
                prev: Union[int, None] = None if not prevs else prevs[-1]
                nexts: List[int] = [
                    x for x in range(max(group) + 1, self.last_exon + 1)
                    if x not in self.all_missing
                ]
                next_: Union[int, None] = None if not nexts else nexts[0]
                # print(f'{prevs=}, {prev=}, {nexts=}, {next_=}, {self.strand=}, {self.first_exon=}, {self.last_exon=}')
                if group_start is None:
                    neighbor: Union[int, None] = prev if self.strand else next_
                    if neighbor is None:
                        self._die(
                            'No upstream neighbor with defined coordinates'
                            f' for group {group} was found'
                        )
                    neighbor_group: List[int] = self._exon2group(neighbor, final_groups)
                    group_start = self._coords_for_defined_group(neighbor_group)[1] + 1 #self.exon_search_spaces[neighbor][3] + 1
                if group_end is None:
                    neighbor: Union[int, None] = next_ if self.strand else prev
                    if neighbor is None:
                        self._die(
                            'No downstream neighbor with defined coordinates'
                            f' for group {group} was found'
                        )
                    neighbor_group: List[int] = self._exon2group(neighbor, final_groups)
                    group_end = self._coords_for_defined_group(neighbor_group)[0] - 1 #self.exon_search_spaces[neighbor][2] - 1
                ## now, try shrinking this coordinates within the better-defined blocks
                self._to_log(
                        f'Shrinking search space for exon group {i}'
                    )
                # print(f'{i=}, {group=}, {self.exon_coords[group[0]]=}')
                # print(f'{group_start=}, {group_stop=}')
                suitable_blocks: List[Tuple[int]] = self.mapper.get_suitable_blocks(
                    group_start, group_end
                )
                group_len: int = sum(
                    [self.annot.exons[x].stop - self.annot.exons[x].start for x in group]
                ) + MIN_INTRON_LENGTH * (len(group) - 1)
                max_inter_coords: Tuple[Union[int, None]] = define_max_space(
                    suitable_blocks, group_len, self.max_space_size
                )
                # print(f'{group_start=}, {group_stop=}, {max_inter_coords=}, {group_len=}, {self.max_space_size=}')
                if max_inter_coords[0] is not None:
                    if (max_inter_coords[1] - max_inter_coords[0]) < (group_end - group_start):
                        self._to_log(
                            f'Previous coords for group {i}: '
                            f'{self.chrom}:{group_start}-{group_end} '
                            f'(size = {group_end-group_start}), '
                            f'shrinked coords: '
                            f'{self.chrom}:{max_inter_coords[0]}-{max_inter_coords[1]} '
                            f'(size = {max_inter_coords[1]-max_inter_coords[0]})'
                        )
                        group_start = max_inter_coords[0]
                        group_end = max_inter_coords[1]
                # self.group_coords[i] = (group_start, group_end)
            self.cesar_exon_grouping.append(group)
            self.group_coords[i] = (group_start, group_end)
                


def exon_belongs_to_group(exon: int, group: List[int]) -> bool:
    return min(group) < exon < max(group)


def cesar_memory_check(ref_lengths: List[int], query_length: int) -> float:
    """
    Estimates memory consumption for CESAR2.0 . Arguments are:
    :block_lengths: is a list of reference (exon) sequence sizes
    :query_length: is a length of the query sequence
    """
    num_states, rlength, extra = 0, 0, 100000
    for block_length in ref_lengths:
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
    transcript: str = ''
    exon: int = 0
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
    s_sites: List[str] = []
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
    ## TODO: All this functionality should be rustified one day
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
            sys.exit(
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
            sys.exit(
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