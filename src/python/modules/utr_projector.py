#!/usr/bin/env python3

"""
Projects reference untranslated regions and adds the annotated blocks 
to the TOGA-produced BED file
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from collections import defaultdict
from dataclasses import dataclass
from extra.bed12ToFraction import make_track
from preprocessing import get_chain
from shared import CommandLineManager, CONTEXT_SETTINGS, intersection
from string_splitter import transcriptwise_subchains ## TODO: Rename the module
from typing import Any, Dict, List, Optional, TextIO, Tuple, Union

import click
import gzip

__author__ = 'Yury Malovichko'
__credits__ = ('Luca Hoppach')
__year__ = '2025'

PLUS: str = '+'

@dataclass
class UtrExon:
    __slots__ = (
        'chrom', 'start', 'end', 'name', 
        'strand', 'upstream', 'utr_type', 
        'adjacent_to_cds'
    )
    chrom: str
    start: Union[int, None]
    end: Union[int, None]
    name: str
    strand: bool
    upstream: bool
    utr_type: str
    adjacent_to_cds: bool

    def coords(self) -> Tuple[int, int]:
        """Returns the (start, end) tuple"""
        return (self.start, self.end)


def bed_data2line(bed_data: List[List[str]]) -> str:
    """
    Formats data record for a BED12 entry back into a BED file string
    """
    out_line: str = ''
    ## iterate over all records
    last_record: int = len(bed_data) - 1
    for i, bed_record in enumerate(bed_data):
        line: str = '\t'.join(bed_record)
        out_line += line
        if i < last_record:
            out_line.append('\n')
    return out_line


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'ref_bed',
    type=click.File('r', lazy=True),
    metavar='REF_BED'
)
@click.argument(
    'query_bed',
    type=click.File('r', lazy=True),
    metavar='QUERY_BED'
)
## TO DO: Add exon meta file support as a mandatory argument
@click.argument(
    'exon_meta',
    type=click.Path(exists=True),
    metavar='EXON_META_FILE'
)
@click.argument(
    'chain_file',
    type=click.Path(exists=True),
    metavar='CHAIN_FILE'
)
@click.option(
    '--relative_threshold',
    '-d',
    type=click.FloatRange(min=0, max=None),
    metavar='FLOAT',
    default=2.0,
    show_default=True,
    help=(
        'Relative reference exon length threshold by which UTR exon\'s projection can be extended; '
        'unaligned exon ends protruding further than this value '
        'are cropped to the last aligned base'
    )
)
@click.option(
    '--absolute_threshold',
    '-s',
    type=click.IntRange(min=0, max=None),
    metavar='INT',
    default=3000,
    show_default=True,
    help=(
        'Absolute length threshold by which UTR exon\'s projection can be extended; '
        'unaligned exon ends protruding further than this value '
        'are cropped to the last aligned base'
    )
)
@click.option(
    '--output',
    '-o',
    type=click.File('w', lazy=True),
    default=sys.stdout,
    show_default=False,
    help='A path to write the results to [default: stdout]'
)
@click.option(
    '--log_file',
    '-l',
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help='A file to log the progress to'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)
@click.option(
    '--debug',
    '-d',
    is_flag=True,
    default=False,
    help='Logs technical data for debugging purposes'
)

class UtrProjector(CommandLineManager):
    __slots__ = (
        'chain_file', 'ref_utrs', 'tr2exon_num', 'query_bed', 'chain2ref_trs',
        'proj2exon_status', 'first_exon_coords', 'last_exon_coords', 
        'abs_threshold', 'rel_threshold', 'output', 'd'
    )

    def __init__(
        self,
        ref_bed: click.File,
        query_bed: click.File,
        exon_meta: click.Path,
        chain_file: click.Path,
        relative_threshold: Optional[float],
        absolute_threshold: Optional[int],
        output: Optional[click.Path],
        log_file: Optional[click.Path],
        verbose: Optional[bool],
        debug: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.d: bool = debug
        self.set_logging()

        self.chain_file: click.Path = chain_file

        self._to_log('Parsing reference annotation file')
        self.ref_utrs: Dict[str, List[Any]] = {}
        self.tr2exon_num: Dict[str, int] = {}
        self.read_ref_bed(ref_bed)

        self._to_log('Parsing query annotation file')
        self.query_bed: Dict[str, List[List[Any]]] = defaultdict(list)
        self.chain2ref_trs: Dict[str, List[str]] = defaultdict(list)
        self.read_query_bed(query_bed)

        self._to_log('Parsing exon meta file')
        self.proj2exon_status: Dict[str, Tuple[bool, bool]] = {}
        self.first_exon_coords: Dict[str, Tuple[str, int, int]] = {}
        self.last_exon_coords: Dict[str, Tuple[str, int, int]] = {}
        self.read_exon_meta(exon_meta)

        self.abs_threshold: int = absolute_threshold
        self.rel_threshold: int = relative_threshold
        self.output: click.File = output

        self.run()

    def run(self) -> None:
        """
        The algorithm so far:
        1) For each chain, sort the UTR exons from the corresponding transcripts;
        2) Run transcriptwise_chains(), get the aligned blocks for each UTR exon
        3) Get aligned coordinates for start and end (or the closest to them) for each block;
        * For 5'-adjacent UTRs, get only the start coordinate;
        * For 3'-adjacent UTR exons, get only the end coordinate
        4) Once all the UTR projections are obtained, add the 
        """
        for chain, trs in self.chain2ref_trs.items():
            self._to_log('Processing chain %s' % chain)
            ## create a transcriptwise storage for projected exons
            tr2exons: Dict[str, List[UtrExon]] = defaultdict(list)
            # print(f'{self.ref_utrs=}')
            # print(f'{chain=}, {trs=}')
            chain_str: str = get_chain(self.chain_file, chain)
            utr_exons: List[Tuple[int, int, str, str]] = [
                x for tr in trs for x in self.ref_utrs.get(tr, [])
            ]
            if not utr_exons:
                self._to_log(
                    'No UTR exons projected via chain %s; skipping to the next chain' % chain,
                    'warning'
                )
                continue
            # print(f'{utr_exons=}')
            utr_exons.sort(key=lambda x: (x.start, x.end))
            utr_starts: List[int] = []
            utr_ends: List[int] = []
            utr_names: List[str] = []
            for exon in utr_exons:
                utr_starts.append(exon.start)
                utr_ends.append(exon.end)
                utr_names.append(exon.name)
            # print(f'{utr_names=}, {utr_starts=}, {utr_ends=}')
            subchains, chain_meta = transcriptwise_subchains(
                chain_str, utr_names, utr_starts, utr_ends
            )
            # print(f'{chain_meta=}')
            # t_chrom: str = chain_meta[0]
            q_chrom: str = chain_meta[1]
            # t_chain_start: int = chain_meta[2]
            # t_chain_stop: int = chain_meta[3]
            t_strand: bool = chain_meta[4]
            q_strand: bool = chain_meta[7]
            codirected: bool = t_strand == q_strand
            # print(f'{t_strand=}, {q_strand=}, {codirected=}')
            for i, exon in enumerate(utr_exons):
                # subchain: List[Tuple[str, Tuple[int]]] = list(subchains[exon.name].items())
                subchain: List[Tuple[str, Tuple[int]]] = subchains[exon.name]
                if not subchain:
                    continue
                else:
                    # if len(subchain) > 0:
                    #     subchain = [subchain[0], subchain[-1]]
                    ## now, get the start and end coordinates
                    ## TODO: Leave only the first and the last block for iterating
                    query_exon: Union[UtrExon, None] = self.project_coords(
                        exon, subchain, q_chrom, q_strand, codirected
                    )
                    ## exon is fully enclosed within a chain gap; skip it 
                    if query_exon is None:
                        continue
                tr: str = exon.name.split('*')[0]
                # print(f'{tr=}, {exon.name=}, {exon.start=}, {exon.end=}, {query_exon=}')
                tr2exons[tr].append(query_exon)

            ## add UTR exons to the output BED file
            for tr, utr_exons in tr2exons.items():
                proj: str = f'{tr}#{chain}'
                proj_data: List[List[str]] = self.query_bed[proj]
                if not utr_exons:
                    self._to_log(
                        'No UTR exons found for projection %s' % proj, 'warning'
                    )
                    bed_line: str = bed_data2line(proj_data)
                    self.output.write(bed_line + '\n')
                    continue
                modified_bed: List[Any] = self.attach_utrs_to_bed(
                    proj_data, utr_exons
                )
                bed_line: str = bed_data2line(modified_bed)
                self.output.write(bed_line + '\n')

    def read_ref_bed(self, file: TextIO) -> None:
        """
        Reads reference BED file, extracting UTR exons from each record 
        """
        for i, line in enumerate(file, start=1):
            line = line.rstrip()
            data: List[str] = line.split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 12:
                self._die('Reference bed file contains less than 12 fields at line {i}')
            tr: str = data[3]
            strand: bool = data[5] == '+'
            thick_start: int = int(data[6])
            thick_end: int = int(data[7])
            utr_exons: List[List[Any]] = make_track(
                line, mode='utr', bed6=True, num2score=True
            )
            exon_num: int = int(data[9])
            exons: List[UtrExon] = []
            # if tr == 'ENST00000377149.5#OR11A1':
            #     print(f'ENST00000377149.5#OR11A1: {exon_num=}')
            for ex in utr_exons:
                # if tr == 'ENST00000377149.5#OR11A1':
                #     print(f'{ex=}')
                if ex[2] <= thick_start:
                    side: str = '5' if strand else '3'
                    adjacent_to_cds: bool = ex[2] == thick_start
                    exon_num -= int(not adjacent_to_cds)
                    upstream: bool = True
                elif ex[1] >= thick_end:
                    side: str = '3' if strand else '5'
                    adjacent_to_cds: bool = ex[1] == thick_end
                    exon_num -= int(not adjacent_to_cds)
                    upstream: bool = False
                name: str = f'{ex[3]}*{ex[4]}'
                utr_exon: UtrExon = UtrExon(
                    ex[0], ex[1], ex[2], name, 
                    strand, upstream, side, adjacent_to_cds
                )
                exons.append(utr_exon)
            # if tr == 'ENST00000377149.5#OR11A1':
            #     print(f'ENST00000377149.5#OR11A1: {exon_num=}')
            self.ref_utrs[tr] = exons
            self.tr2exon_num[tr] = exon_num

    def read_query_bed(self, file: TextIO) -> None:
        """
        Reads query BED file, saving record data and chain numbers
        """
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 12:
                self._die('Reference bed file contains less than 12 fields at line {i}')
            proj: str = data[3]
            proj_components: str = proj.split('#')
            tr: str = '#'.join(proj_components[:-1])
            chains: List[str] = proj_components[-1].split(',')
            for chain in chains:
                if not chain.isdigit():
                    self._die(
                        f'Invalid or missing chain ID at line {i}: {proj}'
                    )
                self.chain2ref_trs[chain].append(tr)
            self.query_bed[proj].append(data)

    def read_exon_meta(self, file: str) -> None:
        """
        Infers the following data from the exon metadata file:
        1) Terminal (first and last) exon presence status;
        2) Terminal exon coordinates
        """
        # proj2max_exon: Dict[str, int] = defaultdict(int)
        proj2first_status: Dict[str, str] = {}
        proj2last_status: Dict[str, str] = {}
        with (gzip.open if file.endswith('.gz') or file.endswith('.gzip') else open)(file, 'r') as h:
            for i, line in enumerate(h):
                data: List[str] = line.rstrip().decode('utf8').split('\t')
                if not data or not data[0]:
                    continue
                if data[0] == 'projection':
                    continue
                proj: str = data[0]
                tr: str = '#'.join(proj.split('#')[:-1])
                if tr not in self.tr2exon_num:
                    continue
                exon: int = int(data[1])
                chrom: str = data[3]
                start: int = int(data[4])
                end: int = int(data[5])
                status: str = data[7]
                if exon == 1:
                    proj2first_status[proj] = status
                    self.first_exon_coords[proj] = (chrom, start, end)
                if exon == self.tr2exon_num[tr]:
                    # proj2max_exon[proj] = exon
                    proj2last_status[proj] = status
                    self.last_exon_coords[proj] = (chrom, start, end)
        self.proj2exon_status = {
            x: (proj2first_status[x], proj2last_status[x]) for x in proj2first_status
        }


    def project_coords(
        self,
        exon: UtrExon, 
        blocks: Dict[str, List[int]],
        q_chrom: str,
        q_strand: bool,
        codirected: bool
    ) -> Union[UtrExon, None]:
        """
        For each coordinate in tuple, gets its counterpart projected via an aligned block.

        Argument `blocks` contains chain block as coordinate tuples in the following format:
        [ref_start, ref_end, query_start, query_end],
        where end is guaranteed to be larger than start (sic)


        """
        if not blocks:
            return None
        print(f'{blocks=}') if self.d else None
        query_coords: List[Union[int, None]] = [None, None]
        ## infer the relative threshold
        ref_size: int = exon.end - exon.start
        rel_threshold: int = int(ref_size * self.rel_threshold)
        ## infer chain boundaries in the reference
        first_block: str = list(blocks.keys())[0]
        last_block: str = list(blocks.keys())[-1]
        chain_start: int = blocks[first_block][0]
        chain_end: int = blocks[last_block][1]
        # chain_start: int = blocks[0][1][0]
        # chain_end: int = blocks[-1][1][1]
        ## handle the marginal cases with exon ends uncovered by the chain
        ## marginal chain blocks in this case are guaranteed to be aligned sequences
        if exon.start < chain_start:
            ## exon is guaranteed to locate upstream to the coding sequence
            ## extend it over the chain boundary is possible
            j: int = 0 if codirected else 1
            offset: int = chain_start - exon.start
            if offset > self.abs_threshold and offset > rel_threshold:
                ## fall back to the first aligned base
                # query_coord: int = blocks[0][1][2]
                query_coord: int = blocks[first_block][2]
            else:
                # query_coord: int = max(0, blocks[first_block][2] + offset * (-1 if codirected else 1))
                if codirected:
                    # query_coord: int = max(0, blocks[0][1][2] - offset)
                    query_coord: int = max(0, blocks[first_block][2] - offset)
                else:
                    # query_coord: int = max(0, blocks[0][1][3] + offset)
                    query_coord: int = max(0, blocks[first_block][3] + offset)
            query_coords[j] = query_coord
        if exon.end > chain_end:
            ## exon is guaranteed to locate downstream to the coding sequence
            ## again, extend it over the chain boundary is possible
            j: int = 1 if codirected else 0
            offset: int = exon.end - chain_end
            if offset > self.abs_threshold and offset > rel_threshold:
                ## fall back to the last aligned base
                # query_coord: int = blocks[-1][1][3]
                query_coord: int = blocks[last_block][3]
            else:
                # query_coord: int = max(0, blocks[last_block][3] + offset * (1 if codirected else -1))
                if codirected:
                    # query_coord: int = max(0, blocks[-1][1][3] + offset)
                    query_coord: int = max(0, blocks[last_block][3] + offset)
                else:
                    # query_coord: int = max(0, blocks[-1][1][2] - offset)
                    query_coord: int = max(0, blocks[last_block][2] - offset)
            query_coords[j] = query_coord
        ## a single UTR exon is highly unlikely to extend over both ends of the chain
        ## raise error if somehow both coordinates have been already define for the exon
        if all(x is not None for x in query_coords):
            self._die(
                'Single UTR exon %s spans over the whole chain' % exon.name
            )
        # print(f'{exon.strand=}, {q_strand=}, {codirected=}')
        ## set a boolean flag indicating whether the block is fully enclosed in a single gap
        enclosed_in_gap: bool = False
        ## iterate over all blocks
        for b_name, block in blocks.items():
            b_ref_start, b_ref_end, b_que_start, b_que_end = block
            ## set the counter for coordinates falling into this block
            ends_in_gap: int = 0
            ## check both coordinates
            for i, ref_coord in enumerate([exon.start, exon.end]):
                ## flip the coordinate if an exon is projected to the opposite strand
                j: int = i if codirected else int(not i)
                if b_ref_start <= ref_coord <= b_ref_end:
                    ## do not infer CDS-adjacent coordinate
                    if exon.adjacent_to_cds:
                        # if codirected:
                        adjacent_in_query: bool = exon.upstream and i == 1 or not exon.upstream and i == 0
                        # else:
                        #     adjacent_in_query: bool = exon.upstream and i == 0 or not exon.upstream and i == 1
                        # if exon.upstream and i == 1 or not exon.upstream and i == 0:
                        if adjacent_in_query:
                            print(f'Exon {exon.name} is adjacent to CDS from the {i} side') if self.d else None
                            ## all's fine, continue
                            continue
                    ## infer coordinate for an unaligned terminus;
                    ## those require special handling
                    if isinstance(b_name, str) and '_' in b_name:
                        print(f'Exon {exon.name} has its coordinate {i} in a gap') if self.d else None
                        print(f'Ref gap: {exon.chrom}:{b_ref_start}-{b_ref_end} | Query gap: {q_chrom}:{b_que_start}-{b_que_end}') if self.d else None
                        ## find the last aligned base; for adjacent exons, make sure 
                        #№ the UTR projection will not overlap with the coding sequence 
                        if exon.upstream:
                            last_aligned_base: int = (
                                b_ref_end if not exon.adjacent_to_cds else  min(b_ref_end, exon.end)
                            )
                        else:
                            last_aligned_base: int = (
                                b_ref_start if not exon.adjacent_to_cds else max(b_ref_start, exon.start)
                            )
                        last_aligned_base: int = b_ref_end if exon.upstream else b_ref_start
                        offset: int = abs(last_aligned_base - (exon.start if i == 0 else exon.end))
                        # print(f'{offset=}')
                        if  offset > self.abs_threshold and offset > rel_threshold:
                            ## offset is too long!
                            ## crop to the closest aligned base
                            query_coord: int = b_que_end if i == 1 else b_que_start
                        else:
                            ## extend the UTR into the wild
                            if i == 1:
                                query_coord: int = b_que_start + offset * (1 if codirected else -1)
                            else:
                                query_coord: int = b_que_end + offset * (-1 if codirected else 1)
                            # query_coord: int = last_aligned_base + offset * (-1 if exon.upstream else -1)
                        ends_in_gap += 1
                    ## terminus is aligned; the procedure is pretty easy:
                    ## move the query block coordinate by an offset in the reference
                    else:
                        print(f'Exon {exon.name} has its coordinate {i} aligned') if self.d else None
                        offset: int = ref_coord - b_ref_start
                        if codirected:
                            query_coord: int = b_que_start + offset
                        else:
                            query_coord: int = b_que_end - offset
                    query_coords[j] = query_coord
            ## if an exon is fully enclosed within a single chain gap,
            ## break the iteration
            if ends_in_gap == 2:
                enclosed_in_gap = True
                break

        ## an exon fully enclosed within a chain gap will not be projected 
        if enclosed_in_gap:
            return None

        ## estimate the strand in query
        query_strand: bool = exon.strand if codirected else not exon.strand
        ## estimate whether the exon is expected to be upstream to CDS in query
        upstream_in_query: bool = (exon.upstream == codirected)

        print(f'Ref coords: {exon.chrom}:{exon.start}-{exon.end} | Query coords: {q_chrom}:{query_coords[0]}-{query_coords[1]}') if self.d else None

        ## sanity checks
        ## check that at least one coordinate was properly projected for the exon
        if all(x is None for x in query_coords):
            self._die(
                'Failed to project either of the two coordinates for exon %s' % exon.name
            )
        ## then, check that the CDS-non-adjacent coordinates were inferred properly
        if query_coords[0] is None:
            ## query start is missing; must be downstream-adjacent in the query, otherwise a clear error
            if not (exon.adjacent_to_cds and not upstream_in_query):
                self._die(
                    'Failed to infer start coordinate for exon %s' % exon.name
                )
        if query_coords[1] is None:
            ## query end is missing; must be upstream-adjacent in the query, otherwise a clear error
            if not (exon.adjacent_to_cds and upstream_in_query):
                self._die(
                    'Failed to infer end coordinate for exon %s' % exon.name
                )
        query_exon: UtrExon = UtrExon(
            q_chrom, query_coords[0], query_coords[1], exon.name, 
            query_strand, upstream_in_query, exon.utr_type, 
            exon.adjacent_to_cds
        )
        return query_exon

    def attach_utrs_to_bed(
        self, projection: List[List[str]], utr_exons: List[UtrExon]
    ) -> List[str]:
        """
        Attach projected exons to the reference BED12 tracks
        """
        output: List[List[str]] = []
        for segment in projection:
            proj: str = segment[3]
            ## assess whether both terminal exons are present in the projection;
            ## if no, no UTRs should be added to the projection
            first_exon_status, last_exon_status = self.proj2exon_status[proj]
            if not first_exon_status and not last_exon_status:
                self._to_log(
                    'Both terminal exons are missing in the query for projection %s; skipping' % proj,
                    'warning'
                )
                output.append(segment)
                continue
            ## define whether this segment harbors any of the terminal exons
            chrom: str = segment[0]
            strand: bool = segment[5] == '+'
            thin_start: int = int(segment[1])
            thin_end: int = int(segment[2])
            first_chrom, first_start, first_end = self.first_exon_coords[proj]
            last_chrom, last_start, last_end = self.last_exon_coords[proj]
            first_is_present: bool = (
                first_exon_status and first_chrom == chrom and intersection(
                    thin_start, thin_end, first_start, first_end
                ) > 0
            )
            last_is_present: bool = (
                last_exon_status and last_chrom == chrom and intersection(
                    thin_start, thin_end, last_start, last_end
                ) > 0
            )
            ## no terminal exons present in this segment -> proceed further
            if not first_is_present and last_is_present:
                output.append(segment)
                continue
            ## showtime!
            thick_start: int = int(segment[6])
            thick_end: int = int(segment[7])
            utrs: List[UtrExon] = sorted(
                [x for x in utr_exons if x.chrom == chrom],
                key=lambda x: x.start if x.start is not None else thick_start
            )
            exon_sizes: List[int] = [
                int(x) for x in segment[10].split(',') if x
            ]
            exon_starts: List[int] = [
                int(x) for x in segment[11].split(',') if x
            ]
            ## check which 
            add_up: bool = strand and first_is_present or not strand and last_is_present
            if not add_up:
                self._to_log(
                    '%i\'-UTR sequences will not be added to projection %s' % ((5 if strand else 3), proj),
                    'warning'
                )
            add_down: bool = strand and last_is_present or not strand and first_is_present
            if not add_down:
                self._to_log(
                    '%i\'-UTR sequences will not be added to projection %s' % ((3 if strand else 5), proj),
                    'warning'
                )
            ## keep track on how many exons have been added from the upstream side
            num_up_added: int = 0
            ## start grafting the exon blocks
            for exon in utrs:
                ## check if this exon will be ever added
                if exon.upstream and not add_up or not exon.upstream and not add_down:
                    ## too bad! exons from this side will not be added
                    continue
                if exon.start is None:
                    overlaps_cds: bool = exon.end <= thick_end
                elif exon.end is None:
                    overlaps_cds: bool = exon.start >= thick_start
                else:
                    overlaps_cds: bool = intersection(exon.start, exon.end, thick_start, thick_end) > 0
                if overlaps_cds:
                    self._to_log(
                        'Projection for UTR exon %s overlaps the annotated coding sequence in %s, skipping' % (
                            exon.name, proj
                        ),
                        'warning'
                    )
                    continue
                ## update the 'thin' coordinates if needed
                if exon.start is not None and exon.start < thin_start:
                    ## starting point updated! update all the exons
                    ## likely will happen only once provided that the UTR exons are properly sorted
                    for i in range(len(exon_starts)):
                        exon_starts[i] = exon_starts[i] + thin_start - exon.start
                    thin_start = exon.start
                if exon.end is not None and exon.end > thin_end:
                    ## update the thin end value, otherwise no fields are affected
                    thin_end = exon.end
                ## now, define the exon start and size
                if exon.end is None:
                    ## means that the exon is upstream-adjacent to CDS
                    ## graft it to the first coding exon 
                    if not exon.adjacent_to_cds:
                        self._die(
                            'Attempting to merge an exon not adjacent to CDS in the reference'
                        )
                    ## remove the first coding exon
                    first_coding_size: int = exon_sizes.pop(num_up_added)
                    _ = exon_starts.pop(num_up_added)
                    ## and calculate the block stats
                    exon_size: int = first_coding_size + thick_start - exon.start
                    exon_start: int = exon.start - thin_start
                elif exon.start is None:
                    ## means that the exon is downstream-adjacent to CDS
                    ## graft it to the last coding exon
                    if not exon.adjacent_to_cds:
                        self._die(
                            'Attempting to merge an exon not adjacent to CDS in the reference'
                        )
                    ## the rest is essentialy the same
                    last_coding_size: int = exon_sizes.pop(-1)
                    exon_start: int = exon_starts.pop(-1)
                    exon_size: int = last_coding_size + exon.end - thick_end
                else:
                    exon_start: int = exon.start - thin_start
                    exon_size: int = exon.end - exon.start
                ## push the size and start  values
                if exon.upstream:
                    exon_sizes.insert(num_up_added, exon_size)
                    exon_starts.insert(num_up_added, exon_start)
                    num_up_added += 1
                else:
                    exon_sizes.append(exon_size)
                    exon_starts.append(exon_start)
            ## collapse the the data fields back into a Bed file line
            if len(exon_starts) != len(exon_sizes):
                print(f'{exon_starts=}')
                print(f'{exon_sizes=}')
                self._die(
                    (
                        'Number of exon start points not equal to number of exon size values '
                        'after UTR grafting for transcript %s'
                    ) % segment[3]
                )
            exon_sizes_line: str = ','.join(map(str, exon_sizes)) + ','
            exon_starts_line: str = ','.join(map(str, exon_starts)) + ','
            exon_num: int = len(exon_starts)
            upd_segment_data: List[str] = [
                chrom, str(thin_start), str(thin_end),
                *segment[3:9], str(exon_num), exon_sizes_line, exon_starts_line
            ]
            output.append(upd_segment_data)
        return output



if __name__ == '__main__':
    UtrProjector()