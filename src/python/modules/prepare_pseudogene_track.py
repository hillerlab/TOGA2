#!/usr/bin/env python3

"""
Prepares a BED12 track for processed pseudogenes found by TOGA2
"""

from cesar_wrapper_constants import PINK
from collections import defaultdict
from shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Any, Dict, List, Optional, Set, TextIO, Tuple

import click
import os
import sys

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__version__ = '2.0'
__email__ = 'yury.malovichko@senckenberg.de'
__credits__ = ('Bogdan Kirilenko', 'Michael Hiller')

TRANSCRIPT: str = 'TRANSCRIPT'

def dfs(
        graph: Dict[Any, List[Any]], 
        node: Any, 
        visited: List[Any],
        component: List[Any]
    ) -> List[Any]:
        """Performs a depth-first search over an adjacency dictionary of nodes"""
        if node not in graph:
            return component
        visited.append(node)
        component.append(node)
        for adj_node in graph[node]:
            if adj_node in visited:
                continue
            component = dfs(graph, adj_node, visited, component)
        return component


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'projection_report',
    type=click.File('r', lazy=True),
    metavar='PROJECTION_FILE'
)
@click.argument(
    'chain_file',
    type=click.Path(exists=True),
    metavar='CHAIN_FILE'
)
@click.option(
    '--output',
    '-o',
    type=click.File('w', lazy=True),
    metavar='FILE',
    default=sys.stdout,
    show_default=False,
    help='A path to save the results to [default: stdout]'
)
@click.option(
    '--log_file',
    '-l',
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help='A file to log the code progress to'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'Controls execution verbosity; '
        'if set, standard output will be used as an additional channel for logging'
    )
)


class PseudogeneTrackBuilder(CommandLineManager):
    __slots__ = ('output', 'chains', 'tr2chain', 'chain2byte', 'chain2coords')

    def __init__(
        self, 
        projection_report: click.File,
        chain_file: click.Path,
        output: Optional[click.File],
        # chain_index_file: Optional[click.File],
        log_file: Optional[click.Path],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging()
        self._to_log('Initializing pseudogene track builder')
        self.output: TextIO = output
        self.chains: Set[str] = set()
        self.tr2chain: List[str, List[str]] = defaultdict(list)
        self.chain2byte: Dict[str, int] = {}
        self.chain2coords: Dict[str, Tuple[str, int, int, str]] = {}
        self.get_pgenes_from_projection_report(projection_report)
        if not self.tr2chain.keys():
            self._to_log('No processed pseudogenes found', 'warning')
            return
        ## if chain index file was provided, parse the chain positions to data
        chain_index_file: str = chain_file.replace('.chain', '.chain_ID_position')
        if not os.path.exists(chain_index_file):
            self._die(f'ERROR: Chain file {chain_file} is not indexed')
        self.parse_chain_index(chain_index_file)
        ## extract chain headers
        self.get_chain_headers(chain_file)
        ## merge overlapping pseudogene projections and output BED lines
        self.prepare_bed_lines()
        

    def get_pgenes_from_projection_report(self, file: TextIO) -> None:
        """Retrieves pseudogene projection data from the TOGA2 projection table"""
        for line in file:
            data: List[str] = line.rstrip().split('\t')
            if not data:
                return
            if len(data) < 5:
                self._die(
                    'ERROR: Projection report contains less than five columns'
                )
            if data[0] == TRANSCRIPT:
                continue
            tr: str = data[0]
            ppgene_col: str = data[4]
            if ppgene_col == '0':
                continue
            pp_chains: List[str] = ppgene_col.split(',')
            self.tr2chain[tr].extend(pp_chains)
            self.chains = self.chains.union(pp_chains)

    def parse_chain_index(self, file: str) -> None:
        """
        Parses the chain index file. By default, expects the file with at least three 
        columns were the second column is a 
        """
        if file is None:
            return
        with open(file, 'r') as h:
            for line in h:
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) < 3:
                    self._die(
                        'ERROR: Chain index contains less than three columns'
                    )
                chain: str = data[0]
                if chain  not in self.chains:
                    continue
                start_byte: int = int(data[1])
                self.chain2byte[chain] = start_byte

    def get_chain_headers(self, file: str) -> None:
        all_chains: List[Tuple[str, int]] = sorted(
            self.chain2byte.items(), key=lambda x: x[1]
        )
        with open(file, 'rb') as h:
            for chain in all_chains:
                h.seek(chain[1])
                header: bytes = b''
                byte: bytes = b''
                while byte != b'\n':
                    byte = h.read(1)
                    header += byte
                utf_header: List[str] = header.decode('utf8').rstrip('\n').split(' ')
                if len(utf_header) != 13:
                    self._die(
                        'ERROR: Chain header contains improper number of fields'
                    )
                chain_num: str = utf_header[12]
                q_chrom: str = utf_header[7]
                q_strand: str = utf_header[9]
                q_start: int = int(utf_header[10])
                q_end: int = int(utf_header[11])
                if q_strand != '+':
                    q_size: int = int(utf_header[8])
                    tmp_start: int = q_start
                    q_start = q_size - q_end
                    q_end = q_size - tmp_start
                self.chain2coords[chain_num] = (q_chrom, q_start, q_end, q_strand)

    def prepare_bed_lines(self) -> None:
        """
        For each transcript, merges overlapping processed pseudogene projections,
        the prepares the BED line for resulting projections
        """
        for tr, chains in self.tr2chain.items():
            ## extract coordinates for all projections for a given transcript,
            ## group them by chromosome
            chrom2chains: Dict[int, List[str]] = defaultdict(list)
            headers: Dict[str, Tuple[str, int, int, str]] = {}
            for chain in chains:
                header: Tuple[str, int, int, str] = self.chain2coords[chain]
                chrom2chains[header[0]].append(chain)
                headers[chain] = header
            for chrom, chrom_chains in chrom2chains.items():
                ## for each chromosome bearing pseudogene projections, 
                ## estumate all overlap groups
                intersecting_projs: Dict[str, List[str]] = defaultdict(list)
                sorted_chains: List[str] = sorted(
                    chrom_chains, key=lambda x: headers[x][1:3]
                )
                for i, chain1 in enumerate(sorted_chains):
                    start1, end1, strand1 = headers[chain1][1:]
                    for chain2 in sorted_chains[i+1:]:
                        start2, end2, strand2 = headers[chain2][1:]
                        if strand1 != strand2:
                            continue
                        if start1 >= end2:
                            continue
                        if start2 >= end1:
                            break
                        intersecting_projs[chain1].append(chain2)
                    if chain1 not in intersecting_projs:
                        intersecting_projs[chain1] = []
                visited: List[str] = []
                for x in intersecting_projs:
                    if x in visited:
                        continue
                    current_component: List[str] = dfs(
                        intersecting_projs, x, visited, []
                    )
                    current_component.sort(key=lambda x: (headers[x][0], headers[x][1]))
                    first_proj: str = current_component[0]
                    last_proj: str = current_component[-1]
                    comp_start: int = headers[first_proj][1]
                    comp_end: int = headers[last_proj][2]
                    comp_strand: str = headers[last_proj][3]
                    name: str = f'{tr}#{",".join(current_component)}'
                    bed_line: List[str] = [
                        chrom, comp_start, comp_end, name, '0', comp_strand,
                        comp_start, comp_start, PINK 
                    ]
                    self.output.write('\t'.join(map(str, bed_line)) + '\n')


if __name__ == '__main__':
    PseudogeneTrackBuilder()