#!/usr/bin/env python3

"""
A module for CESAR job binning based on maximal memory requirements
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from collections import defaultdict, namedtuple
from dataclasses import dataclass
from heapq import heappop, heappush
from math import ceil
from modules.constants import PRE_CLEANUP_LINE
from modules.cesar_wrapper_constants import (
    DEF_BLOSUM_FILE, MIN_PROJ_OVERLAP_THRESHOLD,
    HG38_CANON_U2_ACCEPTOR, HG38_CANON_U2_DONOR,
    FIRST_ACCEPTOR, LAST_DONOR
)
from modules.shared import (
    CommandLineManager, get_upper_dir, intersection, CONTEXT_SETTINGS, SPLIT_JOB_HEADER
)
from pathlib import Path
from shared import get_connected_components
from shutil import which
from typing import Dict, Iterable, List, Optional, TextIO, Tuple, Type, Union

import click
import networkx as nx

__author__  = 'Yury V. Malovichko'
__credits__ = ['Bogdan Kirilenko', 'Michael Hiller']
__year__ = '2024'

TOGA2_ROOT: str = get_upper_dir(__file__, 4)

## define constants
# LOCATION: str = os.path.dirname(os.path.abspath(__file__))
BLOSUM_FILE: str = os.path.join(LOCATION, *DEF_BLOSUM_FILE)
HL_CESAR_PATH: str = os.path.join(
    os.path.sep,
    'projects',
    'hillerlab',
    'genome',
    'src',
    'TOGA_pub',
    'CESAR2.0',
    'cesar'
)
CESAR_WRAPPER_SCRIPT: str = os.path.join(PARENT, 'cesar_exec.py')
REDUNDANT_ENTRY: str = (
    'PROJECTION\t{}\t0\tRedundant projection to the given locus\tREDUNDANT\tN'
)
HEAVY_ENTRY: str = (
    'PROJECTION\t{}\t0\tMaximum memory limit exceeded\tHEAVY\tN'
)
CHIMERIC_ENTRY: str = (
    'PROJECTION\t{}\t0\tPotential chimeric projection\tCHIMERIC\tN'
)

HG38_CANON_U2_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *HG38_CANON_U2_ACCEPTOR)
HG38_CANON_U2_DONOR: str = os.path.join(TOGA2_ROOT, *HG38_CANON_U2_DONOR)
FIRST_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *FIRST_ACCEPTOR)
LAST_DONOR: str = os.path.join(TOGA2_ROOT, *LAST_DONOR)

OK: str = '.ok'
TOUCH: str = 'touch {}'

# @dataclass
# class ProjectionMeta:
#     __slots__ = ('chrom', 'start', 'stop', 'max_mem', 'sum_mem', 'path')
#     chrom: str
#     start: int
#     end: int
#     max_mem: float
#     sum_mem: float
#     path: str

ProjectionMeta: Type = namedtuple(
    'ProjectionMeta',
    ['name', 'chain', 'chrom', 'start', 'end', 'max_mem', 'sum_mem', 'path']
)

def fragmented_projection(chain_id: str) -> bool:
    return ',' in chain_id

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'memory_report',
    metavar='MEMORY_REPORT',
    type=click.File('r', lazy=True)
)
@click.argument(
    'job_directory',
    type=click.Path(exists=False),
    metavar='CESAR_JOB_DIRECTORY'
)
@click.argument(
    'cesar_output_directory',
    type=click.Path(exists=False),
    metavar='CESAR_OUTPUT_DIRECTORY'
)
@click.option(
    '--joblist_file',
    '-jl',
    type=click.Path(exists=False),
    metavar='JOBLIST',
    default=None,
    show_default=False,
    help=(
        'A path to joblist for slurm/Para '
        '[default: CESAR_JOB_DIRECTORY/cesar_joblist]. Note that if jobs are binned, '
        'each memory bin will get its own job list'
    )
)
@click.option(
    '--job_number',
    '-j',
    type=int,
    metavar='INT',
    default=300,
    show_default=True,
    help='A number of cluster jobs to split the commands into'
)
@click.option(
    '--memory_bins',
    '-b',
    type=str,
    metavar='BIN_LIST',
    default=None,
    show_default=True,
    help=(
        'A comma-separated list of memory bin caps, in GB. For each memory bin, '
        'job scheduling will be performed independently. If you want to process '
        'memory-intensive projections as a single cluster call, set the last value '
        'to "big" or enable the --allow_heavy_jobs flag'
    )
)
@click.option(
    '--job_nums_per_bin',
    '-jb',
    type=str,
    metavar='BIN_JOB_NUM_LIST',
    default=None,
    show_default=True,
    help=(
        'A comma-separated list of job numbers per memory bin. '
        'Job jumbers must follow in the same order as memory caps passed '
        'to --memory_bins option'
    )
)
@click.option(
    '--allow_heavy_jobs',
    '-ahj',
    type=bool,
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Aggregate all jobs exceeding the highest memory cap as a single joblist; '
         'if memory bins are provided, duplicates the "big" memory bin behavior'
)
@click.option(
    '--parallel_execution',
    '-p',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, CESAR alignment job partitions share the same output directory, '
        'with output file identity controlled by respective lock files'
    )
)
@click.option(
    '--paralog_list',
    '-pl',
    type=click.File('r', lazy=True),
    metavar='PARALOG_LIST',
    default=None,
    show_default=False,
    help='A single-column file containing known paralogous projections'
)
@click.option(
    '--processed_pseudogene_list',
    '-ppl',
    type=click.File('r', lazy=True),
    metavar='PARALOG_LIST',
    default=None,
    show_default=False,
    help='A single-column file containing known processed pseudogene projections'
)
@click.option(
    '--cesar_binary',
    '-cs',
    type=click.Path(exists=True),
    metavar='CESAR_BINARY',
    default=None,
    show_default=False,
    help='A path to the actual CESAR2.0 binary; default is set for Hiller '
         'lab Delta cluster'
)
@click.option(
    '--matrix',
    '-m',
    type=click.Path(exists=True),
    metavar='BLOSUM_MATRIX_FILE',
    default=BLOSUM_FILE,
    show_default=True,
    help='A file containing the protein alignment matrix'
)
@click.option(
    '--mask_n_terminal_mutations',
    '-m10m',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, masks mutations occurring in the first 10 percents '
         'of query projection length'
)
@click.option(
    '--rescue_missing_stop',
    '-rms',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, scans the downsteam of the search space for inframe stop codons '
         'in case the alignmentn does not end with one'
)
@click.option(
    '--filtered_bed_output',
    '-fbo',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, CESAR wrapper commands will not report missing and/or deleted'
         'exons in their output BED12 files'
)
@click.option(
    '--spliceai_correction_mode',
    '-scm',
    type=click.IntRange(min=0, max=7),
    metavar='MODE_NUM',
    default=0,
    show_default=False,
    help=(
        'Set the mode of SpliceAI-mediated exon boundary correction:\n'
        '0 - no correction [default; equivalent to not providing SpliceAI results directory];\n'
        '1 - use SpliceAI predictions to correct boundaries of missing and deleted exons;\n'
        '2 - correct mutated canonical U2 splice sites;\n'
        '3 - correct all canonical U2 splice sites in the presence of alternatives with higher SpliceAI support;\n'
        '4 - correct all canonical U2 as well as mutated GT-AG U12 splice sites;\n'
        '5 - correct all canonical U2 and mutated and/or unsupported  GT-AG U12 splice sites;\n'
        '6 - correct all canonical U2 and all U12 splice sites;\n'
        '7 - correct all U2 (including known non-canonical sites) and U12 splice sites'
    )
)
@click.option(
    '--min_splice_prob',
    '-msp',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.02,
    show_default=True,
    help='Minimum SpliceAI prediction probability to consider the splice site'
)
@click.option(
    '--splice_prob_margin',
    '-spm',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.02,
    show_default=True,
    help=(
        'For splice sites with SpliceAI support 0<x<min_splice_prob, ignore '
        'alternative sites with support < x + splice_prob_margin'
    )
)
@click.option(
    '--intron_gain_check',
    '-ig',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, performs SpliceAI-guided check for query-specific introns'
    )
)
@click.option(
    '--intron_gain_threshold',
    '-igt',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.9,
    show_default=True,
    help='Minimal intron gain threshold to consider'
)
@click.option(
    '--min_intron_prob_trusted',
    '-mipt',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.1,
    show_default=True,
    help=(
        'Minimal SpliceAI support for query-specific introns supported by the presence of '
        'both extensive alignment gaps and frame-shifting/nonsense mutations'
    )
)
@click.option(
    '--min_intron_prob_supported',
    '-mips',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.2,
    show_default=True,
    help=(
        'Minimal SpliceAI support for query-specific introns supported by the presence of '
        'either extensive alignment gaps and frame-shifting/nonsense mutations'
    )
)
@click.option(
    '--min_intron_prob_unsupported',
    '-mipu',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.8,
    show_default=True,
    help=(
        'Minimal SpliceAI support for query-specific introns not supported by '
        'either extensive alignment gaps and frame-shifting/nonsense mutations'
    )
)
@click.option(
    '--cesar_regular_acceptor',
    '-cra',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_CANON_U2_ACCEPTOR,
    show_default=True,
    help=(
        'Regular acceptor site file for intron gain CESAR alignment. '
        'The script does not evaluate intron class and canonicity '
        'so using canonical U2 profile for reference species is recommended'
    )
)
@click.option(
    '--cesar_regular_donor',
    '-crd',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_CANON_U2_DONOR,
    show_default=True,
    help=(
        'Regular acceptor site file for intron gain CESAR alignment. '
        'The script does not evaluate intron class and canonicity '
        'so using canonical U2 profile for reference species is recommended'
    )
)
@click.option(
    '--cesar_first_acceptor',
    '-cfa',
    type=click.Path(exists=True),
    metavar='PATH',
    default=FIRST_ACCEPTOR,
    show_default=True,
    help='Acceptor site profile for the first exon'
)
@click.option(
    '--cesar_last_donor',
    '-cld',
    type=click.Path(exists=True),
    metavar='PATH',
    default=LAST_DONOR,
    show_default=True,
    help='Donor site profile for the last exon'
)
@click.option(
    '--no_spliceai_correction',
    '-no_sai',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, SpliceAI data will be used exclusively for restoring missing '
         'exons, with no post-CESAR correction'
)
@click.option(
    '--correct_ultrashort_introns',
    '-c_si',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, corrects the introns shorter than 30 bp erroneously introduced '
        'by CESAR by treating them as precise intron deletions and respective '
        'insertion in the acceptor exon'
    )
) ## TODO: Make a default feature
@click.option(
    '--ignore_alternative_frame',
    '-no_alt_frame',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, codons in the alternative reading frame '
        '(=residing between compensated frameshifts) are ignored '
        'when computing sequence intactness features'
    )
)
@click.option(
    '--save_cesar_input',
    '-sci',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, saves CESAR input files to the temporary directory'
)
@click.option(
    '--rejection_report',
    '-rr',
    type=click.Path(exists=False),
    metavar='PATH',
    default=None,
    show_default=False,
    help=(
        'A path to save rejected projections to '
        '[default:PREPROCESS_JOB_DIRECTORY/genes_rejection_reason.tsv]'
    )
)
## benchmarking-related - REMOVE IN THE FINAL VERSION
@click.option(
    '--toga1_compatible',
    '-t1',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'Alignment procedure is fully TOGA1.0-compliant except for exonwise '
        'CESAR alignment; benchmarking feature, do not use in real runs'
    )
)
@click.option(
    '--verbose',
    '-v',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls the execution verbosity'
)


class CesarScheduler(CommandLineManager):
    """
    """
    __slots__ = [
        'memory_report', 'job_directory', 'cesar_output_directory',
        'joblist_file', 'jobnum', 'bins', 'max_bin', 'bin_job_nums', 'bin2jobnum',
        'allow_heavy_jobs', 'parallel_execution', 'paralog_list', 'processed_pseudogene_list',
        'cesar_binary', 'aa_matrix', 'mask_n_terminal_mutations',
        'rescue_missing_stop', 'report_raw_bed',
        'correction_mode', 'min_splice_prob', 'splice_prob_margin',
        'no_sai_correction', 'intron_gain_check', 'min_intron_gain_score',
        # 'min_intron_prob_gapped', 'min_intron_prob_ungapped', 
        'min_intron_prob_trusted', 'min_intron_prob_supported', 'min_intron_prob_unsupported',
        'regular_acceptor', 'regular_donor', 
        'first_acceptor', 'last_donor', 
        'correct_short_introns', 'ignore_alternative_frame',
        'save_cesar_input', 'v',

        'proj2max_mem', 'proj2sum_mem', 'proj2storage', 'rejected_transcripts',
        'heavier_jobs', 'heavy_job_nums', 'heavy_job_max_mem', 'jobs', 'job2mem',
        'joblist_descr', 'rejection_file',
        'toga1'
    ]

    def __init__(
        self,
        memory_report: click.File,
        job_directory: click.Path,
        cesar_output_directory: click.Path,
        joblist_file: Optional[click.Path],
        job_number: Optional[int],
        memory_bins: Optional[str],
        job_nums_per_bin: Optional[str],
        allow_heavy_jobs: Optional[bool],
        parallel_execution: Optional[bool],
        paralog_list: Optional[click.File],
        processed_pseudogene_list: Optional[click.File],
        cesar_binary: Optional[click.Path],
        matrix: Optional[click.Path],
        mask_n_terminal_mutations: Optional[bool],
        rescue_missing_stop: Optional[bool],
        filtered_bed_output: Optional[bool],
        spliceai_correction_mode: Optional[int],
        min_splice_prob: Optional[float],
        splice_prob_margin: Optional[float],
        intron_gain_check: Optional[bool],
        intron_gain_threshold: Optional[float],
        # min_intron_prob_gapped: Optional[float],
        # min_intron_prob_ungapped: Optional[float],
        min_intron_prob_trusted: Optional[float],
        min_intron_prob_supported: Optional[float],
        min_intron_prob_unsupported: Optional[float],
        cesar_regular_acceptor: Optional[click.Path],
        cesar_regular_donor: Optional[click.Path],
        cesar_first_acceptor: Optional[click.Path],
        cesar_last_donor: Optional[click.Path],
        no_spliceai_correction: Optional[bool],
        correct_ultrashort_introns: Optional[bool],
        ignore_alternative_frame: Optional[bool],
        save_cesar_input: Optional[bool],
        rejection_report: Optional[click.Path],
        toga1_compatible: Optional[bool],
        verbose: bool
    ) -> None:
        self.memory_report: click.File = memory_report
        self.job_directory: click.Path = job_directory
        self.cesar_output_directory: click.Path = cesar_output_directory
        self.joblist_file: click.Path = (
            joblist_file if joblist_file is not None else os.path.join(
                self.job_directory, 'joblist'
            )
        )
        self.bins: Union[List[int, str], None] = self.parse_bin_list(memory_bins)
        self.max_bin: Union[int, None] = (
            None if self.bins is None else max(x for x in self.bins if isinstance(x, int))
        )
        self.bin_job_nums: Union[List[int, str], None] = None if (
            self.bins is None or job_nums_per_bin is None
        ) else self.parse_bin_list(job_nums_per_bin)
        self.jobnum: int = max(
            job_number, 0 if self.bin_job_nums is None else sum(self.bin_job_nums)
        )
        self.allow_heavy_jobs: bool = allow_heavy_jobs or 'big' in self.bins
        self.parallel_execution: bool = parallel_execution
        self.paralog_list: Union[List[str], None] = (
            None if paralog_list is None else [x.rstrip() for x in paralog_list.readlines() if x]
        )
        self.processed_pseudogene_list: Union[List[str], None] = (
            None if processed_pseudogene_list is None else [
                x.rstrip() for x in processed_pseudogene_list.readlines() if x
            ]
        )
        if cesar_binary is None:
            cesar_in_path: str = which('cesar')
            if os.path.exists(HL_CESAR_PATH):
                self.cesar_binary: str = HL_CESAR_PATH
            elif cesar_in_path is not None:
                self.cesar_binary: str = cesar_in_path
            else:
                raise FileNotFoundError(
                    'CESAR executable was not provided, '
                    'with no default in the user PATH'
                )
        else:
            self.cesar_binary: str = cesar_binary
        self.aa_matrix: click.Path = matrix
        self.mask_n_terminal_mutations: bool = mask_n_terminal_mutations
        self.rescue_missing_stop: bool = rescue_missing_stop
        self.report_raw_bed: bool = not filtered_bed_output
        self.correction_mode: int = spliceai_correction_mode
        self.min_splice_prob: float = min_splice_prob
        self.splice_prob_margin: float = splice_prob_margin
        self.intron_gain_check: bool = intron_gain_check
        self.min_intron_gain_score: bool = intron_gain_threshold
        # self.min_intron_prob_gapped: float = min_intron_prob_gapped
        # self.min_intron_prob_ungapped: float = min_intron_prob_ungapped
        self.min_intron_prob_trusted: float = min_intron_prob_trusted
        self.min_intron_prob_supported: float = min_intron_prob_supported
        self.min_intron_prob_unsupported: float = min_intron_prob_unsupported
        self.regular_acceptor: click.Path = cesar_regular_acceptor
        self.regular_donor: click.Path = cesar_regular_donor
        self.first_acceptor: click.Path = cesar_first_acceptor
        self.last_donor: click.Path = cesar_last_donor
        self.no_sai_correction: bool = no_spliceai_correction
        self.correct_short_introns: bool = correct_ultrashort_introns
        self.ignore_alternative_frame: bool = ignore_alternative_frame
        self.save_cesar_input: bool = save_cesar_input
        self.toga1: bool = toga1_compatible
        self.v: bool = verbose

        self.proj2max_mem: Dict[str, float] = {}
        self.proj2sum_mem: Dict[str, int] = {}
        self.proj2storage: Dict[str, str] = {}
        self.rejected_transcripts: List[Tuple[str, int]] = []
        # self.heavier_jobs: List[str] = []
        self.heavy_job_nums: List[int] = []
        self.heavy_job_max_mem: int = 0
        self.jobs: Dict[int, List[str]] = defaultdict(list)
        self.job2mem: Dict[int, int] = {}
        self.joblist_descr: str = os.path.join(job_directory, 'joblist_description.tsv')
        self.rejection_file: str = (
            rejection_report if rejection_report is not None else os.path.join(
                job_directory, 'genes_rejection_reason.tsv'
            )
        )

        self.run()

    def run(self) -> None:
        ## initialize all the necessary stuff
        self._mkdir(self.job_directory)
        self._mkdir(self.cesar_output_directory)

        ## parse the memory report
        self.parse_memory_report()

        ## bin the jobs
        self.allocate_job_numbers()
        self.lpt()

        ## format CESAR commands and dump them to job partition files,
        ## then write the joblist file
        self.write_job_files()

        ## if any projections were discarded due to memory consumption reasons,
        ## write them to a separate file
        # if self.heavier_jobs:
        #     self.write_heavier_jobs()

        ## write memory requirements for each joblist
        self.write_joblist_mem_map()

        ## if any projections were discarded for any reasons,
        ## write the data to the rejection report file
        self.rejection_report()

    def parse_bin_list(self, bin_list: str) -> Union[List[Union[int, str]], None]:
        """
        Parses comma-separated list of values. Accepted values are either integer
        numbers or 'big' argument denoting bigger jobs
        """
        if bin_list is None:
            return None
        out_list: List[Union[int, str]] = []
        split_list: List[str] = bin_list.split(',')
        for item in split_list:
            if not item:
                continue
            if item.isdigit():
                out_list.append(int(item))
                continue
            if item == 'big':
                out_list.append(item)
                continue
            raise AttributeError(
                'Invalid value provided in the comma-separated list'
            )
        return out_list

    # def parse_memory_report(self) -> None: ## DONE
    #     """
    #     Given the path to a CESAR preprocessing report,
    #     parses the results producing a storage class instances
    #     """
    #     ## TODO:
    #     ## 1) Save the cumulative RAM requirements
    #     ## 2) Add the same-locus filter by maximal RAM requirements
    #     tr2proj2coords: Dict[str, Dict[str, Tuple[str, int]]] = defaultdict(dict)
    #     for line in self.memory_report.readlines():
    #         data: List[str] = line.rstrip().split('\t')
    #         tr: str = data[0]
    #         chain: str = data[1]
    #         proj: str = f'{tr}.{chain}' ## TODO: For some reasons, I'm still using legacy format in the memory file; switch to storing projection name in a single column and update the code accordingly
    #         max_mem: float = float(data[2].split()[0])
    #         sum_mem: int = ceil(float(data[3].split()[0]) + 0.1)
    #         chrom: str = data[6]
    #         start: int = int(data[7])
    #         stop: int = int(data[8])
    #         max_inter: int = int((stop - start) * 0.3)
    #         path: str = data[-1]
    #         ## check if there are other projections of the same transcripts
    #         ## corresponding to the same locus; if so, check the one with the least
    #         ## memory requirements
    #         if tr not in tr2proj2coords:
    #             self.proj2max_mem[proj] = max_mem
    #             self.proj2sum_mem[proj] = sum_mem
    #             self.proj2storage[proj] = path
    #             tr2proj2coords[tr][chain] = (chrom, start, stop)
    #             print(f'First occurrence; adding {proj}')
    #         else:
    #             to_del: List[str] = []
    #             for _chain in tr2proj2coords[tr]:
    #                 _chrom, _start, _stop = tr2proj2coords[tr][_chain]
    #                 if chrom != _chrom:
    #                     continue
    #                 _max_inter: int = int((_stop - _start) * 0.3)
    #                 inter: int = intersection(start, stop, _start, _stop)
    #                 if inter >= max_inter or inter >= _max_inter:
    #                     _proj: str = f'{tr}.{_chain}'
    #                     _max_mem: float = self.proj2max_mem[_proj]
    #                     if max_mem > _max_mem:
    #                         print(f'{proj} intersects {_proj} but consumes more memory; breaking')
    #                         rej_reason: Tuple[str] = REDUNDANT_ENTRY.format(proj)
    #                         self.rejected_transcripts.append(rej_reason)
    #                         break
    #                     elif max_mem == _max_mem and int(_chain) > (chain):
    #                         rej_reason: Tuple[str] = REDUNDANT_ENTRY.format(proj)
    #                         self.rejected_transcripts.append(rej_reason)
    #                         break
    #                     else:
    #                         print(f'{proj} intersects {_proj} and consumes less memory; removing {_proj}')
    #                         del self.proj2max_mem[_proj]
    #                         del self.proj2sum_mem[_proj]
    #                         del self.proj2storage[_proj]
    #                         to_del.append(_chain)
    #                         rej_reason: str = REDUNDANT_ENTRY.format(_proj)
    #                         self.rejected_transcripts.append(rej_reason)
    #             else:
    #                 self.proj2max_mem[proj] = max_mem
    #                 self.proj2sum_mem[proj] = sum_mem
    #                 self.proj2storage[proj] = path
    #                 tr2proj2coords[tr][chain] = (chrom, start, stop)
    #                 print(f'No intersections among the previous records; adding {proj}')
    #             for chain_to_del in to_del:
    #                 del tr2proj2coords[tr][chain_to_del]
    #         print('*'*30)

    def parse_memory_report(self) -> None:
        """
        Given the path to a CESAR preprocessing report,
        parses the results producing a storage class instances
        """
        ## TODO:
        # tr2proj2coords: Dict[str, Dict[str, Tuple[str, int]]] = defaultdict(dict)
        tr2chrom2graph: Dict[str, Dict[str, nx.Graph]] = defaultdict(dict)
        # print('UAAAAA')
        for line in self.memory_report.readlines():
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if data[0] == 'transcript':
                continue
            tr: str = data[0]
            chain: str = data[1]
            proj: str = f'{tr}#{chain}' ## TODO: For some reasons, I'm still using legacy format in the memory file; switch to storing projection name in a single column and update the code accordingly
            max_mem: float = float(data[2].split()[0])
            sum_mem: int = ceil(float(data[3].split()[0]) + 0.1)
            chrom: str = data[7]
            start: int = int(data[8])
            stop: int = int(data[9])
            max_inter: int = ceil((stop - start) * MIN_PROJ_OVERLAP_THRESHOLD)
            path: str = data[-1]
            entry: ProjectionMeta = ProjectionMeta(
                proj, chain, chrom, start, stop, max_mem, sum_mem, path
            )
            proj_is_ppgene: bool = (
                self.processed_pseudogene_list is not None and 
                proj in self.processed_pseudogene_list
            )
            # print(f'{entry=}')
            # print(f'{tr2chrom2graph=}, {tr=}, {chrom=}')
            if tr not in tr2chrom2graph.keys() or chrom not in tr2chrom2graph[tr].keys():
                # print(f'Creating a new graph for {proj}')
                new_graph: nx.Graph = nx.Graph()
                new_graph.add_node(entry)
                tr2chrom2graph[tr][chrom] = new_graph
            else:
                # print(f'Adding {proj} to an existing graph')
                tr2chrom2graph[tr][chrom].add_node(entry)
                ## do not intersect fragmented and whole projections
                if fragmented_projection(chain):
                    continue
                for node in tr2chrom2graph[tr][chrom].nodes():
                    ## do not intersect orthologs/paralogs and retrocopy projections
                    if self.processed_pseudogene_list is not None:
                        in_is_ppgene: bool = (
                            self.processed_pseudogene_list is not None and 
                            node.name in self.processed_pseudogene_list
                        )
                        if proj_is_ppgene != in_is_ppgene:
                            continue
                    ## if the current projection and a node in the subgraph intersect,
                    ## traverse an edge between them
                    _max_inter: int = ceil((node.end - node.start) * MIN_PROJ_OVERLAP_THRESHOLD)
                    inter: int = intersection(start, stop, node.start, node.end)
                    if inter >= max_inter or inter >= _max_inter:
                        tr2chrom2graph[tr][chrom].add_edge(entry, node)
        ## now, resolve the resulting graph
        for tr, chroms in tr2chrom2graph.items():
            for chrom in chroms:
                proj_graph: nx.Graph = tr2chrom2graph[tr][chrom]
                # debug = any(x.name == 'ENST00000006724.CEACAM7.201279' for x in proj_graph.nodes())
                # if debug:
                #     for ooo in proj_graph.nodes():
                #         print(ooo)
                #     print('-'*30)
                #     for aaa in proj_graph.edges():
                #         print(aaa)
                #     print('-'*30 + '\n')
                # print(proj_graph)
                # print(f'{proj_graph.nodes()=}')
                # print(f'{proj_graph.edges()=}')
                if len(proj_graph.nodes()) == 1:
                    # print('Solitary node:')
                    sole_node: ProjectionMeta = next(iter(proj_graph.nodes()))
                    proj_: str = sole_node.name
                    self.proj2max_mem[proj_] = sole_node.max_mem
                    self.proj2sum_mem[proj_] = sole_node.sum_mem
                    self.proj2storage[proj_] = sole_node.path
                    continue
                components: List[nx.Graph] = get_connected_components(proj_graph)
                for component in components:
                    comp: nx.Graph = component.copy()
                    ## remove potential chimeric projections
                    art_nodes: List[ProjectionMeta] = list(
                        nx.articulation_points(comp)
                    )
                    # if debug and art_nodes:
                    #     print(f'The following nodes are likely chimeric: {[x.name for x in art_nodes]}')
                    for art_node in art_nodes:
                        rej_reason: Tuple[str] = CHIMERIC_ENTRY.format(art_node.name)
                        self.rejected_transcripts.append(rej_reason)
                    comp.remove_nodes_from(art_nodes)
                    subcliques: List[nx.Graph] = get_connected_components(comp)
                    for subclique in subcliques:
                        # print(f'Starting a new subclique: {[x.name for x in subclique.nodes()]}')
                        # min_mem: float = min(x.max_mem for x in subclique.nodes())
                        # min_mem_nodes: List[ProjectionMeta] = [
                        #     x for x in subclique.nodes() if x.max_mem == min_mem
                        # ]
                        # if len(min_mem_nodes) > 1:
                        #     min_mem_nodes.sort(key=lambda x: int(x.chain))
                        # # if debug:
                        # #     for eee in min_mem_nodes:
                        # #         print(f'min_mem_node: {eee}')
                        # best_pick_node: ProjectionMeta = min_mem_nodes[0]
                        # # if debug:
                        # #     print(f'{best_pick_node=}')
                        # # print('-'*30 + '\n')
                        # proj_: str = best_pick_node.name
                        # self.proj2max_mem[proj_] = best_pick_node.max_mem
                        # self.proj2sum_mem[proj_] = best_pick_node.sum_mem
                        # self.proj2storage[proj_] = best_pick_node.path
                        # for redundant_node in subclique.nodes:
                        #     if redundant_node.name == proj_:
                        #         continue
                        #     rej_reason: Tuple[str] = REDUNDANT_ENTRY.format(
                        #         redundant_node.name
                        #     )
                        #     self.rejected_transcripts.append(rej_reason)
                        for node in subclique.nodes():
                            proj_: str = node.name
                            if self.max_bin is None or node.max_mem < self.max_bin or self.allow_heavy_jobs:
                                self.proj2max_mem[proj_] = node.max_mem
                                self.proj2sum_mem[proj_] = node.sum_mem
                                self.proj2storage[proj_] = node.path
                            else:
                                rej_reason: Tuple[str] = REDUNDANT_ENTRY.format(proj_)
                                self.rejected_transcripts.append(rej_reason)
                        # print('-'*30)
                    # print('*'*30)
        # print(f'{self.proj2storage=}')


    def allocate_job_numbers(self) -> None:
        """
        Maps memory bins to joblist sizes
        """
        ## case 1: jobs are not binned anyhow: assign all jobs to a single bin
        if self.bins is None and self.bin_job_nums is None:
            self.bin2jobnum: Dict[int, int] = {0: self.jobnum}
            return
        ## case 2: jobs were binned according to memory caps but job numbers per
        ## each bin were not provided; job numbers will be allocated after projections
        ## are split into memory buckets
        if self.bin_job_nums is None:
            self.bin_job_nums: None = None
            return
        ## case 3: bin caps and binwise job numbers are properly provided;
        ## map bins to job numbers, set placeholders if job number list is shorter
        ## than memory cap list
        self.bin2jobnum: Dict[Union[int, str], int] = dict(zip(self.bins, self.bin_job_nums))
        curr_sum: int = sum(self.bin2jobnum.values())
        for bin in self.bins:
            if bin not in self.bin2jobnum:
                job_num: int = max(self.jobs - curr_sum, 1)
                self.bin2jobnum[bin] = job_num
                curr_sum += job_num

    def lpt(self) -> None: ## TODO: Ask Bogdan and/or Michael about bigmem logic
        """
        A naïve implementation of longest-processing-time-first (LPT)
        scheduling algorithm; given a dictionary of maximum memory requirements per
        projection, returns an iterable of memory-balanced cluster jobs
        """
        memory_buckets: Dict[int, List[Tuple[str, float]]] = defaultdict(list)
        ## Step 1: If memory caps are provided, split projections into memory bins
        for proj, max_mem in self.proj2max_mem.items():
            max_mem = ceil(max_mem + 0.1)
            # mem: int = self.proj2sum_mem[proj]
            if self.bins is not None:
                for bin in self.bins:
                    if bin == 'big' or max_mem <= bin:
                        memory_buckets[bin].append((proj, max_mem))
                        if bin == 'big':
                            self.heavy_job_max_mem = max(self.heavy_job_max_mem, max_mem)
                        break
                else:
                    if 'big' in self.bins or self.allow_heavy_jobs:
                        memory_buckets['big'].append((proj, max_mem))
                        self.heavy_job_max_mem = max(self.heavy_job_max_mem, max_mem)
                    else:
                        rej_reason: Tuple[str] = HEAVY_ENTRY.format(proj)
                        self.rejected_transcripts.append(rej_reason) ## In need of Michael's advice here
            else:
                memory_buckets[0].append((proj, max_mem))

        ## Step 2: If caps are set, calculate the number of jobs corresponding to each bin
        if self.bins is not None:
            job_heap: Dict[int, List[Tuple[float, int]]] = {}
            job_counter: int = 0
            ## if bin-to-job-number mapping was provided, organise heap according
            ## to the defined mapping
            if self.bin2jobnum is not None:
                for bin in self.bins:
                    jobs_per_bin: int = self.bin2jobnum[bin]
                    job_heap[bin] = [
                        (0, i) for i in range(job_counter, job_counter + jobs_per_bin)
                    ]
                    self.job2mem = {
                        **self.job2mem,
                        **{
                            k:(
                                bin if isinstance(bin, int) else 0
                            ) for k in range(job_counter, job_counter + jobs_per_bin)
                        }
                    }
                    if bin == 'big':
                        self.heavy_job_nums.extend(
                            [i for i in range(job_counter, job_counter + jobs_per_bin)]
                        )
                    job_counter += jobs_per_bin
            ## otherwise, allocate job numbers according to the proportion of
            ## respective to
            else:
                all_jobs: int = len(self.proj2max_mem)
                mem_class_proportions: Dict[int, float] = {
                    k:len(v) / all_jobs for k,v in memory_buckets.items()
                }
                available_jobs: int = self.jobnum
                for bin, proportion in sorted(mem_class_proportions.items(), key=lambda x: x[1]):
                    prop_share: float = self.jobnum * proportion
                    jobs_per_bin: int = min(ceil(prop_share), available_jobs)
                    available_jobs -= jobs_per_bin
                    job_heap[bin] = [
                        (0, i) for i in range(job_counter, job_counter + jobs_per_bin)
                    ]
                    self.job2mem = {
                        **self.job2mem,
                        **{
                            k:(bin if isinstance(bin, int) else 0)
                            for k in range(job_counter, job_counter + jobs_per_bin)
                        }
                    }
                    if bin == 'big':
                        self.heavy_job_nums.extend(
                            [i for i in range(job_counter, job_counter + jobs_per_bin)]
                        )
                    job_counter += jobs_per_bin
        else:
            job_heap: Dict[int, List[Tuple[float, int]]] = {
                0: [(0, i) for i in range(self.jobnum)]
            }
            self.job2mem = {k: 0 for k in range(self.jobnum)}
        balanced_jobs: Dict[int, List[str]] = defaultdict(list)

        ## Step 3: For each cap, split projections into jobs in the LPT fashion
        for bin in memory_buckets:
            ## First, ort jobs in each bucket by memory requirements in the descending order
            memory_buckets[bin].sort(key=lambda x: -x[1])
            ## and now, assign each projection to a job
            for proj, mem in memory_buckets[bin]:
                ## get the current least memory-expensive job
                total_mem, jobid = heappop(job_heap[bin])
                ## assign current projection to this job
                self.jobs[jobid].append(proj)
                ## push the updated job back
                heappush(job_heap[bin], (total_mem + mem, jobid))
                ## update maximum memory requirements for job if jobs were not binned
                # if not self.bins:
                self.job2mem[jobid] = max(self.job2mem[jobid], mem)


    def write_job_files(self) -> None: ## DONE?
        """
        """
        # print(f'{self.jobs=}')
        # print(f'{self.job2mem=}')
        for jobid in self.jobs:
            mem: int = self.job2mem[jobid]
            cmds: List[str] = []
            if self.bins:
                if jobid in self.heavy_job_nums:
                    joblist: str = f'{self.joblist_file}_{self.heavy_job_max_mem}GB.txt'
                else:
                    joblist: str = f'{self.joblist_file}_{mem}GB.txt'
            else:
                joblist: str = f'{self.joblist_file}.txt'
            with open(joblist, 'a') as h1:
                cesar_output: str = Path(
                    os.path.join(self.cesar_output_directory, f'batch{jobid}')
                ).absolute()
                for proj in self.jobs[jobid]:
                    proj_name_split: List[str] = proj.split('#')
                    tr: str = '#'.join(proj_name_split[:-1])
                    chain: str = proj_name_split[-1]
                    input_dir: str = Path(
                        self.proj2storage[proj]
                    ).absolute()
                    input_file: str = os.path.join(input_dir, 'exon_storage.hdf5')
                    cmd: str = (
                        f'{CESAR_WRAPPER_SCRIPT} "{tr}" '
                        f'{chain} {input_file} -cs {self.cesar_binary} '
                        f'-scm {self.correction_mode} '
                        f'-msp {self.min_splice_prob} '
                        f'-spm {self.splice_prob_margin} '
                    )
                    if self.toga1:
                        cmd += ' -t1 '
                    if self.aa_matrix:
                        cmd += f' --matrix {self.aa_matrix}'
                    if self.mask_n_terminal_mutations:
                        cmd += ' --mask_terminal_mutations'
                    if self.rescue_missing_stop:
                        cmd += ' --rescue_missing_stop'
                    if not self.report_raw_bed:
                        cmd += ' --filtered_bed_output'
                    if self.paralog_list is not None and proj in self.paralog_list:
                        cmd += ' --paralogous_projection'
                    if self.processed_pseudogene_list is not None and proj in self.processed_pseudogene_list:
                        cmd += ' --processed_pseudogene_projection'
                    if self.no_sai_correction:
                        cmd += ' --no_spliceai_correction'
                    if self.intron_gain_check:
                        cmd += (
                            ' --intron_gain_check '
                            f' --intron_gain_threshold {self.min_intron_gain_score} '
                            # f' --min_intron_prob_gapped {self.min_intron_prob_gapped} '
                            # f' --min_intron_prob_ungapped {self.min_intron_prob_ungapped} '
                            f' --min_intron_prob_trusted {self.min_intron_prob_trusted} '
                            f' --min_intron_prob_supported {self.min_intron_prob_supported} '
                            f' --min_intron_prob_unsupported {self.min_intron_prob_unsupported} '
                            f' -cra {self.regular_acceptor} '
                            f' -crd {self.regular_donor} '
                            f' -cfa {self.first_acceptor} '
                            f' -cld {self.last_donor} '
                        )
                    if self.correct_short_introns:
                        cmd += ' --correct_ultrashort_introns'
                    if self.ignore_alternative_frame:
                        cmd += ' --ignore_alternative_frame'
                    if self.save_cesar_input:
                        cmd += ' --save_cesar_input'
                    if self.parallel_execution:
                        cmd += ' --parallel_job'
                    # cmd += ' -v'
                    cmd += f' -o {cesar_output}'
                    cmds.append(cmd)
                jobfile_dest: str = os.path.join(self.job_directory, f'batch{jobid}.ex')
                with open(jobfile_dest, 'w') as h2:
                    for line in SPLIT_JOB_HEADER:
                        h2.write(line + '\n')
                    h2.write(PRE_CLEANUP_LINE.format(cesar_output) + '\n')
                    for line in cmds:
                        h2.write(line + '\n')
                    ok_file: str = os.path.join(self.cesar_output_directory, f'batch{jobid}', OK)
                    h2.write(TOUCH.format(ok_file) + '\n')
                ## make the partition files executable
                file_mode: bytes = os.stat(jobfile_dest).st_mode
                file_mode |= (file_mode & 0o444) >> 2
                os.chmod(jobfile_dest, file_mode)
                h1.write(jobfile_dest + '\n')

    # def write_heavier_jobs(self) -> None:
    #     """
    #     Write a list of projections discarded due to set memory requirements
    #     to a text file
    #     """
    #     heavy_job_file: str = os.path.join(
    #         self.job_directory, 'heavier_jobs.txt'
    #     )
    #     with open(heavy_job_file, 'w') as h:
    #         for proj in self.heavier_jobs:
    #             h.write(proj + '\n')

    def write_joblist_mem_map(self) -> None: ## DONE
        """
        Write joblist-to-maximum-memory
        """
        with open(self.joblist_descr, 'w') as h:
            if self.bins is not None:
                for bin in self.bins:
                    if bin == 'big':
                        max_mem: int = self.heavy_job_max_mem
                    else:
                        max_mem: int = bin
                    joblist: str = f'{self.joblist_file}_{max_mem}GB.txt'
                    h.write(f'{joblist}\t{max_mem}\n')
            else:
                joblist: str = f'{self.joblist_file}.txt'
                max_mem: int = max(self.job2mem.values())
                h.write(f'{joblist}\t{max_mem}\n')

    def rejection_report(self) -> None:
        if not self.rejected_transcripts:
            return
        with open(self.rejection_file, 'w') as h:
            for line in self.rejected_transcripts:
                h.write(line + '\n')


if __name__ == '__main__':
    CesarScheduler()
