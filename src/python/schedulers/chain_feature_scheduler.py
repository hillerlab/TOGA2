#!/usr/bin/env python3

"""
Schedules chain-projection pair feature extraction jobs
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from _chain_bed_intersect import retrieve_chain_headers, get_bed_coords
from collections import defaultdict
from heapq import heappop, heappush
from modules.shared import CommandLineManager, CONTEXT_SETTINGS, intersection, SPLIT_JOB_HEADER
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import click

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = ('Bogdan M. Kirilenko', 'Michael Hiller')
__all__ = (None)

# FEATURE_EXTRACTION_SCRIPT: str = 'feature_extractor.py'
FEATURE_EXTRACTION_SCRIPT: str = os.path.join(PARENT, 'chain_runner.py')
# print(f'{FEATURE_EXTRACTION_SCRIPT=}')
CMD_STUB: str = f'{FEATURE_EXTRACTION_SCRIPT} {{}} {{}} -i {{}} -o {{}}'
UNCOV_CHROM: str = (
    'TRANSCRIPT\t{}\t0\tNo chain corresponds to reference chromosome\tCHROM_UNALIGNED\tN'
)
UNALIGNED_TR: str = (
    'TRANSCRIPT\t{}\t0\tNo chain corresponds to the transcript\tTRANSCRIPT_UNALIGNED\tN'
)
OK: str = 'ok'
TOUCH: str = 'touch {}'

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'chain_file',
    type=click.Path(exists=True),
    metavar='CHAIN_FILE'
)
@click.argument(
    'bed_file',
    type=click.Path(exists=True),
    metavar='BED_FILE'
)
@click.argument(
    'job_directory',
    type=click.Path(exists=False),
    metavar='JOB_DIR'
)
@click.argument(
    'data_directory',
    type=click.Path(exists=False),
    metavar='DATA_DIR'
)
@click.argument(
    'result_directory',
    type=click.Path(exists=False),
    metavar='RESULT_DIR'
)
@click.option(
    '--job_number',
    '-j',
    type=int,
    default=500, ## NOTE: Default value in TOGA 1.0 is 800,
    show_default=True,
    help='A number of cluster jobs to split the overall command list into'
)
@click.option(
    '--job_list',
    '-jl',
    type=click.Path(exists=False),
    default=None,
    show_default=False,
    help='A job list file for cluster execution [default: JOB_DIR/joblist]'
)
@click.option(
    '--rejection_log',
    '-r',
    type=click.Path(exists=False),
    default=None,
    show_default=False,
    help=(
        'A path to write the rejected transcript data to [default: RESULT_DIR/]'
    )
)
@click.option(
    '--log_name',
    '-ln',
    type=str,
    metavar='STR',
    default=None,
    show_default=True,
    help='Logger name to use; relevant only upon main class import'
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

class ChainFeatureScheduler(CommandLineManager):

    """
    """

    __slots__ = [
        'chain_file', 'bed_file', 'chain_coords', 'bed_coords',
        'job_dir', 'data_dir', 'res_dir',
        'job_number', 'joblist', 'rejection_log',
        'chain2trs', 'rejected_transcripts', 'jobs',
        'log_file'
    ]

    def __init__(
        self,
        chain_file: click.Path,
        bed_file: click.Path,
        job_directory: click.Path,
        data_directory: click.Path,
        result_directory: click.Path,
        job_number: Optional[int],
        job_list: Optional[click.Path],
        rejection_log: Optional[click.Path],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)

        if os.path.islink(chain_file):
            self.chain_file: str = os.readlink(chain_file)
        else:
            self.chain_file: str = Path(chain_file).absolute()
        self.chain_coords: Dict[str, List[Tuple[str, int]]] = retrieve_chain_headers(chain_file)
        self.bed_file: str = Path(bed_file).absolute()
        self.bed_coords: Dict[str, List[Tuple[str, int]]] = get_bed_coords(bed_file)
        self.job_dir: str = Path(job_directory).absolute()
        self.data_dir: str = Path(data_directory).absolute()
        self.res_dir: str = Path(result_directory).absolute()
        self.joblist: str = (
            job_list if job_list is not None else os.path.join(self.job_dir, 'joblist')
        )
        self.job_number: int = job_number
        self.rejection_log: str = (
            rejection_log if rejection_log is not None else (
                os.path.join(self.job_dir, 'rejection_log.tsv')
            )
        )
        self.chain2trs: Dict[str, List[str]] = defaultdict(list)
        self.jobs: Dict[int, List[Tuple[str, str]]] = defaultdict(list)
        self.rejected_transcripts: Set[str] = set()
        self.run()

    def _mkdir(self, d: str) -> None:
        """Safe directory creation function"""
        try:
            os.makedirs(d)
        except FileExistsError:
            pass

    def run(self) -> None:
        """
        """
        ## create the output directories
        self._mkdir(self.job_dir)
        self._mkdir(self.data_dir)
        self._mkdir(self.res_dir)

        ## intersect chain and transcript coordinates, define which transcripts will
        self.intersect_chains_and_transcripts()

        ## balance jobs with LPT
        self.lpt()

        ## write job- and job list files
        self.write_jobs()

        ## if any transcripts were discarded due to any reason, write those as well
        self.rejection_report()

    def intersect_chains_and_transcripts(self) -> None:
        """
        """
        for chrom in self.bed_coords:
            if chrom not in self.chain_coords:
                for tr in self.bed_coords[chrom]:
                    rej_report: str = UNCOV_CHROM.format(tr[0])
                    self.rejected_transcripts.add(rej_report)
                continue
            curr_tr: int = 0
            curr_tr_start: int = 0
            curr_tr_end: int = 0
            aligned_transcripts: Set[str] = set()
            prev_start, prev_stop = 0, 0
            pprev_stop = 0
            first_tr: int = 0
            tr_name = None
            for i, chain in enumerate(self.chain_coords[chrom]):
                chain_id: str = chain[0]
                # print(f'{first_tr=}')
                chain_start: int = chain[1]
                chain_stop: int = chain[2]
                intersects_prev: bool = intersection(
                    prev_start, prev_stop, chain_start, chain_stop
                ) > 0
                # if intersects_prev:
                # tr_to_start: int = (
                #     first_tr if intersects_prev else curr_tr
                # )
                tr_to_start: int = (
                    first_tr if chain_start < pprev_stop else curr_tr
                )
                # if v:
                #     print(f'{chain_id}, {chain_start=}, {chain_stop=}, {prev_start=}, {prev_stop=}, {intersects_prev=}, {first_tr=}, {curr_tr=}, {tr_to_start=}, {tr_name=}')
                pprev_stop = max(pprev_stop, chain_stop)
                # first_tr = None
                found_first_tr: bool = False
                prev_start, prev_stop = chain_start, chain_stop
                for j, transcript in enumerate(
                    self.bed_coords[chrom][tr_to_start:], start=tr_to_start
                ):
                    tr_name = transcript[0]
                    tr_start = transcript[1]
                    tr_stop = transcript[2]
                    if chain_start >= tr_stop: ## chain lies dowstream to this transcript; proceed further
                        if tr_name not in aligned_transcripts:
                            rej_report: str = UNALIGNED_TR.format(tr_name)
                            self.rejected_transcripts.add(rej_report)
                        continue
                    if chain_stop <= tr_start: ## chain lies upstream to the next transcript; the following ones are guaranteed to lie dowstream as well
                        ## since the next chain can start within the last
                        ## intersected locus, do not update curr_tr
                        # if intersection(
                        #     tr_start, tr_stop, curr_tr_start, curr_tr_end
                        # ) > 0 and tr_stop > curr_tr_end:
                        #     print(f'BREAK: Updating current transcript; previous transcript is {curr_tr} ({self.bed_coords[chrom][curr_tr]}); new start is {j} ({self.bed_coords[chrom][j]}); transcript to halt is {tr_name}')
                        #     curr_tr_start = tr_start
                        #     curr_tr_end = tr_stop
                        #     curr_tr = j
                        break
                    ## this is certainly an intersecting pair
                    if intersection(tr_start, tr_stop, curr_tr_start, curr_tr_end) < 0:
                        ## this transcript does not intersect with the previous transcript
                        ## or group of transcripts
                        ## reset pointer and locus coordinates
                        curr_tr_start = tr_start
                        curr_tr_end = tr_stop
                        curr_tr = j
                    else:
                        ## this transcript intersects with the previous transcript
                        ## extend locus coordinates without updating the pointer
                        curr_tr_start = min(curr_tr_start, tr_start)
                        curr_tr_end = max(curr_tr_end, tr_stop)
                    # if tr_stop > curr_tr_end:
                    #     curr_tr_start = tr_start
                    #     curr_tr_end = tr_stop
                    #     curr_tr = j
                    if not found_first_tr:
                        # print('OOGA')
                        found_first_tr = True
                        first_tr = j
                    self.chain2trs[chain_id].append(tr_name)
                    aligned_transcripts.add(tr_name)
            ## add all the unaligned downstream transcripts to the rejection report
            # last_tr: str = self.bed_coords[chrom][curr_tr][0]
            # if last_tr not in aligned_transcripts:
            #     rej_report: str = UNALIGNED_TR.format(last_tr)
            #     self.rejected_transcripts.append(rej_report)
            for tr in self.bed_coords[chrom][curr_tr:]:
                if tr[0] not in aligned_transcripts:
                    rej_report: str = UNALIGNED_TR.format(tr[0])
                    self.rejected_transcripts.add(rej_report)

    def lpt(self) -> None:
        """
        Splits a list of chain-transcript tuples into a defined number of jobs in
        a Longest-Processing-Time-First (LPT) manner. Chain IDs are used as substitutes
        for processing time estimates.
        """
        job_heap: List[Tuple[int, int]] = [(0, i) for i in range(self.job_number)]
        for chain, transcripts in self.chain2trs.items():
            ## get the current fastest job
            chain_id_sum, jobid = heappop(job_heap)
            trs: str = ','.join(transcripts)
            ## assign current chain-transcript tuple to this job
            self.jobs[jobid].append((chain, trs))
            ## convert chain ID into a numeric value
            chain_id: int = int(chain)
            ## push the updated job back
            heappush(job_heap, (chain_id_sum + chain_id, jobid))

    ## TODO: Write chain->transcripts to a separate input file,
    ## organize commands in "one job - one command (I/O operation)" fashion
    def write_jobs(self) -> None:
        """
        Write jobs to jobfiles, along with jobfile paths being written to job list
        """
        with open(self.joblist, 'w') as h1:
            for job_id in self.jobs:
                jobfile: str = os.path.join(self.job_dir, f'batch{job_id}.ex')
                input_file: str = os.path.join(self.data_dir, f'batch{job_id}.tsv')
                with open(input_file, 'w') as h2:
                    for chain, transcripts in self.jobs[job_id]:
                        h2.write(f'{chain}\t{transcripts}\n')
                output_file: str = os.path.join(self.res_dir, f'batch{job_id}')
                # self._mkdir(output_dir)
                with open(jobfile, 'w') as h3:
                    h3.write('\n'.join(SPLIT_JOB_HEADER) + '\n')
                    cmd: str = CMD_STUB.format(
                        self.bed_file,
                        self.chain_file,
                        input_file,
                        output_file
                    )
                    h3.write(cmd + '\n')
                    ok_file: str = f'{output_file}_{OK}'
                    h3.write(TOUCH.format(ok_file) + '\n')
                file_mode: bytes = os.stat(jobfile).st_mode
                file_mode |= (file_mode & 0o444) >> 2
                os.chmod(jobfile, file_mode)
                h1.write(jobfile + '\n')

    def rejection_report(self) -> None:
        """
        """
        if not self.rejected_transcripts:
            return
        with open(self.rejection_log, 'w') as h:
            for rej_report in self.rejected_transcripts:
                h.write(rej_report + '\n')

if __name__ == '__main__':
    ChainFeatureScheduler()
