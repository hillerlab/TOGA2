#!/usr/bin/env python3

"""
Schedules CESAR preprocessing jobs for optimal cluster performance
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.append(PARENT)

from collections import defaultdict
from modules.cesar_wrapper_constants import *
from modules.constants import PRE_CLEANUP_LINE
from heapq import heappop, heappush
from math import ceil
from pathlib import Path
from modules.shared import CommandLineManager, get_upper_dir, CONTEXT_SETTINGS, SPLIT_JOB_HEADER
from shutil import which
from typing import Dict, Iterable, List, Optional, TextIO, Union

import click

TOGA2_ROOT: str = get_upper_dir(__file__, 4)

CESAR_PREPROCESS_SCRIPT: str = os.path.join(PARENT, 'cesar_preprocess.py')
HL_BW2W_PATH: str = os.path.join(
    os.path.sep,
    'projects',
    'hillerlab',
    'genome',
    'bin',
    'x86_64',
    'bigWigToWig'
)
HG38_CANON_U2_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *HG38_CANON_U2_ACCEPTOR)
HG38_CANON_U2_DONOR: str = os.path.join(TOGA2_ROOT, *HG38_CANON_U2_DONOR)
HG38_NON_CANON_U2_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *HG38_NON_CANON_U2_ACCEPTOR)
HG38_NON_CANON_U2_DONOR: str = os.path.join(TOGA2_ROOT, *HG38_NON_CANON_U2_DONOR)
HG38_CANON_U12_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *HG38_CANON_U12_ACCEPTOR)
HG38_CANON_U12_DONOR: str = os.path.join(TOGA2_ROOT, *HG38_CANON_U12_DONOR)
HG38_NON_CANON_U12_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *HG38_NON_CANON_U12_ACCEPTOR)
HG38_NON_CANON_U12_DONOR: str = os.path.join(TOGA2_ROOT, *HG38_NON_CANON_U12_DONOR)
FIRST_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *FIRST_ACCEPTOR)
LAST_DONOR: str = os.path.join(TOGA2_ROOT, *LAST_DONOR)
# HL_COMMON_ACCEPTOR: str = os.path.join(*HL_COMMON_ACCEPTOR)
# HL_COMMON_DONOR: str = os.path.join(*HL_COMMON_DONOR)
# HL_FIRST_ACCEPTOR: str = os.path.join(*HL_FIRST_ACCEPTOR)
# HL_LAST_DONOR: str = os.path.join(*HL_LAST_DONOR)
# HL_EQ_ACCEPTOR: str = os.path.join(TOGA2_ROOT, *HL_EQ_ACCEPTOR)
# HL_EQ_DONOR: str = os.path.join(TOGA2_ROOT, *HL_EQ_DONOR)
LIMIT_EXCEED_REJ: str = (
    'PROJECTION\t{}\t0\tNumber of homologous chains exceeds the set limit ({})\t'
    'CHAIN_LIMIT_EXCEEDED\tN'
)
MULTIPLE_ORTHOLOG_REJ: str = (
    'TRANSCRIPT\t{}\t0\tMultiple orthologs detected\tMULTIPLE_ORTHOLOGY\tN'
)
NO_CHAINS_REJ: str = (
    'TRANSCRIPT\t{}\t0\tNo covering chains detected\tNO_CHAINS\tN'
)
ZERO_ORTHOLOGY_REJ: str = (
    'TRANSCRIPT\t{}\t0\tNo orthologous chains detected\tZERO_ORTHOLOGY\tN'
)
OK: str = '.ok'
TOUCH: str = 'touch {}'

class PreprocessingScheduler(CommandLineManager):
    __slots__ = [
        'job_directory', 'preprocessing_directory', 'chain_map', 'ref_annotation',
        'ref', 'query', 'chain_file', 'ref_chrom_sizes', 'query_chrom_sizes',
        'joblist_file', 'segments',
        'fragmented_projections', 'jobnum', 'projections_per_command',
        'max_chain_number', 'orthologs_only', 'one2one_only', 'parallel_execution',
        'disable_spanning_chains', 'no_inference', 'memory_limit',
        'max_space_size', 'extrapolation_modifier', 'minimal_covered_fraction',
        'exon_locus_flank',
        'cesar_canon_u2_acceptor', 'cesar_canon_u2_donor',
        'cesar_non_canon_u2_acceptor', 'cesar_non_canon_u2_donor',
        'cesar_canon_u12_acceptor', 'cesar_canon_u12_donor',
        'cesar_non_canon_u12_acceptor', 'cesar_non_canon_u12_donor',
        'cesar_first_acceptor', 'cesar_last_donor', 'separate_site_treat',
        'assembly_gap_size', 'u12', 'spliceai_dir',
        'bigwig2wig_binary', 'min_splice_prob', 'annotate_ppgenes',
        'chain2trs', 'transcript_spans', 'job2cmds',
        'paralog_list', 'paralog_report', 'ppgene_list', 'ppgene_report',
        'rejected_transcripts', 'rejection_report',
        'toga1', 'toga1_plus_cesar',
        'v'
    ]

    def __init__(
        self,
        chain_map: click.Path,
        ref_annotation: click.Path,
        ref: click.Path,
        query: click.Path,
        chain_file: click.Path,
        ref_chrom_sizes: click.Path,
        query_chrom_sizes: click.Path,
        job_directory: click.Path,
        preprocessing_directory: click.Path,
        joblist_file: Optional[Union[click.Path, None]] = None,
        segments: Optional[Union[click.Path, None]] = None,
        fragmented_projections: Optional[Union[click.Path, None]] = None,
        job_number: Optional[int] = 300,
        projections_per_command: Optional[int] = 50,
        max_chain_number: Optional[int] = 100,
        orthologs_only: Optional[bool] = False,
        one2one_only: Optional[bool] = False,
        parallel_execution: Optional[bool] = False,
        disable_spanning_chains: Optional[bool] = False,
        no_inference: Optional[bool] = False,
        memory_limit: Optional[Union[float, None]] = None,
        max_space_size: Optional[int] = 500000,
        extrapolation_modifier: Optional[float] = 1.2,
        minimal_covered_fraction: Optional[float] = 0.0,
        exon_locus_flank: Optional[int] = FLANK_SPACE,
        cesar_canon_u2_acceptor: Optional[click.Path] = HG38_CANON_U2_ACCEPTOR,
        cesar_canon_u2_donor: Optional[click.Path] = HG38_CANON_U2_DONOR,
        cesar_non_canon_u2_acceptor: Optional[click.Path] = HG38_NON_CANON_U2_ACCEPTOR,
        cesar_non_canon_u2_donor: Optional[click.Path] = HG38_NON_CANON_U2_DONOR,
        cesar_canon_u12_acceptor: Optional[click.Path] = HG38_CANON_U12_ACCEPTOR,
        cesar_canon_u12_donor: Optional[click.Path] = HG38_CANON_U12_DONOR,
        cesar_non_canon_u12_acceptor: Optional[click.Path] = HG38_NON_CANON_U12_ACCEPTOR,
        cesar_non_canon_u12_donor: Optional[click.Path] = HG38_NON_CANON_U12_DONOR,
        cesar_first_acceptor: Optional[click.Path] = FIRST_ACCEPTOR,
        cesar_last_donor: Optional[click.Path] = LAST_DONOR,
        separate_splice_site_treatment: Optional[bool] = False,
        assembly_gap_size: Optional[int] = MIN_ASMBL_GAP_SIZE,
        u12: Optional[Union[click.Path, None]] = None,
        spliceai_dir: Optional[Union[click.Path, None]] = None,
        bigwig2wig_binary: Optional[Union[click.Path, None]] = None,
        min_splice_prob: Optional[bool] = 0.5,
        paralog_report: Optional[Union[click.File, None]] = None,
        annotate_processed_pseudogenes: Optional[bool] = False,
        processed_pseudogene_report: Optional[Union[click.File, None]] = None,
        rejection_report: Optional[Union[click.File, None]] = None,
        toga1_compatible: Optional[bool] = False,
        toga1_plus_corrected_cesar: Optional[bool] = False,
        log_name: Optional[Union[str, None]] = None,
        verbose: Optional[bool] = False
    ) -> None:
        """
        """
        self.v: bool = verbose
        self.set_logging(log_name)

        self.job_directory: click.Path = job_directory
        self.preprocessing_directory: click.Path = preprocessing_directory
        self.chain_map: click.Path = chain_map
        self.ref_annotation: click.Path = ref_annotation
        self.ref: click.Path = ref
        self.query: click.Path = query
        self.chain_file: click.Path = chain_file
        self.ref_chrom_sizes: click.Path = ref_chrom_sizes
        self.query_chrom_sizes: click.Path = query_chrom_sizes
        self.joblist_file: click.Path = (
            joblist_file if joblist_file is not None else
            os.path.join(self.job_directory, 'joblist_preprocess')
        )
        self.segments: Union[click.Path, None] = segments
        self.fragmented_projections: Union[click.Path, None] = fragmented_projections
        self.jobnum: int = job_number
        self.projections_per_command: bool = projections_per_command
        self.max_chain_number: int = max_chain_number
        self.orthologs_only: bool = orthologs_only
        self.one2one_only: bool = one2one_only
        self.parallel_execution: bool = parallel_execution
        self.disable_spanning_chains: bool = disable_spanning_chains
        self.no_inference: bool = no_inference
        self.memory_limit: bool = memory_limit
        self.max_space_size: int = max_space_size
        self.extrapolation_modifier: float = max(0.0, extrapolation_modifier)
        self.minimal_covered_fraction: float = max(0.0, minimal_covered_fraction)
        self.exon_locus_flank: int = exon_locus_flank
        self.cesar_canon_u2_acceptor: click.Path = cesar_canon_u2_acceptor
        self.cesar_canon_u2_donor: click.Path = cesar_canon_u2_donor
        self.cesar_non_canon_u2_acceptor: click.Path = cesar_non_canon_u2_acceptor
        self.cesar_non_canon_u2_donor: click.Path = cesar_non_canon_u2_donor
        self.cesar_canon_u12_acceptor: click.Path = cesar_canon_u12_acceptor
        self.cesar_canon_u12_donor: click.Path = cesar_canon_u12_donor
        self.cesar_non_canon_u12_acceptor: click.Path = cesar_non_canon_u12_acceptor
        self.cesar_non_canon_u12_donor: click.Path = cesar_non_canon_u12_donor
        self.cesar_first_acceptor: click.Path = cesar_first_acceptor
        self.cesar_last_donor: click.Path = cesar_last_donor
        self.separate_site_treat: bool = separate_splice_site_treatment

        self.assembly_gap_size: int = assembly_gap_size
        self.u12: Union[click.Path, None] = u12
        self.spliceai_dir: Union[click.Path, None] = spliceai_dir
        if bigwig2wig_binary is None:
            bw2w_in_path: str = which('bigWigToWig')
            if bw2w_in_path is not None:
                self.bigwig2wig_binary: click.Path = Path(bw2w_in_path).absolute()
            else:
                self.bigwig2wig_binary: click.Path = HL_BW2W_PATH
        else:
            self.bigwig2wig_binary: click.Path = bigwig2wig_binary
        self.min_splice_prob: float = max(0.0, min(min_splice_prob, 1.0))
        self.annotate_ppgenes: bool = annotate_processed_pseudogenes
        self.toga1: bool = toga1_compatible
        self.toga1_plus_cesar: bool = toga1_plus_corrected_cesar

        # self.tr2orth: Dict[str, List[str]] = defaultdict(list)
        # self.tr2par: Dict[str, List[str]] = defaultdict(list)
        # self.tr2trans: Dict[str, List[str]] = defaultdict(list)
        # self.tr2pgenes: Dict[str, List[str]] = defaultdict(list)
        self.chain2trs: Dict[str, List[str]] = defaultdict(list)
        self.transcript_spans: Dict[str, int] = {}
        self.job2cmds: List[int, Tuple[str, str]] = defaultdict(list)
        self.paralog_list: List[str] = []
        self.paralog_report: str = (
            paralog_report if paralog_report is not None else os.path.join(
                self.job_directory, 'paralogous_projections_to_align.txt'
            )
        )
        self.ppgene_list: List[str] = []
        self.ppgene_report: str = (
            processed_pseudogene_report if processed_pseudogene_report is not None else os.path.join(
                self.job_directory, 'processed_pseudogene_projections_to_align.txt'
            )
        )
        self.rejected_transcripts: List[Tuple[str]] = []
        self.rejection_report: str = (
            rejection_report if rejection_report is not None else os.path.join(
                self.job_directory, 'genes_rejection_reason.tsv'
            )
        )

        self.run()

    def run(self) -> None:
        """
        """
        ## create output directory
        self._mkdir(self.job_directory)
        ## upload the necessary data
        self.parse_mapper_file()
        if self.fragmented_projections:
            self.parse_fragment_file()
        self.infer_transcript_spans()
        ## split projections into jobwise bins
        self.split_commands_into_jobs()
        ## write job files and the job list
        self.write_job_files()
        ## if any transcripts were explicitly discarded, write a rejection report
        self.write_rejection_report()
        ## write the list of paralogs if any were processed
        self.write_paralog_report()
        ## presto!

    def _add_chain2trs(
        self, tr: str, chains: List[str], paralogs: bool = False, ppgenes: bool = False
    ) -> None:
        """
        For each chain in an iterable of chains, adds the transcript to a value
        in the self.chain2trs dictionary
        """
        if len(chains) > self.max_chain_number and not ppgenes:
            chains = sorted(chains, key=lambda x: int(x))
            relevant, dropped = chains[:self.max_chain_number], chains[self.max_chain_number:]
            self._to_log(
                f'Number of chains for transcript {tr} exceeds the set '
                f'chain number limit {self.max_chain_number}; '
                'dropping the excessive chains',
                'warning'
            )
        else:
            relevant, dropped = chains, []
        for chain in relevant:
            self.chain2trs[chain].append(tr)
            if paralogs:
                self.paralog_list.append(f'{tr}#{chain}')
            if ppgenes:
                self.ppgene_list.append(f'{tr}#{chain}')
        for chain in dropped:
            self.rejected_transcripts.append(
                LIMIT_EXCEED_REJ.format(f'{tr}#{chain}', self.max_chain_number)
            )

    def parse_mapper_file(self) -> None:
        """
        """
        with open(self.chain_map, 'r') as h:
            for line in h.readlines():
                data: List[str] = line.rstrip().split('\t')
                if len(data) != 5:
                    continue
                if data[0] == 'TRANSCRIPT':
                    continue
                tr: str = data[0]
                orth: List[str] = sorted(data[1].split(',')) if data[1] != '0' else []
                if self.one2one_only and (len(orth) > 1 or tr in self.tr2orth):
                    self.rejected_transcripts.append(
                        MULTIPLE_ORTHOLOG_REJ.format(tr) ## TODO: Clarify loss status with Michael
                    )
                    continue
                # self.tr2orth[tr].extend(orth)
                if self.orthologs_only:
                    if not orth:
                        self.rejected_transcripts.append(
                            ZERO_ORTHOLOGY_REJ.format(tr) ## Clarify loss status with Michael
                        )
                    continue
                par: List[str] = sorted(data[2].split(',')) if data[2] != '0' else []
                # self.tr2par[tr].extend(par)
                spanning: List[str] = (
                    sorted(data[3].split(',')) if data[3] != '0' #and not self.disable_spanning_chains 
                    else []
                )
                ppgenes: List[str] = (
                    sorted(data[4].split(',')) if data[4] != '0' and self.annotate_ppgenes else []
                )
                if not orth and not par and not spanning and not ppgenes:
                    self.rejected_transcripts.append(
                        NO_CHAINS_REJ.format(tr) ## Clarify loss status with Michael
                    )
                ## processed pseudogenes, if they are considered, 
                ## are added regardless of what other projections are present
                if ppgenes:
                    self._add_chain2trs(tr, ppgenes, ppgenes=True)
                if orth:
                    self._add_chain2trs(tr, orth)
                    continue
                if spanning:
                    # self.tr2spanning[tr].extend(spanning)
                    self._add_chain2trs(tr, spanning)
                    continue
                # else:
                #     self._add_chain2trs(tr, par)
                self._to_log(
                    f'No orthologs or spanning chains found for transcript {tr}; '
                    'processing paralogs instead',
                    'warning'
                )
                self._add_chain2trs(tr, par, paralogs=True)


    def parse_fragment_file(self) -> None:
        """
        """
        with open(self.fragmented_projections, 'r') as h:
            for line in h.readlines():
                data: List[str] = line.rstrip().split('\t')
                if not data or not data[0]:
                    continue
                if data[0] == 'transcript':
                    continue
                # self.tr2ort[data[0]].extend(data[1])
                tr: str = data[0]
                chains: str = data[1]
                self._add_chain2trs(tr, [chains])
                chains_split: List[str] = [x for x in chains.split(',') if x]
                for chain in chains_split:
                    self.chain2trs[chain].remove(tr)

    def infer_transcript_spans(self) -> None:
        """
        From the annotation file, selects the line corresponding to the focal
        transcripts and store its contents in the AnnotationEntry object
        """
        self._to_log('Uploading transcript data')
        with open(self.ref_annotation, 'r') as h:
            for line in h.readlines():
                data: List[str] = line.rstrip().split('\t')
                if len(data) < 8:
                    self._die('Reference annotation BED file contains insufficient data')
                transcript: str = data[3]
                cds_start: int = int(data[6])
                cds_stop: int = int(data[7])
                span: int = cds_stop - cds_start
                # self._echo(f'Coding region span for transcript {transcript}: {span} nt')
                self.transcript_spans[transcript] = span
            self._to_log(
                'Transcript span data for all reference transcripts '
                'have been successfully uploaded'
            )

    def split_commands_into_jobs(self) -> None:
        """
        An LPT-based job splitting procedure. Since preprocessing step speed
        is mostly defined by the number of chain blocks, reference transcript
        span lengths are used as proxies for memory requirements/processing time.
        At the first step, projections are aggregated chainwise to form input for
        preprocessing commands. If the number of transcripts for a given chain
        exceeds the set cap (self.proj_per_cmd attribute), projection list is
        then split into multiple commands, with memory requirements equalized via
        LPT according to cumulative transcript span.
        Then, resulting chain-transcript tuples are split into a given number of
        jobs controlled by self.job_number attribute via a similar LPT procedure.
        The resulting buckets are furthe rused to generate cesar_preprocess.py commands.
        """
        ## Step 1: Prepare input chain-transcripts pairs; split projections for
        ## highly syntenic chains
        input_tuples: List[Tuple[str, List[str]]] = []
        for chain, trs in self.chain2trs.items():
            if not trs:
                continue
            if len(trs) <= self.projections_per_command:
                input_tuples.append((chain, trs))
                continue
            cmd_num: int = ceil(len(trs) / self.projections_per_command)
            cmd_heap: List[Tuple[int, int]] = [(0, i) for i in range(cmd_num)]
            cmd2tr: List[int, List[str]] = {i: [] for i in range(cmd_num)}
            for tr in trs:
                total_span, cmd = heappop(cmd_heap)
                span: int = self.transcript_spans[tr]
                cmd2tr[cmd].append(tr)
                heappush(cmd_heap, (total_span + span, cmd))
            for cmd in cmd2tr:
                input_tuples.append((chain, cmd2tr[cmd]))
        ## Step 2: Split chain-transcript tuples into cluster jobs with LPT
        job_heap: List[Tuple[int, int]] = [(0, i) for i in range(self.jobnum)]
        for chain, trs in input_tuples:
            total_span, jobid = heappop(job_heap)
            self.job2cmds[jobid].append((chain, ','.join(trs)))
            span: int = sum(self.transcript_spans[tr] for tr in trs)
            heappush(job_heap, (total_span + span, jobid))

    def write_job_files(self) -> None:
        """
        Writes job files and a jo
        """
        # print(f'{self.ppgene_list=}')
        with open(self.joblist_file, 'w') as h1:
            for jobid, inputs in self.job2cmds.items():
                job_file: str = Path(
                    os.path.join(self.job_directory, f'batch{jobid}.ex')
                ).absolute()
                prepr_output: str = Path(
                    os.path.join(self.preprocessing_directory, f'batch{jobid}')
                ).absolute()
                with open(job_file, 'w') as h2:
                    h2.write('\n'.join(SPLIT_JOB_HEADER) + '\n')
                    h2.write(PRE_CLEANUP_LINE.format(prepr_output) + '\n')
                    for chain, trs in inputs:
                        cmd: str = (
                            f'{CESAR_PREPROCESS_SCRIPT} "{trs}" {chain} {self.ref_annotation} '
                            f'{self.ref} {self.query} {self.chain_file} '
                            f'{self.ref_chrom_sizes} {self.query_chrom_sizes} '
                            f' --cesar_canon_u2_acceptor {self.cesar_canon_u2_acceptor}'
                            f' --cesar_canon_u2_donor {self.cesar_canon_u2_donor}'
                            f' --cesar_non_canon_u2_acceptor {self.cesar_non_canon_u2_acceptor}'
                            f' --cesar_non_canon_u2_donor {self.cesar_non_canon_u2_donor}'
                            f' --cesar_canon_u12_acceptor {self.cesar_canon_u12_acceptor}'
                            f' --cesar_canon_u12_donor {self.cesar_canon_u12_donor}'
                            f' --cesar_non_canon_u12_acceptor {self.cesar_non_canon_u12_acceptor}'
                            f' --cesar_non_canon_u12_donor {self.cesar_non_canon_u12_donor}'
                            f' --cesar_first_acceptor {self.cesar_first_acceptor}'
                            f' --cesar_last_donor {self.cesar_last_donor}'
                            f' --assembly_gap_size {self.assembly_gap_size}'
                            f' --minimal_covered_fraction {self.minimal_covered_fraction}'
                            f' --max_space_size {self.max_space_size}'
                            f' --extrapolation_modifier {self.extrapolation_modifier}'
                            f' --exon_locus_flank {self.exon_locus_flank}'
                        )
                        if self.toga1 and not self.toga1_plus_cesar:
                            cmd += ' -t1 '
                        if self.toga1_plus_cesar:
                            cmd += ' -t1c '
                        if self.segments is not None:
                            cmd += f'  --segments {self.segments}'
                        if self.separate_site_treat:
                            cmd += ' --separate_splice_site_treatment '
                        if self.memory_limit is not None:
                            cmd += f' --memory_limit {self.memory_limit}'
                        if self.parallel_execution:
                            cmd += ' --parallel_job'
                        if self.disable_spanning_chains:
                            cmd += ' --disable_spanning_chains'
                        if self.no_inference:
                            cmd += ' --no_inference'
                        if self.u12 is not None:
                            cmd += f' --u12 {self.u12}'
                        if self.annotate_ppgenes:
                            # print(f'{trs=}')
                            ppgenes_in_batch: List[str] = [x for x in trs.split(',') if f'{x}#{chain}' in self.ppgene_list]
                            if ppgenes_in_batch:
                                ppgene_str: str = ','.join(ppgenes_in_batch)
                                cmd += f' --processed_pseudogene_list {ppgene_str}'
                        if self.spliceai_dir is not None:
                            cmd += f' --spliceai_dir {self.spliceai_dir}'
                            cmd += f' --bigwig2wig_binary {self.bigwig2wig_binary}'
                            cmd += f' --min_splice_prob {self.min_splice_prob}'
                        cmd += f' --output {prepr_output}'# -v'
                        h2.write(cmd + '\n')
                        ok_file: str = os.path.join(prepr_output, OK)
                    h2.write(TOUCH.format(ok_file) + '\n')
                file_mode: bytes = os.stat(job_file).st_mode
                file_mode |= (file_mode & 0o444) >> 2
                os.chmod(job_file, file_mode)
                h1.write(f'{job_file}\n')

    def write_paralog_report(self) -> None:
        """
        Writes paralogous and processed pseudogene projections permitted by the scheduler 
        to the respective single-column text files
        """
        if self.paralog_list:
            with open(self.paralog_report, 'w') as h:
                for par in self.paralog_list:
                    h.write(par + '\n')
        if self.ppgene_list:
            with open(self.ppgene_report, 'w') as h:
                for ppgene in self.ppgene_list:
                    h.write(ppgene + '\n')

    def write_rejection_report(self) -> None:
        """
        Records transcripts which were discarded prior to preprocessing step
        """
        if not self.rejected_transcripts:
            return
        with open(self.rejection_report, 'w') as h:
            for line in self.rejected_transcripts:
                h.write(line + '\n')

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'chain_map',
    type=click.Path(exists=True),
    metavar='TRANSCRIPT_CHAIN_MAP'
)
@click.argument(
    'ref_annotation',
    type=click.Path(exists=True),
    metavar='REF_ANNOT_BED'
)
@click.argument(
    'ref',
    type=click.Path(exists=True),
    metavar='REF_GENOME'
)
@click.argument(
    'query',
    type=click.Path(exists=True),
    metavar='QUERY_GENOME'
)
@click.argument(
    'chain_file',
    type=click.Path(exists=True),
    metavar='CHAIN_FILE'
)
@click.argument(
    'ref_chrom_sizes',
    type=click.Path(exists=True),
    metavar='REF_CHROM_SIZE_FILE'
)
@click.argument(
    'query_chrom_sizes',
    type=click.Path(exists=True),
    metavar='QUERY_CHROM_SIZE_FILE'
)
@click.argument(
    'job_directory',
    type=click.Path(exists=False),
    metavar='PREPROCESS_JOB_DIRECTORY'
)
@click.argument(
    'preprocessing_directory',
    type=click.Path(exists=False),
    metavar='PREPROCESS_OUTPUT_DIRECTORY'
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
        '[default: PREPROCESSING_JOB_DIRECTORY/cesar_joblist]. '
    )
)
@click.option(
    '--segments',
    '-s',
    type=click.Path(exists=True),
    metavar='BED_FILE',
    default=None,
    show_default=True,
    help='A BED12 file containing segments with external evidence'
)
@click.option(
    '--fragmented_projections',
    '-f',
    type=click.Path(exists=True),
    metavar='TSV',
    default=None,
    show_default=True,
    help='A two-column file containing fragmented (multi-chain) projections. '
         'By default, TOGA saves these data to tmp/gene_fragments.tsv'
)

@click.option(
    '--job_number',
    '-j',
    type=int,
    metavar='INT',
    default=300,
    show_default=True,
    help='A number of cluster jobs to split the projection list into'
)
@click.option(
    '--projections_per_command',
    '-p',
    type=int,
    metavar='INT',
    default=50,
    show_default=True,
    help='A maximum number of projections per each preprocessing command. '
         'Chains projecting the number of transcripts exceeding this number '
         'will be split into multiple commands with balanced memory requirements. '
)
@click.option(
    '--max_chain_number',
    '-mc',
    type=int,
    default=100,
    show_default=True,
    help=(
        'Maximum number of chains per transcript; entries with '
        'number of homologous (usually paralogous) chains exceeding this limit '
        'will be discarded'
    )
)
@click.option(
    '--orthologs_only',
    '-r',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, only orthologous projections are considered'
)
@click.option(
    '--one2one_only',
    '-o2o',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, only transcript with a single orthologous projection '
         'are considered'
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
    '--disable_spanning_chains',
    '-nospan',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=False,
    help='If set, ignores spanning chains '
         '(i.e. chains with no coding blocks corresponding to transcript exons)'
)
@click.option(
    '--no_inference',
    '-nf',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, disables extrapolating missing exons\' location. Temporary feature'
)
@click.option(
    '--memory_limit',
    '-ml',
    type=float,
    metavar='FLOAT',
    default=None,
    show_default=True,
    help='Upper memory limit for CESAR jobs. If limit is exceeded, the program '
         'terminates with zero exit status'
)
@click.option(
    '--max_space_size',
    '-mss',
    type=int,
    metavar='INT',
    default=500000,
    show_default=True,
    help='Maximum search space size used for locus shrinking for missing exons, bps'
)
@click.option(
    '--extrapolation_modifier',
    '-em',
    type=float,
    metavar='FLOAT',
    default=1.2,
    show_default=True,
    help='Multiply extrapolated extension by this value'
)
@click.option(
    '--minimal_covered_fraction',
    '-mincov',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.0,
    help='Minimal fraction of reference CDS to be covered by alignment data. '
         'Projections covering less than this portion will be discarded.'
)
@click.option(
    '--exon_locus_flank',
    '-ef',
    type=int,
    metavar='INT',
    default=FLANK_SPACE,
    show_default=True,
    help='Flank size to extend the estimated exon loci by'
)
@click.option(
    '--cesar_canon_u2_acceptor',
    '-cca',
    type=str,
    metavar='REL_PATH',
    default=HG38_CANON_U2_ACCEPTOR,
    show_default=True,
    help='A path to canonical (GT/GC-AG) U2 acceptor profile'
)
@click.option(
    '--cesar_canon_u2_donor',
    '-ccd',
    type=str,
    metavar='REL_PATH',
    default=HG38_CANON_U2_DONOR,
    show_default=True,
    help='A path to canonical (GT/GC-AG) U2 donor profile'
)
@click.option(
    '--cesar_non_canon_u2_acceptor',
    '-cnca',
    type=str,
    metavar='REL_PATH',
    default=HG38_NON_CANON_U2_ACCEPTOR,
    show_default=True,
    help='A path to non-canonical (non GT/GC-AG) U2 acceptor profile'
)
@click.option(
    '--cesar_non_canon_u2_donor',
    '-cncd',
    type=str,
    metavar='REL_PATH',
    default=HG38_NON_CANON_U2_DONOR,
    show_default=True,
    help='A path to non-canonical (non GT/GC-AG) U2 donor profile'
)
@click.option(
    '--cesar_canon_u12_acceptor',
    '-cua',
    type=str,
    metavar='REL_PATH',
    default=HG38_CANON_U12_ACCEPTOR,
    show_default=True,
    help='A path to canonical (GT-AG) U12 exon acceptor profile'
)
@click.option(
    '--cesar_canon_u12_donor',
    '-cud',
    type=str,
    metavar='REL_PATH',
    default=HG38_CANON_U12_DONOR,
    show_default=True,
    help='A path to canonical (GT-AG) U12  donor profile'
)
@click.option(
    '--cesar_non_canon_u12_acceptor',
    '-cnua',
    type=str,
    metavar='REL_PATH',
    default=HG38_NON_CANON_U12_ACCEPTOR,
    show_default=True,
    help='A path to non-canonical (non-GT-AG) U12 exon acceptor profile'
)
@click.option(
    '--cesar_non_canon_u12_donor',
    '-cnud',
    type=str,
    metavar='REL_PATH',
    default=HG38_NON_CANON_U12_DONOR,
    show_default=True,
    help=(
        'A path to non-canonical (non-GT-AG) U12 exon donor profile'
    )
)
@click.option(
    '--cesar_first_acceptor',
    '-cfa',
    type=str,
    metavar='REL_PATH',
    default=FIRST_ACCEPTOR,
    show_default=True,
    help='A (relative to CESAR2 location) path to first exon acceptor profile'
)
@click.option(
    '--cesar_last_donor',
    '-cld',
    type=str,
    metavar='REL_PATH',
    default=LAST_DONOR,
    show_default=True,
    help='A (relative to CESAR2 location) path to last exon donor profile'
)
@click.option(
    '--separate_splice_site_treatment',
    '-ssst',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, donor and acceptor intron splice sites are treated '
        'as (non-)canonical indepent of each other'
    )
)
@click.option(
    '--assembly_gap_size',
    '-gs',
    type=int,
    metavar='INT',
    default=MIN_ASMBL_GAP_SIZE,
    show_default=True,
    help='Minimum number of consecutive N symbols to be considered an assembly gap'
)
@click.option(
    '--u12',
    '-u12',
    type=click.Path(exists=True),
    metavar='U12_FILE',
    default=None,
    show_default=True,
    help=(
        'A three-column tab-separated file containing information on the '
        'non-canonical splice sites'
    )
)
@click.option(
    '--spliceai_dir',
    '-sai',
    type=click.Path(exists=True),
    metavar='SPLICEAI_OUT_DIR',
    help='A path to the SpliceAI pipeline output directory'
)
@click.option(
    '--bigwig2wig_binary',
    '-bw2w',
    type=click.Path(exists=True),
    metavar='BIGWIG2WIG_BINARY',
    default=None,
    help=(
        'A path to the UCSC bigWigToWig binary; if none is provided or found '
        'in PATH, default is set to Hiller Lab Delta cluster utility location'
    )
)
@click.option(
    '--min_splice_prob',
    '-msp',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.5,
    show_default=True,
    help='Minimum SpliceAI prediction probability to consider the splice site'
)
@click.option(
    '--paralog_report',
    '-pr',
    type=click.Path(exists=False),
    metavar='PATH',
    default=None,
    show_default=False,
    help=(
        'A path to save analysed paralogous projections to '
        '[default:PREPROCESS_JOB_DIRECTORY/paralogous_projections_to_align.txt]'
    )
)
@click.option(
    '--annotate_processed_pseudogenes',
    '-pp',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, spanning chains (i.e., chains with alignment gap corresponding '
        'to the projected transcript) are not considered for CESAR alignment'
    )
)
@click.option(
    '--processed_pseudogene_report',
    '-ppr',
    type=click.Path(exists=False),
    metavar='PATH',
    default=None,
    show_default=False,
    help=(
        'A path to save analysed processed pseudogene projections to '
        '[default:PREPROCESS_JOB_DIRECTORY/processed_pseudogene_projections_to_align.txt]'
    )
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
    '--toga1_plus_corrected_cesar',
    '-t1c',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'Alignment procedure is fully TOGA1.0-compliant except for exonwise '
        'CESAR alignment and corrected CESAR-related bugs; '
        'benchmarking feature, do not use in real runs'
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

def main(**kwargs) -> None:
    PreprocessingScheduler(**kwargs)

if __name__ == '__main__':
    main()
    # PreprocessingScheduler()
