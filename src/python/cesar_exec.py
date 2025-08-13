#!/usr/bin/env python3

"""
Actual CESAR wrapper for TOGA 1.5+
"""

from collections import defaultdict
from datetime import datetime
from filelock import FileLock # remove if not needed
# from GLP_values import * # TODO: make sure that paths to local packages are provided correctly
from math import ceil
from modules.cesar_wrapper_constants import *
from modules.cesar_wrapper_executables import *
from modules.constants import Constants
from modules.intron_gain_check import *
from modules.preprocessing import cesar_memory_check
from modules.processed_segment import *
from modules.shared import (
    CommandLineManager, CONTEXT_SETTINGS, dir_name_by_date, get_upper_dir, hex_code
)
# from twobitreader import TwoBitFile
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

import click
import h5py
import logging
# import numpy as np
import os

# logging.basicConfig(level=logging.INFO)

__author__  = 'Yury V. Malovichko'
__credits__ = ['Bogdan Kirilenko', 'Michael Hiller']
__year__ = '2024'



## create or update the needed constants
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
LOCATION: str = get_upper_dir(__file__, 3)
BLOSUM_FILE: str = os.path.join(LOCATION, *DEF_BLOSUM_FILE)
UNALIGNED_REJ: str = (
    'PROJECTION\t{}\t{}\tNo exons were aligned\tZERO_ALIGNED\t{}'
)

HG38_CANON_U2_ACCEPTOR: str = os.path.abspath(os.path.join(LOCATION, *HG38_CANON_U2_ACCEPTOR))
FIRST_ACCEPTOR: str = os.path.abspath(os.path.join(LOCATION, *FIRST_ACCEPTOR))
HG38_CANON_U2_DONOR: str = os.path.abspath(os.path.join(LOCATION, *HG38_CANON_U2_DONOR))
LAST_DONOR: str = os.path.abspath(os.path.join(LOCATION, *LAST_DONOR))

# CESAR_FIRST_ACC = '/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/scripts/TOGA2.0/CESAR2.0/extra/tables/human/firstCodon_profile.txt'
# CESAR_REG_ACC = '/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/scripts/TOGA2.0/CESAR2.0/extra/tables/human/acc_profile.txt'
# CESAR_LAST_DONOR = '/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/scripts/TOGA2.0/CESAR2.0/extra/tables/human/lastCodon_profile.txt'
# CESAR_REG_DONOR = '/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/scripts/TOGA2.0/CESAR2.0/extra/tables/human/do_profile.txt'

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'transcript',
    type=str,
    metavar='TRANSCRIPT'
)
@click.argument(
    'chain',
    type=str,
    metavar='CHAIN_ID'
)
@click.argument(
    'data',
    type=click.Path(exists=True),
    metavar='INPUT_FILE'
)
@click.option(
    '--parallel_job',
    '-p',
    is_flag=True,
    default=False,
    show_default=True,
    help='Indicates that the process is a part of parallel job batch; '
         'output file names are standardized and provided with lock.files'
)
@click.option(
    '--cesar_binary',
    '-cs',
    type=click.Path(exists=True),
    metavar='CESAR_BINARY',
    default=HL_CESAR_PATH,
    show_default=False,
    help='A path to the actual CESAR2.0 binary; default is set for Hiller '
         'lab Delta cluster'
)
@click.option(
    '--matrix',
    '-m',
    type=click.File(lazy=False),
    metavar='BLOSUM_MATRIX_FILE',
    default=BLOSUM_FILE,
    show_default=True,
    help='A file containing the protein alignment matrix'
)
@click.option(
    '--mask_terminal_mutations',
    '-m10m',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, masks mutations occurring in first and '
         'last 10 percents of query length'
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
    '--rescue_missing_start',
    '-rma',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, scans the upstream of the search space for inframe start codons '
         'in case the alignmentn does not start with one'
)
@click.option(
    '--rescue_missing_stop',
    '-rmo',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, scans the downsteam of the search space for inframe stop codons '
         'in case the alignmentn does not end with one'
)
@click.option(
    '--paralogous_projection',
    '-pg',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, specifies that the analysed projection is a paralogous one'
)
@click.option(
    '--processed_pseudogene_projection',
    '-pp',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, specifies that the analysed projection is a processed_pseudogene'
)
@click.option(
    '--filtered_bed_output',
    '-fbo',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, does not report missing and/or deleted exons '
         'in the output BED12 file'
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
    default=0.5,
    show_default=True,
    help='Minimal intron gain threshold to consider'
)
# @click.option(
#     '--min_intron_prob_gapped',
#     '-mipg',
#     type=click.FloatRange(min=0.0, max=1.0),
#     metavar='FLOAT',
#     default=0.1,
#     show_default=True,
#     help=(
#         'Minimal SpliceAI support for both sites in the query-specific introns '
#         'if intron location is supported by alignment gaps.'
#     )
# )
# @click.option(
#     '--min_intron_prob_ungapped',
#     '-mipu',
#     type=click.FloatRange(min=0.0, max=1.0),
#     metavar='FLOAT',
#     default=0.4,
#     show_default=True,
#     help=(
#         'Minimal SpliceAI support for both sites in the query-specific introns '
#         'if intron location is not supported by alignment gaps. WARNING: If a value provided is less '
#         'than --min_intron_prob_gapped, it is set equal to --min_intron_prob_gapped instead'
#     )
# )
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
    '--max_intron_number',
    '-mnn',
    type=click.IntRange(min=1, max=None),
    metavar='INT',
    default=4, 
    show_default=True,
    help=(
        'Maximum gained intron number per exon. Highly recommended to increase this beyond 5-6'
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
    help=(
        'If set, SpliceAI data will be used exclusively for restoring missing '
        'exons, with no post-CESAR correction'
    )
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
)
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
    '--output',
    '-o',
    type=click.Path(exists=False),
    metavar='OUT_DIR',
    default=dir_name_by_date('segment_reconstruction_'),
    show_default=False,
    help=(
        'A directory to write the output to '
        '[default: segment_reconstruction_{run date}'   
    )
)
@click.option(
    '--tmp',
    '-tmp',
    type=click.Path(exists=False),
    metavar='TMP_DIR',
    default=None,
    show_default=False,
    help='A directory to write the CESAR input to [default: OUT_DIR/tmp]'
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

class CesarExecutor(CommandLineManager):
    """
    The actual CESAR wrapper for TOGA 1.5+ pipeline. Given the transcript-chain
    pair and a path to the HDF5 file storing the CESAR input, runs CESAR,
    processes the results, and performs mutation check.\n
    NOTE: As a part of the TOGA 1.5+ pipeline, this script operates on the
    cesar_preprocess.py output. If you want to run CESAR step manually for a given
    transcript-chain pair, you most likely need the _CESAR_wrapper_toga2.0.py
    script.
    """
    __slots__ = [
        'transcript', 'chain', 'projection_name', 'data', 'p_job', 'cesar_binary', 'aa_matrix',
        'mask_terminal_mutations', 'assembly_gap_size', 'rescue_missing_start',
        'rescue_missing_stop', 'is_paralog', 'is_processed_pseudogene', 'report_raw_bed',
        'min_splice_prob', 'splice_prob_margin', 'intron_gain_check',
        'min_intron_gain_score', 'max_intron_number', 
        'min_intron_prob_gapped', 'min_intron_prob_ungapped',
        'min_intron_prob_trusted', 'min_intron_prob_supported', 'min_intron_prob_unsupported',
        'regular_acceptor', 'regular_donor', 'first_acceptor', 'last_donor',
        'no_sai_correction', 'sai_correction_mode',
        'correct_short_introns', 'ignore_alternative_frame', 
        'save_cesar_input', 'output', 'tmp', 'v', 'max_exon_num',
        'chains', 'fragmented', 'u12_sites', 'raw_cesar_input', 'cesar_results',
        'rejected_projections', 'acceptor_flanks', 'donor_flanks', 'cesar_err',
        'bed_path', 'bed_lock', 'filt_bed_path', 'filt_bed_lock',
        'browser_bed_lock', 'browser_bed_file',
        'id_stub', 'id_lock', 'exon_meta_stub', 'exon_meta_lock',
        'cesar_res_stub', 'cesar_res_lock', 'mutation_stub', 'mutation_lock',
        'rejection_file', 'rejection_lock', 'codon_stub', 'codon_lock',
        'prot_stub', 'prot_lock', 'exon_fa_stub', 'exon_fa_lock',
        'splice_site_stub', 'splice_site_lock', 'log_file',
        'intron_evidence_stub', 'intron_evidence_lock', 'intron2evidence',
        'splice_site_shifts', 'splice_shift_lock',
        'cds_fasta', 'cds_fasta_lock',
        'selenocysteine_codons', 'selenocysteine_lock',
        'max_mem',
        'toga1'
    ]

    def __init__(
        self,
        transcript: str,
        chain: str,
        data: click.Path,
        parallel_job: Optional[bool],
        cesar_binary: Optional[click.Path],
        matrix: Optional[click.File],
        mask_terminal_mutations: Optional[bool],
        assembly_gap_size: Optional[int],
        rescue_missing_start: Optional[bool],
        rescue_missing_stop: Optional[bool],
        paralogous_projection: Optional[bool],
        processed_pseudogene_projection: Optional[bool],
        filtered_bed_output: Optional[bool],
        spliceai_correction_mode: Optional[int],
        min_splice_prob: Optional[float],
        splice_prob_margin: Optional[float],
        intron_gain_check: Optional[bool],
        intron_gain_threshold: Optional[bool],
        # min_intron_prob_gapped: Optional[float],
        # min_intron_prob_ungapped: Optional[float],
        min_intron_prob_trusted: Optional[bool],
        min_intron_prob_supported: Optional[bool],
        min_intron_prob_unsupported: Optional[bool],
        max_intron_number: Optional[int],
        cesar_regular_acceptor: Optional[click.Path],
        cesar_regular_donor: Optional[click.Path],
        cesar_first_acceptor: Optional[click.Path],
        cesar_last_donor: Optional[click.Path],
        no_spliceai_correction: Optional[bool],
        correct_ultrashort_introns: Optional[bool],
        ignore_alternative_frame: Optional[bool],
        save_cesar_input: Optional[bool],
        output: Optional[click.Path],
        tmp: Optional[click.Path],
        toga1_compatible: Optional[bool],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.output: str = output
        if not os.path.isdir(self.output):
            self._echo('Creating output directory')
            self._mkdir(self.output)
        self.log_file: str = os.path.join(self.output, 'log.txt')
        self.set_logging(f'{__name__}_{hex_code()}')
        
        self._to_log('Initializing the CESAR wrapper')
        self.transcript: str = transcript
        self._to_log(f'Input transcript is {self.transcript}')
        self.chain: str = chain
        self.projection_name: str = f'{self.transcript}#{self.chain}'
        self.chains: List[str] = self.chain.split(',')
        self.fragmented: bool = len(self.chains) > 1
        self._to_log(
            f'Input chain is {chain}; expecting complete transcript'
            if len(self.chains) == 1 else
            f'Multiple chains were provided ({chain}); expecting fragmented transcript'
        )

        self.u12_sites: Dict[int, Set[int]] = {}

        self.data: click.Path = data
        self.p_job: bool = parallel_job
        self.cesar_binary: str = cesar_binary
        self.mask_terminal_mutations: bool = mask_terminal_mutations
        self.assembly_gap_size: int = assembly_gap_size
        self.rescue_missing_start: bool = rescue_missing_start
        self.rescue_missing_stop: bool = rescue_missing_stop
        self.is_paralog: bool = paralogous_projection
        self.is_processed_pseudogene: bool = processed_pseudogene_projection
        self.report_raw_bed: bool = not filtered_bed_output
        self.no_sai_correction: bool = no_spliceai_correction
        self.sai_correction_mode: int = spliceai_correction_mode
        self.min_splice_prob: float = min_splice_prob
        self.splice_prob_margin: float = splice_prob_margin
        self.intron_gain_check: bool = intron_gain_check
        self.min_intron_gain_score: float = intron_gain_threshold
        # self.min_intron_prob_gapped: float = min_intron_prob_gapped
        # self.min_intron_prob_ungapped: float = min_intron_prob_ungapped
        self.min_intron_prob_trusted: bool = min_intron_prob_trusted
        self.min_intron_prob_supported: bool = min_intron_prob_supported
        self.min_intron_prob_unsupported: bool = min_intron_prob_unsupported
        self.max_intron_number: int = max_intron_number
        self.regular_acceptor: click.Path = os.path.abspath(cesar_regular_acceptor)
        self.regular_donor: click.Path = os.path.abspath(cesar_regular_donor)
        self.first_acceptor: click.Path = os.path.abspath(cesar_first_acceptor)
        self.last_donor: click.Path = os.path.abspath(cesar_last_donor)
        self.correct_short_introns: bool = correct_ultrashort_introns
        self.ignore_alternative_frame: bool = ignore_alternative_frame
        self.save_cesar_input: bool = save_cesar_input
        self.tmp: str = os.path.join(output, 'tmp') if not tmp else tmp
        self.toga1: bool = toga1_compatible

        ## parse the protein matrix
        self._to_log('Parsing the protein score matrix')
        self.aa_matrix: Dict[str, Dict[str, int]] = make_matrix(matrix)

        self.max_exon_num: int = 1
        self.cesar_results: Dict[int, Dict[int, List[CesarExonEntry]]] = defaultdict(
            lambda: defaultdict(list)
        )
        self.rejected_projections: List[str] = []
        self.acceptor_flanks: Dict[str, int] = {}
        self.donor_flanks: Dict[str, int] = {}
        self.intron2evidence: Dict[int, List[IntronMeta]] = {}
        self.max_mem: int = 0

        self.cesar_err: str = 'CESAR died with the following error'

        self.bed_path: str = os.path.join(self.output, 'query_annotation.with_discarded_exons.bed')
        self.bed_lock: str = self.bed_path + '.lock'
        self.filt_bed_path: str = os.path.join(self.output, 'query_annotation.bed')
        self.filt_bed_lock: str = self.filt_bed_path + '.lock'
        self.browser_bed_file: str = os.path.join(self.output, 'query_annotation.for_browser.bed')
        self.browser_bed_lock: str = self.browser_bed_file + '.lock'
        self.rejection_file: str = os.path.join(self.output, 'genes_rejection_reason.tsv')
        self.rejection_lock: str = self.rejection_file + '.lock'
        # if self.p_job:
        self.id_stub: str = os.path.join(self.output, 'transcript_meta.tsv')
        self.id_lock: str = self.id_stub + '.lock'
        self.exon_meta_stub: str = os.path.join(self.output, 'exon_meta.tsv')
        self.exon_meta_lock: str = self.exon_meta_stub + '.lock'
        self.cesar_res_stub: str = os.path.join(self.output, 'cesar_results.cesout')
        self.cesar_res_lock: str = self.cesar_res_stub + '.lock'
        self.mutation_stub: str = os.path.join(self.output, 'inactivating_mutations.tsv')
        self.mutation_lock: str = self.mutation_stub + '.lock'
        # self.memory_lock: str = self.memory_file + '.lock'
        self.codon_stub: str = os.path.join(self.output, 'codon_aln.fa')
        self.codon_lock: str = self.codon_stub + '.lock'
        self.prot_stub: str = os.path.join(self.output, 'protein_aln.fa')
        self.prot_lock: str = self.prot_stub + '.lock'
        self.exon_fa_stub: str = os.path.join(self.output, 'exon_aln.fa')
        self.exon_fa_lock: str = self.exon_fa_stub + '.lock'
        self.splice_site_stub: str = os.path.join(self.output, 'splice_sites.tsv')
        self.splice_site_lock: str = self.splice_site_stub + '.lock'
        self.intron_evidence_stub: str = os.path.join(self.output, 'gained_intron_summary.tsv')
        self.intron_evidence_lock: str = self.intron_evidence_stub + '.lock'
        self.splice_site_shifts: str = os.path.join(self.output, 'splice_site_shifts.tsv')
        self.splice_shift_lock: str = self.splice_site_shifts + '.lock'
        self.cds_fasta: str = os.path.join(self.output, 'nucleotide.fa')
        self.cds_fasta_lock: str = os.path.join(self.cds_fasta + '.lock')
        self.selenocysteine_codons: str = os.path.join(self.output, 'selenocysteine_codons.tsv')
        self.selenocysteine_lock: str = self.selenocysteine_codons + '.lock'
        # else:
        #     self.id_stub: str = os.path.join(self.output, 'transcript_meta_')
        #     self.exon_meta_stub: str = os.path.join(self.output, 'exon_meta_')
        #     self.cesar_res_stub: str = os.path.join(self.output, 'cesar_results_')
        #     self.mutation_stub: str = os.path.join(self.output, 'mutation_check_')
        #     self.codon_stub: str = os.path.join(self.output, 'codon_aln_')
        #     self.prot_stub: str = os.path.join(self.output, 'protein_aln_')
        #     self.exon_fa_stub: str = os.path.join(self.output, 'exon_aln_')
        #     self.splice_site_stub: str = os.path.join(self.output, 'splice_sites_')

        self.run()

    # def _to_log(self, msg: str) -> None:
    #     """Report a line to standard output if verbosity is enabled"""
    #     click.echo(msg) if self.v else None
    #
    # def _die(self, msg: str) -> None:
    #     """Error-exit with a given message"""
    #     click.echo(msg)
    #     exit(1)
    #
    # def _mkdir(self, d: str) -> None:
    #     """Safe directory creation function"""
    #     try:
    #         os.makedirs(d)
    #     except FileExistsError:
    #         pass
    #
    # def _exec(self, cmd: str, err_msg: str, input_: bytes = None) -> None:
    #     """Runs subprocesses, handles the exceptions"""
    #     pr: subprocess.Pipe = subprocess.Popen(
    #         cmd, shell=True,
    #         stdin=subprocess.PIPE,
    #         stdout=subprocess.PIPE,
    #         stderr=subprocess.PIPE
    #     )
    #     stdout, stderr = pr.communicate(input=input_)
    #     rc: int = pr.returncode
    #     if rc != 0:
    #         print(stderr)
    #         msg: str = f'{err_msg}:\n{stderr.decode("utf8")}'
    #         self._die(msg)
    #     return stdout.decode("utf8")

    # def set_logging(self, name: str = __name__) -> None:
    #     """
    #     Sets up logging system for a TogaMain instance
    #     """
    #     print('Setting logger')
    #     self.logger: logging.Logger = logging.getLogger(name)
    #     file_handler: logging.FileHandler = logging.FileHandler(
    #         self.log_file, mode='a', encoding=Constants.UTF8
    #     )
    #     file_handler.setFormatter(Constants.FORMATTER)
    #     self.logger.addHandler(file_handler)
    #     if self.v:
    #         console_handler: logging.StreamHandler = logging.StreamHandler()
    #         console_handler.setFormatter(Constants.FORMATTER)
    #         self.logger.addHandler(console_handler)
    #     self.logger.propagate = False
    #     print('Logger is set')

    def run(self) -> None:
        ## Prepare the output structures
        self._to_log('Starting CESAR aligning module')
        self._to_log('Checking if output directory exists')
        self._to_log('Checking if temporary directory exists')
        if not os.path.isdir(self.tmp):# and self.save_cesar_input:
            self._to_log(f'Creating temporary directory at {self.tmp}')
            self._mkdir(self.tmp)
        self._to_log('Output and temporary directory were successfully created')

        ## read the input data
        try:
            key: str = f'{self.transcript}|{self.chain}'
            with h5py.File(self.data, 'r') as f:
                dataset: Iterable[Iterable[str]] = f[key][:]
                self.raw_cesar_input: Dict[int, Dict[int, CesarInput]] = a_table_to_cesar_group(
                    dataset,
                )
        except KeyError:
            self._die(
                'Transcript-chain pair was not found in the source HDF5 file'
            )

        groups_finished: int = 0
        segments_finished: int = 0
        for tr in self.raw_cesar_input:
            self.acceptor_flanks[tr] = None
            self.donor_flanks[tr] = None
            for g, group in self.raw_cesar_input[tr].items():
                self.max_exon_num = max(self.max_exon_num, max(group.exons))
                bed_group: List[CesarExonEntry] = []
                # input_: str = group.input_str
                # exon_seqs: List[str] = [
                #     x for x in input_.split('####')[0].split('\n') if x and x[0] != '>'
                # ]
                exon_nums: List[int] = group.exons
                self.u12_sites = {**self.u12_sites, **group.u12_data}
                if self.acceptor_flanks[tr] is None:
                    self.acceptor_flanks[tr] = group.acceptor_flank
                elif self.donor_flanks[tr] != group.acceptor_flank:
                    self._die(
                        'Inconsistent acceptor flank sizes across the groups '
                        f'for segment {tr}'
                    )
                if self.donor_flanks[tr] is None:
                    self.donor_flanks[tr] = group.donor_flank
                elif self.donor_flanks[tr] != group.donor_flank:
                    self._die(
                        'Inconsistent donor flank sizes across the groups '
                        f'for segment {tr}'
                    )
                if not group.will_be_aligned:
                    self._to_log(
                        f'Group {g} for projection {self.projection_name} was excluded from alignment'
                    )
                    result_dummy: RawCesarOutput = RawCesarOutput(
                        group.chain,
                        group.chrom,
                        group.start,
                        group.stop,
                        group.strand,
                        ' ' + ' '.join(group.exon_seqs) + ' ',
                        ' ' + group.query_seq.replace('|', ' ') + ' ',
                        group.exons,
                        {x: (None, None) for x in group.exons}, # group.expected_coords,
                        {x: (None, None) for x in group.exons}, # group.search_space_coords,
                        group.spliceai_data,
                        group.gap_located_exons,
                        group.out_of_chain_exons,
                        not group.will_be_aligned,
                        group.intersects_asmbl_gaps,
                        {} ## for intron coordinates
                    )
                    self.cesar_results[tr][g].append(result_dummy)
                    continue

                ## format CESAR input
                exon_nums: str = (
                    ','.join(map(str, group.exons)) if len(group.exons) <= 10
                    else f'{group.exons[0]}-{group.exons[-1]}'
                )
                cesar_input_file: str = os.path.join(
                    self.tmp, f'{self.projection_name}_exons{exon_nums}_{tr}.cesin'
                ) if self.save_cesar_input else None
                space_name: str = f'{self.projection_name} {group.chrom}:{group.start}-{group.stop}'
                cesar_input: Union[str, None] = dump_for_cesar( ## TODO: Consider turning it into a main()'s method
                    group.exon_headers,
                    group.exon_seqs,
                    [space_name],
                    [group.query_seq.replace('|', '')],
                    cesar_input_file
                )
                if cesar_input is not None:
                    cesar_input = cesar_input.encode()

                ## create a CESAR command
                mem: float = ceil(group.memory + 0.1)
                self.max_mem = max(mem, self.max_mem)
                cesar_cmd: str = (
                    f'{self.cesar_binary} '
                    f'{cesar_input_file if self.save_cesar_input else "/dev/stdin"} '
                    f'--max-memory {mem}'
                )
                # if len(exon_group) == 1 and 1 in exon_group:
                if 1 in group.exons and not self.toga1:
                    cesar_cmd += ' -f'
                # if len(exon_group) == 1 and self.annot_entry.exon_number in exon_group:
                if group.contains_last_exon and not self.toga1:
                    cesar_cmd += ' -l'

                ## run CESAR
                self._to_log(
                    f'Aligning exon group {g} for projection {self.projection_name}'
                )
                cesar_out: str = self._exec(cesar_cmd, self.cesar_err, cesar_input)

                ## parse the CESAR results and store the CesarExonEntry objects
                cesar_lines = cesar_out.split('\n')
                result: RawCesarOutput = RawCesarOutput(
                    group.chain,
                    group.chrom,
                    group.start,
                    group.stop,
                    group.strand,
                    cesar_lines[1],
                    cesar_lines[3],
                    group.exons,
                    group.expected_coords,
                    group.search_space_coords,
                    group.spliceai_data,
                    group.gap_located_exons,
                    group.out_of_chain_exons,
                    not group.will_be_aligned,
                    group.intersects_asmbl_gaps,
                    {} ## for intron coordinates
                )
                self.cesar_results[tr][g].append(result)
                groups_finished += 1
            segments_finished += 1 if groups_finished else 0
            if not groups_finished:
                status: str = self._classify_unaligned_projection(tr)
                self.rejected_projections.append(UNALIGNED_REJ.format(tr, status))

        ## write rejected projections data to the file
        if self.rejected_projections:
            if self.p_job:
                with FileLock(self.rejection_lock, timeout=5):
                    with open(self.rejection_file, 'a', buffering=1) as h:
                        for line in self.rejected_projections:
                            h.write(line + '\n')
                        h.flush()
                        os.fsync(h.fileno())
            else:
                with open(self.rejection_file, 'w') as h:
                    for line in self.rejected_projections:
                        h.write(line + '\n')

        ## once all CESAR jobs are run, process the results
        for tr, cesar_results in self.cesar_results.items():
            ## for each segment, construct all the possible configurations of
            ## alternative exons/exon groups
            for alt_num, alt_segment in enumerate(ddfs(cesar_results, 0, []), 1):
                alt_segment = sorted(alt_segment, key=lambda x: x.exons[0])
                if self.intron_gain_check:
                    self._to_log('Checking segment %s for intron gain events' % alt_num)
                    uiu = IntronGainChecker(
                        alt_segment,
                        self.min_intron_gain_score,
                        # self.min_intron_prob_gapped,
                        # self.min_intron_prob_ungapped,
                        self.min_intron_prob_trusted,
                        self.min_intron_prob_supported,
                        self.min_intron_prob_unsupported,
                        self.max_intron_number,
                        self.logger
                    )
                    uiu.record_intron_gains()
                    updated_segment = uiu.return_updated_cesar_output()
                    upd_segment: List[RawCesarOutput] = []
                    direction: bool = all(
                        alt_segment[x-1].exons[-1] < alt_segment[x].exons[0]
                        for x in range(1, len(alt_segment))
                    ) if len(alt_segment) > 1 else True
                    for uup, upd_portion in enumerate(updated_segment):
                        if isinstance(upd_portion, RawCesarOutput):
                            # print(f'{upd_portion.exons=}, {alt_segment[uup].exons=}')
                            upd_segment.append(upd_portion)
                        else:
                            elements: List[int] = sorted(upd_portion.ref_seqs.keys())
                            combined_ref_seq: str = ''
                            combined_query_seq: str = ''
                            intron_coordinates: Dict[int, List[Tuple[int, int]]] = {}
                            
                            for el in elements:
                                if upd_portion.needs_realigment[el]:
                                    ex: int = upd_portion.el2exons[el][0]
                                    up_ref, up_query = upd_portion.upstream_seqs[el]
                                    combined_ref_seq += up_ref
                                    combined_query_seq += up_query
                                    layout2id: Dict[int, float] = {}
                                    layout2seqs: Dict[int, Tuple[str, str]] = {}
                                    layout2intron_coords: Dict[int, List[Tuple[int, int]]] = {}
                                    # print(f'{len(upd_portion.exon_seqs)=}, {len(elements)=}')
                                    for layout in range(len(upd_portion.ref_seqs[el])):
                                        # print(f'{len(upd_portion.exon_seqs[el])=}, {layout=}')
                                        extra_cesar_input: str = ''
                                        subexon_num: int = len(upd_portion.exon_seqs[el][layout]) - 1
                                        # extra_cesar_headers: List[str] =
                                        extra_headers: List[str] = []
                                        extra_seqs: List[str] = []
                                        for s, sub in enumerate(upd_portion.exon_seqs[el][layout]):
                                            acc_profile = (
                                                self.first_acceptor if ex == 1 and s == 0 
                                                else self.regular_acceptor
                                            )
                                            donor_profile = (
                                                self.last_donor if ex == self.max_exon_num and s == subexon_num 
                                                else self.regular_donor
                                            )
                                            extra_headers.append(f'{s}\t{acc_profile}\t{donor_profile}')
                                            extra_seqs.append(sub)
                                        extra_query_name: str = f'>{self.transcript}'
                                        extra_cesar_query: str = upd_portion.ref_seqs[el][layout]
                                        extra_input_file: str = os.path.join(
                                            self.tmp, f'{self.transcript}_intron_gain_exon{ex}.cesin'
                                        ) if self.save_cesar_input else None
                                        extra_cesar_input: Union[str, None] = dump_for_cesar(
                                            extra_headers,
                                            extra_seqs,
                                            [extra_query_name],
                                            [extra_cesar_query],
                                            extra_input_file
                                        )
                                        if extra_cesar_input is not None:
                                            extra_cesar_input = extra_cesar_input.encode('utf8')
                                        # mem = 3 ## sic!
                                        raw_extra_mem: float = cesar_memory_check(
                                            [len(x) for x in extra_seqs],
                                            len(extra_cesar_query)
                                        )
                                        adj_extra_mem: int = ceil(raw_extra_mem + 0.1)
                                        # print(f'{raw_extra_mem=}, {adj_extra_mem=}, {self.max_mem=}')
                                        if adj_extra_mem > self.max_mem:
                                            self._die(
                                                (
                                                    'ERROR: CESAR2 realignment job for intron annotation in exon %i '
                                                    'requires more memory (%i) than the most memory-consuming regular CESAR 2'
                                                    'job (%i GB)'
                                                ) % (ex, adj_extra_mem, self.max_mem)
                                            )

                                        extra_cesar_cmd: str = (
                                            f'{self.cesar_binary} '
                                            f'{extra_input_file if self.save_cesar_input else "/dev/stdin"} '
                                            f'--max-memory {adj_extra_mem}'
                                        )
                                        if ex == 1:
                                            cesar_cmd += ' -f'
                                        # if len(exon_group) == 1 and self.annot_entry.exon_number in exon_group:
                                        if ex == self.max_exon_num:
                                            cesar_cmd += ' -l'
                                        # print(f'{extra_cesar_input=}')
                                        cesar_out: str = self._exec(extra_cesar_cmd, self.cesar_err, extra_cesar_input)

                                        ## parse the CESAR results and store the CesarExonEntry objects
                                        cesar_lines = cesar_out.split('\n')
                                        # print(f'{cesar_lines=}')
                                        ## first, check for frame integrity
                                        if frameshift_in_cesar_data(cesar_lines[1], cesar_lines[3]):
                                            continue
                                        # print(f'{upd_portion.introns[el][layout]=}')
                                        corr_ref, corr_query, subexon_coords, frameshifting_aln = introduce_intron_sequences(
                                            cesar_lines[3], cesar_lines[1], upd_portion.introns[el][layout]
                                            # *upd_portion.exon2extra_seq[ex]
                                        )
                                        if frameshifting_aln:
                                            continue
                                        prev_seq_len: int = sum(len(x.query) for x in upd_segment) + len(combined_ref_seq)
                                        upd_subexon_coords: List[Tuple[int, int]] = []
                                        for subexon in subexon_coords:
                                            # print(f'{subexon=}')
                                            upd_subexon: Tuple[int, int] = (
                                                subexon[0] + prev_seq_len,# + offset,
                                                subexon[1] + prev_seq_len# + offset
                                            )
                                            # print(f'{upd_subexon=}')
                                            upd_subexon_coords.append(upd_subexon)
                                        # print(f'{subexon_coords=}, {upd_subexon_coords=}')
                                        # print(f'{len(extra_cesar_query)=}, {len(corr_ref)=}')
                                        # print(f'{extra_cesar_query=}')
                                        # print(f'{corr_ref=}')
                                        layout2intron_coords[layout] = upd_subexon_coords
                                        ref_left_phase, ref_right_phase = upd_portion.ref_phases[el]
                                        corr_ref_len: int = len(corr_ref)
                                        # print(f'{ref_left_phase=}, {ref_right_phase=}')
                                        corr_ref = (
                                            corr_ref[:ref_left_phase].lower() +
                                            corr_ref[ref_left_phase:corr_ref_len-ref_right_phase] +
                                            corr_ref[corr_ref_len - ref_right_phase:].lower()
                                        )
                                        # intron_coordinates[ex] = [
                                        #     (x[0] + len(combined_ref_seq), x[1] + len(combined_ref_seq))
                                        #     for x in intron_coords
                                        # ]
                                        # print(f'{corr_ref=}')
                                        # print(f'{corr_query=}')
                                        identity: float = fast_seq_id(corr_ref, corr_query)
                                        # print(f'{ex=}, {identity=}, {upd_portion.orig_id[ex]=}')
                                        if identity >= upd_portion.orig_id[ex] - 0.1:
                                            layout2id[layout] = identity
                                            layout2seqs[layout] = (corr_ref, corr_query)
                                            # self._to_log(f'Evidence for the winning layout: {upd_portion.intron_evidence[ex][layout]}')
                                            self.intron2evidence[ex] = upd_portion.intron_evidence[ex][layout]
                                    # print(f'{layout2id=}', upd_portion.orig_id[ex])
                                    if layout2id:
                                        best_cand: int = max(layout2id.items(), key=lambda x: x[1])[0]
                                        corr_ref, corr_query = layout2seqs[best_cand]
                                        combined_ref_seq += corr_ref
                                        combined_query_seq += corr_query
                                        intron_coordinates[ex] = layout2intron_coords[best_cand]
                                    else:
                                        combined_ref_seq += upd_portion.orig_ref[ex]
                                        combined_query_seq += upd_portion.orig_query[ex]
                                    down_ref, down_query = upd_portion.downstream_seqs[el]
                                    combined_ref_seq += down_ref
                                    combined_query_seq += down_query
                                else:
                                    combined_ref_seq += upd_portion.ref_seqs[el]
                                    combined_query_seq += upd_portion.query_seqs[el]
                            # print(f'{combined_ref_seq=}')
                            # print(f'{combined_query_seq=}')
                            # print(f'{upd_portion.subexon_coords=}')
                            # print(f'{intron_coordinates=}')
                            # print(f'{len(alt_segment)=}, {len(updated_segment)=}')
                            orig_num: int = (
                                uup if direction else len(alt_segment) - uup - 1
                            )
                            orig_raw_output: RawCesarOutput = alt_segment[orig_num]
                            # print(f'{len(orig_raw_output.reference)=}, {len(orig_raw_output.query)=}, {len(orig_raw_output.reference.strip(" "))=}, {len(orig_raw_output.query.strip(" "))=}, {len(combined_ref_seq)=}, {len(combined_query_seq)=}')
                            # print(f'{orig_raw_output.exons=}, {upd_portion.el2exons=}')
                            # print(f'{orig_raw_output.reference=}\n{orig_raw_output.query=}\n{combined_ref_seq=}\n{combined_query_seq=}')
                            # print(f'{orig_raw_output.chrom}\t{orig_raw_output.start}\t{orig_raw_output.stop}')
                            upd_raw_output: RawCesarOutput = RawCesarOutput(
                                orig_raw_output.chain,
                                orig_raw_output.chrom,
                                orig_raw_output.start,
                                orig_raw_output.stop,
                                orig_raw_output.strand,
                                combined_ref_seq,
                                combined_query_seq,
                                orig_raw_output.exons,
                                orig_raw_output.exon_expected_loci,
                                orig_raw_output.exon_search_spaces,
                                orig_raw_output.spliceai_sites,
                                orig_raw_output.gap_located_exons,
                                orig_raw_output.out_of_chain_exons,
                                orig_raw_output.was_not_aligned,
                                orig_raw_output.assembly_gap,
                                # upd_portion.subexon_coords
                                intron_coordinates
                            )
                            upd_segment.append(upd_raw_output)
                    # print(f'{len(alt_segment)=}, {len(updated_segment)=}')
                    alt_segment = upd_segment
                # return
                self._to_log('Processing segment data')
                out_file_postfix: str = f'segment{tr}_config{alt_num}'
                processed_segment: ProcessedSegment = ProcessedSegment(
                    self.transcript,
                    tr,
                    alt_num,
                    alt_segment,
                    self.aa_matrix,
                    # ref_intron_lengths,
                    self.u12_sites,
                    self.mask_terminal_mutations,
                    self.rescue_missing_start,
                    self.rescue_missing_stop,
                    self.assembly_gap_size,
                    self.sai_correction_mode,
                    self.min_splice_prob,
                    self.splice_prob_margin,
                    self.acceptor_flanks[tr],
                    self.donor_flanks[tr],
                    self.correct_short_introns and not self.toga1,
                    self.ignore_alternative_frame,
                    self.is_paralog,
                    self.is_processed_pseudogene,
                    self.logger,
                    self.v
                )
                ## run the processor
                processed_segment.run()
                ## save the results
                if self.p_job:
                    with FileLock(self.bed_lock, timeout=5):
                        with open(self.bed_path, 'a', buffering=1) as h:
                            self._to_log(f'Writing raw BED data for {self.projection_name}')
                            for subsegment in processed_segment.bed12():
                                h.write(subsegment + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.filt_bed_lock, timeout=5):
                        with open(self.filt_bed_path, 'a', buffering=1) as h:
                            self._to_log(f'Writing final BED data for {self.projection_name}')
                            for subsegment in processed_segment.bed12(raw=False):
                                h.write(subsegment + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.browser_bed_lock, timeout=5):
                        with open(self.browser_bed_file, 'a', buffering=1) as h:
                            self._to_log(f'Writing UCSC browser format BED for {self.projection_name}')
                            for subsegment in processed_segment.bed12(raw=False, browser=True):
                                h.write(subsegment + '\n')
                    with FileLock(self.cds_fasta_lock, timeout=5):
                        with open(self.cds_fasta, 'a') as h:
                            self._to_log(
                                f'Writing nucleotide sequence FASTA for {self.projection_name}'
                            )
                            h.write(processed_segment.cds_nuc() + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    # with FileLock(self.cesar_res_lock, timeout=5):
                    #     with open(self.cesar_res_stub, 'a', buffering=1) as h:
                    #         self._to_log(f'Writing CESOUT data for {self.transcript}#{self.chain}')
                    #         h.write(processed_segment.bdb() + '\n')
                    #         h.flush()
                    #         os.fsync(h.fileno())
                    with FileLock(self.id_lock, timeout=5):
                        with open(self.id_stub, 'a+', buffering=1) as h:
                            self._to_log(f'Writing transcript metadata for {self.projection_name}')
                            h.write(
                                processed_segment.transcript_meta() + '\n'
                            )
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.exon_meta_lock, timeout=5):
                        with open(self.exon_meta_stub, 'a+', buffering=1) as h:
                            self._to_log(f'Writing exon metadata for {self.projection_name}')
                            for line in processed_segment.exon_meta():
                                h.write(line + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.mutation_lock, timeout=5):
                        with open(self.mutation_stub, 'a+', buffering=1) as h:
                            self._to_log(f'Writing mutation data for {self.projection_name}')
                            h.write(processed_segment.mutation_file())
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.codon_lock, timeout=5):
                        with open(self.codon_stub, 'a+', buffering=1) as h:
                            self._to_log(f'Writing codon alignment for {self.projection_name}')
                            h.write(processed_segment.codon_fasta() + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.prot_lock, timeout=5):
                        with open(self.prot_stub, 'a+', buffering=1) as h:
                            self._to_log(f'Writing protein alignment for {self.projection_name}')
                            h.write(processed_segment.aa_fasta() + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.exon_fa_lock, timeout=5):
                        with open(self.exon_fa_stub, 'a+', buffering=1) as h:
                            self._to_log(f'Writing exon alignment for {self.projection_name}')
                            h.write(processed_segment.exon_fasta() + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.splice_site_lock, timeout=5):
                        with open(self.splice_site_stub, 'a+', buffering=1) as h:
                            self._to_log(
                                f'Writing splice site dinucleotide data for {self.projection_name}'
                            )
                            shifts: str = processed_segment.splice_site_table()
                            if shifts:
                                h.write(shifts + '\n')
                                h.flush()
                                os.fsync(h.fileno())
                    with FileLock(self.intron_evidence_lock, timeout=5):
                        with open(self.intron_evidence_stub, 'a+', buffering=1) as h:
                            self._to_log(
                                f'Writing gained intron evidence data for {self.projection_name}'
                            )
                            for ex, introns in self.intron2evidence.items():
                                for intron in introns:
                                    out_line: str = f'{self.projection_name}\t{ex}\t{intron.__repr__()}'
                                    h.write(out_line + '\n')
                                    h.flush()
                                    os.fsync(h.fileno())
                    with FileLock(self.splice_shift_lock, timeout=5):
                        self._to_log(
                            f'Writing splice site shift data for {self.projection_name}'
                        )
                        with open(self.splice_site_shifts, 'a', buffering=1) as h:
                            h.write(processed_segment.splice_site_shifts() + '\n')
                            h.flush()
                            os.fsync(h.fileno())
                    with FileLock(self.selenocysteine_lock, timeout=5):
                        self._to_log(
                            'Writing selenocysteine codon data for %s' % self.projection_name
                        )
                        seleno_table: str = processed_segment.selenocysteine_codon_table()
                        with open(self.selenocysteine_codons, 'a') as h:
                            if seleno_table:
                                h.write(seleno_table + '\n')
                    continue
                else:
                    self._to_log(f'Writing raw BED data for {self.projection_name}')
                    with open(self.bed_path, 'a') as h:
                        for subsegment in processed_segment.bed12():
                            h.write(subsegment + '\n')
                    self._to_log(f'Writing final BED data for {self.projection_name}')
                    with open(self.filt_bed_path, 'a') as h:
                        for subsegment in processed_segment.bed12(raw=False):
                            h.write(subsegment + '\n')
                    with open(self.browser_bed_file, 'a') as h:
                        for subsegment in processed_segment.bed12(raw=False, browser=True):
                            h.write(subsegment + '\n')
                    self._to_log(
                        f'Writing nucleotide sequence FASTA for {self.projection_name}'
                    )
                    with open(self.cds_fasta, 'a') as h:
                        h.write(processed_segment.cds_nuc() + '\n')
                    # with open(self.cesar_res_stub, 'a') as h:
                    #     h.write(processed_segment.bdb())
                    self._to_log(f'Writing transcript metadata for {self.projection_name}')
                    with open(self.id_stub, 'a') as h:
                            h.write(processed_segment.transcript_meta() + '\n')
                    self._to_log(f'Writing exon metadata for {self.projection_name}')
                    with open(self.exon_meta_stub, 'a') as h:
                        for line in processed_segment.exon_meta():
                            h.write(line + '\n')
                    self._to_log(f'Writing mutation data for {self.projection_name}')
                    with open(self.mutation_stub, 'a') as h:
                        h.write(processed_segment.mutation_file())
                    self._to_log(f'Writing codon alignment for {self.projection_name}')
                    with open(self.codon_stub, 'a') as h:
                        h.write(processed_segment.codon_fasta() + '\n')
                    self._to_log(f'Writing protein alignment for {self.projection_name}')
                    with open(self.prot_stub, 'a') as h:
                        h.write(processed_segment.aa_fasta() + '\n')
                    self._to_log(f'Writing exon alignment for {self.projection_name}')
                    with open(self.exon_fa_stub, 'a') as h:
                        h.write(processed_segment.exon_fasta() + '\n')
                    self._to_log(
                        f'Writing splice site dinucleotide data for {self.projection_name}'
                    )
                    with open(self.splice_site_stub, 'a') as h:
                        h.write(processed_segment.splice_site_table() + '\n')
                    self._to_log(
                        f'Writing gained intron evidence data for {self.projection_name}'
                    )
                    with open(self.intron_evidence_stub, 'a', buffering=1) as h:
                        for ex, introns in self.intron2evidence.items():
                            for intron in introns:
                                out_line: str = f'{self.projection_name}\t{ex}\t{intron.__repr__()}'
                                h.write(out_line + '\n')
                    self._to_log(
                        f'Writing splice site shift data for {self.projection_name}'
                    )
                    shifts: str = processed_segment.splice_site_shifts()
                    with open(self.splice_site_shifts, 'a', buffering=1) as h:
                        if shifts:
                            h.write(shifts + '\n')
                    self._to_log(
                        'Writing selenocysteine codon data for %s' % self.projection_name
                    )
                    seleno_table: str = processed_segment.selenocysteine_codon_table()
                    with open(self.selenocysteine_codons, 'a') as h:
                        if seleno_table:
                            h.write(seleno_table + '\n')
                    continue

if __name__ == '__main__':
    CesarExecutor()
