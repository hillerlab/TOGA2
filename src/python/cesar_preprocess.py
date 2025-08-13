#!/usr/bin/env python3

"""
A chainwise CESAR data preprocessing module. Given the chain number and list of
projections, assesses the number of projections and CESAR exon groups,
prepares the CESAR input and dumps it to an HDF5 file
"""

__author__  = 'Yury V. Malovichko'
__credits__ = ['Bogdan Kirilenko', 'Michael Hiller']
__year__ = '2024'

from collections import Counter, defaultdict
from filelock import FileLock # remove if not needed
from numpy import array, str_
from math import ceil
from modules.cesar_wrapper_constants import (
    ACCEPTOR_SITE, ACCEPTOR_SITE_U12,
    DONOR_SITE, DONOR_SITE_U12,
    EXTRA_FLANK, FLANK_SPACE,
    HL_BW2W_PATH,
    HG38_CANON_U2_ACCEPTOR, HG38_CANON_U2_DONOR,
    HG38_NON_CANON_U2_ACCEPTOR, HG38_NON_CANON_U2_DONOR,
    HG38_CANON_U12_ACCEPTOR, HG38_CANON_U12_DONOR,
    HG38_NON_CANON_U12_ACCEPTOR, HG38_NON_CANON_U12_DONOR,
    FIRST_ACCEPTOR, LAST_DONOR, 
    LARGE_EXON_UNALIGNED_WARNING, MEM_UNALIGNED_WARNING,
    MIN_ASMBL_GAP_SIZE, MIN_REF_LEN_PERC,
    SHORT_SPACE_UNALGNED_WARNING, 
    SINGLE_EXON_MIN_REF_LEN_PERC,
    SPLICEAI_PROCESS_SCRIPT, U2, U12
)
from modules.preprocessing import (
    AnnotationEntry, Exon, ExonDict, Exon2BlockMapper, 
    ProjectionCoverageData, ProjectionGroup,
    Segment, cesar_memory_check, 
    get_chain, find_gaps, intersect_exons_to_blocks,
    parse_extr_exon_fasta, prepare_exons
)
from pathlib import Path
from string_splitter import transcriptwise_subchains ## TODO: Rename the module
from modules.shared import (
    CommandLineManager, CONTEXT_SETTINGS, dir_name_by_date, intersection, 
    hex_code, parts, reverse_complement
)
from typing import Dict, Iterable, List, Optional, Set, TextIO, Tuple, Union

import click
import h5py
import os

## create or update the needed constants
LOCATION: str = os.path.dirname(os.path.abspath(__file__))
HG38_CANON_U2_ACCEPTOR: str = os.path.join(LOCATION, *HG38_CANON_U2_ACCEPTOR)
HG38_CANON_U2_DONOR: str = os.path.join(LOCATION, *HG38_CANON_U2_DONOR)
HG38_NON_CANON_U2_ACCEPTOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U2_ACCEPTOR)
HG38_NON_CANON_U2_DONOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U2_DONOR)
HG38_CANON_U12_ACCEPTOR: str = os.path.join(LOCATION, *HG38_CANON_U12_ACCEPTOR)
HG38_CANON_U12_DONOR: str = os.path.join(LOCATION, *HG38_CANON_U12_DONOR)
HG38_NON_CANON_U12_ACCEPTOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U12_ACCEPTOR)
HG38_NON_CANON_U12_DONOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U12_DONOR)
FIRST_ACCEPTOR: str = os.path.join(LOCATION, *FIRST_ACCEPTOR)
LAST_DONOR: str = os.path.join(LOCATION, *LAST_DONOR)
# HL_COMMON_ACCEPTOR: str = os.path.join(*HL_COMMON_ACCEPTOR)
# HL_COMMON_DONOR: str = os.path.join(*HL_COMMON_DONOR)
# HL_FIRST_ACCEPTOR: str = os.path.join(*HL_FIRST_ACCEPTOR)
# HL_LAST_DONOR: str = os.path.join(*HL_LAST_DONOR)
# HL_EQ_ACCEPTOR: str = os.path.join(LOCATION, *HL_EQ_ACCEPTOR)
# HL_EQ_DONOR: str = os.path.join(LOCATION, *HL_EQ_DONOR)
PREPROCESSING_REJ: str = 'PROJECTION\t{}\t{}\t{}\t{}\t{}'
SPANNING_CHAIN_REASON: str = 'PROJECTION\t{}\t0\tSpanning chain\tSPANNING\t{}'
LOST: str = 'L'
MISSING: str = 'M'

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'transcripts',
    type=str,
    metavar='TRANSCRIPTS'
)
@click.argument(
    'chains',
    type=str,
    metavar='CHAIN(S)'
)
@click.argument(
    'ref_annot',
    type=click.File('r', lazy=True),
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
@click.option(
    '--segments',
    '-s',
    type=click.File('r', lazy=True),
    metavar='BED_FILE',
    default=None,
    show_default=True,
    help='A BED12 file containing segments with external evidence'
)
@click.option(
    '--parallel_job',
    '-p',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'Indicates that the process is a part of parallel job batch; '
        'output file names are standardized and provided with lock.files'
    )
)
@click.option(
    '--disable_spanning_chains',
    '-nospan',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=False,
    help=(
        'If set, enables spanning chains '
        '(i.e. chains with no coding blocks corresponding to transcript exons)'
    )
)
@click.option(
    '--processed_pseudogene_list',
    '-pplist',
    type=str,
    metavar='PPGENE_LIST',
    default=None,
    show_default=True,
    help=(
        'A comma-separated list of transcript whose projections through this chain '
        'lead to retroposed genes'
    )
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
    help=(
        'Upper memory limit for CESAR jobs. If limit is exceeded, the program '
        'terminates with zero exit status'
    )
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
    type=float,
    metavar='FLOAT',
    default=0.0,
    help=(
        'Minimal fraction of reference CDS to be covered by alignment data. '
        'Projections covering less than this portion will be discarded.'
    )
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
    help='A path to first exon acceptor profile'
)
@click.option(
    '--cesar_last_donor',
    '-cld',
    type=str,
    metavar='REL_PATH',
    default=LAST_DONOR,
    show_default=True,
    help='A path to last exon donor profile'
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
    help='An HDF5 file containing data on U12 and non-canonical U2 introns' ## TODO: Expand the description
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
    default=HL_BW2W_PATH,
    help=(
        'A path to the UCSC bigWigToWig binary; default is set to Hiller Lab '
        'Delta cluster utility location'
    )
)
@click.option(
    '--min_splice_prob',
    '-msp',
    type=float,
    metavar='FLOAT',
    default=0.5,
    show_default=True,
    help='Minimum SpliceAI prediction probability to consider the splice site'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(exists=False),
    metavar='OUT_DIR',
    default=dir_name_by_date('segment_reconstruction_'),
    show_default=False,
    help='A directory to write the output to '
    '[default: cesar_preprocessing_{run date}'
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
    '--verbose',
    '-v',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls the execution verbosity'
)


class CesarPreprocessor(CommandLineManager):
    __slots__ = [
        'transcripts', 'chains', 'annot_entries', 'ref', 'query', 'chain_file',
        'ref_chrom_sizes', 'query_chrom_sizes',
        'annot_entries', 'segments', 'u12_sites',
        'disable_spanning_chains', 'ppgene_list',
        'p_job', 'no_inference', 'memory_limit',
        'max_space_size', 'extrapolation_modifier', 'min_cov_portion',
        'exon_locus_flank',
        'canon_u2_acceptor', 'canon_u2_donor',
        'non_canon_u2_acceptor', 'non_canon_u2_donor',
        'canon_u12_acceptor', 'canon_u12_donor',
        'non_canon_u12_acceptor', 'non_canon_u12_donor',
        'cesar_first_acceptor', 'cesar_last_donor',
        'separate_site_treat',
        'assembly_gap_size', 'spliceai_dir',
        'bigwig2wig_binary', 'min_splice_prob', 'output',
        'fragmented', 'chainmappers', 'groups', 'universally_missing', 'found_exons',
        'query_2bit',
        'prepared_exons', 'exon_names_for_cesar', 'missing2chain', 'chain2coords',
        'copy2groups', 'copy2proj', 'rejected_transcripts', 'spanning_chains',
        'rejection_reasons', 'unaligned_exons', 'max_est_ram', 'cumulative_ram',
        'largest_ss', 'largest_exon', 'transcripts_finished', 'exon_table', 'tr2cov',
        'memory_file', 'rejection_file', 'spanning_chain_file', 'exon_storage',
        'memory_lock', 'rejection_lock', 'spanning_chain_lock','exon_storage_lock', 'bw2w_err',
        'log_file',

        'toga1', 'toga1_plus_cesar',

        'v'
    ]

    def __init__(
        self,
        transcripts: str,
        chains: str,
        ref_annot: click.Path,
        ref: click.Path,
        query: click.Path,
        chain_file: click.Path,
        ref_chrom_sizes: click.Path,
        query_chrom_sizes: click.Path,
        segments: Optional[click.File],
        parallel_job: Optional[bool],
        disable_spanning_chains: Optional[bool],
        processed_pseudogene_list: Optional[str],
        no_inference: Optional[bool],
        memory_limit: Optional[float],
        max_space_size: Optional[int],
        extrapolation_modifier: Optional[float],
        minimal_covered_fraction: Optional[float],
        exon_locus_flank: Optional[int],
        cesar_canon_u2_acceptor: Optional[str],
        cesar_canon_u2_donor: Optional[str],
        cesar_non_canon_u2_acceptor: Optional[str],
        cesar_non_canon_u2_donor: Optional[str],
        cesar_canon_u12_acceptor: Optional[str],
        cesar_canon_u12_donor: Optional[str],
        cesar_non_canon_u12_acceptor: Optional[str],
        cesar_non_canon_u12_donor: Optional[str],
        cesar_first_acceptor: Optional[str],
        cesar_last_donor: Optional[str],
        separate_splice_site_treatment: Optional[bool],
        assembly_gap_size: Optional[int],
        u12: Optional[Union[click.Path, None]],
        spliceai_dir: Optional[click.Path],
        bigwig2wig_binary: Optional[click.Path],
        min_splice_prob: Optional[float],
        output: Optional[click.Path],
        toga1_compatible: Optional[bool],
        toga1_plus_corrected_cesar: Optional[bool],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self._echo('Initializing CESAR preprocessing module')
        self._echo('Checking if output directory exists')
        self.output: str = output if output else dir_name_by_date('cesar_preprocessing_')
        if not os.path.isdir(self.output):
            self._to_log('Creating output directory')
            self._mkdir(self.output)
        self.log_file: str = os.path.join(self.output, 'log.txt')
        self.set_logging(f'{__name__}_{hex_code()}')
        self.chains: Tuple[str] = tuple(chains.split(','))
        self._to_log(
            f'Input chain{" is" if len(self.chains) == 1 else "a are"}: {chains}'
        )
        self.fragmented: bool = len(self.chains) > 1
        if self.fragmented:
            self._to_log(
                'Multiple chain IDs were provided; expecting all '
                'projections to be fragmented'
            )
        else:
            self._to_log(
                'Single chain ID was provided; expecting all '
                'projections to be contained within the chain'
            )
        self.transcripts: Tuple[str] = tuple(transcripts.split(','))
        self._to_log(
            f'Input transcript{" is" if len(self.transcripts) == 1 else "s are"}: '
            f'{transcripts}'
        )
        if len(self.chains) > 1 and len(self.transcripts) > 1:
            self._to_log(
                'WARNING: multiple transcripts were provided for fragmented '
                'assembly; make sure it is intended behavior'
            )
        self.ref: str = ref
        self.query: str = query
        self.annot_entries: Dict[str, AnnotationEntry] = {}
        self._to_log('Parsing reference annotation BED file')
        self.parse_annotation_bed(ref_annot)
        self.chain_file: click.Path = chain_file
        self.segments: Dict[int, List[Segment]] = defaultdict(dict)
        self.parse_segments(segments)
        self.ref_chrom_sizes: Dict[str, int] = {}
        self.parse_chrom_size_file(ref_chrom_sizes, ref=True)
        self.query_chrom_sizes: Dict[str, int] = {}
        self.parse_chrom_size_file(query_chrom_sizes, ref=False)
        ## if U12 data were provided, parse them
        # self.u12_sites: Dict[int, Set[int]] = defaultdict(lambda: defaultdict(set))
        self.u12_sites: Dict[str, Dict[int, Dict[str, Tuple[str]]]] = {
            x: {y: {'donor': (), 'acceptor': ()} for y in self.annot_entries[x].exons}
            for x in self.transcripts
        }
        self.separate_site_treat: bool = separate_splice_site_treatment
        if u12:
            self._to_log('Parsing the U12 intron file')
            # for tr in self.transcripts:
            #     self.u12_sites[tr] = parseU12(u12, tr)
            self.parse_u12_file(u12)
        else:
            self._to_log('No external U12 data were provided')

        self.disable_spanning_chains: bool = disable_spanning_chains
        self.ppgene_list: str = (
            [] if processed_pseudogene_list is None else 
            [x for x in processed_pseudogene_list.split(',') if x]
        )
        self.p_job: bool = parallel_job
        self.no_inference: bool = no_inference
        self.memory_limit: float = memory_limit
        self.max_space_size: int = max_space_size
        self.extrapolation_modifier: float = extrapolation_modifier
        self.min_cov_portion: float = minimal_covered_fraction
        self.exon_locus_flank: int = exon_locus_flank
        self.canon_u2_acceptor: str = cesar_canon_u2_acceptor
        self.canon_u2_donor: str = cesar_canon_u2_donor
        self.non_canon_u2_acceptor: str = cesar_non_canon_u2_acceptor
        self.non_canon_u2_donor: str = cesar_non_canon_u2_donor
        self.canon_u12_acceptor: str = cesar_canon_u12_acceptor
        self.canon_u12_donor: str = cesar_canon_u12_donor
        self.non_canon_u12_acceptor: str = cesar_non_canon_u12_acceptor
        self.non_canon_u12_donor: str = cesar_non_canon_u12_donor
        self.cesar_first_acceptor: str = cesar_first_acceptor
        self.cesar_last_donor: str = cesar_last_donor
        self.assembly_gap_size: int = assembly_gap_size
        self.bigwig2wig_binary: str = bigwig2wig_binary
        self.spliceai_dir: str = spliceai_dir
        self.min_splice_prob: float = (
            1.0 if min_splice_prob > 1 else 0.0 if min_splice_prob < 0 else min_splice_prob
        )

        self.toga1: bool = toga1_compatible
        self.toga1_plus_cesar: bool = toga1_plus_corrected_cesar

        self.chainmappers: Dict[str, Dict[str, Exon2BlockMapper]] = defaultdict(dict)
        self.groups: Dict[str, Dict[int, Tuple[int]]] = defaultdict(dict)
        self.universally_missing: Dict[str, Dict[int, Set[int]]] = defaultdict(
            lambda: defaultdict(set)
        )
        self.found_exons: Dict[str, Dict[int, Set[int]]] = defaultdict(
            lambda: defaultdict(set)
        )
        self.prepared_exons: Dict[str, Dict[int, str]] = {}
        self.exon_names_for_cesar: Dict[str, Dict[int, str]] = defaultdict(dict)
        self.missing2chain: Dict[int, Dict[int, str]] = defaultdict(dict)
        self.copy2groups: Dict[str, Dict[int, Set[int]]] = defaultdict(
            lambda: defaultdict(set)
        )
        self.copy2proj: Dict[str, Dict[int, List[ProjectionGroup]]] = defaultdict(
            lambda: defaultdict(dict)
        )
        self.rejected_transcripts: List[str] = []
        self.rejection_reasons: Dict[str, Dict[int, List[str]]] = defaultdict(
            lambda: defaultdict(list)
        )
        self.unaligned_exons: Dict[str, Dict[int, Dict[str]]] = defaultdict(
            lambda: defaultdict(dict)
        )
        self.spanning_chains: List[Tuple[str, int]] = []
        self.max_est_ram: Dict[str, float] = {}
        self.cumulative_ram: Dict[str, float] = defaultdict(float)
        self.largest_ss: Dict[str, int] = {}
        self.largest_exon: Dict[str, str] = {}
        self.tr2cov: Dict[str, float] = {}
        self.transcripts_finished: int = 0
        self.exon_table: Dict[List[str]] = defaultdict(list)
        self.chain2coords: Dict[str, Tuple[int, int]] = {}

        self.memory_file: str = os.path.join(self.output, 'max_memory_requirements.tsv')
        self.rejection_file: str = os.path.join(self.output, 'genes_rejection_reason.tsv')
        self.spanning_chain_file: str = os.path.join(self.output, 'spanning_chains_ref_coords.tsv')
        self.exon_storage: str = os.path.join(output, 'exon_storage.hdf5')
        if self.p_job:
            self.memory_lock: str = self.memory_file + '.lock'
            self.rejection_lock: str = self.rejection_file + '.lock'
            self.spanning_chain_lock: str = self.spanning_chain_file + '.lock'
            self.exon_storage_lock: str = self.exon_storage + '.lock'

        self.bw2w_err: str = 'ERROR: bigWigToBed run failed with the following error'

        self.run()

    def run(self) -> None:
        """
        """
        ## Part I: Prepare output directory structure
        self._to_log('Checking if temporary directory exists')
        self._to_log('Output and temporary directory were successfully created')

        ## Part II: Extract the necessary data
        ## for each chain-transcript pair, get the projection subchain mapper
        self._to_log('Retrieving chain data')
        self.get_chain_mappers()

        ## if a fragmented projection is provided, reconcile the shared exons
        if self.fragmented:
            self.resolved_shared_exons()

        ## get exon sequences
        self._to_log(
            'Preparing reference exons for CESAR alignment'
        )
        self.get_exon_sequences()

        ## Part III: Combine chain and segment data to define projections; adjust
        ## projection and exon coordinates
        self.generate_projection_groups()

        ## for fragmented transcripts, projection groups should be sorted by exon numbers
        ## (zero is a placeholder for fragmented copies)
        if self.fragmented:
            for tr in self.transcripts:
                self.copy2proj[tr][0] = {
                    x:v for x, v in enumerate(
                        sorted(self.copy2proj[tr][0].values(), key=lambda p: p.first_exon)
                    )
                }

        ## at this point, each group has defined boundaries, defined exon loci,
        ## assessed evidence for exon coordinates, and complete grouping for CESAR;
        ## finally, we can proceed to CESAR alignment per each group
        self.prepare_cesar_input()

        ## if no projection for any transcript was found suitable for alignment,
        ## exit without throwing an error
        if not self.transcripts_finished:
            self._to_log('EXIT: All projections were discarded for certain reasons', 'warning')
        else:
            ## all set! to round it all out, dump the CESAR input data
            if self.p_job:
                with FileLock(self.exon_storage_lock, timeout=.5):
                    self.save_results_to_hdf5()
            else:
                self.save_results_to_hdf5()

            ## as well as memory requirements data
            if self.p_job:
                with FileLock(self.memory_lock, timeout=50):
                    self.save_cesar_guide_data()
            else:
                self.save_cesar_guide_data()

        ## write the rejected projection data
        if self.p_job:
            with FileLock(self.rejection_lock, timeout=50):
                self.rejection_report()
        else:
            self.rejection_report()

        ## write the spanning chain coorindates in the reference
        if self.p_job:
            with FileLock(self.spanning_chain_lock, timeout=50):
                self.spanning_chain_report()
        else:
            self.spanning_chain_report()


    def parse_annotation_bed(self, ref_annot: TextIO) -> None:
        """
        From the annotation file, selects the line corresponding to the focal
        transcripts and store its contents in the AnnotationEntry object
        """
        for line in ref_annot.readlines():
            data: List[str] = line.rstrip().split('\t')
            transcript: str = data[3]
            if transcript not in self.transcripts:
                continue
            chrom: str = data[0]
            start: int = int(data[6])
            stop: int = int(data[7])
            strand: str = data[5] == '+'
            exon_num: int = int(data[9])
            exon_sizes: List[int] = list(map(int, data[10].split(',')[:-1]))
            exon_starts: List[int] = list(map(int, data[11].split(',')[:-1]))
            exons: ExonDict[int, Exon] = ExonDict()
            for i in range(exon_num):
                exon_start: int = start + exon_starts[i]
                exon_stop: int = exon_start + exon_sizes[i]
                e_num: int = i + 1 if strand else exon_num - i
                exon: Exon = Exon(
                    e_num,
                    exon_start,
                    exon_stop
                )
                exons[e_num] = exon
            self.annot_entries[transcript] = AnnotationEntry(
                transcript, chrom, start, stop, strand, exon_num, exons
            )
            self._to_log(
                f'Annotation data on {transcript} are successfully uploaded'
            )
            if len(self.annot_entries) == len(self.transcripts):
                return
        if not self.annot_entries:
            self._die(
                'None of the stated transcripts are in the provided '
                'annotation file'
            )

    def parse_u12_file(self, file: str) -> None:
        """
        Parses the HDF5 storage containing U12/non-canonical U2 splice site data
        """
        with h5py.File(file, 'r') as h:
            for tr in self.transcripts:
                if tr not in h:
                    # self.u12_sites[tr] = {}
                    continue
                data: Iterable[bytes] = h[tr][:]
                for intron in data:
                    intron_num: int = int(intron[0].decode('utf8'))
                    intron_class: str = intron[1].decode('utf8')
                    donor_site: str = intron[2].decode('utf8')
                    don_: str = donor_site.lower()
                    acc_site: str = intron[3].decode('utf8')
                    acc_: str = acc_site.lower()
                    if intron_class == U2:
                        divergent_donor: bool = (
                            don_ not in DONOR_SITE and don_ != 'nn'
                        )
                        divergent_acc: bool = (
                            acc_ not in ACCEPTOR_SITE and acc_ != 'nn'
                        )
                    else:
                        divergent_donor: bool = (
                            don_ != DONOR_SITE_U12 and don_ != 'nn'
                        )
                        divergent_acc: bool = (
                            acc_ != ACCEPTOR_SITE_U12 and acc_ != 'nn'
                        )
                    regular_intron: int = int(not (divergent_donor or divergent_acc))
                    regular_donor: int = (
                        int(not divergent_donor) if self.separate_site_treat
                        else regular_intron
                    )
                    regular_acc: int = (
                        int(not divergent_acc) if self.separate_site_treat
                        else regular_intron
                    )
                    self.u12_sites[tr][intron_num]['donor'] = (
                        intron_class, donor_site, regular_donor
                    )
                    self.u12_sites[tr][intron_num+1]['acceptor'] = (
                        intron_class, acc_site, regular_acc
                    )

    def parse_chrom_size_file(self, chrom_file: str, ref: bool = True) -> None:
        """
        Parses a two-column file containing chromosome sizes in nucleotide numbers
        """
        with open(chrom_file, 'r') as h:
            for line in h.readlines():
                data: List[str] = line.rstrip().split('\t')
                chrom: str = data[0]
                size: int = int(data[1])
                if ref:
                    self.ref_chrom_sizes[chrom] = size
                else:
                    self.query_chrom_sizes[chrom] = size

    def parse_segments(self, segments: TextIO) -> None:
        """
        Parses the segment BED12 file
        """
        if segments is None:
            self._to_log('No segment file was provided')
            return
        for segnum, line in enumerate(segments.readlines()):
            data: List[str] = line.strip().split('\t')
            if len(data) != 12:
                self._die(
                    f'Segment BED file contains number of columns not equal to 12 at line {segnum}'
                )
            chrom: str = data[0]
            start: int = int(data[1])
            stop: int = int(data[2])
            strand: str = data[5] == '+'
            id_split: List[str] = data[3].split('|')
            source_proj: str = id_split[1]
            source_trans: str = '#'.join(source_proj.split('#')[:-1])
            if source_trans not in self.transripts:
                continue
            chain: str = source_proj.split('#')[-1]
            if chain not in self.chains:
                continue
            exon_nums: List[int] = list(map(int, id_split[2].split(':')[1].split(',')))
            min_exon: int = min(exon_nums)
            max_exon: int = max(exon_nums)
            exon_sizes: List[int] = list(map(int, data[10].split(',')[:-1]))
            exon_starts: List[int] = list(map(int, data[11].split(',')[:-1]))
            exons: ExonDict[int, Exon] = ExonDict()
            for i in range(len(exon_nums)):
                exon_start: int = start + exon_starts[i]
                exon_stop: int = exon_start + exon_sizes[i]
                exon: Exon = Exon(
                    exon_nums[i],
                    exon_start,
                    exon_stop
                )
                exons[exon_nums[i]] = exon
            self.segments[source_trans] = Segment(
                source_proj, chrom,
                start, stop, strand,
                min_exon, max_exon,
                exons
            )
        if not self.segments:
            self._to_log('WARNING: empty segment file was provided')

    def extract_query_sequence(self, chrom: str, start: int, stop: int, strand: bool) -> str:
        """
        """
        cmd: str = (
            f'twoBitToFa -seq={chrom} -start={start} -end={stop} '
            f'{self.query} stdout'
        )
        res: str = self._exec(cmd, 'ERROR: Query search space extraction failed')
        seq: str = ''.join([x for x in res.split('\n') if x and x[0] != '>'])
        if not strand:
            seq = reverse_complement(seq)
        return seq

    def get_chain_mappers(self) -> None:
        """
        For each chain provided, extract the subchain and convert it to a
        Exon2BlockMapper object
        """
        for chain_id in self.chains:
            self._to_log(f'Retrieving chain data for chain {chain_id}')
            chain_str: str = get_chain(self.chain_file, chain_id)
            # chain_meta: List[str] = chain_str.split('\n')[0].split()
            # print(chain_meta)
            sorted_transcripts: List[str] = sorted(
                self.annot_entries.keys(),
                key=lambda x: (
                    self.annot_entries[x].start, self.annot_entries[x].stop
                )# if self.annot_entries[x].strand else (
                #     self.annot_entries[x].stop, self.annot_entries[x].start
                # )
            )
            transcript_starts: List[int] = [
                min(self.annot_entries[x].start, self.annot_entries[x].stop)
                for x in sorted_transcripts
            ]
            transcript_stops: List[int] = [
                max(self.annot_entries[x].start, self.annot_entries[x].stop)
                for x in sorted_transcripts
            ]
            subchains, chain_meta = transcriptwise_subchains(
                chain_str, sorted_transcripts, transcript_starts, transcript_stops
            )
            # print(f'{subchains=}, {transcript_starts=}, {transcript_stops=}, {chain_meta=}')
            # chain_header = chain_str.split("\n")[0]
            # print(f'{chain_header=}')
            # u = {x:len(y) for x,y in subchains.items()}
            t_chrom: str = chain_meta[0]
            q_chrom: str = chain_meta[1]
            t_chain_start: int = chain_meta[2]
            t_chain_stop: int = chain_meta[3]
            t_strand: bool = chain_meta[4]
            q_chain_start: int = chain_meta[5]
            q_chain_stop: int = chain_meta[6]
            q_strand: bool = chain_meta[7]
            q_chrom_size: int = self.query_chrom_sizes[q_chrom]
            for tr in self.annot_entries:
                self._to_log('Retrieving subchain for transcript %s' % tr)
                codirected: bool = t_strand == q_strand
                subchain_oriented: Dict[str, Tuple[int]] = subchains[tr]
                self._to_log(f'Creating a mapper object for {tr}#{chain_id}')
                mapper: Exon2BlockMapper = intersect_exons_to_blocks(
                    self.annot_entries[tr].exons,
                    subchain_oriented,
                    chain_id,
                    t_chain_start,
                    t_chain_stop,
                    t_chrom,
                    q_chrom,
                    q_chrom_size,
                    t_strand,
                    q_strand,
                    codirected,
                    self.exon_locus_flank,
                    self.exon_locus_flank,
                    EXTRA_FLANK, ## TODO: make this number variable
                    self.logger,
                    self.v
                )
                if mapper.spanning_chain:
                    self.spanning_chains.append(
                        (f'{tr}#{chain_id}', mapper.tchrom, mapper.tstart, mapper.tstop)
                    )
                    if self.disable_spanning_chains:
                        start, stop = sorted((mapper.qstart, mapper.qstop))
                        space_seq: str = self.extract_query_sequence(
                            mapper.qchrom, start, stop, mapper.qstrand
                        )
                        space_gaps: bool = find_gaps(space_seq, self.assembly_gap_size)
                        status: List[str] = MISSING if space_gaps else LOST
                        rej_line: str = SPANNING_CHAIN_REASON.format(f'{tr}#{chain_id}', status)
                        self.rejected_transcripts.append(rej_line)
                        continue
                self.chainmappers[tr][chain_id] = mapper
                self.chain2coords[chain_id] = (q_chain_start, q_chain_stop)

    def resolved_shared_exons(self) -> None:
        """
        For fragmented projections, estimate exons shared between fragments,
        and assign them to the one projection with better coverage
        """
        self._to_log('Attributing exons shared by more than one projection fragment')
        for tr, chains in self.chainmappers.items():
            if len(chains) == 1:
                continue
            mappers: List[Exon2BlockMapper] = [self.chainmappers[tr][x] for x in chains]
            mappers.sort(key=lambda x: (min(x.e2c), max(x.e2c)))
            for i in range(len(mappers) - 1):
                prev_mapper: Exon2BlockMapper = mappers[i]
                next_mapper: Exon2BlockMapper = mappers[i+1]
                shared_exons: List[int] = [
                    x for x in prev_mapper.e2c if x in next_mapper.e2c
                ]
                # print(f'{prev_mapper.e2c=}, {prev_mapper.missing=}')
                # print(f'{next_mapper.e2c=}, {next_mapper.missing=}')
                if not shared_exons:
                    continue
                prev_sole_exon: bool = len(prev_mapper.e2c) == 1
                next_sole_exon: bool = len(next_mapper.e2c) == 1
                prev_exon_assigned_to_next: bool = False
                for x in sorted(shared_exons):
                    if prev_exon_assigned_to_next:
                        del prev_mapper.e2c[x]
                        prev_mapper.missing.add(x)
                        continue
                    cov_in_prev: int = prev_mapper.init_cov.get(x, 0)
                    cov_in_next: int = next_mapper.init_cov.get(x, 0)
                    leave_in_prev: bool = (
                        cov_in_prev > cov_in_next and not next_sole_exon
                    ) or prev_sole_exon
                    if leave_in_prev:
                        del next_mapper.e2c[x]
                        next_mapper.missing.add(x)
                    else:
                        prev_exon_assigned_to_next = True
                        del prev_mapper.e2c[x]
                        prev_mapper.missing.add(x)

    def update_splice_sites(self, tr: str, sites: List[str]) -> None: ## TODO: Needs revision
        """Updates the splice site dictionary"""
        sites_split: List[List[str]] = parts(sites[1:-1], 2)
        for exon, intron in enumerate(sites_split, 1):
            don, acc = intron
            don_, acc_ = map(lambda x: x.lower(), intron)
            divergent_donor: bool = don_ not in DONOR_SITE and don_ != 'nn'
            divergent_acc: bool = acc_ not in ACCEPTOR_SITE and acc_ != 'nn'
            regular_intron: int = int(not (divergent_donor or divergent_acc))
            if not self.u12_sites[tr][exon]['donor']:
                if not self.toga1 and not self.toga1_plus_cesar:
                    self.u12_sites[tr][exon]['donor'] = ('U2', don_, regular_intron)
                else:
                    if regular_intron:
                        self.u12_sites[tr][exon]['donor'] = ('U2', don_, regular_intron)
                    else:
                        self.u12_sites[tr][exon]['donor'] = ('U12', don_, regular_intron)
            if not self.u12_sites[tr][exon + 1]['acceptor']:
                if not self.toga1 and not self.toga1_plus_cesar:
                    self.u12_sites[tr][exon + 1]['acceptor'] = ('U2', acc_, regular_intron)
                else:
                    if regular_intron:
                        self.u12_sites[tr][exon + 1]['acceptor'] = ('U2', acc_, regular_intron)
                    else: ## TODO: REVERT ONCE FINISHED!!!
                        self.u12_sites[tr][exon + 1]['acceptor'] = ('U12', acc_, regular_intron)

    def get_exon_sequences(self) -> None:
        """
        Retrieves exon sequences and splice sites, then formats exon sequences
        and header suitable for CESAR input
        """
        # ref_chroms: Set[str] = {x.chrom for x in self.annot_entries.values()}
        # stderr.write(f'{ref_chroms=}\n')
        # init = datetime.now()
        # ref_2bit: Dict[str, str] = {
        #     k:v[:] for k,v in TwoBitFile(self.ref).items() if k in ref_chroms
        # }
        # twobit_point = datetime.now()
        # stderr.write(f'Reference genome upload time: {str((twobit_point - init).total_seconds())}\n')
        # for tr, entry in self.annot_entries.items():
        #     init = datetime.now()
        #     chrom: str = self.annot_entries[tr].chrom
        #     exon_coordinates, exon_sequences, s_sites, exon_flanks = get_exons(
        #         entry, ref_2bit[chrom], mask_stops=True ## True stands for inframe stop codon masking
        #     )
        #     self.u12_sites[tr] = infer_u12_sites(s_sites, self.u12_sites[tr])
        #     # exon_sequences, sec_codons = check_ref_exons(exon_sequences, True) ## True stands for inframe stop codon masking
        #     # self.prepared_exons[tr] = prepare_exons_for_cesar(exon_sequences)
        #     self.prepared_exons[tr] = exon_sequences
        #     exon_prep_time = datetime.now()
        #     stderr.write(f'Exon preparation time for tr {tr}: {str((exon_prep_time - init).total_seconds())}\n')
        #     init = datetime.now()
        #     self.exon_names_for_cesar[tr] = add_exon_headers(
        #         sorted(entry.exons.keys()),
        #         tr,
        #         self.u12_sites[tr],
        #         self.cesar_common_acceptor,
        #         self.cesar_common_donor,
        #         self.cesar_first_acceptor,
        #         self.cesar_last_donor,
        #         self.cesar_u12_acceptor,
        #         self.cesar_u12_donor
        #     )
        #     header_point = datetime.now()
        #     stderr.write(f'Header addition time for tr {tr}: {str((header_point - init).total_seconds())}\n')
        bed_file: str = ''
        s_site_flanks: Dict[str, Dict[int, Tuple[int]]] = defaultdict(dict)
        for tr, entry in self.annot_entries.items():
            # init2 = datetime.now()
            chrom: str = entry.chrom
            chrom_size: int = self.ref_chrom_sizes[chrom]
            strand: str = '+' if entry.strand else '-'
            for e_num, exon in entry.exons.items():
                start, stop = sorted((exon.start, exon.stop))
                adj_start: int = max(start-2, 0)
                adj_stop: int = min(stop+2, chrom_size)
                l_flank: int = start - adj_start
                r_flank: int = adj_stop - stop
                s_site_flanks[tr][e_num] = (
                    (l_flank, r_flank) if entry.strand else (r_flank, l_flank)
                )
                bed_line: str = (
                    f'{chrom}\t{adj_start}\t{adj_stop}\t{tr}|{e_num}\t0\t{strand}'
                )
                bed_file += bed_line + '\n'
            # tr_extr_point = datetime.now()
            # stderr.write(f'Exon extraction time for transcript {tr}: {str(tr_extr_point - init2)}\n')
        bed_file = bed_file.encode('utf-8')
        res: str = self._run_2bit2fa(self.ref, bed_file)
        parsed_fasta: Dict[str, Dict[str, str]] = parse_extr_exon_fasta(res)
        for tr in self.transcripts:
            exon_subset: Dict[int, str] = parsed_fasta[tr]
            flanks: Dict[int, int] = s_site_flanks[tr]
            prepared_exons, splice_sites = prepare_exons(
                exon_subset, flanks, mask_stops=True
            )
            # self.u12_sites[tr] = infer_u12_sites(splice_sites, self.u12_sites[tr]) ## TODO: Obsolete
            self.update_splice_sites(tr, splice_sites)
            self.prepared_exons[tr] = prepared_exons
            # self.exon_names_for_cesar[tr] = add_exon_headers( ## TODO: Move to a separate method
            #         sorted(self.prepared_exons[tr]),
            #         tr,
            #         self.u12_sites[tr],
            #         self.cesar_common_acceptor,
            #         self.cesar_common_donor,
            #         self.cesar_first_acceptor,
            #         self.cesar_last_donor,
            #         self.cesar_u12_acceptor,
            #         self.cesar_u12_donor
            #     )
        self._to_log('Adding headers to CESAR exon input')
        self.add_exon_headers()

    def get_min_max_exon(
        self,
        chain: Exon2BlockMapper,
        segment: Segment,
        tr: str
    ) -> None:
        """
        Return minimum and maximum exon numbers for exons covered by the
        provided chain-segment group
        """
        mapper_exons: List[int] = [
            x for x in chain.e2c if x not in chain.missing
        ] if chain else []
        segment_exons: List[int] = [
            x for x in segment.exons
        ] if segment is not None else []
        all_exons: List[int] = mapper_exons + segment_exons
        if not all_exons:
            self._to_log('Chain-segment group contains no exons', 'warning')
            return (None, None)
            # exit(0)
        # all_exons = [x for x in all_exons if x in self.universally_missing]
        # print(f'{min(all_exons)=}, {max(all_exons)=}')
        return min(all_exons), max(all_exons)

    def get_group_boundaries(
        self,
        min_exon: int,
        max_exon: int,
        chain_id: str,
        tr: str,
        group: int
    ) -> Tuple[int, int]:
        """
        Defines which exons are expected to be found within the provided group.
        1) If single chain identifier was provided or provided group does not
        contain a chain mapper, all exons are expected to be found.
        2) Otherwise, extend exon boundaries to the last exon absent in all the groups
        """
        # print(f'{self.universally_missing[tr]=}, {self.missing2chain[tr]=}, {max(self.universally_missing[tr]) + 1=}')
        if self.no_inference: ## TODO: Think how properly handle these cases in 'fragmented mode'
            return min_exon, max_exon
        if not self.fragmented or not chain_id:
            return 1, self.annot_entries[tr].exon_number
        if not self.universally_missing[tr][group]:
            self._to_log(f'No universally missing exons for transcript {tr}')
            return min_exon, max_exon
        # print(f'{self.found_exons[tr][group]=}')
        for emin in range(min_exon, max_exon + 1):
            if emin not in self.found_exons[tr][group]:
                break
        min_exon = emin
        for emax in range(max_exon, min_exon - 1, -1):
            if emax not in self.found_exons[tr][group]:
                break
        max_exon = emax
        # print(f'Inside get_group_boundaries, before inference: {min_exon=}, {max_exon=}')
        if min_exon == max_exon and min_exon in self.found_exons[tr][group]:
            return -1, -1
        f, l = 0, 0
        for f in range(min_exon - 1, min(self.universally_missing[tr][group]) - 1, -1):
            if f not in self.universally_missing[tr][group]:
                f += 1
                break
            elif f in self.missing2chain[tr].keys():
                f += 1
                break
            if f in self.universally_missing[tr][group]:
                self.missing2chain[tr][f] = chain_id
        # print(f'{f=}')
        f = f if f else min_exon
        for l in range(max_exon + 1, max(self.universally_missing[tr][group]) + 1):
            if l not in self.universally_missing[tr][group]:
                l -= 1
                break
            elif l in self.missing2chain[tr].keys():
                l -= 1
                break
            if l in self.universally_missing[tr][group]:
                self.missing2chain[tr][l] = chain_id
        # print(f'{l=}')
        l = l if l else max_exon
        # print(f'Inside get_group_boundaries, after inference: {f=}, {l=}')
        return f, l

    def correct_projection_boundaries(
        self,
        proj: Tuple[int, int],
        tr: str,
        group: int
    ) -> Tuple[int]:
        """
        Corrects search boundaries for a given projection
        """
        c_id, s_id = proj
        chain: Exon2BlockMapper = self.chainmappers[tr].get(c_id, None)
        segment: Segment = self.segments[tr].get(s_id, None)
        min_exon, max_exon = self.get_min_max_exon(chain, segment, tr)
        # print(f'{proj=}, {min_exon=}, {max_exon=}')
        if min_exon is None:
            return None
        chrom: str = chain.qchrom if chain else segment.chrom
        chrom_len: int = self.query_chrom_sizes[chrom]
        first_exon, last_exon = self.get_group_boundaries(
            min_exon, max_exon, c_id, tr, group
        )
        # print(f'{proj=}, {min_exon=}, {max_exon=}, {first_exon=}, {last_exon=}')
        if first_exon == -1:
            return None
        for found_ex in range(first_exon, last_exon + 1):
            self.found_exons[tr][group].add(found_ex)
        # print(f'{proj=}, {first_exon=}, {last_exon=}, {min_exon=}, {max_exon=}')
        strand: bool = (
            self.annot_entries[tr].strand == chain.qstrand if chain else segment.strand
        )
        starts, stops = [], []
        if chain:
            starts.append(chain.e2c[min_exon if strand else max_exon][0])
            stops.append(chain.e2c[max_exon if strand else min_exon][1])
        if segment:
            starts.append(
                (segment.exons.min() if strand else segment.exons.max()).start
            )
            stops.append(
                (segment.exons.max() if strand else segment.exons.min()).stop
            )
        init_start: int = min(starts)
        init_stop: int = max(stops)
        upd_start: int = init_start
        upd_stop: int = init_stop
        if min_exon > first_exon and tr not in self.ppgene_list:
            self._to_log('Adjusting the coordinate from the minimum exon side')
            if self.annot_entries[tr].strand:
                min_side_flank: int = (
                    self.annot_entries[tr].exons[min_exon].start -
                    self.annot_entries[tr].exons[first_exon].start
                )
            else:
                min_side_flank: int = (
                    self.annot_entries[tr].exons[first_exon].stop -
                    self.annot_entries[tr].exons[min_exon].stop
                )
            self._to_log(f'Upstream flank size before adjustment: {min_side_flank}')
            min_side_flank = int(min_side_flank * self.extrapolation_modifier)
            self._to_log(f'Upstream flank size before adjustment: {min_side_flank}')
            self._to_log(f'Unadjusted coordinate for the sequence start is: {upd_start if strand else upd_stop}')
            if strand:
                upd_start -= min_side_flank
            else:
                upd_stop += min_side_flank
            self._to_log(f'Adjusted coordinate for the sequence start is: {upd_start if strand else upd_stop}')
        if max_exon < last_exon and tr not in self.ppgene_list:
            self._to_log('Adjusting the coordinate from the maximum exon side')
            if self.annot_entries[tr].strand:
                max_side_flank: int = (
                    self.annot_entries[tr].exons[last_exon].stop -
                    self.annot_entries[tr].exons[max_exon].stop
                )
            else:
                max_side_flank: int = (
                    self.annot_entries[tr].exons[max_exon].start -
                    self.annot_entries[tr].exons[last_exon].start
                )
            max_side_flank = int(max_side_flank * self.extrapolation_modifier)
            if strand:
                upd_stop += max_side_flank
            else:
                upd_start -= max_side_flank
        upd_start = max(0, upd_start)
        upd_stop = max(0, min(upd_stop, chrom_len))
        qstart: int = chain.qstart if chain else init_start
        qstop: int = chain.qstop if chain else init_stop
        return ProjectionCoverageData(
            min_exon, max_exon, first_exon, last_exon, upd_start, upd_stop, 
            init_start, init_stop, qstart, qstop
        )

    def generate_projection_groups(self) -> None:
        """
        For each input reference transcript, resolves chain-inferred and segment
        data to define query projection loci and CESAR input structure.

        Each reference transcript can potentially correspond more than one query
        transcript, which, in its turn, can comprise more than one locus in the
        query genome assembly. The hierarchy of grouping objects goes as follows:
        1. Transcript: reference transcript projected via stated chain(s);
        2. Copy: a single orthologous transcript in the query, occupying one
           (if one chain id was provided) or more (in case of fragmented projections)
           loci in the query assembly;
        3. Group: a combination of chain fragment and/or LASTZ segments describing
           orthologous exon location in a single locus. Usually one copy corresponds
           to one group, except for the case of fragmented projections.

        Groups serve as elementary units of CESAR input organization. Once evidence
        sources are defined for each group, the method identifies search loci for
        each exon belonging to the group and organizes them for further CESAR
        alignment
        """
        ## Loop 1: For each reference transcript, define its projections
        ## (combinations of chain and segment data corresponding to a single query transcript)
        # spliceai_bed_input: Dict[bool, str] = {True: '', False: ''}
        # spliceai_bed_input: Dict[bool, List[Tuple[str, int, int, str]]] = {True: [], False: []}
        spliceai_bed_input: Dict[bool, Dict[Tuple[str, int, int], Tuple[str, int, int]]] = {
            True: {}, False: {}
        }
        for tr in self.transcripts:
            # print(f'BEFORE: {self.groups[tr]=}')
            mappers_grouped: Dict[str, bool] = {
                x:False for x in self.chainmappers[tr]
            }
            segments_grouped: Dict[int, bool] = {
                x:False for x in self.segments[tr]
            }
            group_strand: Dict[str, bool] = {}
            group_marginal_exons: Dict[int, Tuple[int, int]] = {}
            group_num: int = 1
            copy_num: int = 1
            for i, mapper in self.chainmappers[tr].items():
                for j, segment in self.segments[tr].items():
                    ## segment and mapper are considered corresponding if they
                    ## intersect by at least 1bp
                    if intersection(
                        mapper.qstart, mapper.qstop, segment.start, segment.stop
                    ) > 0 and mapper.qchrom == segment.chrom:
                        ## if a chain intersects more than one segment
                        ## and multiple chains were provided, segment data are
                        ## disregarded for this chain
                        ## TODO: Likely to be reassessed; see LiftOn manuscript for details
                        if mappers_grouped[i] and self.fragmented:
                            mappers_grouped[i] = False
                            segments_grouped[j-1] = False
                            group_num -= 1
                            del self.groups[tr][group_num]
                            copy_num  -= 1
                            del self.copy2groups[tr][copy_num]
                            break
                        ## otherwise, create a group of intersecting chain and segment
                        self.groups[tr][group_num] = (i, j)
                        self.copy2groups[tr][group_num].add(copy_num)
                        group_num += 1
                        copy_num += 1
                        mappers_grouped[i] = True
                        segments_grouped[j] = True
                ## chain mappers not intersecting any segment form their own groups
                if not mappers_grouped[i]:
                    self.groups[tr][group_num] = (i, None)
                    if self.fragmented: ## attribute it to the shared fragmented projection
                        self.copy2groups[tr][0].add(group_num)
                    else: ## make a separate projection out of it
                        self.copy2groups[tr][copy_num].add(group_num)
                        copy_num += 1
                    group_num += 1
            ## segments not intersecting any chain mapper form their own groups
            for j, segment in self.segments[tr].items():
                if segments_grouped[j]:
                    continue
                self.groups[tr][group_num] = (None, j)
                group_strand[group_num] = segment.strand
                group_marginal_exons[group_num] = (
                    min(segment.exons), max(segment.exons)
                )
                self.copy2groups[tr][copy_num].add(group_num)
                group_num += 1
                copy_num += 1

            # print(f'AFTER: {self.groups[tr]=}')

            ## track exons missing from both mapper and segment in each group
            for copy, copy_groups in self.copy2groups[tr].items():
                # for g, group in self.groups[tr].items():
                for g in copy_groups:
                    group: Tuple[str, int] = self.groups[tr][g]
                    # if any(map(lambda x: x is None, group)):
                    #     continue
                    # print(f'{group[0]=}')
                    if group[0] is not None:
                        exons_not_in_chain: Set[int] = (
                            self.chainmappers[tr][group[0]].missing
                        )
                    else:
                        exons_not_in_chain: Set[int] = {
                            x for x in range(1, self.annot_entries[tr].exon_number + 1)
                        }
                    if group[1] is not None:
                        exons_not_in_segment: Set[int] = {
                            x for x in range(1, self.annot_entries[tr].exon_number + 1)
                            if x not in self.segments[group[1]].exons.keys()
                        }
                    else:
                        exons_not_in_segment: Set[int] = {
                            x for x in range(1, self.annot_entries[tr].exon_number + 1)
                        }
                    missing_in_projection: Set[int] = exons_not_in_chain.intersection(
                        exons_not_in_segment
                    )
                    # print(f'{missing_in_projection=}')
                    if copy in self.universally_missing[tr] and missing_in_projection:
                        self.universally_missing[tr][copy] = self.universally_missing[tr][copy].intersection(
                            missing_in_projection
                        )
                    elif missing_in_projection:
                        self.universally_missing[tr][copy] = missing_in_projection
                    # print(f'{self.universally_missing[tr][copy]=}')

            ## finally, prepare ProjectionGroup objects containing all the
            ## necessary information for CESAR run
            for copy, copy_groups in self.copy2groups[tr].items():
                valid: bool = True
                for g in copy_groups:
                    group: Tuple[str, int] = self.groups[tr][g]
                    chain: Exon2BlockMapper = self.chainmappers[tr].get(group[0], None)
                    # print(f'{group[0]=}')
                    # print(f'{chain.e2c=}')
                    # print(f'{chain.gap_located=}')
                    # print(f'{chain.missing=}')
                    segment: Segment = self.segments[tr].get(group[1], None)
                    covered_space: ProjectionCoverageData = self.correct_projection_boundaries(group, tr, copy)
                    # print(f'{covered_space=}')
                    ## a botched group with no defined exons has been encountered;
                    ## discard the current copy
                    if covered_space is None:
                        if chain is not None and chain.spanning_chain:
                            # self._die('Chain-segment group contains no exons')
                            # blocks: List[Tuple[int, int, int, int]] = sorted(
                            #     chain.blocks.values(), key=lambda x: x[0]
                            # )
                            # query_start: int = blocks[0][2]
                            # query_end: int = blocks[-1][3]
                            query_start, query_end = sorted((chain.qstart, chain.qstop))
                            if query_end - query_start:
                                seq_in_query: str = self.extract_query_sequence(
                                    chain.qchrom, query_start, query_end, chain.qstrand
                                )
                                # print(f'{seq_in_query=}')
                                self._to_log(
                                    'Checking space gaps in the likely search space of unaligned projection'
                                )
                                space_gaps: bool = find_gaps(seq_in_query, self.assembly_gap_size)
                                print(f'{query_start=}, {query_end=}, {space_gaps=}, {chain.tstart=}, {chain.tstop=}, {chain.qstart=}, {chain.qstop=}')
                                status: str = MISSING if space_gaps else LOST
                            else:
                                self._to_log(
                                    'Zero span in query; marking projection as missing', 'warning'
                                )
                                status: str = LOST
                            # self.spanning_chains.append(
                            #     (f'{tr}#{chain.chainid}', chain.tchrom, chain.tstart, chain.tstop)
                            # )
                        else:
                            status: str = LOST
                        rej_info: Tuple[str, int] = (
                            tr, '0', status, 'No aligned exons found'
                        )
                        self.rejected_transcripts.append(rej_info)
                        # exons_harbored: List[int] = list(chain.e2c.keys())#all_projection_exons((chain, segment))
                        # unique_exons_lost: bool = any(
                        #     x for x in exons_harbored if x not in self.found_exons[tr][copy]
                        # )
                        # # if unique_exons_lost:
                        # #     self.rejection_reasons[tr][copy].append('NO_EXONS_FOUND')
                        # #     rej_info: Tuple[str, int] = (
                        # #         tr, copy, LOST, 'No aligned exons found'
                        # #     )
                        #     self.rejected_transcripts.append(rej_info)
                        valid = False
                        break
                        # else:
                        #     continue
                    # splice_sites: Dict[str, Dict[int, float]] = {'donor': {}, 'acceptor': {}}
                    if self.spliceai_dir:
                        chrom: str = chain.qchrom if chain else segment.chrom
                        strand: bool = (
                            chain.qstrand == self.annot_entries[tr].strand if chain else
                            segment.strand
                        )
                        spliceai_key: Tuple[str, int, int] = (tr, copy, g)
                        spliceai_bed_input[strand][spliceai_key] = (
                            chrom,
                            max(0, covered_space.start - self.exon_locus_flank),
                            min(covered_space.stop + self.exon_locus_flank, self.query_chrom_sizes[chrom]),
                        )
                    ## extract exons lying outside of the chain span
                    out_of_chain_exons: Set[int] = (
                        chain.out_of_chain if chain is not None else set()
                    )
                    ## the rest of the procedure is held within a specialized class instance
                    chrom: str = chain.qchrom if chain is not None else segment.chrom
                    chrom_size: int = self.query_chrom_sizes[chrom]
                    proj: ProjectionGroup = ProjectionGroup(
                        chain,
                        self.annot_entries[tr],
                        group[0] if group[0] else '-1',
                        covered_space.start,
                        covered_space.stop,
                        covered_space.init_start,
                        covered_space.init_stop,
                        covered_space.qstart,
                        covered_space.qstop,
                        chrom_size,
                        covered_space.first_exon,
                        covered_space.last_exon,
                        self.exon_locus_flank,
                        self.exon_locus_flank,
                        out_of_chain_exons,
                        self.max_space_size,
                        self.logger
                    )
                    self.copy2proj[tr][copy][g] = proj
                if not valid:
                    self._to_log(
                        'Copy %i for transcript %s has been discarded' % (copy, tr), 'warning'
                    )
                    if tr in self.copy2proj or copy in self.copy2proj[tr]:
                        del self.copy2proj[tr][copy]
                        if not self.copy2proj[tr]:
                            del self.copy2proj[tr]
                            for strand, keys in spliceai_bed_input.items():
                                items_to_del: List[Tuple[str, int, int]] = []
                                for key in keys:
                                    if key[0] == tr:
                                        items_to_del.append(key)
                                for key in items_to_del:
                                    del spliceai_bed_input[strand][key]
        if self.spliceai_dir:
            self._to_log('Extracting valid SpliceAI predictions')
            self.__extract_splice_sites(spliceai_bed_input)
            self._to_log('SpliceAI data successfully extracted')

    def prepare_cesar_input(self) -> None:
        """
        """
        for tr in self.transcripts:
            segments_finished: int = 0
            for copy, projections in self.copy2proj[tr].items():
                ## check if projection has an overall insufficient coverage in the query
                ## TODO: move to generate_projection_groups()
                if self._insufficient_coverage(projections, tr):
                    msg: str = ( ## TODO: TOOOOOOOO UGLY PUT AS A METHOD FUNCTION
                        'Spanning chain' if
                        any(map(lambda x: x.mapper and x.mapper.spanning_chain, projections.values())) else
                        'Insufficient sequence coverage'
                    )
                    self.rejection_reasons[tr][copy].append('INSUFFICIENT_SEQ_COVERAGE')
                    rej_line: Tuple[str, int] = (tr, copy, LOST, msg)
                    self.rejected_transcripts.append(rej_line)
                    continue
                # if self.disable_spanning_chains:
                # discard_spanning_chain: bool = any(x.spanning_chain for x in projections.values()) and self.disable_spanning_chains
                # if discard_spanning_chain:
                #     self._to_log(
                #         (
                #             'Spanning projection %s will not be aligned due to settings; '
                #             'assessing loss status'
                #         ) % tr,
                #         'warning'
                #     )
                #     start, stop = sorted((projections[0].mapper.qstart, projections[0].mapper.qstop))
                #     space_seq: str = self.extract_query_sequence(
                #         projections[0].chrom, start, stop, projections.strand
                #     )
                #     space_gaps: bool = find_gaps(space_seq, self.assembly_gap_size)
                #     status: List[str] = 'M' if space_gaps else 'L'
                #     chain: str = ','.join(self.chains)
                #     proj_name: str = f'{tr}.{chain}'
                #     rej_line: str = SPANNING_CHAIN_REASON.format(proj_name, status)
                #     self.rejected_transcripts.append(rej_line)
                #     continue
                    
                numerated_exon_groups: Dict[int, int] = {}
                prim_num: int = 0
                groups_finished: int = 0
                for p, proj in projections.items():
                    # print(f'{tr=}, {proj.cesar_exon_grouping=}, {proj.group_coords=}')
                    for e, exon_group in enumerate(proj.cesar_exon_grouping):
                        exp_coords: Dict[int, Tuple[int]] = {}
                        exon_search_spaces: Dict[int, Tuple[int]] = {}
                        gap_located: Set[int] = {
                            x for x in exon_group if x in proj.gap_located
                        }
                        primary_group_num: int = prim_num
                        prim_num += 1
                        exon_group = sorted(exon_group)
                        for ex in exon_group:
                            exp_coords[ex] = (
                                proj.exon_expected_loci[ex][2:]
                                if ex in proj.exon_expected_loci else (None, None)
                            )
                            exon_search_spaces[ex] = (
                                proj.exon_search_spaces[ex][2:]
                                if ex in proj.exon_search_spaces else (None, None)
                            )
                        # group_start: int = proj.exon_coords[exon_group[0]][2]
                        # group_stop: int = proj.exon_coords[exon_group[-1]][3]
                        group_start, group_stop = proj.group_coords[e]
                        rejection_key: Tuple[int] = tuple(exon_group)
                        if group_start > group_stop:
                            self._to_log(
                                f'Group {e} comprising of exons '
                                f'{",".join(map(str, exon_group))} has shifted exons '
                                f'(start={group_start}, end={group_stop}); marking the group as unaligned',
                                'warning'
                            )
                            group_start = group_stop
                        numerated_exon_groups[primary_group_num] = exon_group
                        contains_last_exon: bool = self.annot_entries[tr].exon_number in exon_group
                        out_of_chain_exons: Set[int] = {
                            x for x in exon_group if x in proj.out_of_chain_exons
                        }
                        chrom_len: int = self.query_chrom_sizes[proj.chrom]
                        # print(f'{group_start=}, {group_stop=}, {proj.missing=}, {proj.gap_located=}')
                        # if any(x not in proj.missing and x not in proj.gap_located for x in exon_group):
                        #     self._to_log(f'Adding exon flanks to search space of group {exon_group}')
                        #     group_start = group_start - self.exon_locus_flank
                        #     group_stop = group_stop + self.exon_locus_flank
                        #     self._to_log(f'Updated coordinates for exon group {exon_group} are: {proj.chrom}:{group_start}-{group_stop}')
                        # else:
                        #     self._to_log(f'Exon group {exon_group} contains missing or gap-located exons, coordinates are: {group_start}-{group_stop}')
                        group_start = max(group_start, 0)
                        group_stop = min(group_stop, chrom_len)
                        exon_seqs: List[str] = [self.prepared_exons[tr][x] for x in exon_group]
                        # exon_nums: str = (
                        #     ','.join(map(str, exon_group)) if len(exon_group) <= 10
                        #     else f'{exon_group[0]}-{exon_group[-1]}'
                        # )

                        ## extract query sequence
                        if group_stop - group_start > 1:
                            space_name: str = f'{proj.chrom}:{group_start}-{group_stop}'
                            space_seq: str = self.extract_query_sequence(
                                proj.chrom, group_start, group_stop, proj.strand
                            )
                            self._to_log(f'Checking space gaps in the search space of exon group {exon_group}')
                            space_gaps: bool = find_gaps(space_seq, self.assembly_gap_size)

                            ## calculate the expected memory requirements
                            mem: float = cesar_memory_check(
                                [len(x) for x in exon_seqs],
                                len(space_seq),
                            )
                        else:
                            space_name: str = f'{proj.chrom}:None-None'
                            space_seq: str = ''
                            space_gaps: bool = False
                            mem: float = 0.0
                        self._to_log(f'Memory requirements for group {e} is {mem} GB')

                        ## check if the current group exceeds memory/length caps
                        exceeds_mem_limit: bool = (
                            self.memory_limit is not None and mem > self.memory_limit
                        )
                        exceeds_len_limit: bool = len(space_seq) > self.max_space_size
                        ## check if query sequence is too short to accommodate all the exons
                        query_too_short: bool = (
                            sum((len(x) for x in exon_seqs)) * MIN_REF_LEN_PERC > len(space_seq)
                            if len(exon_seqs) > 1 else
                            len(exon_seqs[0]) * SINGLE_EXON_MIN_REF_LEN_PERC > len(space_seq)
                        )
                        ## extra large exons missing from the chain oftentimes feature 
                        ## lost genes and/or off-target alignments; they can significantly
                        ## slow down intron gain search and are therefore recommended to skip
                        large_exons_are_missing: bool = any(
                            self.annot_entries[tr].exons[x].length() > 2000 ## TODO: Make this threshold adjustable?
                            for x in out_of_chain_exons
                        )
                        ## a workaround for CESAR bug relating to single-base queries
                        single_bp_query: bool = len(space_seq) <= 2
                        valid_for_alignment: bool = not (
                            exceeds_mem_limit or exceeds_len_limit or
                            query_too_short or large_exons_are_missing or 
                            single_bp_query
                        )
                        ## if exon group is not to be aligned due to any reason,
                        ## create a CESAR alignment substitute object
                        if not valid_for_alignment:
                            if exceeds_len_limit or exceeds_mem_limit:
                                unaligned_warn_msg: str = MEM_UNALIGNED_WARNING.format(
                                    e, ",".join(map(str, exon_group)), mem, len(space_seq)
                                )
                            elif large_exons_are_missing:
                                large_and_unaligned: str = ','.join(
                                    map(
                                        str, 
                                        (
                                            x for x in out_of_chain_exons if 
                                            self.annot_entries[tr].exons[x].length() > 2000
                                        )
                                    )
                                )
                                unaligned_warn_msg: str = LARGE_EXON_UNALIGNED_WARNING.format(
                                    e, ",".join(map(str, exon_group)), large_and_unaligned
                                )
                            elif query_too_short or single_bp_query:
                                unaligned_warn_msg: str = SHORT_SPACE_UNALGNED_WARNING.format(
                                    e, ",".join(map(str, exon_group)), len(space_seq)
                                )
                            self._to_log(unaligned_warn_msg, 'warning')
                            ## CESAR output requires a separate storage format
                            ## TODO: Dump dummy alignments in a separate HDF5 storage?
                            space_seq: str = '|'.join(
                                (''.join('-' for bp in range(len(self.prepared_exons[tr][ex]))) for ex in exon_group)
                            )
                            exp_coords: Dict[int, Tuple[None]] = {
                                x: (None, None) for x in exon_group
                            }
                            spliceai_sites: Dict[str, Dict[str, int]] = {'donor': {}, 'acceptor': {}}
                            spliceai_donor_line: str = ''
                            spliceai_acc_line: str = ''
                            if exceeds_mem_limit:
                                rej_reason: str = (
                                    'EXCEEDS_MEMORY+GAP' if space_gaps else 'EXCEEDS_MEMORY'
                                )
                            elif exceeds_len_limit:
                                rej_reason: str = (
                                    'EXCEEDS_SPACE+GAP' if space_gaps else 'EXCEEDS_SPACE'
                                )
                            else:
                                rej_reason: str = (
                                    'INSUFFICIENT_SEARCH_SPACE+GAP' if space_gaps else
                                    'INSUFFICIENT_SEARCH_SPACE'
                                )
                            self.rejection_reasons[tr][copy].append(rej_reason)
                            self.unaligned_exons[tr][copy][rejection_key] = rej_reason
                        else:
                            self.max_est_ram[tr] = max(self.max_est_ram.get(tr, 0.0), mem)
                            if self.max_est_ram[tr] == mem:
                                self.largest_ss[tr] = len(space_seq) ## TODO: record these data copywise?
                                self.largest_exon[tr] = ','.join(
                                    str(len(self.prepared_exons[tr][x])) for x in exon_group
                                ) ## TODO: record these data copywise?
                            self.cumulative_ram[tr] += ceil(mem + 0.1)
                            # spliceai_sites: Dict[str, List[int]] = {
                            #     'donor': {
                            #         x:y for x, y in proj.spliceai_sites['donor'].items() if
                            #         group_start <= x <= group_stop
                            #     },
                            #     'acceptor': {
                            #         x:y for x, y in proj.spliceai_sites['acceptor'].items() if
                            #         group_start <= x <= group_stop
                            #     }
                            # }
                            spliceai_donor_line: str = ','.join(
                                f'{x}:{y}' for x, y in proj.spliceai_sites['donor'].items() if
                                group_start <= x <= group_stop
                            )
                            spliceai_acc_line: str = ','.join(
                                f'{x}:{y}' for x, y in proj.spliceai_sites['acceptor'].items() if
                                group_start <= x <= group_stop
                            )
                            groups_finished += 1
                        ## format exon headers and sequences for HDF5 storage
                        exon_names: str = '|'.join([self.exon_names_for_cesar[tr][x] for x in exon_group])
                        exon_seqs: str = '|'.join(self.prepared_exons[tr][x] for x in exon_group)
                        exon_splice_data: str = self._prepare_splice_data_line(tr, exon_group) ## TODO: For upcoming U12 update

                        ex_nums: str = ','.join(map(str, exon_group))
                        exp_coord_line: str = ','.join(
                            f'{x}:{exp_coords[x][0]}-{exp_coords[x][1]}'
                            for x in exon_group
                        )
                        search_space_line: str = ','.join(
                            f'{x}:{exon_search_spaces[x][0]}-{exon_search_spaces[x][1]}'
                            for x in exon_group
                        )
                        gap_located_line: str = ','.join(
                            str(x) for x in sorted(gap_located)
                        )
                        out_of_chain_line: str = ','.join(
                            str(x) for x in sorted(out_of_chain_exons)
                        )
                        # print(f'{tr=}')
                        # print(f'{exon_group=}')
                        # print(f'{exp_coord_line=}')
                        # print(f'{search_space_line=}')
                        # print(f'{group_start=}, {group_stop=}')

                        hdf5_line: Tuple[str] = (
                            str(copy),
                            str(p),
                            str(primary_group_num),
                            proj.chain,
                            proj.chrom,
                            str(group_start),
                            str(group_stop),
                            '+' if proj.strand else '-',
                            ex_nums,
                            exon_names,
                            exon_seqs,
                            space_seq,
                            exon_splice_data,
                            # acc_u12_status,
                            # donor_u12_status,
                            spliceai_acc_line,
                            spliceai_donor_line,
                            exp_coord_line,
                            search_space_line,
                            gap_located_line,
                            out_of_chain_line,
                            # str(int(len(space_gaps) > 0)),
                            str(int(space_gaps)),
                            str(int(contains_last_exon)),
                            str(self.exon_locus_flank),
                            str(self.exon_locus_flank),
                            str(int(valid_for_alignment)),
                            mem
                        )
                        self.exon_table[tr].append(hdf5_line)

                if not groups_finished:
                    status: str = self._classify_unaligned_projection(tr, copy)
                    rej_line: Tuple[int, str] = (
                        tr, copy, status, 'No exons were aligned'
                    )
                    self.rejected_transcripts.append(rej_line)
                else:
                    segments_finished += 1
            if segments_finished:
                self.transcripts_finished += 1

    def save_results_to_hdf5(self) -> None:
        """
        """
        with h5py.File(self.exon_storage, 'a') as f:
            for tr in self.transcripts:
                ## transcripts for which no copies reached the memory estimation step
                ## will not be aligned anyway
                if tr not in self.max_est_ram:
                    continue
                ## prepare a NumPy array to put into the HDF5 file
                ds_name: str = f'{tr}|{",".join(self.chains)}'
                exon_data: Iterable[Iterable[str]] = array(self.exon_table[tr])
                try:
                    ds = f.create_dataset(
                        ds_name,
                        shape=exon_data.shape,
                        dtype=h5py.string_dtype(encoding='utf-8'),
                        chunks=True
                    )
                    ds[:] = exon_data.astype(str_)
                except ValueError:
                    del f[ds_name]
                    ds = f.create_dataset(
                        ds_name,
                        shape=exon_data.shape,
                        dtype=h5py.string_dtype(encoding='utf-8'),
                        chunks=True
                    )
                    ds[:] = exon_data.astype(str_)

    def save_cesar_guide_data(self) -> None:
        """
        """
        with open(self.memory_file, 'a', buffering=1) as h:
            for tr in self.transcripts:
                ## transcripts for which no copies reached the memory estimation step
                ## will not be aligned anyway
                if tr not in self.max_est_ram:
                    continue
                if (
                    len(self.copy2proj[tr]) > 1
                ):
                    chrom, start, stop, init_start, init_stop, qstart, qstop = (
                        'NA', 'NA', 'NA', 'NA',  'NA', 'NA', 'NA'
                    )
                else:
                    sole_proj: ProjectionGroup = next(
                        y for x in self.copy2proj[tr].values() for y in x.values()
                    )
                    chrom: str = sole_proj.chrom
                    start: int = sole_proj.start
                    stop: int = sole_proj.stop
                    init_start: int = sole_proj.init_start
                    init_stop: int = sole_proj.init_stop
                    # qstart: int = sole_proj.qstart
                    # qstop: int = sole_proj.qstop
                    qstart, qstop = self.chain2coords[sole_proj.chain]
                location: str = Path(self.output).absolute()
                mem_line: str = (
                    f'{tr}\t{",".join(self.chains)}\t{self.max_est_ram[tr]} GB\t'
                    f'{self.cumulative_ram[tr]} GB\t{self.largest_exon[tr]}\t'
                    f'{self.largest_ss[tr]}\t{self.tr2cov[tr]}\t{chrom}\t{start}\t{stop}\t'
                    f'{init_start}\t{init_stop}\t{qstart}\t{qstop}\t{location}'
                )
                h.write(mem_line + '\n')

    def rejection_report(self) -> None:
        """
        Writes a tab-separated report on discarded projections
        """
        with open(self.rejection_file, 'a') as h:
            for x in self.rejected_transcripts:
                if isinstance(x, str):
                    h.write(x + '\n')
                    continue
                tr, proj, status, msg = x
                reason_col: List[str] = sorted(
                    [
                        f'{k}:{v}' for k,v in
                        Counter(self.rejection_reasons[tr][proj]).items()
                    ]
                )
                report: str = PREPROCESSING_REJ.format(
                    f'{tr}#{",".join(self.chains)}',
                    proj,
                    msg,
                    ';'.join(reason_col),
                    status
                )
                h.write(report + '\n')

    def spanning_chain_report(self) -> None:
        """
        Writes the reference coordinates for spanning chain projections
        """
        if not self.spanning_chains:
            return
        with open(self.spanning_chain_file, 'a') as h:
            for proj in self.spanning_chains:
                h.write('\t'.join(map(str, proj)) + '\n')

    def _insufficient_coverage(
        self,
        projections: List[ProjectionGroup],
        tr: str
    ) -> bool:
        """
        Estimates if a projection group defines loci for less than a fraction of
        total transcript length specified by self.min_cov_portion
        """
        total_proj_len: int = 0
        spanning: bool = False
        min_cov_portion: int = int(
            sum(
                [self.annot_entries[tr].exons[x].length() for x in self.annot_entries[tr].exons]
            ) * min(max(0, self.min_cov_portion), 1)
        )
        for p, proj in projections.items():
            exon_length: int = sum(
                [
                    self.annot_entries[tr].exons[x].length()
                    for x in self.annot_entries[tr].exons
                    if x not in proj.missing
                ]
            )
            total_proj_len += exon_length
            spanning =  False if not proj.mapper else (spanning or proj.mapper.spanning_chain) ## TODO: adjust for accommodating segment data
        self.tr2cov[tr] = total_proj_len / sum(
            x.length() for x in self.annot_entries[tr].exons.values()
        )
        spanning_allowed = spanning and not self.disable_spanning_chains
        return (total_proj_len < min_cov_portion) and not spanning_allowed

    def __extract_splice_sites(
        self, bed_dict: Dict[bool, Dict[Tuple[str, int, int], Tuple[str, int, int]]]
    ) -> None:
        """
        Extracts splice sites from the SpliceAI results for all groups, with input
        data organized as BED6 file
        """
        for strand in bed_dict:
            _strand: str = 'Plus' if strand else 'Minus'
            donor_file: str = os.path.join(
                self.spliceai_dir, f'spliceAiDonor{_strand}.bw'
            )
            acc_file: str = os.path.join(
                self.spliceai_dir, f'spliceAiAcceptor{_strand}.bw'
            )
            self.__run_bw2w(acc_file, bed_dict[strand], acceptor=True)
            self.__run_bw2w(donor_file, bed_dict[strand], acceptor=False)

    def _run_2bit2fa(self, file: str, bed: str) -> str:
        """
        """
        cmd = f"twoBitToFa -bed=/dev/stdin {file} stdout"
        res_str: str =  self._exec(cmd, self.bw2w_err, bed)
        return res_str

    def __run_bw2w(
        # self, file: str, bed_file: str, acceptor: bool
        self,
        file: str,
        # bed_file: Tuple[str, int, int, str],
        bed_dict: Dict[Tuple[str, int, int], Tuple[str, int, int]],
        acceptor: bool
    ) -> Dict[str, Dict[int, Dict[int, List[int]]]]:
        """
        Runs bigWigToWig and extracts splice site data for a given file and
        a BED-formatted input
        """

        bed_file: str = '\n'.join(
            map(lambda x: '\t'.join(map(str, x)), bed_dict.values())
        )
        process_line: str = SPLICEAI_PROCESS_SCRIPT#.format(self.min_splice_prob)
        cmd: str = (
            f"{self.bigwig2wig_binary} -bed=/dev/stdin -header {file} stdout | "
            f"{process_line}"
        )
        res_str: str =  self._exec(cmd, self.bw2w_err, bed_file.encode('utf8'))
        key: str = 'acceptor' if acceptor else 'donor'

        # print(f'{self.copy2proj=}')
        for line in res_str.split('\n'):
            # print(f'{line=}')
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            chrom_: str = data[0]
            pos: int = int(data[1])
            prob: float = float(data[2])
            for (tr, c, g), (chrom, start, stop) in bed_dict.items():
                if chrom == chrom_ and start <= pos <= stop:
                    # print(f'{chrom=}, {start=}, {stop=}, {tr=}, {c=}, {g=}')
                    self.copy2proj[tr][c][g].spliceai_sites[key][pos] = prob

    def _classify_unaligned_projection(self, tr: str, copy: int) -> str:
        """
        For a transcript with none exons being aligned, defines its loss status
        """
        cds_length: int = sum(x.length() for x in self.annot_entries[tr].exons.values())
        lost_portion: int = 0
        missing_portion: int = 0
        for exon_group, reason in self.unaligned_exons[tr][copy].items():
            missing: bool = '+GAP' in reason
            for exon in exon_group:
                exon_length: int = self.annot_entries[tr].exons[exon].length()
                if exon_length / cds_length > 0.4:
                    return MISSING if missing else LOST
                if missing:
                    missing_portion += exon_length
                else:
                    lost_portion += exon_length
        return MISSING if missing_portion > lost_portion else LOST

    def _prepare_splice_data_line(self, tr: str, exons: List[int]) -> None:
        """
        Converts data on splice site status into a string for further dumping
        """
        exon_num: int = self.annot_entries[tr].exon_number
        out_string: str = ''
        for i, exon in enumerate(exons):
            if exon != 1:
                acc_data: Tuple[str, int] = self.u12_sites[tr][exon]['acceptor']
                acc_class, acc_site = acc_data[:2]
                # prev_intron: Tuple[str] = self.u12_sites[tr][exon-1]
                # acc_class: str = prev_intron[0]
                # acc_site: str = prev_intron[2]
            else:
                acc_class, acc_site = None, None
            if exon != exon_num:
                donor_data: Tuple[str, int] = self.u12_sites[tr][exon]['donor']
                donor_class, donor_site = donor_data[:2]
                # next_intron: Tuple[str] = self.u12_sites[exon-1]
                # donor_class: str = next_intron[0]
                # donor_site: str = next_intron[1]
            else:
                donor_class, donor_site = None, None
            out_string += f'{exon}:{acc_class},{acc_site},{donor_class},{donor_site}'
            if i != len(exons) - 1:
                out_string += ';'

        # for intron in range(1, exon_num):
        #     if intron not in self.u12_sites[tr]:
        #         continue
        #     intron_class, donor_dinuc, acc_dinuc = self.u12_sites[intron]
        #     out_string += f'{intron}:{intorn_class},{donor_dinuc},{acc_dinuc}'
        #     if intron < exon_num - 1:
        #         out_string += ';'
        return out_string

    def add_exon_headers(self) -> None:
        """
        """
        for tr in self.transcripts:
            exon_nums: List[int] = sorted(self.prepared_exons[tr])
            last_exon: int = self.annot_entries[tr].exon_number
            for num in exon_nums:
                header: str = f'>{tr}_exon{num} '
                if num == 1:
                    header += f'\t{self.cesar_first_acceptor}'
                elif self.u12_sites[tr][num]['acceptor']:  ## TODO: Add an auxiliary method
                    if self.u12_sites[tr][num]['acceptor'][0] == U12:
                        if self.u12_sites[tr][num]['acceptor'][2]:
                            header += f'\t{self.canon_u12_acceptor}'
                        else:
                            header += f'\t{self.non_canon_u12_acceptor}'
                    else:
                        if self.u12_sites[tr][num]['acceptor'][2]:
                            header += f'\t{self.canon_u2_acceptor}'
                        else:
                            header += f'\t{self.non_canon_u2_acceptor}'
                else:
                    header += f'\t{self.canon_u2_acceptor}'

                if num == last_exon:
                    header += f'\t{self.cesar_last_donor}'
                elif self.u12_sites[tr][num]['donor']:
                    if self.u12_sites[tr][num]['donor'][0] == U12:
                        if self.u12_sites[tr][num]['donor'][2]:
                            header += f'\t{self.canon_u12_donor}'
                        else:
                            header += f'\t{self.non_canon_u12_donor}'
                    else:
                        if self.u12_sites[tr][num]['donor'][2]:
                            header += f'\t{self.canon_u2_donor}'
                        else:
                            header += f'\t{self.non_canon_u2_donor}'
                else:
                    header += f'\t{self.canon_u2_donor}'
                if num == 1 and num == last_exon and self.toga1:
                    header = f'>{tr}_exon{num}'
                self.exon_names_for_cesar[tr][num] = header



if __name__ == '__main__':
    CesarPreprocessor()
