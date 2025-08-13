#!/usr/bin/env python3

"""
Master script for TOGA2
"""

from src.python.modules.cesar_wrapper_constants import (
    DEF_BLOSUM_FILE, MIN_ASMBL_GAP_SIZE,
    HG38_CANON_U2_ACCEPTOR, HG38_CANON_U2_DONOR,
    HG38_NON_CANON_U2_ACCEPTOR, HG38_NON_CANON_U2_DONOR,
    HG38_CANON_U12_ACCEPTOR, HG38_CANON_U12_DONOR,
    HG38_NON_CANON_U12_ACCEPTOR, HG38_NON_CANON_U12_DONOR,
    FIRST_ACCEPTOR, LAST_DONOR
)
from src.python.modules.constants import TOGA2_EPILOG, Constants
from src.python.modules.shared import CONTEXT_SETTINGS, PrettyGroup
from src.python.modules.toga_main import TogaMain
from typing import Optional

import click
import logging
import os

__author__ = 'Yury V. Malovichko'
__version__ = '2.0.3'
__year__ = '2024'
__credits__ = ('Bogdan M. Kirilenko', 'Michael Hiller')

logging.basicConfig(level=logging.INFO)

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
BLOSUM_FILE: str = os.path.join(LOCATION, *DEF_BLOSUM_FILE)


input_options: PrettyGroup = PrettyGroup(
    'Additional input'
)
control_flow_options: PrettyGroup = PrettyGroup(
    'Pipeline', help='Control flow settings'
)
extraction_options: PrettyGroup = PrettyGroup(
    'Feature extraction', help='Settings for the feature extraction step'
)
class_options: PrettyGroup = PrettyGroup(
    'Classification', help='Projection classification settings'
)
gene_select_options: PrettyGroup = PrettyGroup(
    'Query gene selection', help='Controls orthology/completeness classes of query projections to annotate'
)
prepr_options: PrettyGroup = PrettyGroup(
    'Preprocessing', help='Data preprocessing for CESAR alignment'
)
cesar_options: PrettyGroup = PrettyGroup(
    'CESAR alignment', help='Exon alignment & gene annotation with CESAR'
)
parallel_options: PrettyGroup = PrettyGroup(
    'Parallel execution', help='Execution parameters for the pipeline\'s parallel steps'
)
spliceai_options: PrettyGroup = PrettyGroup(
    'SpliceAI use', help=(
        'SpliceAI use for exon annotation, splice site correction, and intron gain search'
    )
)
annot_options: PrettyGroup = PrettyGroup(
    'Annotation', help='Post-CESAR gene annotation & mutation check settings'
)
loss_options: PrettyGroup = PrettyGroup(
    'Gene loss', help='Gene conservation/loss classification settings'
)
orth_options: PrettyGroup = PrettyGroup(
    'Orthology resolution', help=(
        'Orthology resolution settings, including the gene tree-based orthology refinement'
    )
)
browser_options: PrettyGroup = PrettyGroup(
    'UCSC browser', help='UCSC genome browser report parameters'
)
utr_options: PrettyGroup = PrettyGroup(
    'UTR annotation', help='Settings for the UTR annotation module'
)
legacy_and_experimental: PrettyGroup = PrettyGroup(
    'Legacy & experimental features'
)
verbosity_options: PrettyGroup = PrettyGroup(
    'Verbosity', help='Verbosity & notifications controls'
)
out_options: PrettyGroup = PrettyGroup('Output')
misc_options: PrettyGroup = PrettyGroup('Miscellaneous')


@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
def toga2() -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.     .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM  
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM  
         MM     MM        MM MM             ,M  `MM          MM    MM  
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM  
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM  
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    TOGA2 - Tool for Ortholog Inference from Genome Alignment
    """

@toga2.command(
    context_settings=CONTEXT_SETTINGS, 
    no_args_is_help=True, 
    epilog=TOGA2_EPILOG,
    short_help='Run TOGA2 pipeline with command line arguments'
)
@click.argument(
    'ref_2bit',
    type=click.Path(exists=True),
    metavar='REF_2BIT'
)
@click.argument(
    'query_2bit',
    type=click.Path(exists=True),
    metavar='QUERY_2BIT'
)
@click.argument(
    'chain_file',
    type=click.Path(exists=True),
    metavar='ALIGNMENT_CHAINS'
)
@click.argument(
    'ref_annotation',
    type=click.Path(exists=True),
    metavar='REF_ANNOTATION_BED'
)
@control_flow_options.option(
    '--resume_from',
    '-res',
    type=click.Choice(Constants.RESUME_OPTIONS, case_sensitive=False),
    metavar='STEP',
    default='all',
    show_default=True,
    help=(
        """If you have an unfinished run and want to resume with the same results, 
        select the step from the following list:\b
        all: a placeholder for full starts the pipeline from the very beginning;\b\n
        setup: input data filtering, indexing, and format conversion;\b\n
        feature_extraction: projection and chain feature extraction for projection classification;\b\n
        classification: projection classification in terms of orthology;\b\n
        preprocessing: projection data preprocessing for further alignment; parallel step;\b\n
        aggregate_preprocessing_res: data aggregation and summary across independent preprocessing batches;\b\n
        alignment: CESAR alignment, mutation check, and loss inference for each projection; parallel step;\b\n
        aggregate_cesar_res: data aggregation and summary across independent alignment batches;\b\n
        loss_summary: gene loss data summary;\b\n
        orthology: orthology relationship resolution; if "-st" flag is set, gene tree orthology batches are run at this step;\b\n
        summarize_trees: if "-st" flag was set, individual gene tree batch results are summarized and added to the original orthology data at this step;\b\n
        finalize: gene renaming and output Bed file filtering;\b\n
        ucsc_report: BigBed file preparation\n"""
    )
)
@control_flow_options.option(
    '--halt_at',
    '-halt',
    type=click.Choice(Constants.RESUME_OPTIONS, case_sensitive=False),
    metavar='STEP',
    default='all',
    show_default=True,
    help=(
        'Halts the pipeline at the selected step (see above). Option "all" '
        'implies running the pipeline to the last step'
    )
)
@control_flow_options.option(
    '--selected_feature_batches',
    '-feat_b',
    type=str,
    metavar='COMMA_SEPARATED_LIST',
    default=None,
    show_default=True,
    help=(
        'A comma-separated list of batch numbers for the feature extraction step to be run. '
        'Valid only if --resume_from is set to "feature_extraction" or lower, and '
        '--legacy_chain_feature_extraction option was set in both interrupted and scheduled runs'
    )
)
@control_flow_options.option(
    '--selected_preprocessing_batches',
    '-prep_b',
    type=str,
    metavar='COMMA_SEPARATED_LIST',
    default=None,
    show_default=True,
    help=(
        'A comma-separated list of batch numbers for the CESAR preprocessing step to be run. '
        'Valid only if --resume_from is set to "preprocessing" or lower'
    )
)
@control_flow_options.option(
    '--selected_alignment_batches',
    '-aln_b',
    type=str,
    metavar='COMMA_SEPARATED_LIST',
    default=None,
    show_default=True,
    help=(
        'A comma-separated list of batch numbers for the CESAR alignment step to be run. '
        'Valid only if --resume_from is set to "alignment" or lower'
    )
)
@control_flow_options.option(
    '--no_utr_annotation',
    '-no_utr',
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, UTR sequences are not added to the final annotation file'
)
@input_options.option(
    '--isoform_file',
    '-i',
    type=click.Path(exists=True),
    metavar='ISOFORMS_FILE',
    default=None,
    show_default=True,
    help='A path to a two-column tab-separated file containing gene-to-isoform mapping'
)
@input_options.option(
    '--u12_file',
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
@input_options.option(
    '--spliceai_dir',
    '-sai',
    type=click.Path(exists=True),
    metavar='SPLICEAI_OUT_DIR',
    help='A path to the SpliceAI pipeline output directory'
)
@extraction_options.option(
    '--min_chain_score',
    '-mcs',
    type=click.IntRange(min=0, max=None),
    metavar='INT',
    default=5000,
    show_default=True,
    help=(
        'Minimal score for chains to be considered for classification. Setting '
        'this value to zero disables chain filtering'
    )
)
@extraction_options.option(
    '--min_orthologous_chain_score',
    '-minscore',
    type=click.IntRange(min=0, max=None),
    metavar='INT',
    default=15000,
    show_default=True,
    help=(
        'Minimal score for chains to be potentially classified as orthologous. Chains with '
        'min_chain_score <= X < min_orthologous_chain_score are discarded unless they are '
        'classified as retrogenes/processed pseudogenes'
    )
)
@extraction_options.option(
    '--feature_jobs',
    '-fj',
    type=int,
    metavar='INT',
    default=100,
    show_default=True,
    help='A number of jobs to projection feature extraction commands into'
)
@class_options.option(
    '--orthology_threshold',
    '-ot',
    type=float,
    metavar='FLOAT',
    default=0.5,
    show_default=True,
    help='Probability threshold for considering projections as orthologous'
)
@class_options.option(
    '--single_exon_model',
    '-se_model',
    type=click.Path(exists=True),
    metavar='PATH',
    default=os.path.join('models', 'se_model.dat'),
    show_default=True,
    help=(
        'A path to a orthology classification model for single exon reference transcripts'
    )
)
@class_options.option(
    '--multi_exon_model',
    '-me_model',
    type=click.Path(exists=True),
    metavar='PATH',
    default=os.path.join('models', 'me_model.dat'),
    show_default=True,
    help=(
        'A path to a orthology classification model for multi-exon reference transcripts'
    )
)
@class_options.option(
    '--use_long_distance_model',
    '-use_ld',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, applies extra classifier for distantly related species; '
        'relevant at molecular distances >=1sps'
    )
)
@class_options.option(
    '--long_distance_model',
    '-ld_model',
    type=click.Path(exists=True),
    metavar='PATH',
    default=os.path.join('models', 'se_model.dat'),
    show_default=True,
    help=(
        'A path to a refinement classification model for distantly related '
        'reference-query species pairs'
    )
)
@gene_select_options.option(
    '--disable_fragment_assembly',
    '-no_f',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, does not attempt to recover fragmented projections '
        'from individual chains'
    )
)
@gene_select_options.option(
    '--orthologs_only',
    '-r',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, only orthologous projections are considered'
)
@gene_select_options.option(
    '--one2ones_only',
    '-o2o',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, only transcript with a single orthologous projection '
        'are considered'
    )
)
@gene_select_options.option(
    '--enable_spanning_chains',
    '-nospan',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, spanning chains (i.e., chains with alignment gap corresponding '
        'to the projected transcript) are considered for CESAR alignment; otherwise '
        'spanning chains are used only to discriminate between Lost and Missing projections'
    )
)
@gene_select_options.option(
    '--annotate_processed_pseudogenes',
    '-pp',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, processed pseudogene projections are considered for CESAR alignment; '
        'otherwise, a separate BED9 track "processed_pseudogenes.bed" is added to output directory'
    )
)
@prepr_options.option(
    '--preprocessing_jobs',
    '-pj',
    type=int,
    metavar='INT',
    default=300,
    show_default=True,
    help='A number of jobs to partition CESAR preprocessing commands into'
)
@prepr_options.option(
    '--max_chains_per_transcript',
    '-mc',
    type=int,
    metavar='INT',
    default=100,
    show_default=True,
    help=(
        'A maximum number of chains to project each transcript through. If the '
        'number of projections exceeds the given value, only the first N chains '
        'will be considered'
    )
)
@prepr_options.option(
    '--max_search_space_size',
    '-mss',
    type=int,
    metavar='INT (BP)',
    default=1_000_000,
    show_default=True,
    help=(
        'Query sequence length limit for CESAR jobs. Projections in which any '
        'exon group aligns to sequence beyond this limit will be discarded '
        'from alignment step'
    )
)
@prepr_options.option(
    '--extrapolation_modifier',
    '-em',
    type=float,
    metavar='FLOAT',
    default=1.2,
    show_default=True,
    help=(
        'Modifier by which counterparts of missing reference sequence are modified '
        'when extrapolating ambiguous projection termini'
    )
)
@prepr_options.option(
    '--minimal_covered_fraction',
    '-mincov',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.0,
    show_default=True,
    help=(
        'Minimal fraction of reference CDS sequence covered by aligned chain '
        'blocks to consider the projection for CESAR alignment. If the '
        '--enable_spanning_chains flag is set, this parameter is automatically '
        'set to 0.0'
    )
)
@prepr_options.option(
    '--exon_locus_flank',
    '-ef',
    type=int,
    metavar='INT (BP)',
    default=100,
    show_default=True,
    help=(
        'Size of flank to be added to each defined exon locus from both sides. '
    )
)
@prepr_options.option(
    '--assembly_gap_size',
    '-gs',
    type=int,
    metavar='INT',
    default=MIN_ASMBL_GAP_SIZE,
    show_default=True,
    help=(
        'Minimum number of consecutive N symbols to be considered an assembly gap'
    )
)
@cesar_options.option(
    '--cesar_binary',
    '-cs',
    type=click.Path(exists=True),
    metavar='CESAR_BINARY',
    default=None,
    show_default=False,
    help=(
        'A path to the actual CESAR2.0 binary; if not provided, will look for one '
        'in the PATH, otherwise defaulting to the CESAR2.0 instance '
        'in the current directory'
    )
)
@cesar_options.option(
    '--memory_bins',
    '-b',
    type=str,
    metavar='BIN_LIST',
    default='3,5,7,10,15',
    show_default=True,
    help=(
        'A comma-separated list of memory bin caps, in GB. For each memory bin, '
        'job scheduling will be performed independently. If you want to process '
        'memory-intensive projections as a single cluster call, set the last value '
        'to "big" or enable the --allow_heavy_jobs flag'
    )
)
@cesar_options.option(
    '--job_nums_per_bin',
    '-jb',
    type=str,
    metavar='BIN_JOB_NUM_LIST',
    default='500,20,20,10,5',
    show_default=True,
    help=(
        'A comma-separated list of job numbers per memory bin. '
        'Job jumbers must follow in the same order as memory caps passed '
        'to --memory_bins option'
    )
)
@cesar_options.option(
    '--allow_heavy_jobs',
    '-ahj',
    type=bool,
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'Aggregate all jobs exceeding the highest memory cap as a single joblist; '
        'if memory bins are provided, duplicates the "big" memory bin behavior'
    )
)
@cesar_options.option(
    '--cesar_memory_limit',
    '-ml',
    type=float,
    metavar='FLOAT (GB)',
    default=15.0,
    show_default=True,
    help=(
        'Upper memory limit for CESAR jobs. Projections in which any exon group '
        'requires memory beyond this limit will be discarded before the alignment step. '
        'Value of 0 denotes unlimited memory'
    )
)
@cesar_options.option(
    '--cesar_canon_u2_acceptor',
    '-cca',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_CANON_U2_ACCEPTOR,
    show_default=True,
    help='A path to canonical (GT/GC-AG) U2 acceptor profile'
)
@cesar_options.option(
    '--cesar_canon_u2_donor',
    '-ccd',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_CANON_U2_DONOR,
    show_default=True,
    help='A path to canonical (GT/GC-AG) U2 donor profile'
)
@cesar_options.option(
    '--cesar_non_canon_u2_acceptor',
    '-cnca',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_NON_CANON_U2_ACCEPTOR,
    show_default=True,
    help='A path to non-canonical (non GT/GC-AG) U2 acceptor profile'
)
@cesar_options.option(
    '--cesar_non_canon_u2_donor',
    '-cncd',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_NON_CANON_U2_DONOR,
    show_default=True,
    help='A path to non-canonical (non GT/GC-AG) U2 donor profile'
)
@cesar_options.option(
    '--cesar_canon_u12_acceptor',
    '-cua',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_CANON_U12_ACCEPTOR,
    show_default=True,
    help='A path to canonical (GT-AG) U12 exon acceptor profile'
)
@cesar_options.option(
    '--cesar_canon_u12_donor',
    '-cud',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_CANON_U12_DONOR,
    show_default=True,
    help='A path to canonical (GT-AG) U12  donor profile'
)
@cesar_options.option(
    '--cesar_non_canon_u12_acceptor',
    '-cnua',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_NON_CANON_U12_ACCEPTOR,
    show_default=True,
    help='A path to non-canonical (non-GT-AG) U12 exon acceptor profile'
)
@cesar_options.option(
    '--cesar_non_canon_u12_donor',
    '-cnud',
    type=click.Path(exists=True),
    metavar='PATH',
    default=HG38_NON_CANON_U12_DONOR,
    show_default=True,
    help=(
        'A path to non-canonical (non-GT-AG) U12 exon donor profile'
    )
)
@cesar_options.option(
    '--cesar_first_acceptor',
    '-cfa',
    type=click.Path(exists=True),
    metavar='PATH',
    default=FIRST_ACCEPTOR,
    show_default=True,
    help='A (relative to CESAR2 location) path to first exon acceptor profile'
)
@cesar_options.option(
    '--cesar_last_donor',
    '-cld',
    type=click.Path(exists=True),
    metavar='PATH',
    default=LAST_DONOR,
    show_default=True,
    help='A (relative to CESAR2 location) path to last exon donor profile'
)
@cesar_options.option(
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
@spliceai_options.option(
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
@spliceai_options.option(
    '--spliceai_correction_mode',
    '-scm',
    type=click.IntRange(min=0, max=7),
    metavar='MODE_NUM',
    default=0,
    show_default=False,
    help=(
        """Set the mode of SpliceAI-mediated exon boundary correction:\b\n
        0 - no correction [default; equivalent to not providing SpliceAI data directory];\b\n
        1 - use SpliceAI predictions to correct boundaries of missing and deleted exons;\b\n
        2 - correct mutated canonical U2 splice sites;\b\n
        3 - correct all canonical U2 splice sites in the presence of alternatives with higher SpliceAI support;\b\n
        4 - correct all canonical U2 as well as mutated GT-AG U12 splice sites;\b\n
        5 - correct all canonical U2 and mutated and/or unsupported  GT-AG U12 splice sites;\b\n
        6 - correct all canonical U2 and all U12 splice sites;\b\n
        7 - correct all U2 (including known non-canonical sites) and U12 splice sites.\b\n
        Current recommendation (at least, for vertebrate queries) is 3"""
    )
)
@spliceai_options.option(
    '--min_splice_prob',
    '-msp',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.02,
    show_default=True,
    help='Minimum SpliceAI prediction probability to consider the splice site'
)
@spliceai_options.option(
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
@spliceai_options.option(
    '--intron_gain_check',
    '-ig',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, performs SpliceAI-guided check for query-specific introns'
    )
)
@spliceai_options.option(
    '--intron_gain_threshold',
    '-igt',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.5,
    show_default=True,
    help='Minimal intron gain threshold to consider'
)
@spliceai_options.option(
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
@spliceai_options.option(
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
@spliceai_options.option(
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
@spliceai_options.option(
    '-max_intron_number',
    '-mnn',
    type=click.IntRange(min=1, max=None),
    metavar='INT',
    default=4, 
    show_default=True,
    help=(
        'Maximum gained intron number per exon. Highly recommended to increase this beyond 5-6'
    )
)
@annot_options.option(
    '--matrix',
    '-m',
    type=click.Path(exists=True),
    metavar='BLOSUM_MATRIX_FILE',
    default=BLOSUM_FILE,
    show_default=True,
    help='A file containing the protein alignment matrix'
)
@annot_options.option(
    '--mask_n_terminal_mutations',
    '-m10m',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, masks mutations occurring in the first 10 percents '
        'of query projection length regardless of alternative start codon presence'
    )
)
@annot_options.option(
    '--disable_missing_stop_search',
    '-rmo',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, TOGA2 does not search for alternative stop codons '
        'downstream to CESAR alignment boundary if the original alignment '
        'does not end with one'
    )
)
@loss_options.option(
    '--accepted_loss_symbols',
    '-l',
    type=str,
    metavar='LOSS_SYMBOLS',
    default=','.join(Constants.DEFAULT_LOSS_SYMBOLS),
    show_default=True,
    help=(
        'A comma-separated list of loss status symbols; only projections of '
        'respective statuses will be considered when creating a connection graph. '
        'Keyword ALL lets all possible statuses in.'
    )
)
@orth_options.option(
    '--skip_gene_trees',
    '-st',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, skips gene-tree based resolution step for convoluted many2many clades'
    )
)
@orth_options.option(
    '--use_raxml',
    '-raxml',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'Use raxmlHPC-PTHREADS-AVX instead of IqTree2 for tree inference'
    )
)
@orth_options.option(
    '--max_clique_size',
    '-mqs',
    type=int,
    metavar='INT',
    default=50,
    show_default=True,
    help=(
        'A maximum number of sequences in many:many cliques to be resolved with '
        'gene trees'
    )
)
@orth_options.option(
    '--orthology_jobs',
    '-oj',
    type=int,
    metavar='INT',
    default=100,
    show_default=True,
    help=(
        'A number of jobs to split orthology fine resolution commands into'
    )
)
@orth_options.option(
    '--prank_binary',
    '-pb',
    type=click.Path(exists=True),
    metavar='PRANK_BINARY',
    default=None,
    show_default=True,
    help=(
        'A path to the PRANK executable to be used at fine resolution step. '
        'If not provided, the program will try to infer its location from the PATH'
    )
)
@orth_options.option(
    '--tree_binary',
    '-rb',
    type=click.Path(exists=True),
    metavar='TREE_BINARY',
    default=None,
    show_default=True,
    help=(
        'A path to the IqTree2/raxmlHPC-PTHREADS-AVX executable to be used at fine resolution step. '
        'If not provided, the program will try to infer its location from the PATH.'
    )
)
@orth_options.option(
    '--tree_cpus',
    '-rc',
    type=int,
    metavar='INT',
    default=1,
    show_default=True,
    help='A maximum number of CPUs to run IqTree2/RAxML with'
)
@utr_options.option(
    '--utr_abs_threshold',
    '-utr_abs',
    type=click.IntRange(min=0),
    default=3000,
    show_default=True,
    help='Absolute threshold by which the projected UTR block/exon can exceed the reference counterpart'
)
@utr_options.option(
    '--utr_rel_threshold',
    '-utr_rel',
    type=click.FloatRange(min=0.0),
    default=2.5,
    show_default=True,
    help=(
        'Relative (to reference exon length) threshold by which the projected UTR block/exon '
        'can exceed the reference counterpart'
    )
)
@utr_options.option(
    '--no_utr_boundary_extrapolation',
    '-no_utr_extra',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, UTR block/exon boundaries will be projected by the last aligned UTR position, '
        'without extrapolating the terminal unaligned regions'
    )
)
@utr_options.option(
    '--no_adjacent_utr_extra',
    '-no_adj_utr',
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, unaligned CDS-adjacent UTR sequences will not be extrapolated'

)
@utr_options.option(
    '--fixed_adjacent_utr_extra',
    '-fixed_adj_utr',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, unaligned CDS-adjacent UTR sequences will be annotated as block of fixed length '
        'equal to --utr_abs_threshold value'
    )
)
@browser_options.option(
    '--link_file',
    '-lf',
    type=click.Path(exists=True),
    metavar='FILE',
    default=None,
    show_default=False,
    help=(
        'A path to the two-column file containing HTML-formatted hyperlinks '
        'to external sources on reference transcripts'
    )
)
@browser_options.option(
    '--ucsc_prefix',
    '-up',
    type=str,
    metavar='PREFIX',
    default=Constants.DEFAULT_UCSC_PREFIX,
    show_default=True,
    help='A prefix to use in the output file names'
)
@browser_options.option(
    '--bedtobigbed_binary',
    type=click.Path(exists=True),
    metavar='BEDTOBIGBED_PATH',
    default=None,
    show_default=True,
    help=(
        'A path to UCSC bedToBigBed executable; if not set, the executable with this name '
        'will be sought for in PATH'
    )
)
@parallel_options.option(
    '--parallel_strategy',
    '-c',
    type=click.Choice(Constants.ALL_PARALLEL_EXECS),
    metavar='PARALLEL_EXECUTOR',
    default='local',
    show_default=True,
    help=(
        """Specify the HPC strategy. By default, TOGA2 uses Nextflow to handle 
        parallel processes and supports, at least in theory, all Nextflow executors. 
        Please consult the full list of options at Nextflow help page. Note that setting 
        executor to "local" will parallel the processes over the local machine CPUs.\b\n
        If you want to use Parasol as parallel process manager, set this option to "para"\b
        If you want to implement a fully custom parallel manager strategy, modify the CustomStrategy class 
        and set this option to "custom" or contact the TOGA2 team
        """
    )
)
@parallel_options.option(
    '--nextflow_exec_script',
    type=click.Path(exists=True),
    metavar='NEXTFLOW_SCRIPT_PATH',
    default=None,
    show_default=True,
    help=(
        'A path to a user-defined Nextflow script used for parallel process execution. '
        'If not specified, TOGA2 will use a minimal boilerplate script instead. '
        'Ignored if Parasol or custom HPC strategy were specified for parallel steps'
    )
)
@parallel_options.option(
    '--max_number_of_retries',
    type=click.IntRange(min=1),
    metavar='INT',
    default=3,
    show_default=True,
    help=(
        'Maximum number of retries per parallel job before reporting job failure. '
        'Ignored if Parasol or custom HPC strategy were specified for parallel steps'
    )
)
@parallel_options.option(
    '--nextflow_config_dir',
    '-nc',
    type=click.Path(exists=True),
    metavar='NEXTFLOW_CONFIG_DIR',
    default=None,
    show_default=True,
    help=(
        'A path to a directory containing user-defined Nextflow configuration files; '
        'for detauls, please see nextflow_config_files/readme.txt'
    )
)
@parallel_options.option(
    '--max_parallel_time',
    '-max_t',
    type=click.IntRange(min=1),
    metavar='HOURS',
    default=24,
    show_default=True,
    help='Maximum time duration (in hours) for Nextflow parallel processes'
)
@parallel_options.option(
    '--cluster_queue_name',
    '-q',
    type=str,
    metavar='QUEUE_NAME',
    default='medium',
    show_default=True,
    help=(
        'Cluster partition/queue name used. Default value assumes that name '
        '"batch" is available on your machine. Please consult your cluster '
        'administrator for available and recommended queues'
    )
)
@parallel_options.option(
    '--keep_nextflow_log',
    '-knf',
    is_flag=True,
    default=False,
    show_default=True,
    help='Do not delete Nextflow logs after successful TOGA execution'
)
@parallel_options.option(
    '--ignore_crashed_parallel_batches',
    '-ignore',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, proceeds through parallel steps even if some batches failed. '
        'Failed batches are further added to the "failed_batches_<project_name>.tsv file stepwise. '
        'Note that the results of the failed batches will be missing from the final output'
    )
)
@legacy_and_experimental.option(
    '--legacy_chain_feature_extraction',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, legacy Python implementation will be used for chain feature extraction. '
        'NOTE: Legacy feature extraction is a parallel step, and a time-consuming one'
    )
)
@legacy_and_experimental.option(
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
@legacy_and_experimental.option(
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
@legacy_and_experimental.option(
    '--account_for_alternative_frame',
    '-alt_frame',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, codons in the alternative reading frame '
        '(=residing between compensated frameshifts) are considered '
        'when computing sequence intactness features'
    )
)
@out_options.option(
    '--output',
    '-o',
    type=click.Path(exists=False),
    metavar='PATH',
    default=None,
    show_default=False,
    help='A directory to store results into [default: toga2_run_<date_time>]'
)
@out_options.option(
    '--keep_temporary_files',
    '-k',
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, temporary directory (tmp) is left intact after execution is complete'
)
@verbosity_options.option(
    '--verbose',
    '-v',
    is_flag=True,
    default=False,
    help='Control logging verbosity'
)
@verbosity_options.option(
    '--email',
    type=str,
    metavar='EMAIL_ADDRESS',
    default=None,
    show_default=True,
    help=(
        'A valid e-mail address to send notifications to. If provided, TOGA2 '
        'will notify the user on pipeline crash, successful completeion, and '
        'potential problems after certain pipeline steps'
    )
)
@verbosity_options.option(
    '--mailx_binary',
    type=click.Path(exists=True),
    metavar='MAILX_PATH',
    default=None,
    show_default=True,
    help=(
        'A path to mailx executable; if not set, the executable with this name will be sought for in $PATH'
    )
)
@misc_options.option(
    '--fatotwobit_binary',
    type=click.Path(exists=True),
    metavar='FATOTWOBIT_PATH',
    default=None,
    show_default=True,
    help=(
        'A path to UCSC faToTwoBit executable; if not set, the executable with this name '
        'will be sought for in $PATH'
    )
)
@misc_options.option(
    '--twobittofa_binary',
    type=click.Path(exists=True),
    metavar='TWOBITTOFA_PATH',
    default=None,
    show_default=True,
    help=(
        'A path to UCSC twoBitToFa executable; if not set, the executable with this name '
        'will be sought for in $PATH'
    )
)

def run(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.     .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM  
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM  
         MM     MM        MM MM             ,M  `MM          MM    MM  
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM  
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM  
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    run - Run TOGA2 main pipeline

    \b
    Mandatory arguments are:
    * REF_2BIT - a path to reference genome assembly in .2bit format;
    * QUERY_2BIT - a path to query genome assembly, also in .2bit format;
    * ALIGNMENT_CHAINS - a path to genome alignment chains, with REF_2BIT as reference and QUERY_2BIT as query. Can be compressed in .gzip format;
    * REF_ANNOTATION_BED - protein-coding gene annotation for the reference genome, in Bed12 format
    """
    TogaMain(**kwargs)


@toga2.command(
    context_settings=CONTEXT_SETTINGS, 
    no_args_is_help=True, 
    short_help='Run TOGA2 pipeline with a configuration file'
)
@click.argument(
    'config_file',
    type=click.File('r', lazy=True),
    metavar='CONFIG_FILE',
)
@click.option(
    '--override',
    type=str,
    metavar='SETINGS_LIST',
    default=None,
    show_default=True,
    help=(
        'Additional settings for TOGA2 listed in double quotation marks. '
        'Settings provided this way will supersede those listed in the configuration file'
    )
)

def from_config(config_file: click.File, override: Optional[str]) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.     .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM  
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM  
         MM     MM        MM MM             ,M  `MM          MM    MM  
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM  
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM  
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    from_config - Run TOGA2 pipeline with a predefined configuration file

    \b
    NOTE: This mode is currently under development and has not been tested yet.
    """
    pass


@toga2.command(
    context_settings=CONTEXT_SETTINGS, 
    no_args_is_help=True, 
    short_help="List example commands for 'run' mode"
)

def cookbook() -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.     .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM  
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM  
         MM     MM        MM MM             ,M  `MM          MM    MM  
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM  
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM  
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    cookbook - A detailed list of TOGA2 example commands & best practices.
    NOTE: This mode is currently under development, with the list of commands being gradually expanded
    """
    pass


@toga2.command(
    context_settings=CONTEXT_SETTINGS, 
    no_args_is_help=True, 
    short_help='Test TOGA2 with companion dataset'
)
def test() -> None:
    pass


if __name__ == '__main__':
    toga2()