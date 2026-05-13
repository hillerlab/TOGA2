#!/usr/bin/env python3

"""
Master script for TOGA2
"""

import logging
import os
import sys
from typing import Any, Dict, List, Optional

import click

from __version__ import __version__
from src.python.modules.cesar_wrapper_constants import (
    DEF_BLOSUM_FILE,
    EQUIPROBABLE_ACCEPTOR,
    FIRST_ACCEPTOR,
    HG38_CANON_U2_ACCEPTOR,
    HG38_CANON_U2_DONOR,
    HG38_CANON_U12_ACCEPTOR,
    HG38_CANON_U12_DONOR,
    HG38_NON_CANON_U2_ACCEPTOR,
    HG38_NON_CANON_U2_DONOR,
    HG38_NON_CANON_U12_ACCEPTOR,
    HG38_NON_CANON_U12_DONOR,
    LAST_DONOR,
    MIN_ASMBL_GAP_SIZE,
)
from src.python.modules.codon_alignment import ALIGNERS_TO_USE, PRANK
from src.python.modules.constants import TOGA2_EPILOG, Constants
from src.python.modules.input_producer import (
    DEFAULT_MEMORY_LIMIT,
    MIN_INTRON_LENGTH_FOR_CLASSIFICATION,
    MIN_INTRON_LENGTH_FOR_PROFILES,
)
from src.python.modules.shared import (
    CONTEXT_SETTINGS,
    DependentOption,
    PrettyGroup,
)

__author__ = "Yury V. Malovichko"
__year__ = "2025"
__credits__ = ("Bogdan M. Kirilenko", "Michael Hiller")

logging.basicConfig(level=logging.INFO)

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
BIN: str = os.path.join(LOCATION, "bin")
POSTOGA_DIR: str = os.path.join(LOCATION, "postoga")
TEST_DIR: str = os.path.join(LOCATION, "test_input")
CESAR_PROFILE_PATH: str = os.path.join(LOCATION, "supply", "CESAR2.0", "profiles", "human")

HG38_CANON_U2_ACCEPTOR: str = os.path.join(LOCATION, *HG38_CANON_U2_ACCEPTOR)
HG38_CANON_U2_DONOR: str = os.path.join(LOCATION, *HG38_CANON_U2_DONOR)
HG38_NON_CANON_U2_ACCEPTOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U2_ACCEPTOR)
HG38_NON_CANON_U2_DONOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U2_DONOR)
HG38_CANON_U12_ACCEPTOR: str = os.path.join(LOCATION, *HG38_CANON_U12_ACCEPTOR)
HG38_CANON_U12_DONOR: str = os.path.join(LOCATION, *HG38_CANON_U12_DONOR)
HG38_NON_CANON_U12_ACCEPTOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U12_ACCEPTOR)
HG38_NON_CANON_U12_DONOR: str = os.path.join(LOCATION, *HG38_NON_CANON_U12_DONOR)
EQUIPROBABLE_ACCEPTOR: str = os.path.join(LOCATION, *EQUIPROBABLE_ACCEPTOR)
FIRST_ACCEPTOR: str = os.path.join(LOCATION, *FIRST_ACCEPTOR)
LAST_DONOR: str = os.path.join(LOCATION, *LAST_DONOR)
# HL_COMMON_ACCEPTOR: str = os.path.join(*HL_COMMON_ACCEPTOR)
# HL_COMMON_DONOR: str = os.path.join(*HL_COMMON_DONOR)
# HL_FIRST_ACCEPTOR: str = os.path.join(*HL_FIRST_ACCEPTOR)
# HL_LAST_DONOR: str = os.path.join(*HL_LAST_DONOR)
# HL_EQ_ACCEPTOR: str = os.path.join(LOCATION, *HL_EQ_ACCEPTOR)
# HL_EQ_DONOR: str = os.path.join(LOCATION, *HL_EQ_DONOR)
BLOSUM_FILE: str = os.path.join(LOCATION, *DEF_BLOSUM_FILE)

## Option groups
mandatory: PrettyGroup = PrettyGroup("Mandatory arguments for this mode")
input_options: PrettyGroup = PrettyGroup("General input")
bundle_options: PrettyGroup = PrettyGroup(
    "Annotation bundle",
    help="Input configuration through input directory structure template",
)
control_flow_options: PrettyGroup = PrettyGroup(
    "Pipeline", help="Control flow settings"
)
extraction_options: PrettyGroup = PrettyGroup(
    "Feature extraction",
    help="Settings for the feature extraction step",
)
class_options: PrettyGroup = PrettyGroup(
    "Classification",
    help="Projection classification settings",
)
gene_select_options: PrettyGroup = PrettyGroup(
    "Query gene selection",
    help="Controls orthology/completeness classes of query projections to annotate",
)
prepr_options: PrettyGroup = PrettyGroup(
    "Preprocessing",
    help="Data preprocessing for CESAR alignment",
)
cesar_options: PrettyGroup = PrettyGroup(
    "CESAR alignment",
    help="Exon alignment & gene annotation with CESAR",
)
parallel_options: PrettyGroup = PrettyGroup(
    "Parallel execution",
    help="Execution parameters for the pipeline's parallel steps",
)
container_options: PrettyGroup = PrettyGroup(
    "Container execution",
    help="Container execution options",
)
spliceai_options: PrettyGroup = PrettyGroup(
    "SpliceAI use",
    help="SpliceAI use for exon annotation, splice site correction, and intron gain search",
)
spliceai_run_options: PrettyGroup = PrettyGroup(
    "SpliceAI settings", help=("SpliceAI and SpliceAI wrapper options")
)
annot_options: PrettyGroup = PrettyGroup(
    "Annotation", help="Post-CESAR gene annotation & mutation check settings"
)
loss_options: PrettyGroup = PrettyGroup(
    "Gene loss", help="Gene conservation/loss classification settings"
)
orth_options: PrettyGroup = PrettyGroup(
    "Orthology resolution",
    help="Orthology resolution settings, including the gene tree-based orthology refinement",
)
browser_options: PrettyGroup = PrettyGroup(
    "UCSC browser",
    help="UCSC genome browser report parameters",
)
utr_options: PrettyGroup = PrettyGroup(
    "UTR annotation",
    help="Settings for the UTR annotation module",
)
legacy_and_experimental: PrettyGroup = PrettyGroup("Legacy & experimental features")
verbosity_options: PrettyGroup = PrettyGroup(
    "Verbosity",
    help="Verbosity & notifications controls",
)
binary_options: PrettyGroup = PrettyGroup(
    "Executables",
    help="Auxiliary executables & third party software",
)
aligner_options: PrettyGroup = PrettyGroup(
    "Aligner options",
    help="Nucleotide aligner selection & settings",
)
intronic_options: PrettyGroup = PrettyGroup(
    "IntronIC options",
    help="IntronIC execution settings & control",
)
integration_options: PrettyGroup = PrettyGroup(
    "Integration options",
    help="Parameters for meta-gene inference and isoform selection",
)
out_options: PrettyGroup = PrettyGroup("Output")
misc_options: PrettyGroup = PrettyGroup("Miscellaneous")
qc_options: PrettyGroup = PrettyGroup("Quality Control")


@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.version_option(__version__, "--version", "-V", prog_name="TOGA2")
def toga2() -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
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
    short_help="Run TOGA2 pipeline with command line arguments",
)
@input_options.option(
    "--ref_2bit",
    type=click.Path(exists=True),
    metavar="REF_2BIT",
    cls=DependentOption,
    required=True,
    not_required_if=["input_directory"],
    help="Reference genome assembly file, in .2bit format",
)
@input_options.option(
    "--query_2bit",
    type=click.Path(exists=True),
    metavar="QUERY_2BIT",
    cls=DependentOption,
    required=True,
    not_required_if=["input_directory"],
    help="Query genome assembly file, in .2bit format",
)
@input_options.option(
    "--chain_file",
    type=click.Path(exists=True),
    metavar="ALIGNMENT_CHAINS",
    cls=DependentOption,
    required=True,
    not_required_if=["input_directory"],
    help="""Genome alignment chains, with REF_2BIT as reference and QUERY_2BIT as query. 
TOGA2 annotates query genome by projecting reference transcripts through the chains 
contained in this file. Can be compressed in .gzip format""",
)
@input_options.option(
    "--ref_annotation",
    type=click.Path(exists=True),
    metavar="REF_ANNOTATION_BED",
    cls=DependentOption,
    required=True,
    not_required_if=["input_directory"],
    help="""Reference annotation file, in BED12 format. TOGA2 annotates transcripts in the 
query genome by projecting reference transcripts contained in this file""",
)
@input_options.option(
    "--isoform_file",
    "-i",
    type=click.Path(exists=True),
    metavar="ISOFORMS_FILE",
    cls=DependentOption,
    competes_with=["no_isoform_file"],
    required_mutex=True,
    not_required_if=["input_directory"],
    help="A path to a two-column tab-separated file containing gene-to-isoform mapping",
)
@input_options.option(
    "--no_isoform_file",
    is_flag=True,
    default=False,
    cls=DependentOption,
    competes_with=["isoform_file"],
    required_mutex=True,
    not_required_if=["input_directory"],
    help="""A flag indicating that TOGA2 will be used without an isoform file. 
Highly discouraged unless assigning transcript to genes is complicated for your reference 
or you want each transcript in REF_ANNOTATION_BED to be treated as a separate gene""",
)
@input_options.option(
    "--u12_file",
    "-u12",
    type=click.Path(exists=True),
    metavar="U12_FILE",
    cls=DependentOption,
    competes_with=["no_u12_file"],
    required_mutex=True,
    not_required_if=["input_directory"],
    help="""A three-column tab-separated file containing information on the 
non-canonical splice sites""",
)
@input_options.option(
    "--no_u12_file",
    is_flag=True,
    default=False,
    cls=DependentOption,
    competes_with=["u12_file"],
    required_mutex=True,
    not_required_if=["input_directory"],
    help="""A flag indicating that TOGA2 will be used without a U12/non-canonical U2 classification file. 
Highly discouraged unless you really do not want to discriminate between GT/GC-AG U2 and other 
intron classes in in your annotation""",
)
@input_options.option(
    "--spliceai_dir",
    "-sai",
    type=click.Path(exists=True),
    metavar="SPLICEAI_OUT_DIR",
    cls=DependentOption,
    competes_with=["no_spliceai"],
    required_mutex=True,
    not_required_if=["input_directory"],
    help="""A path to the SpliceAI predictions directory produced by TOGA2 `run_spliceai` mode. 
These data are used for improved exon annotation in the query""",
)
@input_options.option(
    "--no_spliceai",
    is_flag=True,
    default=False,
    cls=DependentOption,
    competes_with=["spliceai_dir"],
    required_mutex=True,
    not_required_if=["input_directory"],
    help="""A flag indicating that TOGA2 will be used without SpliceAI predictions for exon annotation. 
Highly discouraged unless you cannot obtain SpliceAI annotation for your query genome, 
suspect that SpliceAI performs suboptimally for your query species, or have concerns about 
TOGA2 performance speed and memory consumption""",
)
@bundle_options.option(
    "--input_directory",
    type=click.Path(exists=True),
    metavar="INPUT_DIR",
    cls=DependentOption,
    requires=["ref_name", "query_name"],
    default=None,
    show_default=True,
    help="""A single input directory containing all the necessary for reference (genome assembly 
in 2bit format, reference annotation, reference isoforms, and U12 files), query (genome assembly 
in 2bit format, SpliceAI annotation), and genome alignment (a single chain file). Note that:\n
    \ta) you can override existing and complement missing files in the directory with respective flags 
from the "General input" section\n
    \tb) you can provide flag placeholders for isoform file, U12 intron file, and SpliceAI directories, 
if you don't have the respective data\n
See Manual for more details on input directory structure and formatting""",
)
@bundle_options.option(
    "--ref_name",
    type=str,
    metavar="REF_NAME",
    cls=DependentOption,
    # requires=["input_directory", "query_name"],
    default=None,
    show_default=True,
    help="""Reference organism/assembly name used in input files' names. 
See Manual for more details on input directory formatting""",
)
@bundle_options.option(
    "--query_name",
    type=str,
    metavar="QUERY_NAME",
    cls=DependentOption,
    # requires=["input_directory", "ref_name"],
    default=None,
    show_default=True,
    help="""Query organism/assembly name used in input files' names. 
See Manual for more details on input directory formatting""",
)
@control_flow_options.option(
    "--resume_from",
    "-res",
    type=click.Choice(Constants.RESUME_OPTIONS, case_sensitive=False),
    metavar="STEP",
    default="all",
    show_default=True,
    help="""If you have an unfinished run and want to resume with the same results,
        select the step from the following list:\b
        all: a placeholder for full starts the pipeline from the very beginning;\b\n
        setup: input data filtering, indexing, and format conversion;\b\n
        feature_extraction: projection and chain feature extraction for projection classification;\b\n
        classification: projection classification in terms of orthology;\b\n
        preprocessing: projection data preprocessing for further alignment; parallel step;\b\n
        aggregate_preprocessing_res: data aggregation and summary across independent preprocessing batches;\b\n
        alignment: CESAR alignment, mutation check, and loss inference for each projection; parallel step;\b\n
        aggregate_cesar_res: data aggregation and summary across independent alignment batches;\b\n
        gene_inference: annotate query genes based on the transcript-level annotation results;\b\n
        loss_summary: gene loss data summary;\b\n
        orthology: orthology relationship resolution; if "-st" flag is set, gene tree orthology batches are run at this step;\b\n
        summarize_trees: if "-st" flag was set, individual gene tree batch results are summarized and added to the original orthology data at this step;\b\n
        finalize: gene renaming and output Bed file filtering;\b\n
        ucsc_report: BigBed file preparation\n""",
)
@control_flow_options.option(
    "--halt_at",
    "-halt",
    type=click.Choice(Constants.RESUME_OPTIONS, case_sensitive=False),
    metavar="STEP",
    default="all",
    show_default=True,
    help="""Halts the pipeline at the selected step (see above). Option "all" 
implies running the pipeline to the last step""",
)
@control_flow_options.option(
    "--selected_feature_batches",
    "-feat_b",
    type=str,
    metavar="COMMA_SEPARATED_LIST",
    default=None,
    show_default=True,
    help="""A comma-separated list of batch numbers for the feature extraction step to be run. 
Valid only if --resume_from is set to "feature_extraction" or lower, and 
"--legacy_chain_feature_extraction option was set in both interrupted and scheduled runs""",
)
@control_flow_options.option(
    "--selected_preprocessing_batches",
    "-prep_b",
    type=str,
    metavar="COMMA_SEPARATED_LIST",
    default=None,
    show_default=True,
    help="""A comma-separated list of batch numbers for the CESAR preprocessing step to be run. 
Valid only if --resume_from is set to "preprocessing" or lower""",
)
@control_flow_options.option(
    "--selected_alignment_batches",
    "-aln_b",
    type=str,
    metavar="COMMA_SEPARATED_LIST",
    default=None,
    show_default=True,
    help="""A comma-separated list of batch numbers for the CESAR alignment step to be run. 
Valid only if --resume_from is set to "alignment" or lower""",
)
@control_flow_options.option(
    "--no_utr_annotation",
    "-no_utr",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, UTR sequences are not added to the final annotation file",
)
@extraction_options.option(
    "--min_chain_score",
    "-mcs",
    type=click.IntRange(min=0, max=None),
    metavar="INT",
    default=15000,
    show_default=True,
    help="""Minimal score for chains to be considered for classification. Setting 
this value to zero disables chain filtering""",
)
@extraction_options.option(
    "--min_orthologous_chain_score",
    "-minscore",
    type=click.IntRange(min=0, max=None),
    metavar="INT",
    default=15000,
    show_default=True,
    help="""Minimal score for chains to be potentially classified as orthologous. Chains with 
min_chain_score <= X < min_orthologous_chain_score are discarded unless they are 
classified as retrogenes/processed pseudogenes""",
)
@extraction_options.option(
    "--feature_jobs",
    "-fj",
    type=int,
    metavar="INT",
    default=100,
    show_default=True,
    help="A number of jobs to projection feature extraction commands into",
)
@class_options.option(
    "--orthology_threshold",
    "-ot",
    type=float,
    metavar="FLOAT",
    default=Constants.DEFAULT_ORTH_THRESHOLD,
    show_default=True,
    help="Probability threshold for considering projections as orthologous",
)
@class_options.option(
    "--single_exon_model",
    "-se_model",
    type=click.Path(exists=True),
    metavar="PATH",
    default=os.path.join(LOCATION, "models", "se_model.dat"),
    show_default=True,
    help="A path to a orthology classification model for single exon reference transcripts",
)
@class_options.option(
    "--multi_exon_model",
    "-me_model",
    type=click.Path(exists=True),
    metavar="PATH",
    default=os.path.join(LOCATION, "models", "me_model.dat"),
    show_default=True,
    help="A path to a orthology classification model for multi-exon reference transcripts",
)
@class_options.option(
    "--use_long_distance_model",
    "-use_ld",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, applies extra classifier for distantly related species; 
relevant at molecular distances >=1sps""",
)
@class_options.option(
    "--long_distance_model",
    "-ld_model",
    type=click.Path(exists=True),
    metavar="PATH",
    default=os.path.join(LOCATION, "models", "ld_model.dat"),
    show_default=True,
    help="""A path to a refinement classification model for distantly related 
reference-query species pairs""",
)
@gene_select_options.option(
    "--disable_fragment_assembly",
    "-no_f",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, does not attempt to recover fragmented projections 
from individual chains""",
)
@gene_select_options.option(
    "--orthologs_only",
    "-o_only",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, only orthologous projections are considered",
)
@gene_select_options.option(
    "--one2ones_only",
    "-o2o",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, only transcript with a single orthologous projection are considered",
)
@gene_select_options.option(
    "--paralogs_over_spanning",
    "-p_first",
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        "If set, paralogous projections take priority over spanning-chain projections "
        "when determining the chains to project the transcript through"
    ),
)
@gene_select_options.option(
    "--enable_spanning_chains",
    "-nospan",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, spanning chains (i.e., chains with alignment gap corresponding 
to the projected transcript) are considered for CESAR alignment; otherwise 
spanning chains are used only to discriminate between Lost and Missing projections""",
)
@gene_select_options.option(
    "--annotate_processed_pseudogenes",
    "-pp",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, processed pseudogene projections are considered for CESAR alignment; 
otherwise, a separate BED9 track "processed_pseudogenes.bed" is added to output directory""",
)
@prepr_options.option(
    "--preprocessing_jobs",
    "-pj",
    type=int,
    metavar="INT",
    default=300,
    show_default=True,
    help="A number of jobs to partition CESAR preprocessing commands into",
)
@prepr_options.option(
    "--max_chains_per_transcript",
    "-mc",
    type=int,
    metavar="INT",
    default=100,
    show_default=True,
    help="""A maximum number of chains to project each transcript through. If the 
number of projections exceeds the given value, only the first N chains will be considered""",
)
@prepr_options.option(
    "--max_search_space_size",
    "-mss",
    type=int,
    metavar="INT (BP)",
    default=1_000_000,
    show_default=True,
    help="""Query sequence length limit for CESAR jobs. Projections in which any 
exon group aligns to sequence beyond this limit will be discarded "
from alignment step""",
)
@prepr_options.option(
    "--extrapolation_modifier",
    "-em",
    type=float,
    metavar="FLOAT",
    default=1.2,
    show_default=True,
    help="""Modifier by which counterparts of missing reference sequence are modified 
when extrapolating ambiguous projection termini""",
)
@prepr_options.option(
    "--minimal_covered_fraction",
    "-mincov",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.0,
    show_default=True,
    help="""Minimal fraction of reference CDS sequence covered by aligned chain 
blocks to consider the projection for CESAR alignment. If the 
--enable_spanning_chains flag is set, this parameter is automatically 
set to 0.0""",
)
@prepr_options.option(
    "--exon_locus_flank",
    "-ef",
    type=int,
    metavar="INT (BP)",
    default=100,
    show_default=True,
    help="Size of flank to be added to each defined exon locus from both sides. ",
)
@prepr_options.option(
    "--assembly_gap_size",
    "-gs",
    type=int,
    metavar="INT",
    default=MIN_ASMBL_GAP_SIZE,
    show_default=True,
    help="Minimum number of consecutive N symbols to be considered an assembly gap",
)
@cesar_options.option(
    "--cesar_binary",
    "-cs",
    type=click.Path(exists=True),
    metavar="CESAR_BINARY",
    default=None,
    show_default=False,
    help="""A path to the actual CESAR2.0 binary; if not provided, will look for one 
in the PATH, otherwise defaulting to the CESAR2.0 instance 
in the current directory""",
)
@cesar_options.option(
    "--memory_bins",
    "-b",
    type=str,
    metavar="BIN_LIST",
    default="3,5,7,10,15",
    show_default=True,
    help="""A comma-separated list of memory bin caps, in GB. For each memory bin, 
job scheduling will be performed independently. If you want to process 
memory-intensive projections as a single cluster call, set the last value 
to "big" or enable the --allow_heavy_jobs flag""",
)
@cesar_options.option(
    "--job_nums_per_bin",
    "-jb",
    type=str,
    metavar="BIN_JOB_NUM_LIST",
    default="500,20,20,10,5",
    show_default=True,
    help="""A comma-separated list of job numbers per memory bin. 
Job jumbers must follow in the same order as memory caps passed 
to --memory_bins option""",
)
@cesar_options.option(
    "--allow_heavy_jobs",
    "-ahj",
    type=bool,
    metavar="FLAG",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Aggregate all jobs exceeding the highest memory cap as a single joblist; 
if memory bins are provided, duplicates the "big" memory bin behavior""",
)
@cesar_options.option(
    "--cesar_memory_limit",
    "-ml",
    type=float,
    metavar="FLOAT (GB)",
    default=15.0,
    show_default=True,
    help="""Upper memory limit for CESAR jobs. Projections in which any exon group 
requires memory beyond this limit will be discarded before the alignment step. 
Value of 0 denotes unlimited memory""",
)
@cesar_options.option(
    "--cesar_profile_dir",
    "-prof",
    type=click.Path(exists=True),
    metavar="PATH",
    default=CESAR_PROFILE_PATH,
    show_default=True,
    help="""A path to a single directory containing the full set of CESAR2 profiles.
The following files are expected:\n
\tcanon_U2_donor.tsv - for canonical (GT/GC-AG) U2 donors;\n
\tcanon_U2_acceptor.tsv - for canonical (GT/GC-AG) U2 acceptors;\n
\tnonCanon_U2_donor.tsv - for non-canonical U2 donors;\n
\tnonCanon_U2_acceptor.tsv - for non-canonical U2 acceptors;\n
\tcanon_U12_donor.tsv - for canonical (GT-AG) U12 donors;\n
\tcanon_U12_acceptor.tsv - for canonical (GT-AG) U12 acceptors'\n
\tnonCanon_U12_donor.tsv - for non-canonical U12 donors;\n
\tnonCanon_U12_acceptor.tsv - for non-canonical U12 acceptors\n
Uniform profiles, equiprobable_donor.tsv and equiprobable_acceptor.tsv, generated 
in the `prepare-input` mode, are not expected by default but are highly recommended 
to keep.\n
NOTE: the individual CESAR2 profile arguments above override the expected profiles in the CESAR2 
profile directory. You can also provide an incomplete profile directory and specify the missing 
profiles with the respective arguments, but TOGA2 missing any of the eight profile files will 
result in an error.""",
)
@cesar_options.option(
    "--cesar_canon_u2_acceptor",
    "-cca",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # HG38_CANON_U2_ACCEPTOR,
    show_default=True,
    help="A path to canonical (GT/GC-AG) U2 acceptor profile",
)
@cesar_options.option(
    "--cesar_canon_u2_donor",
    "-ccd",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # HG38_CANON_U2_DONOR,
    show_default=True,
    help="A path to canonical (GT/GC-AG) U2 donor profile",
)
@cesar_options.option(
    "--cesar_non_canon_u2_acceptor",
    "-cnca",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # HG38_NON_CANON_U2_ACCEPTOR,
    show_default=True,
    help="A path to non-canonical (non GT/GC-AG) U2 acceptor profile",
)
@cesar_options.option(
    "--cesar_non_canon_u2_donor",
    "-cncd",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # HG38_NON_CANON_U2_DONOR,
    show_default=True,
    help="A path to non-canonical (non GT/GC-AG) U2 donor profile",
)
@cesar_options.option(
    "--cesar_canon_u12_acceptor",
    "-cua",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # HG38_CANON_U12_ACCEPTOR,
    show_default=True,
    help="A path to canonical (GT-AG) U12 exon acceptor profile",
)
@cesar_options.option(
    "--cesar_canon_u12_donor",
    "-cud",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # HG38_CANON_U12_DONOR,
    show_default=True,
    help="A path to canonical (GT-AG) U12  donor profile",
)
@cesar_options.option(
    "--cesar_non_canon_u12_acceptor",
    "-cnua",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # EQUIPROBABLE_ACCEPTOR,
    show_default=True,
    help="A path to non-canonical (non-GT-AG) U12 exon acceptor profile",
)
@cesar_options.option(
    "--cesar_non_canon_u12_donor",
    "-cnud",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,  # HG38_NON_CANON_U12_DONOR,
    show_default=True,
    help="A path to non-canonical (non-GT-AG) U12 exon donor profile",
)
@cesar_options.option(
    "--cesar_first_acceptor",
    "-cfa",
    type=click.Path(exists=True),
    metavar="PATH",
    default=FIRST_ACCEPTOR,
    show_default=True,
    help="A (relative to CESAR2 location) path to first exon acceptor profile",
)
@cesar_options.option(
    "--cesar_last_donor",
    "-cld",
    type=click.Path(exists=True),
    metavar="PATH",
    default=LAST_DONOR,
    show_default=True,
    help="A (relative to CESAR2 location) path to last exon donor profile",
)
@cesar_options.option(
    "--joint_splice_site_treatment",
    "-j_sst",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, donor and acceptor intron splice sites are treated 
in conjunction in terms of their `canonicity`""",
)
@spliceai_options.option(
    "--spliceai_correction_mode",
    "-scm",
    type=click.IntRange(min=0, max=7),
    metavar="MODE_NUM",
    default=3,
    show_default=False,
    help="""Set the mode of SpliceAI-mediated exon boundary correction:\b\n
0 - no correction [default; equivalent to not providing SpliceAI data directory];\b\n
1 - use SpliceAI predictions to correct boundaries of missing and deleted exons;\b\n
2 - correct mutated canonical U2 splice sites;\b\n
3 - correct all canonical U2 splice sites in the presence of alternatives with higher SpliceAI support;\b\n
4 - correct all canonical U2 as well as mutated GT-AG U12 splice sites;\b\n
5 - correct all canonical U2 and mutated and/or unsupported  GT-AG U12 splice sites;\b\n
6 - correct all canonical U2 and all U12 splice sites;\b\n
7 - correct all U2 (including known non-canonical sites) and U12 splice sites.\b\n
Current recommendation (at least, for vertebrate queries) is 3""",
)
@spliceai_options.option(
    "--min_splice_prob",
    "-msp",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.02,
    show_default=True,
    help="Minimum SpliceAI prediction probability to consider the splice site",
)
@spliceai_options.option(
    "--splice_prob_margin",
    "-spm",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.02,
    show_default=True,
    help="""For splice sites with SpliceAI support 0<x<min_splice_prob, ignore 
alternative sites with support < x + splice_prob_margin""",
)
@spliceai_options.option(
    "--intron_gain_check",
    "-ig",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, performs SpliceAI-guided check for query-specific introns",
)
@spliceai_options.option(
    "--intron_gain_threshold",
    "-igt",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.5,
    show_default=True,
    help="Minimal intron gain threshold to consider",
)
@spliceai_options.option(
    "--min_intron_prob_trusted",
    "-mipt",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.1,
    show_default=True,
    help="""Minimal SpliceAI support for query-specific introns supported by the presence of 
both extensive alignment gaps and frame-shifting/nonsense mutations""",
)
@spliceai_options.option(
    "--min_intron_prob_supported",
    "-mips",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.2,
    show_default=True,
    help="""Minimal SpliceAI support for query-specific introns supported by the presence of 
either extensive alignment gaps and frame-shifting/nonsense mutations""",
)
@spliceai_options.option(
    "--min_intron_prob_unsupported",
    "-mipu",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.8,
    show_default=True,
    help="""Minimal SpliceAI support for query-specific introns not supported by 
either extensive alignment gaps and frame-shifting/nonsense mutations""",
)
@spliceai_options.option(
    "-max_intron_number",
    "-mnn",
    type=click.IntRange(min=1, max=None),
    metavar="INT",
    default=4,
    show_default=True,
    help="""Maximum gained intron number per exon. 
Highly recommended not to increase this beyond 5-6""",
)
@annot_options.option(
    "--matrix",
    "-m",
    type=click.Path(exists=True),
    metavar="BLOSUM_MATRIX_FILE",
    default=BLOSUM_FILE,
    show_default=True,
    help="A file containing the protein alignment matrix",
)
@annot_options.option(
    "--mask_n_terminal_mutations",
    "-m10m",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, masks mutations occurring in the first 10 percents 
of query projection length regardless of alternative start codon presence""",
)
@annot_options.option(
    "--disable_missing_stop_search",
    "-rmo",
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, TOGA2 does not search for alternative stop codons 
downstream to CESAR alignment boundary if the original alignment 
does not end with one""",
)
@loss_options.option(
    "--accepted_loss_symbols",
    "-l",
    type=str,
    metavar="LOSS_SYMBOLS",
    default=",".join(Constants.DEFAULT_LOSS_SYMBOLS),
    show_default=True,
    help="""A comma-separated list of loss status symbols; only projections of 
respective statuses will be considered for orthology resolution. 
Supported symbols are: %s. Keyword ALL lets all possible statuses in."""
    % ",".join(Constants.ALL_LOSS_SYMBOLS),
)
@orth_options.option(
    "--skip_gene_trees",
    "-st",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, skips gene-tree based resolution step for convoluted many2many clades",
)
@orth_options.option(
    "--use_raxml",
    "-raxml",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use raxmlHPC-PTHREADS-AVX instead of IqTree2 for tree inference",
)
@orth_options.option(
    "--max_clique_size",
    "-mqs",
    type=int,
    metavar="INT",
    default=50,
    show_default=True,
    help="A maximum number of sequences in many:many cliques to be resolved with gene trees",
)
@orth_options.option(
    "--orthology_jobs",
    "-oj",
    type=int,
    metavar="INT",
    default=100,
    show_default=True,
    help="A number of jobs to split orthology fine resolution commands into",
)
@orth_options.option(
    "--prank_binary",
    "-pb",
    type=click.Path(exists=True),
    metavar="PRANK_BINARY",
    default=None,
    show_default=True,
    help="""A path to the PRANK executable to be used at fine resolution step. 
If not provided, the program will try to infer its location from the PATH""",
)
@orth_options.option(
    "--tree_binary",
    "-rb",
    type=click.Path(exists=True),
    metavar="TREE_BINARY",
    default=None,
    show_default=True,
    help="""A path to the IqTree2/raxmlHPC-PTHREADS-AVX executable to be used at fine resolution step. 
If not provided, the program will try to infer its location from the PATH.""",
)
@orth_options.option(
    "--tree_cpus",
    "-rc",
    type=int,
    metavar="INT",
    default=1,
    show_default=True,
    help="A maximum number of CPUs to run IqTree2/RAxML with",
)
@utr_options.option(
    "--utr_abs_threshold",
    "-utr_abs",
    type=click.IntRange(min=0),
    default=3000,
    show_default=True,
    help="Absolute threshold by which the projected UTR block/exon can exceed the reference counterpart",
)
@utr_options.option(
    "--utr_rel_threshold",
    "-utr_rel",
    type=click.FloatRange(min=0.0),
    default=2.5,
    show_default=True,
    help="""Relative (to reference exon length) threshold by which the projected UTR block/exon 
can exceed the reference counterpart""",
)
@utr_options.option(
    "--no_utr_boundary_extrapolation",
    "-no_utr_extra",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, UTR block/exon boundaries will be projected by the last aligned UTR position, 
without extrapolating the terminal unaligned regions""",
)
@utr_options.option(
    "--no_adjacent_utr_extra",
    "-no_adj_utr",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, unaligned CDS-adjacent UTR sequences will not be extrapolated",
)
@utr_options.option(
    "--fixed_adjacent_utr_extra",
    "-fixed_adj_utr",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, unaligned CDS-adjacent UTR sequences will be annotated as block of fixed length 
equal to --utr_abs_threshold value""",
)
@browser_options.option(
    "--link_file",
    "-lf",
    type=click.Path(exists=True),
    metavar="FILE",
    default=None,
    show_default=False,
    help="""A path to the two-column file containing HTML-formatted hyperlinks 
to external sources on reference transcripts""",
)
@browser_options.option(
    "--ucsc_prefix",
    "-up",
    type=str,
    metavar="PREFIX",
    default=Constants.DEFAULT_UCSC_PREFIX,
    show_default=True,
    help="A prefix to use in the output file names",
)
@parallel_options.option(
    "--parallel_strategy",
    "-c",
    type=click.Choice(Constants.ALL_PARALLEL_EXECS),
    metavar="PARALLEL_EXECUTOR",
    default="local",
    show_default=True,
    help="""Specify the HPC strategy. By default, TOGA2 uses Nextflow to handle
parallel processes and supports, at least in theory, all Nextflow executors.
Please consult the full list of options at Nextflow help page. Note that setting
executor to "local" will parallel the processes over the local machine CPUs.\b\n
If you want to use Parasol as parallel process manager, set this option to "para"\b
If you want to implement a fully custom parallel manager strategy, modify the CustomStrategy class
and set this option to "custom" or contact the TOGA2 team""",
)
@parallel_options.option(
    "--nextflow_exec_script",
    type=click.Path(exists=True),
    metavar="NEXTFLOW_SCRIPT_PATH",
    default=None,
    show_default=True,
    help="""A path to a user-defined Nextflow script used for parallel process execution. 
If not specified, TOGA2 will use a minimal boilerplate script instead. 
Ignored if Parasol or custom HPC strategy were specified for parallel steps""",
)
@parallel_options.option(
    "--max_number_of_retries",
    type=click.IntRange(min=1),
    metavar="INT",
    default=3,
    show_default=True,
    help="""Maximum number of retries per parallel job before reporting job failure. 
Ignored if Parasol or custom HPC strategy were specified for parallel steps""",
)
@parallel_options.option(
    "--nextflow_config_dir",
    "-nc",
    type=click.Path(exists=True),
    metavar="NEXTFLOW_CONFIG_DIR",
    default=None,
    show_default=True,
    help="""A path to a directory containing user-defined Nextflow configuration files; 
for detauls, please see nextflow_config_files/readme.txt""",
)
@parallel_options.option(
    "--max_parallel_time",
    "-max_t",
    type=click.IntRange(min=1),
    metavar="HOURS",
    default=24,
    show_default=True,
    help="Maximum time duration (in hours) for Nextflow parallel processes",
)
@parallel_options.option(
    "--cluster_queue_name",
    "-q",
    type=str,
    metavar="QUEUE_NAME",
    default="batch",
    show_default=True,
    help="""Cluster partition/queue name used. Default value assumes that name 
"batch" is available on your machine. Please consult your cluster 
administrator for available and recommended queues""",
)
@parallel_options.option(
    "--keep_nextflow_log",
    "-knf",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not delete Nextflow logs after successful TOGA execution",
)
@parallel_options.option(
    "--ignore_crashed_parallel_batches",
    "-ignore",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, proceeds through parallel steps even if some batches failed. 
Failed batches are further added to the "failed_batches_<project_name>.tsv" file stepwise. 
Note that the results of the failed batches will be missing from the final output""",
)
@container_options.option(
    "--container_image",
    type=click.Path(exists=True),
    default=None,
    show_default=True,
    help="""A path to the executable TOGA2 container image. 
All the parallel step scripts will be executed by invoking this container. """,
)
@container_options.option(
    "--container_executor",
    type=str,
    default="apptainer",
    show_default=True,
    help='A name for container executor engine. WARNING: Currently only "apptainer" is supported',
)
@container_options.option(
    "--bindings",
    type=str,
    metavar="STRING",
    default=None,
    show_default=True,
    help="""A list of directory mounts to provide to the container instances at parallel steps. 
Binginds should be provided as expected by the container executor engine and wrapped in double 
quotes, e.g. "/tmp,/src/,~/:/home\"""",
)
@legacy_and_experimental.option(
    "--legacy_chain_feature_extraction",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, legacy Python implementation will be used for chain feature extraction. 
NOTE: Legacy feature extraction is a parallel step, and a time-consuming one""",
)
@legacy_and_experimental.option(
    "--toga1_compatible",
    "-t1",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Alignment procedure is fully TOGA1.0-compliant except for exonwise 
CESAR alignment; benchmarking feature, do not use in real runs""",
)
@legacy_and_experimental.option(
    "--toga1_plus_corrected_cesar",
    "-t1c",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Alignment procedure is fully TOGA1.0-compliant except for exonwise 
CESAR alignment and corrected CESAR-related bugs; 
benchmarking feature, do not use in real runs""",
)
@legacy_and_experimental.option(
    "--account_for_alternative_frame",
    "-alt_frame",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, codons in the alternative reading frame 
(=residing between compensated frameshifts) are considered 
when computing sequence intactness features""",
)
@out_options.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    metavar="PATH",
    default=None,
    show_default=False,
    help="A directory to store results into [default: toga2_run_<date_time>]",
)
@out_options.option(
    "--project_name",
    "-name",
    type=str,
    metavar="PROJECT_NAME",
    default="TOGA2",
    show_default=True,
    help="""A name for the current TOGA2 project. This name will be used 
as a prefix followed by run start date and time to name current runs' 
log and metadata files""",
)
@out_options.option(
    "--project_arg_format",
    "-argf",
    type=click.Choice(Constants.CONFIG_FORMATS),
    default=Constants.YAML,
    show_default=True,
    help="Project argument file format",
)
@out_options.option(
    "--keep_temporary_files",
    "-k",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, temporary directory (tmp) is left intact after execution is complete",
)
@verbosity_options.option(
    "--verbose", "-v", is_flag=True, default=False, help="Control logging verbosity"
)
@verbosity_options.option(
    "--debug",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Increase logging verbosity by logging debugging messages. 
Automatically enables \"--verbose\" flag""",
)
@verbosity_options.option(
    "--email",
    type=str,
    metavar="EMAIL_ADDRESS",
    default=None,
    show_default=True,
    help="""A valid e-mail address to send notifications to. If provided, TOGA2 "
will notify the user on pipeline crash, successful completeion, and "
potential problems after certain pipeline steps""",
)
@verbosity_options.option(
    "--mailx_binary",
    type=click.Path(exists=True),
    metavar="MAILX_PATH",
    default=None,
    show_default=True,
    help="""A path to mailx executable; if not set, the executable with 
this name will be sought for in $PATH""",
)
@binary_options.option(
    "--fatotwobit_binary",
    type=click.Path(exists=True),
    metavar="FATOTWOBIT_PATH",
    default=os.path.join(BIN, "faToTwoBit"),
    show_default=True,
    help="A path to UCSC faToTwoBit executable",
)
@binary_options.option(
    "--twobittofa_binary",
    type=click.Path(exists=True),
    metavar="TWOBITTOFA_PATH",
    default=os.path.join(BIN, "twoBitToFa"),
    show_default=True,
    help="A path to UCSC twoBitToFa executable",
)
@binary_options.option(
    "--bigwig2wig_binary",
    "-bw2w",
    type=click.Path(exists=True),
    metavar="BIGWIG2WIG_BINARY",
    default=os.path.join(BIN, "bigWigToWig"),
    help="A path to the UCSC bigWigToWig binary",
)
@binary_options.option(
    "--bedtobigbed_binary",
    type=click.Path(exists=True),
    metavar="BEDTOBIGBED_PATH",
    default=None,
    show_default=True,
    help="A path to UCSC bedToBigBed executable",
)
@binary_options.option(
    "--ixixx_binary",
    type=click.Path(exists=True),
    metavar="IXIXX_PATH",
    default=None,
    show_default=True,
    help="A path to UCSC ixIxx executable",
)
def run(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    run - Run TOGA2 main pipeline
    """
    cmd: str = " ".join(sys.argv)
    from src.python.modules.toga_main import TogaMain

    TogaMain(**kwargs, cmd=cmd)


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Run TOGA2 pipeline with a configuration file",
)
@click.option(
    "--config_file",
    "-c",
    type=click.File("r", lazy=True),
    metavar="CONFIG_FILE",
    required=True,
    default=None,
    show_default=True,
    help="A path to a config file",
)
@click.option(
    "--config_format",
    "-f",
    type=click.Choice(Constants.CONFIG_FORMATS),
    default="tsv",
    show_default=True,
    help=("Input config file format"),
)
@click.option(
    "--override",
    type=str,
    metavar="SETINGS_LIST",
    default=None,
    show_default=True,
    help="""Additional settings for TOGA2 listed in double quotation marks. 
Settings provided this way will supersede those listed in the configuration file. 
WARNING: Currenly does not accept argument short (single-dash) names""",
)
def from_config(
    config_file: click.File, config_format: Optional[str], override: Optional[str]
) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    from_config - Run TOGA2 pipeline with a predefined configuration file

    \b
    The only positional argument is a path to a two-column configuration file where first column corresponds to command line argument full names and the second column stores respective values.
    Configuration file boilerplate is provided at supply/project_args.tsv . Alternatively, you can use and/or modify logs/project_args_${proj_name}.tsv for any successful run.

    \b
    NOTE: This mode is currently under development.
    """
    from src.python.modules.toga_configured import Toga2ConfiguredLauncher
    from src.python.modules.toga_main import TogaMain

    args: List[str] = Toga2ConfiguredLauncher(
        config_file, config_format=config_format, override=override
    ).run()
    TogaMain(**args)


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Prepare reference annotation files for TOGA2 input",
)
@mandatory.option(
    "--ref_2bit",
    type=click.Path(exists=True),
    metavar="REF_2BIT",
    cls=DependentOption,
    required=True,
)
@mandatory.option(
    "--ref_annot",
    type=click.Path(exists=True),
    metavar="REF_ANNOTATION_BED",
    cls=DependentOption,
    required=True,
)
@input_options.option(
    "--ref_isoforms",
    "-i",
    type=click.Path(exists=True),
    metavar="ISOFORMS_FILE",
    default=None,
    show_default=True,
    help="""A path to a two-column tab-separated file containing gene-to-isoform mapping. 
The contents will be also checked for consistency, with transcripts missing 
gene mapping further removed from the annotation file""",
)
@input_options.option(
    "--disable_transcript_filtering",
    "-no_filter",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, copies the input annotation file to the output directory as-is, 
with no filtering. Highly discouraged unless you want to use compromised transcripts 
for intron classificaiton/CESAR2 profile generation""",
)
@input_options.option(
    "--contigs",
    type=str,
    metavar="CONTIG_LIST",
    default=None,
    show_default=True,
    help="""A comma-separated list of contig (scaffold, chromosome, etc.) names 
to restrict the fitlered annotation to. Transcripts located in other contigs 
will be excluded from the annotation.""",
)
@input_options.option(
    "--excluded_contigs",
    type=str,
    metavar="CONTIG_LIST",
    default=None,
    show_default=True,
    help="""A comma-separated list of deprecated contig (scaffold, chromosome, etc.) names 
to exclude from the filtered annotation. Transcripts located in this contigs 
will be excluded from the annotation. Contigs appearing in both --contigs 
and --excluded_contigs are treated as excluded.""",
)
@control_flow_options.option(
    "--disable_intron_classification",
    "-no_intronic",
    is_flag=True,
    default=False,
    help="""If set, stops the procedure before intron classification step. Recommended if 
you already have U12 intron file list or do not have access to intronIC. 
NOTE: Setting this flag will also disable CESAR2 profile generation.""",
)
@control_flow_options.option(
    "--disable_cesar_profiles",
    "-no_cesar",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, stops the procedure before reference-specific CESAR2 profile generation. 
TOGA2 comes with a set of trusted CESAR2 profiles for selected mammalian and avian 
references; generating custom profiles is recommended if your reference belongs to a 
distant clade with highly divergent intron structure. 
NOTE: The step is skipped automatically if --disable_intron_classification flag is set.""",
)
@intronic_options.option(
    "--intronic_binary",
    type=click.Path(exists=True),
    metavar="INTRONIC_PATH",
    default=None,
    show_default=True,
    help="""A path to intronIC binary. If not set, will check 
for executable intronIC instance in $PATH""",
)
@intronic_options.option(
    "--intronic_cores",
    type=click.IntRange(min=1),
    metavar="INT",
    default=1,
    show_default=True,
    help="Number of CPUs to run intronIC with",
)
@intronic_options.option(
    "--min_intron_length_intronic",
    type=click.IntRange(min=1),
    metavar="INT",
    default=MIN_INTRON_LENGTH_FOR_CLASSIFICATION,
    show_default=True,
    help="Minimal intron length for intronIC to consider for classification",
)
@binary_options.option(
    "--twobittofa_binary",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,
    show_default=True,
    help="""A path to UCSC twoBitToFa executable; if not set, the executable with this name 
will be sought for in $PATH""",
)
@binary_options.option(
    "--fatotwobit_binary",
    type=click.Path(exists=True),
    metavar="PATH",
    default=None,
    show_default=True,
    help="""A path to UCSC faToTwoBit executable; if not set, the executable with this name 
will be sought for in $PATH""",
)
@cesar_options.option(
    "--min_intron_length_cesar",
    "-cesar_min_l",
    type=click.IntRange(min=1),
    metavar="INT",
    default=MIN_INTRON_LENGTH_FOR_PROFILES,
    show_default=True,
    help="Minimal intron length to consider for CESAR2 profile generation",
)
@out_options.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    metavar="PATH",
    default=None,
    show_default=False,
    help="A path to save the results to [default: TOGA2_ref_annotation_<hex_code>]",
)
@out_options.option(
    "--prefix",
    type=str,
    metavar="STR",
    default=None,
    show_default=True,
    help="A prefix to prepend to the output file names.",
)
@out_options.option(
    "--keep_temporary",
    "-k",
    is_flag=True,
    default=False,
    show_default=True,
    help="If set, does not remove the temporary files directory",
)
def prepare_input(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    prepare-input: Check mandatory reference annotation files & generate additional TOGA2 input


    Input preparation comprises four consecutive steps:\n
    \t1) Reference annotation file is checked for format consistency; the code crashed if
    file format deviates from the Bed12 standard, and transcripts not suitable for use
    with TOGA2 (non-coding, frameshifted, or containing deprecated symbols in the name field) are filtered out
    and reported in the rejection log. Input transcripts can be further restricted to specific contigs (scaffolds, chromosomes, etc.),
    or certain contigs can be specifically excluded from the annotation.\n
    \t2) If isoforms file is provided, genes whose all transcripts were filtered out at the previous step are further
    removed from the isoforms file. Likewise, transcripts not mapped to any gene in the isoforms file are removed from the annotation.\n
    \t3) Unless disabled, intrones in the filtered annotation are classified into U2 and U12 spliceosomal classes
    with intronIC. The resulting file can be further used with TOGA2 --u12_file/-u12 argument to improve exon annotation.\n
    \t4) Unless disabled, and provided that introns were classified with intronIC, reference-specific CESAR2 HMM profiles are generated.
    TOGA2 provides CESAR2 profiles for selected mammalian and avian references used for the original companion dataset generation; these profiles were prepared with custom
    training data and further manually adjusted for improved performance, and acquired results suggest that human-specific profiles
    perform sufficiently well with highly distant clades (e.g., insects). However, if your reference of choice is known to have peculiar intron content
    (ultra-short intron prevalence, highly divergent splice site dinucleotide profiles), this step is highly recommended to go through.\n

    \b
    Mandatory arguments are:
    \t*REF_2BIT is a reference genome file in .2bit format. The same genome file is expected to be further used for TOGA2 runs;
    \t*REF_ANNOT_FILE is user-provided reference annotation file in Bed12 format. Each entry is expected to be a single reference species protein-coding transctipt.

    \b
    For further details on input files, intron classification adopted by TOGA2, and CESAR2 profiles, consult TOGA2 paper, GitHub Wiki, or TOGA2 cookbook reference.
    """
    from src.python.modules.input_producer import InputProducer

    InputProducer(**kwargs)


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Generate SpliceAI predictions for query assembly",
)
@mandatory.option(
    "--query_2bit",
    type=click.Path(exists=True),
    metavar="QUERY_2BIT",
    cls=DependentOption,
    required=True,
)
@spliceai_run_options.option(
    "--chunk_size",
    "-c",
    type=click.IntRange(min=1),
    metavar="INT",
    default=6000000,
    show_default=True,
    help="""Sequence chunk size for parallel SpliceAi annotation jobs, in bp. 
Each contig (scaffold, chromosome, etc.) in the genome file is split into chunks 
of the specified size to facilitate parallel computation""",
)
@spliceai_run_options.option(
    "--flank_size",
    "-f",
    type=click.IntRange(min=1),
    metavar="INT",
    default=50000,
    show_default=True,
    help="""Sequence chunk flank size, in bp. Each sequence chunk is furher added this many 
nucleotides from each side overlapping neighbouring chunks to mititgate 
potential boundary effects""",
)
@spliceai_run_options.option(
    "--min_contig_size",
    "-m",
    type=click.IntRange(min=1),
    metavar="INT",
    default=500,
    show_default=True,
    help="Minimal contig (scaffold, chromosome, etc.) size fo consider for SpliceAI annotation",
)
@spliceai_run_options.option(
    "--round_to",
    type=click.IntRange(min=1),
    metavar="INT",
    default=4,
    show_default=True,
    help="Number of decimal digits to round SpliceAI predicted probabilities to",
)
@spliceai_run_options.option(
    "--min_prob",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.001,
    show_default=True,
    help="""Minimal SpliceAI-predicted probability to consider. Values lower than this 
will not be reported in the output bigWig files""",
)
@parallel_options.option(
    "--job_number",
    "-j",
    type=click.IntRange(min=1),
    metavar="INT",
    default=500,
    show_default=True,
    help="Number of parallel jobs to split the task into",
)
@parallel_options.option(
    "--parallel_strategy",
    "-s",
    type=click.Choice(Constants.ALL_PARALLEL_EXECS),
    metavar="PARALLEL_EXECUTOR",
    default="local",
    show_default=True,
    help="""Specify the HPC strategy. By default, TOGA2 uses Nextflow to handle
parallel processes and supports, at least in theory, all Nextflow executors.
Please consult the full list of options at Nextflow help page. Note that setting
executor to "local" will parallel the processes over the local machine CPUs.\b\n
If you want to use Parasol as parallel process manager, set this option to "para"\b
If you want to implement a fully custom parallel manager strategy, modify the CustomStrategy class
and set this option to "custom" or contact the TOGA2 team
""",
)
@parallel_options.option(
    "--nextflow_exec_script",
    type=click.Path(exists=True),
    metavar="NEXTFLOW_SCRIPT_PATH",
    default=None,
    show_default=True,
    help="""A path to a user-defined Nextflow script used for parallel process execution. 
If not specified, TOGA2 will use a minimal boilerplate script instead. 
Ignored if Parasol or custom HPC strategy were specified for parallel steps""",
)
@parallel_options.option(
    "--max_number_of_retries",
    type=click.IntRange(min=1),
    metavar="INT",
    default=3,
    show_default=True,
    help="""Maximum number of retries per parallel job before reporting job failure. 
Ignored if Parasol or custom HPC strategy were specified for parallel steps""",
)
@parallel_options.option(
    "--nextflow_config_file",
    "-nc",
    type=click.Path(exists=True),
    metavar="NEXTFLOW_CONFIG_DIR",
    default=None,
    show_default=True,
    help="""A path to a custom Nextflow configuration file. Settings in this file 
will override those provided via command line""",
)
@parallel_options.option(
    "--max_parallel_time",
    "-max_t",
    type=click.IntRange(min=1),
    metavar="HOURS",
    default=24,
    show_default=True,
    help="Maximum time duration (in hours) for Nextflow parallel processes",
)
@parallel_options.option(
    "--cluster_queue_name",
    "-q",
    type=str,
    metavar="QUEUE_NAME",
    default="batch",
    show_default=True,
    help="""Cluster partition/queue name used. Default value assumes that name 
"batch" is available on your machine. Please consult your cluster 
administrator for available and recommended queues""",
)
@parallel_options.option(
    "--memory_limit",
    type=click.IntRange(min=1),
    metavar="GB",
    default=DEFAULT_MEMORY_LIMIT,
    show_default=True,
    help="Memory limit for parallel SpliceAI jobs, in GB",
)
@control_flow_options.option(
    "--resume_from",
    "-res",
    type=click.Choice(
        ("all", "prepare", "schedule", "run", "aggregate"), case_sensitive=False
    ),
    metavar="STEP",
    default="all",
    show_default=True,
    help="""If you have an unfinished run and want to resume with the same results,
        select the step from the following list:\b
        all: a placeholder for full starts the pipeline from the very beginning;\b\n
        prepare: process the input data for SpliceAI annotation;\n
        schedule: schedule parallel jobs\n
        run: run parallel jobs for SpliceAI annotation\n
        aggregate: aggregate the results of the parallel step\n""",
)
@control_flow_options.option(
    "--halt_at",
    "-halt",
    type=click.Choice(
        ("all", "prepare", "schedule", "run", "aggregate"), case_sensitive=False
    ),
    metavar="STEP",
    default="all",
    show_default=True,
    help="""Halts the pipeline at the selected step (see above). Option "all" 
implies running the pipeline to the last step""",
)
@binary_options.option(
    "--twobittofa_binary",
    type=click.Path(exists=True),
    metavar="PATH",
    default=os.path.join(BIN, "twoBitToFa"),
    show_default=True,
    help="A path to UCSC twoBitToFa binary",
)
@binary_options.option(
    "--fatotwobit_binary",
    type=click.Path(exists=True),
    metavar="PATH",
    default=os.path.join(BIN, "faToTwoBit"),
    show_default=True,
    help="A path to UCSC faToTwoBit binary",
)
@binary_options.option(
    "--wigtobigwig_binary",
    type=click.Path(exists=True),
    metavar="PATH",
    default=os.path.join(BIN, "wigToBigWig"),
    show_default=True,
    help="A path to UCSC wigToBigWigs binary",
)
@out_options.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    metavar="PATH",
    default=None,
    show_default=True,
    help="A path to save the results to [default: spliceai_<date_time>]",
)
@out_options.option(
    "--project_name",
    "-name",
    type=str,
    metavar="PROJECT_NAME",
    default="TOGA2",
    show_default=True,
    help="""A name for the current TOGA2 project. This name will be used 
as a prefix followed by run start date and time to name current runs' 
log and metadata files""",
)
@out_options.option(
    "--keep_temporary_files", "-k", is_flag=True, default=False, show_default=True
)
@verbosity_options.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls execution verbosity; if set, progress log will be repeated to stdout",
)
def spliceai(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    spliceai - Predict putative splice sites in the query assembly with SpliceAI\n
    If you notice any substantial differences from the expected results, pleae contact TOGA2 developer team.\n
    \n
    TOGA2 uses SpliceAI predictions for the query genome to improve exon annotation and record
    unique evolutionary events, such as distant splice site shifts and intron gains, in the query.
    To learn more on how TOGA2 uses SpliceAI predictions, consult `toga2.py cookbook` or GitHub Wiki page.\n
    \n
    NOTE: TOGA2 does not invoke SpliceAI during runtime and relies on predictions provided beforehand.
    If you want to improve your TOGA2 annotation results with SpliceAI data,
    please run this mode for query genome or consult GitHub Wiki page for alternative solutions.
    """
    from src.python.modules.spliceai_manager import SpliceAiManager

    SpliceAiManager(**kwargs)


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Merge complementing TOGA2 results for the same reference and query",
)
def merge(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    merge - Merge complementing TOGA2 results for the same reference and query
    WARNING: This mode is currently under development
    """


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="""Prepare an integrated TOGA2 annotation 
by combining annotation with different references""",
)
@mandatory.option(
    "--ref_data",
    type=click.Path(exists=True),
    metavar="INPUT_JSON",
    cls=DependentOption,
    required=True,
    help="Input map, in JSON format",
)
@out_options.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help="A path to a directory to save the results to",
)
@integration_options.option(
    "--accepted_statuses",
    "-l",
    type=str,
    metavar="ACCEPTED_LOSS_STATUSES",
    default=",".join(Constants.DEFAULT_LOSS_SYMBOLS),
    show_default=True,
    help="""A comma-separated list of preferential loss status symbols. 
Projections corresponding to these loss statuses will have preference 
for query gene inference and final annotation content unless no projection 
of these statuses occurs in a query locus. 
Supported symbols are: %s. Keyword ALL lets all possible statuses in."""
    % ",".join(Constants.ALL_LOSS_SYMBOLS),
)
@integration_options.option(
    "--paralog_rel_novelty_threshold",
    "-parrel",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.3,
    show_default=True,
    help=(
        "Minimal fraction of exon length not overlapping with orthologous items for paralogous "
        "sequences to be retained"
    ),
)
@integration_options.option(
    "--paralog_abs_novelty_threshold",
    "-parabs",
    type=click.IntRange(min=1, max=None),
    metavar="INT",
    default=15,
    show_default=True,
    help=(
        "Minimal number of exon bases not overlapping with orthologous items for paralogous "
        "sequences to be retained"
    ),
)
@integration_options.option(
    "--lost_rel_novelty_threshold",
    "-lostrel",
    type=click.FloatRange(min=0.0, max=1.0),
    metavar="FLOAT",
    default=0.3,
    show_default=True,
    help=(
        "Minimal fraction of exon length not overlapping with functionally intact "
        "ortholog for lost orthologs to be retained"
    ),
)
@integration_options.option(
    "--lost_abs_novelty_threshold",
    "-lostabs",
    type=click.IntRange(min=1, max=None),
    metavar="INT",
    default=15,
    show_default=True,
    help=(
        "Minimal number of exon bases not overlapping with functionally intact "
        "ortholog for lost orthologs to be retained"
    ),
)
@browser_options.option(
    "--prefix",
    "-p",
    type=str,
    metavar="UCSC_PREFIX",
    default="HLTOGA2combined",
    show_default=True,
    help="A prefix to use in UCSC browser file names for integrated annotation",
)
@browser_options.option(
    "--skip_ucsc",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, skips the integrated UCSC BigBed preparation step altogether, 
even if the reference-wise UCSC files were provided""",
)
@browser_options.option(
    "--chrom_sizes",
    type=click.Path(exists=True),
    metavar="CHROM_SIZES_FILE",
    default=None,
    show_default=True,
    help="""A path to two-column, tab-separated file containing query chromosome 
(contig, scaffold, etc.) sizes, in bp. Required if UCSC BigBed file 
preparation is requested""",
)
@binary_options.option(
    "--bigbedtobed_binary",
    type=click.Path(exists=True),
    default=None,
    show_default=True,
    help="""A path to UCSC bigBedToBed binary. 
If none provided, the code will look for the binary in bin/ directory 
and then for an available executable in $PATH""",
)
@binary_options.option(
    "--bedtobigbed_binary",
    type=click.Path(exists=True),
    metavar="BEDTOBIGBED_PATH",
    default=None,
    show_default=True,
    help="""A path to UCSC bedToBigBed executable. 
If none provided, the code will look for the binary in bin/ directory 
and then for an available executable in $PATH""",
)
@binary_options.option(
    "--ixixx_binary",
    type=click.Path(exists=True),
    metavar="IXIXX_PATH",
    default=None,
    show_default=True,
    help="""A path to UCSC ixIxx executable. 
If none provided, the code will look for the binary in bin/ directory 
and then for an available executable in $PATH""",
)
@misc_options.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls execution verbosity",
)
def integrate(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    integrate - Prepare an integrated TOGA2 annotation by combining annotation with different references.

    \n\n
    Integrated annotation is prepared from multiple independent TOGA2 runs for a single query.
    Since you are most likely interested in integrating annotation from multiple references, each of the
    individual input instances is referenced here and in the code as reference for simplcity, and
    the TOGA2 annotation produced with a given reference and the query of choice are called a run.\n
    \n
    To simplify the otherwise convoluted input structure, TOGA2 expects you to provide a run:{files}
    mapping in JSON format as a single mandatory argument. Each run must have the following mandatory input files:\n
    \t* "query_bed": a reference query file for this run ("query_annotation.bed" or "query_annotation.with_utrs.bed");\n
    \t* "exon_meta": a query exon metafata file for this run ("meta/exon_meta.tsv(.gz)")\n
    The following files are not mandatory but facilitate proper gene integration and annotation refinement:\n
    \t* "ucsc_bigbed": a BigBed file produced by TOGA2 for visualising this run's results in UCSC Browser.
    Used for both preparing the BigBed file for the combined annotation and extracting nucleotide and protein sequences
    if the respective arguments are not provided ("ucsc_browser_files/${your_prefix}.bb");\n
    \t* "protein_file": a FASTA file of query protein sequences for this run; if not provided,
    sequences will be extracted from the "ucsc_bigbed" file ("protein.fa(.gz)");\n
    \t* "nucleotide_file": a FASTA file of query nucleotide sequences for this run;
    if not provided, sequences will be extracted from the the "ucsc_bigbed" file ("nucleotide.fa(.gz)");\n
    \t* "reference_isoforms": a two-column file with gene-to-transcript mapping for the reference genome. NOTE:
    this is a REFERENCE isoforms file, potentially used as TOGA2 input, not the "query_genes.tsv"
    file produced by TOGA2 for this run.\n
    Finally, for each run you can provide "priority" specifying which runs' results are retained in case
    of a tie. If not set, priority is inferred based on the order the runs appear in the JSON file.\n
    \n
    An example JSON file looks as follows:\n
    {\n
        "hg38": {\n
            "query_bed": "hg38/TOGA2/vs_my_query/query_annotation.bed",\n
            "exon_meta": "hg38/TOGA2/vs_my_query/meta/exon_meta.tsv.gz",\n
            "reference_isoforms": "hg38/TOGA2/currentAnnotation/hg38.toga.isoforms.tsv",\n
            "protein_file": "hg38/TOGA2/vs_my_query/protein.fa.gz",\n
            "nucleotide_file": "hg38/TOGA2/vs_my_query/nucleotide.fa.gz",\n
            "ucsc_bigbed": "hg38/TOGA2/vs_my_query/ucsc_browser_files/HLTOGAannotVsHg38.bb",\n
            "priority": 1\n
        },\n
        "mm10": {\n
            "query_bed": "mm10/TOGA2/vs_my_query/query_annotation.bed",\n
            "exon_meta": "mm10/TOGA2/vs_my_query/meta/exon_meta.tsv.gz",\n
            "reference_isoforms": "mm10/TOGA2/currentAnnotation/mm10.toga.isoforms.tsv",\n
            "protein_file": "mm10/TOGA2/vs_my_query/protein.fa.gz",\n
            "nucleotide_file": "mm10/TOGA2/vs_my_query/nucleotide.fa.gz",\n
            "ucsc_bigbed": "mm10/TOGA2/vs_my_query/ucsc_browser_files/HLTOGAannotVsHg38.bb",\n
            "priority": 2\n
        },\n
    }
    """
    from src.python.modules.integrate import AnnotationIntegrator

    AnnotationIntegrator(**kwargs).run()


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Align orthologous sequences from multiple TOGA2 results",
)
@mandatory.option(
    "--input_dirs",
    type=click.File("r", lazy=True),
    metavar="INPUT_DIRS_FILE",
    cls=DependentOption,
    required=True,
)
@mandatory.option(
    "--transcript_id",
    type=str,
    metavar="TRANSCRIPT_ID",
    cls=DependentOption,
    required=True,
)
@input_options.option(
    "--exon_numbers",
    "-en",
    type=str,
    metavar="EXON_NUMS",
    default=None,
    show_default=True,
    help=(
        "Exon numbers to restrict the analysis to. "
        "If not set, all exons for the projection will be used"
    ),
)
@input_options.option(
    "--reference_exons",
    "-re",
    type=click.Path(exists=True),
    metavar="REF_EXON_STORAGE",
    default=None,
    show_default=True,
    help=("A path to 2bit storage with reference exon sequences."),
)
@input_options.option(
    "--reference_name",
    "-asref",
    type=str,
    metavar="REFERENCE_NAME",
    default=None,
    show_default=True,
    help=(
        "Name of the reference assembly. "
        "If --reference_exons were provided, this name will be applied to assembly "
        "from which reference exon sequences come. Otherwise, the provided name is sought for among "
        "the reference assemblies, and its exon sequences are used to restore reference exon phases."
    ),
)
@input_options.option(
    "--tree",
    "-t",
    type=click.Path(exists=True),
    metavar="TREE_FILE",
    default=None,
    show_default=True,
    help="""A path to the tree file to pass to the alignment command.
Relevant if PRANK aligner is selected""",
)
@loss_options.option(
    "--accepted_loss_status",
    "-l",
    type=str,
    metavar="LOSS_STATUS_LIST",
    default=None,
    show_default=True,
    help=(
        "A comma-separated list of loss statuses to consider. If set, the script will parse loss_summary.py "
        "to check whether the found orthologs comply with this loss status. If not set, all projections listed "
        "in the orthology classification file are considered"
    ),
)
@aligner_options.option(
    "--aligner",
    "-a",
    type=click.Choice(ALIGNERS_TO_USE, case_sensitive=False),
    metavar="ALIGNER_NAME",
    default=PRANK,
    show_default=True,
    help=(
        "Aligner program to use. Options are: %s. Case-insensitive"
        % ",".join(ALIGNERS_TO_USE)
    ),
)
@aligner_options.option(
    "--aligner_exe",
    type=click.Path(exists=True),
    metavar="CALLER_NAME",
    default=None,
    show_default=True,
    help=(
        "A path to your aligner of choice. If not set, the path will be inferred from user's PATH"
    ),
)
@aligner_options.option(
    "--phase_split_codons",
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        "If set, codons split between two exons are added in full to the shorter exon. "
        "For query sequences, split codons attributed to Missing or Deleted exons "
        "are fully removed from the sequence. Use this feature if you want to produce codon alignment; "
        "automatically enabled if macse2 set as alignment tool"
    )
)
@aligner_options.option(
    "--confidence_threshold",
    type=click.IntRange(min=0, max=9),
    default=6,
    show_default=True,
    help=(
        "If MUSCLE is set to be aligner of choice, bases with letter confidence "
        "below this value will be replaced with gaps in the final alignment"
    ),
)
@aligner_options.option(
    "--show_ancestors",
    "-anc",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Coerces PRANK to output ancestral sequence reconstruction. 
Does not work with other aligners.""",
)
@aligner_options.option(
    "--seed",
    type=str,
    metavar="SEED",
    default="12345",
    show_default=True,
    help="Fixed random seed value for PRANK",
)
@aligner_options.option(
    "--muscle_threads",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="BETA: Maximum number of threads for MUSCLE alignment jobs to use",
)
@out_options.option(
    "--output",
    "-o",
    type=click.File("w", lazy=True),
    metavar="OUTPUT_FILE",
    default=sys.stdout,
    show_default=False,
    help="A path to write the results to [default: stdout]",
)
@out_options.option(
    "--add_projection_names",
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        "If set, projection names are appended to the query identifiers "
        "in the FASTA file headers (NOTE: Current applies to the final "
        "nucleotide alignment file only)"
    ),
)
@out_options.option(
    "--amino_acids_output",
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help="""If set and MACSE is selected as aligner of choice, 
saves amino acid sequence alignment to the specified file""",
)
@out_options.option(
    "--path_to_ancestor_files",
    "-anc_path",
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help="""BETA: A path to an output directory containin exonwise 
ancestral sequence reconstructions""",
)
@out_options.option(
    "--confidence_scores",
    type=click.File("w", lazy=True),
    default=sys.stdout,
    show_default=False,
    help="""BETA: A path to the output file to write the column confidence scores to. Valid only 
if aligner program is set to MUSCLE""",
)
@out_options.option(
    "--tmp_dir",
    type=click.Path(exists=True),
    metavar="TMP_DIR",
    default=os.getcwd(),
    show_default=False,
    help=("A directory to store temporary files in [default: current directory]"),
)
@out_options.option(
    "--keep_tmp",
    is_flag=True,
    default=False,
    show_default=True,
    help="""If set, all temporary files (inlcuding tmp_dir, if it had not existed before the run) 
will be kept""",
)
@binary_options.option(
    "--twobit2fa",
    type=click.Path(exists=True),
    default=None,
    show_default=True,
    help="""A path to UCSC twoBitToFa executable. If not set, 
the executable will be sought for in PATH""",
)
@misc_options.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls logging verbosity",
)
def sequence_alignment(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    sequence_alignment - Align orthologous sequences from multiple TOGA2 queries.
    Spiritual successor of `extract_codon_alignment.py` companion script from TOGA1 suite.
    """
    from src.python.modules.codon_alignment import CodonAligner

    CodonAligner(**kwargs)


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Run postprocessing analysis with Postoga",
)
@click.option(
    "-td", "--togadir", required=True, type=str, help="Path to TOGA results directory"
)
@click.option(
    "-bc",
    "--by-orthology-status",
    "orthology_status",
    required=False,
    type=str,
    help="Include certain orthology classes (FI, I, PI, UL, M, PM, L, UL)",
)
@click.option(
    "-br",
    "--by-orthology-class",
    "orthology_class",
    required=False,
    type=str,
    help="Include certain orthology relationships (o2o, o2m, m2m, m2m, o2z)",
)
@click.option(
    "-bs",
    "--by-orthology-score",
    "orthology_score",
    required=False,
    type=float,
    help="Preserve orthology scores ≥ threshold (0.0–1.0)",
)
@click.option(
    "-to",
    "--to",
    type=click.Choice(["gtf", "gff", "bed"]),
    default="gtf",
    show_default=True,
    help="Conversion format for .bed file",
)
@click.option(
    "-tg",
    "--target",
    "bed_type",
    type=click.Choice(["bed", "utr"]),
    default="utr",
    show_default=True,
    help="Which .bed input file to use",
)
@click.option(
    "-bp",
    "--by-paralog-score",
    "min_paralog_score",
    required=False,
    type=float,
    help="Preserve transcripts with paralog projection probabilities ≥ score",
)
@click.option(
    "-w",
    "--with-isoforms",
    required=False,
    type=str,
    default=None,
    show_default=True,
    help="Path to custom isoform table",
)
@click.option(
    "-o",
    "--outdir",
    required=False,
    type=str,
    default=None,
    show_default=True,
    help="Path to posTOGA output directory",
)
@click.option(
    "-ext",
    "--extract",
    type=click.Choice(["query", "reference"]),
    required=False,
    default=None,
    help="Extract sequences from filtered projections",
)
@click.option(
    "-ot", "--only-table", is_flag=True, help="Only produce the toga.table file"
)
@click.option(
    "-oc",
    "--only-convert",
    is_flag=True,
    help="Only convert the toga.table file to gtf/gff",
)
@click.option(
    "-L",
    "--level",
    "log_level",
    type=click.Choice(["debug", "info", "warn", "off"]),
    default="info",
    show_default=True,
    help="Logging verbosity",
)
@click.option(
    "-d", "--depure", is_flag=True, help="Remove any trace of other postoga runs/files"
)
def postoga(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.     .g8\"""bgd     db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    postoga - Invoke Postoga for in-depth analysis of TOGA2 results
    """
    from postoga.run import TogaDir

    TogaDir(**kwargs).run()


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="List example commands for 'run' mode",
)
def cookbook() -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    cookbook - A detailed list of TOGA2 example commands & best practices.
    WARNING: This mode is currently under development, with the list of commands being gradually expanded
    """
    pass


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Generate a summary for a finished run from its config file",
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    metavar="CONFIG_FILE",
    required=True,
    default=None,
    show_default=True,
    help=(
        "A path to a TOGA2 run's configuration file (logs/project_args_<>). "
        "The respective output directory is expected to exist and contain "
        "all the default output files."
    ),
)
@click.option(
    "--config_format",
    "-f",
    type=click.Choice(Constants.CONFIG_FORMATS),
    default="tsv",
    show_default=True,
    help=("Input config file format"),
)

def summary(
    config: click.Path,
    config_format: Optional[str] = "tsv",
) -> None:
    from src.python.modules.results_checks import LogParserForSummary, SummaryStat

    kwargs: Dict[str, Any] = LogParserForSummary(config, config_format).extract_settings()
    print(SummaryStat(**kwargs).summary())


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=False,
    short_help="Test TOGA2 with companion dataset",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    metavar="PATH",
    default=Constants.DEFAULT_OUTPUT_DIR,
    show_default=True,
    help="A path to store the results in",
)
def test(output: Optional[click.Path]) -> None:
    # from src.python.modules.toga_configured import Toga2ConfiguredLauncher
    # from src.python.modules.toga_main import TogaMain
    # config_file: str = Constants.DEFAULT_CONFIG
    # override: str = f'--output {output} -v'
    # with open(config_file, 'r') as h:
    #     args: List[str] = Toga2ConfiguredLauncher(h, override=override).run()
    #     TogaMain(**args)
    from src.python.modules.defaults import DEFAULT_ARGS
    from src.python.modules.toga_main import TogaMain

    DEFAULT_ARGS["output"] = output
    DEFAULT_ARGS["no_spliceai"] = True
    cmd: str = " ".join(sys.argv)
    # DEFAULT_ARGS["no_u12_file"]
    # DEFAULT_ARGS["no_isoform_file"] = True
    TogaMain(**DEFAULT_ARGS, cmd=cmd)


@toga2.command(
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    short_help="Build orthogroups from TOGA2 pairwise orthology annotations",
)
@mandatory.option(
    "--toga_dir",
    "-t",
    type=click.Path(exists=True),
    metavar="DIR",
    cls=DependentOption,
    required=True,
    help="Directory with per-species TOGA2 output subdirectories",
)
@mandatory.option(
    "--species_list",
    "-s",
    "species_list_path",
    type=click.Path(exists=True),
    metavar="FILE",
    cls=DependentOption,
    required=True,
    help="Newline-separated list of species names",
)
@mandatory.option(
    "--transcripts_bed",
    "-b",
    type=click.Path(exists=True),
    metavar="FILE",
    cls=DependentOption,
    required=True,
    help="Reference transcript BED file (used to filter out sex-chromosome transcripts)",
)
@mandatory.option(
    "--isoforms",
    "-i",
    type=click.Path(exists=True),
    metavar="FILE",
    cls=DependentOption,
    required=True,
    help="TOGA2 isoforms file (gene-to-transcript mapping)",
)
@out_options.option(
    "--output",
    "-o",
    "out_dir",
    type=click.Path(exists=False),
    metavar="DIR",
    required=True,
    help="Output directory",
)
@input_options.option(
    "--panther",
    type=click.Path(exists=True),
    metavar="FILE",
    default=None,
    show_default=True,
    help="PANTHER database flat file; if provided, enables PANTHER-guided family merging",
)
@control_flow_options.option(
    "--force",
    "-f",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite output files if they already exist",
)
@control_flow_options.option(
    "--one-to-one",
    is_flag=True,
    default=False,
    show_default=True,
    help="Only write list of one-to-one orthologs",
)
@misc_options.option(
    "--include_ul",
    "-ul",
    is_flag=True,
    default=False,
    show_default=True,
    help="Include UL (Uncertain Loss) transcripts when filtering by loss status",
)
@qc_options.option(
    "--span_z",
    type=float,
    metavar="FLOAT",
    default=3.0,
    show_default=True,
    help="SpanZ threshold for spanning-rate outlier detection (default: 3.0)",
)
@qc_options.option(
    "--fam_z",
    type=float,
    metavar="FLOAT",
    default=3.0,
    show_default=True,
    help="FamZ threshold for copy-number family outlier detection (default: 3.0)",
)
@qc_options.option(
    "--no_qc",
    is_flag=True,
    default=False,
    show_default=True,
    help="Skip species QC diagnostics",
)
@verbosity_options.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Print per-species processing stats (enables DEBUG logging)",
)
def orthogroups(**kwargs) -> None:
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    orthogroups - Build orthogroups from TOGA2 pairwise orthology annotations.

    \b
    Constructs orthogroups by modelling (reference_gene - query_gene) orthology
    relationships as a graph and extracting connected components with Union-Find.
    Optionally merges with PANTHER family assignments for broader family grouping.

    \b
    Output files (written to --output directory):
      orthogroups_matrix.tsv  — copy-number count table (CAFE5-compatible)
      orthogroups_map.tsv     — full orthogroup membership per family
      one2one.lst             - list of one-to-one orthologs across all species (optional)
    """
    from src.python.modules.toga2orthogroups import run as _run

    _log = logging.getLogger("src.python.modules.toga2orthogroups")
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(message)s"))
    _log.addHandler(_handler)
    _log.propagate = False
    if kwargs.pop("verbose", False):
        _log.setLevel(logging.DEBUG)
    _run(**kwargs)


if __name__ == "__main__":
    toga2()
