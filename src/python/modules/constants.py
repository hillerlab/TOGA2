"""Class holding all project-wide constants."""

from logging import Formatter
from typing import Dict, List, Tuple

import os


__author__ = 'Yury V. Malovichko'
__credits__ = ('Bogdan M. Kirilenko')


class Constants:
    LOCATION = os.path.dirname(__file__)

    BINARIES_TO_CHECK: Dict[str, str] = {
        'bedtobigbed_binary': 'bedToBigBed',
        'bigwig2wig_binary': 'bigWigToWig',
        'fatotwobit_binary': 'faToTwoBit',
        'twobittofa_binary': 'twoBitToFa',
        'mailx_binary': 'mailx'
    }

    SPLICEAI_FILES: Tuple[str] = (
        'spliceAiDonorPlus.bw', 'spliceAiAcceptorPlus.bw',
        'spliceAiDonorMinus.bw', 'spliceAiAcceptorMinus.bw'
    )

    MEM_FILE: str = 'preprocessing_report'

    PROJECTION_OUTPUT: Tuple[str] = (
        'query_annotation.with_discarded_exons.bed', 'query_annotation.bed'
    )
    CESAR_REJECTION_LOG: str = 'genes_rejection_reason.tsv'
    CDS_FASTA: str = 'nucleotide.fa'
    CODON_ALN: str = 'codon_aln.fa'
    EXON_ALN: str = 'exon_aln.fa'
    EXON_META: str = 'exon_meta.tsv'
    GAINED_INTRONS: str = 'gained_intron_summary.tsv'
    MUTATIONS: str = 'inactivating_mutations.tsv'
    QUERY_BED_CLEAN: str = 'query_annotation.bed'
    QUERY_BED_RAW: str = 'query_annotation.with_discarded_exons.bed'
    PROT_ALN: str = 'protein_aln.fa'
    SELENO_CODONS: str = 'selenocysteine_codons.tsv'
    SPLICE_SITES: str = 'splice_sites.tsv'
    SPLICE_SITE_SHIFTS: str = 'splice_site_shifts.tsv'
    TRANSCRIPT_META: str = 'transcript_meta.tsv'
    UCSC_STUB: str = 'query_annotation.for_browser.bed'
    CESAR_OUT_FILES: Tuple[str] = [
        'query_annotation.with_discarded_exons.bed', 'query_annotation.bed', 
        'query_annotation.for_browser.bed', 'transcript_meta.tsv',
        'exon_meta.tsv', 'inactivating_mutations.tsv', 'exon_aln.fa', 'codon_aln.fa',
        'protein_aln.fa', 'nucleotide.fa',
        'splice_sites.tsv', 'genes_rejection_reason.tsv',
        'gained_intron_summary.tsv', 'splice_site_shifts.tsv',
        'selenocysteine_codons.tsv'
    ]
    CESAR_FILE_TO_DEST: Dict[str, str] = {
        CESAR_REJECTION_LOG: 'alignment_rejection_log',
        CDS_FASTA: 'cds_fasta',
        CODON_ALN: 'codon_fasta',
        EXON_ALN: 'exon_fasta',
        EXON_META: 'query_exon_meta',
        GAINED_INTRONS: 'gained_intron_summary',
        MUTATIONS: 'mutation_report',
        PROT_ALN: 'aa_fasta',
        QUERY_BED_CLEAN: 'query_annotation_filt',
        QUERY_BED_RAW: 'query_annotation_raw',
        SELENO_CODONS: 'selenocysteine_codons',
        SPLICE_SITES: 'splice_sites',
        SPLICE_SITE_SHIFTS: 'splice_site_shifts',
        TRANSCRIPT_META: 'transcript_meta',
        UCSC_STUB: 'aggr_ucsc_stub',
    }
    FINAL_UCSC_FILES: Tuple[str] = (
        '{}.bb', '{}.ix', '{}.ixx'
    ) ## TODO: Devise file naming
    SCORE_CORRECTION_CMD: str = (
        "awk -F'\t' 'BEGIN{{OFS=\"\t\"}}{{$5=1000; print $0}}' {} > {}"
    )
    ALL_LOSS_SYMBOLS: Tuple[str] = ('FI', 'I', 'PI', 'UL', 'M', 'L', 'PG', 'PP', 'N')
    DEFAULT_LOSS_SYMBOLS: Tuple[str] = ('FI', 'I', 'PI', 'UL')

    FORMATTER: Formatter = Formatter(
        '[{asctime}][{filename}] - {levelname}: {message}',
        style='{',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    RESUME_OPTIONS: List[str] = [
        'all', 'setup', 'feature_extraction', 'classification', 
        'preprocessing', 'aggregate_preprocessing_res',
        'alignment', 'aggregate_cesar_res',
        'loss_summary', 'orthology', 'summarize_trees',
        'finalize', 'ucsc_report'
    ]

    RESUME_ORDER: Dict[str, int] = {x:i for i,x in enumerate(RESUME_OPTIONS)}
    CESAR_AGGREGATION_RANK: int = RESUME_ORDER['aggregate_cesar_res']
    OK_FILE: str = '.ok'

    CLEANUP_TARGETS: Dict[str, Tuple[str]] = {
        'all': ('meta', 'nextflow_dir', 'tmp', 'ucsc_dir'),
        'setup': ('meta', 'nextflow_dir', 'tmp', 'ucsc_dir'),
        'feature_extraction': (
            'feature_job_dir', 'feature_data_dir', 'feature_res_dir', 
            'feature_rejection_log', 'feature_table'
        ),
        'classification': (
            'classification_dir', 'pred_scores', 'tr2chain_classes'
        ),
        'preprocessing': (
            'preprocessing_job_dir', 'preprocessing_res_dir',
            'paralog_report', 'processed_pseudogene_report', 'preprocessing_rejection_log',
            'spanning_chain_coords',
            'preprocessing_report', 'fragmented_projection_list'
        ),
        'aggregate_preprocessing_res': (
            'preprocessing_rejection_log', 'preprocessing_report', 'spanning_chain_coords'
        ),
        'alignment': (
            'alignment_job_dir', 'alignment_res_dir'
        ),
        'aggregate_cesar_res': (
            'annot_dir', 'vis_input_dir',
            'query_annotation_raw', 'query_annotation_filt',
            'transcript_meta', 'query_exon_meta', 
            'splice_sites', 'mutation_report',
            'aa_fasta', 'cds_fasta', 'codon_fasta', 
            'exon_fasta', 'exon_2bit',
            'codon_gzip', 'exon_gzip', 'prot_gzip', 'cds_gzip',
            'splice_sites_gzip', 'exon_meta_gzip', 'transcript_meta_gzip',
            'alignment_rejection_log', 'gained_intron_summary', 
            'splice_site_shifts', 'selenocysteine_codons'
        ),
        'loss_summary': (
            'final_rejection_log', 'gene_loss_summary', 'pseudogene_annotation'
        ),
        'orthology': (
            'query_genes_raw', 'query_genes_bed_raw', 'orthology_resolution_dir', 
            'discarded_overextended_projections', 
            'rejected_by_graph', 'weak_ortholog_names',
            'orthology_job_dir', 'orthology_input_dir', 'orthology_res_dir',
            'resolved_leaves_file', 'orth_resolution_report', 'one2zero_genes',
            'aa_hdf5', 'orth_resolution_raw'
        ),
        'summarize_trees': (
            'resolved_leaves_file', 'orth_resolution_report', 'one2zero_genes'
        ),
        'finalize': (
            'query_annotation_final', 'query_annotation_with_utrs', 'processed_pseudogene_annotation',
            'finalized_output_dir', 'query_genes', 'query_genes_bed',
            # 'cds_gzip', 'codon_gzip', 'exon_gzip', 'prot_gzip',
            # 'splice_sites_gzip', 'exon_meta_gzip', 'transcript_meta_gzip'
        ),
        'ucsc_report': (
            'ucsc_dir',
        )
    }

    FILE2HEADER: Dict[str, str] = {
        'selenocysteine_codons': 'SELENO_HEADER',
        'gained_intron_summary': 'GAINED_INTRON_HEADER',
        'query_exon_meta': 'EXON_META_HEADER',
        'preprocessing_report': 'MEM_FILE_HEADER',
        'mutation_report': 'MUT_FILE_HEADER',
        'transcript_meta': 'TRANSCRIPT_META_HEADER',
        'splice_sites': 'SPLICE_SITE_HEADER',
        'splice_site_shifts': 'SPLICE_SHIFT_HEADER',
        'gene_loss_summary': 'LOSS_FILE_HEADER',
        'final_rejection_log': 'REJ_LOG_HEADER',
        'fragmented_projection_list': 'FRAGM_PROJ_HEADER',
        'resolved_leaves_file': 'RESOLVED_LEAVES_HEADER',
        'spanning_chain_coords': 'SPANNING_CHAIN_HEADER'
    }

    FILES_TO_GZIP: Tuple[str] = (
        'aa_fasta', 'cds_fasta', 'codon_fasta', 'exon_fasta',
        'transcript_meta', 'query_exon_meta', 'splice_sites'
    )


    U12_FILE_COLS = 3
    U12_AD_FIELD = {'A', 'D'}
    ISOFORMS_FILE_COLS: int = 2
    NF_DIR_NAME = 'nextflow_logs'
    NEXTFLOW = 'nextflow'
    CESAR_PUSH_INTERVAL: int = 30  # CESAR jobs push interval
    ITER_DURATION: int = 60  # CESAR jobs check interval
    UTF8: str = 'utf-8'

    NEXTFLOW_SUPPORTED_EXECS: Tuple[str] = (
        'awsbatch', 'azurebatch', 'bridge', 'flux',
        'google-batch', 'condor', 'hq', 'k8s',
        'local', 'lsf', 'moab', 'nqsii', 'oar',
        'pbs', 'pbspro', 'sge', 'slurm'
    )
    ALL_PARALLEL_EXECS: Tuple[str] = (*NEXTFLOW_SUPPORTED_EXECS, 'para', 'custom')
    UNIQUE_CONFIGS: Dict[str, str] = {
        'preprocessing': 'preprocessing.nf', 
        'orthology': 'orthology.nf' 
    }
    ALN_CONFIG: str = 'alignment_{}.nf'

    NEXTFLOW_STUB: str = """#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.joblist = 'NONE'  // file containing jobs

if (params.joblist == "NONE"){{
    println("Usage: nextflow execute_joblist.nf  --joblist [joblist file] -c [config file]")
    System.exit(2);
}}

lines = Channel.fromPath(params.joblist).splitText()

process execute_jobs {{

    errorStrategy 'retry'
    maxRetries {}

    input:
    val line

    // one line represents an independent command
    script:
    \"\"\"
    ${{line}}
    \"\"\"
}}

workflow {{
    execute_jobs(lines)
}}"""

    NUM_CESAR_MEM_PRECOMP_JOBS = 500
    PARA_STRATEGIES = ("nextflow", "para", "custom")  # TODO: add snakemake

    TEMP_CHAIN_CLASS = "temp_chain_trans_class"
    MODULES_DIR = "modules"
    RUNNING = "RUNNING"
    CRASHED = "CRASHED"
    TEMP = "temp"

    SETUP: str = 'set -eu; set -o pipefail;'

    ## projection classification features and thresholds
    SE_MODEL_FEATURES: List[str] = ['gl_exo', 'flank_cov', 'exon_perc', 'synt_log']
    ME_MODEL_FEATURES: List[str] = ['gl_exo', 'loc_exo', 'flank_cov', 'synt_log', 'intr_perc']
    LD_MODEL_FEATURES: List[str] = [
        'gl_exo', 'flank_cov', 'exon_perc', 'synt_log', 'loc_exo', 'intr_perc', 'score', 'single_exon'
    ]
    PP_FEATURES: List[str] = ['clipped_exon_qlen', 'clipped_intr_cover']
    PP_CLIPPED_EXON_QLEN: float = 0.3
    PP_CLIPPED_INTRON_QLEN: float = 0.1

    # from CESAR_wrapper.py #
    FRAGMENT_CHAIN_ID = -1
    ORTH_LOC_LINE_SUFFIX = "#ORTHLOC"
    UNDEF_REGION = "None:0-0"

    # Sequence related #
    ATG_CODON = "ATG"
    XXX_CODON = "XXX"
    GAP_CODON = "---"
    NNN_CODON = "NNN"
    STOP_CODONS = {"TAG", "TGA", "TAA"}

    ACCEPTOR_SITE = ("ag",)
    DONOR_SITE = (
        "gt",
        "gc",
    )

    DEFAULT_UCSC_PREFIX: str = 'HLTOGAannot'

    MAILX_TEMPLATE: str = 'echo -e "{}" | {} -s "{}" {}'
    SUCCESS_EMAIL_HEADER: str = '{} - Success'
    SUCCESS_EMAIL: str = (
        'This is an automated notification on TOGA2 project {} hosted at directory {} having finished successfully'
    )
    PARTIAL_RUN_NOTE: str = '\nThe run finished before the "{}" as requested by the user'

    SANITY_CHECK_HEADER: str = ''
    SANITY_CHECK_EMAIL: str = ''

    CRASH_HEADER: str = '{} - Crashed!'
    CRASH_EMAIL: str = """
    This is an automated notification on TOGA2 project {} hosted at directory {} having crashed with the following error:\n{}
    """


class ConstColors:
    BLUE = "0,0,200"
    LIGHT_BLUE = "0,200,255"
    LIGHT_RED = "255,50,50"
    SALMON = "255,160,120"
    GREY = "130,130,130"
    BROWN = "159,129,112"
    BLACK = "10,10,10"


class InactMutClassesConst:
    MISS_EXON = "Missing exon"
    DEL_EXON = "Deleted exon"
    DEL_MISS = {MISS_EXON, DEL_EXON}
    COMPENSATION = "COMPENSATION"
    SSM = "SSM"
    # (ag)acceptor-EXON-donor(gt)
    SSM_D = "SSMD"  # Donor, right, GT,GC
    SSM_A = "SSMA"  # Acceptor, left, AG

    START_MISSING = "START_MISSING"
    ATG = "ATG"
    FS_DEL = "FS_DEL"
    FS_INS = "FS_INS"
    BIG_DEL = "BIG_DEL"
    BIG_INS = "BIG_INS"
    STOP = "STOP"

    STOPS = {"TAG", "TAA", "TGA"}
    D_M = {"D", "M"}
    LEFT_SPLICE_CORR = ("ag",)  # acceptor
    RIGHT_SPLICE_CORR = (
        "gt",
        "gc",
    )  # donor
    LEFT_SSID = 0
    RIGHT_SSID = 1
    ACCEPTOR = 0
    DONOR = 1

    BIG_INDEL_SIZE = 50
    SAFE_EXON_DEL_SIZE = 40  # actually 39
    FIRST_LAST_DEL_SIZE = 20
    BIG_EXON_THR = BIG_INDEL_SIZE * 5


class Headers:
    EXON_META_HEADER: str = '\t'.join(
        (
            'projection', 'exon', 'chain', 
            'chrom', 'start', 'end', 'strand',
            'exon_presence', 'was_aligned', 'alignment_class',
            'start_from_cesar', 'end_from_cesar', 
            'acceptor_support', 'acceptor_prob',
            'donor_support', 'donor_prob',
            'expected_locus', 'found_in_expected_locus',
            'assembly_gaps', 'spanned_by_chain', 
            'nuc_id%', 'blosum_id%'
        )
    ) + '\n'
    FRAGM_PROJ_HEADER: str = '\t'.join(('transcript', 'chains')) + '\n'
    GAINED_INTRON_HEADER: str = '\t'.join(
        (
            'projection', 'exon', 'big_gap_support', 'mutation_support', 'acceptor_prob', 'donor_prob'
        )
    ) + '\n'
    LOSS_FILE_HEADER: str = '\t'.join(('level', 'entry', 'status')) + '\n'
    MEM_FILE_HEADER: str = '\t'.join(
        (
            'transcript', 'chain', 'max_mem', 'sum_mem', 
            'largest_target', 'largest_query', 
            'chrom', 'locus_start', 'locus_end', 
            'init_start', 'init_end', 'chain_start', 'chain_end',
            'batch_path'
        )
    ) + '\n'
    MUT_FILE_HEADER: str = '\t'.join(
        (
            'projection', 'exon', 'ref_codon', 'triplet',
            'chrom', 'start', 'end', 'type', 'description',
            'is_masked', 'masking_reason', 'mut_id'
        )
    ) + '\n'
    ORTHOLOGY_TABLE_HEADER: str = 't_gene\tt_transcript\tq_gene\tq_transcript\torthology_class\n'
    QUERY_GENE_HEADER: str = '\t'.join(
        ('query_gene', 'projection')
    ) + '\n'
    REJ_LOG_HEADER: str = '\t'.join(
        ('level', 'item', 'segment', 'rejection_reason', 'rej_id', 'loss_status')
    ) + '\n'
    RESOLVED_LEAVES_HEADER: str = '\t'.join(('reference', 'query')) + '\n'
    SELENO_HEADER: str = '\t'.join(
        (
            'projection', 'exon', 'codon_num', 'chrom', 'start', 'end', 'query_codon'
        )
    ) + '\n'
    SPANNING_CHAIN_HEADER: str = '\t'.join(
        ('projection', 'ref_chrom', 'ref_start', 'ref_end')
    ) + '\n'
    SPLICE_SITE_HEADER: str = '\t'.join(
        (
            'projection', 'exon', 'acceptor', 'donor'
        )
    ) + '\n'
    SPLICE_SHIFT_HEADER: str = '\t'.join(
        (
            'projection', 'exon', 'exon_coords', 'strand',
            'site', 'shift', 'intron_type', 'spliceai_prob', 
            'dinucleotide'
        )
    ) + '\n'
    TRANSCRIPT_META_HEADER: str = '\t'.join(
        (
            'projection', 'loss_status', 'nuc_id%', 'blosum_id%', 
            'longest_intact_fraction', 'longest_nondeleted_fraction',
            'total_intact_fraction', 'middle_80%_intact', 'middle_80%_present'
        )
    ) + '\n'
    TREE_SUMMARY_HEADER: str = '\t'.join(
        (
            'batch', 'clique', '#genes', '#resolved_pairs', 'model'
        )
    ) + '\n'


# Standalone constants #

TOGA2_EPILOG: str = """\b

For detailed explanation of TOGA2 options and example commands, run 'toga2.py cookbook'
"""

BEST_PRACTICES: str = """\bExamples:

    \b
    Minimal functionality command; run TOGA2 with ref.2bit as reference and query.2bit as query, saving the results to toga2_run_%H:%M_%d.%m.%y
    \ttoga2.py ref.2bit query.2bit chains.chain.gz ref_annotation.bed\n

    \b
    Save the results to ./toga_results directory:
    \ttoga2.py ref.2bit query.2bit chains.chain.gz ref_annotation.bed -o toga_results \n

    \b
    Provide reference gene-to-transcript table to summarize orthology results at the gene level:
    \ttoga2.py ref.2bit query.2bit chains.chain.gz ref_annotation.bed -i ref_isoforms.tsv\n

    \b
    Terminate TOGA2 run before the gene loss summary step:
    \ttoga2.py ref.2bit query.2bit chains.chain.gz ref_annotation.bed --halt_at loss_summary\n

    \b
    Resume a halted/failed TOGA2 run stored at ./failed_run directory starting from the alignment step:
    \ttoga2.py ref.2bit query.2bit chains.chain.gz ref_annotation.bed --resume_from alignment\n

    \b
    HERE GOES SPLICEAI COMMAND EXAMPLE

    \b
    HERE GOES INTRON CHECK EXAMPLE

    \b
    HERE GOES 

    """

PRE_CLEANUP_LINE: str = 'rm -rf {}/*'
IQTREE_ACCEPTED_MODELS: str = ','.join(
    (
        'JTT', 'WAG', 'JTTDCMut', 'Q.LG', 'Q.pfam', 'Q.pfam_gb', 
        'Q.mammal'#, 'Q.bird', 'Q.insect', 'Q.plant', 'Q.yeast'
    )
)
PHYLO_NOT_FOUND: str = '{} was not found in PATH, with no defaults'


COMPLEMENT_BASE: Dict[str, str] = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "G",
    "n": "n",
}


GENETIC_CODE: Dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "---": "-",
    "NNN": "X",
}


## Slots for command line managers
TOGA2_SLOTS: Tuple[str] = (
    'ref_2bit', 'query_2bit', 'chain_file', 'ref_annotation',
    'resume_from', 'halt_at',
    'selected_feature_batches', 'selected_preprocessing_batches', 
    'selected_alignment_batches', 'skip_utr',
    'isoform_file', 'u12_file', 'spliceai_dir', 
    'min_chain_score', 'min_orth_chain_score',
    'feature_job_num', 'orthology_threshold', 
    'se_model', 'me_model', 'use_ld_model', 'ld_model',
    'disable_fragment_assembly', 'orthologs_only', 'one2ones_only',
    'enable_spanning_chains', 'annotate_ppgenes',
    'preprocessing_job_num', 'max_chains_per_transcript', 'cesar_memory_limit',
    'max_search_space_size', 'extrapolation_modifier', 'minimal_covered_fraction',
    'exon_locus_flank', 'assembly_gap_size', 
    'bigwig2wig_binary','bedtobigbed_binary', 'fatotwobit_binary', 'twobittofa_binary',
    'min_splice_prob', 'splice_prob_margin',
    'intron_gain_check', 'min_intron_gain_score', 
    'min_intron_prob_gapped', 'min_intron_prob_ungapped',
    'min_intron_prob_trusted', 'min_intron_prob_supported', 'min_intron_prob_unsupported',
    'cesar_binary', 'cesar_memory_bins', 'job_nums_per_bin', 'allow_heavy_jobs',
    'matrix_file', 'mask_terminal_mutations', 'leave_missing_stop', 'consider_alt_frame',

    'spliceai_correction_mode',

    'cesar_canon_u2_acceptor', 'cesar_canon_u2_donor',
    'cesar_non_canon_u2_acceptor', 'cesar_non_canon_u2_donor',
    'cesar_canon_u12_acceptor', 'cesar_canon_u12_donor',
    'cesar_non_canon_u12_acceptor', 'cesar_non_canon_u12_donor',
    'cesar_first_acceptor', 'cesar_last_donor', 'separate_site_treat',

    'accepted_loss_symbols', 'skip_tree_resolver', 'max_clique_size', 'use_raxml',
    'orth_job_num', 'prank_bin', 'tree_bin', 'tree_cpus',

    'utr_abs_threshold', 'utr_rel_threshold', 'no_utr_extrapolation',
    'no_adjacent_utrs', 'fixed_adjacent_utrs',

    'ref_link_file', 'ucsc_prefix',


    'parallel_strategy', 'nextflow_exec_script', 
    'max_number_of_retries', 'nextflow_config_dir',
    'max_parallel_time', 
    'keep_nextflow_log',
    'output', 'keep_tmp', 
    'v', 'email', 'mailx_binary',

    'toga1', 'toga1_plus_cesar',

    'tmp', 'logs', 'meta', 'ucsc_dir', 
    'nextflow_dir', 'arg_file', 'log_file', 'failed_batches_file',

    'project_name', 'local_executor', 'logger', 
    'cluster_queue_name', 'parallel_process_names', 'ignore_crashed_parallel_batches',
    'legacy_chain_feature_extraction',

    'input_data',
    'bed_file_copy', 'ref_cds_unfilt', 'cds_bed_file', 'prefiltered_transcripts',
    'chain_file_copy', 'chain_index', 'chain_index_txt',
    'se_model', 'me_model', 'ld_model', 
    'full_bed_hdf5', 'cds_bed_hdf5', 'u12_hdf5',
    'ref_contig_size_file', 'query_contig_size_file',
    'feature_table', 'feature_rejection_log',
    'tr2chain_classes', 'pred_scores', 'class_rejection_log',
    'fragmented_projection_list', 'weak_ortholog_names',
    'discarded_overextended_projections', 'rejected_by_graph',
    'preprocessing_report', 'missing_transcripts', 'resolved_leaves_file',
    'unresolved_clades_file', 'temporary_orth_report',
    'preprocessing_rejection_log', 'spanning_chain_coords',
    'processed_pseudogene_report', 'paralog_report',
    'cesar_job_list_summary', 'redundancy_rejection_log', 
    'alignment_rejection_log', 'gene_inference_rejection_log',
    'redundant_paralogs', 'redundant_ppgenes', 'discarded_proj_bed', 'orth_resolution_raw',

    'transcript_meta', 'query_annotation_raw', 'query_annotation_filt', 
    'final_rejection_log', 'gene_loss_summary',
    'query_exon_meta', 'tree_summary_table',
    'aa_fasta', 'cds_fasta', 'codon_fasta', 
    'exon_fasta', 'exon_2bit',
    'query_genes_raw', 'query_genes_bed_raw', 
    'query_genes', 'query_genes_bed',
    'one2zero_genes', 'orth_resolution_report',
    'mutation_report', 'splice_sites', 
    'gained_intron_summary', 'splice_site_shifts', 'selenocysteine_codons',
    'pseudogene_annotation', 'aa_hdf5',
    'cds_gzip', 'codon_gzip', 'exon_gzip', 'prot_gzip',
    'splice_sites_gzip', 'exon_meta_gzip', 'transcript_meta_gzip',
    'query_annotation_final', 'query_annotation_with_utrs',
    'processed_pseudogene_annotation',

    'feature_job_dir', 'feature_data_dir', 'feature_res_dir',
    'preprocessing_job_dir', 'preprocessing_res_dir',
    'alignment_job_dir', 'alignment_res_dir', 'orthology_resolution_dir',
    'orthology_job_dir', 'orthology_input_dir', 'orthology_res_dir',
    'orthology_results_dir',
    'rejection_dir', 'classification_dir', 'vis_input_dir', 'finalized_output_dir', 
    'aggr_ucsc_stub', 'all_deprecated_projs',
    'annot_dir',

    'feature_extraction_joblist', 'cesar_preprocess_joblist', 'cesar_align_joblist',
    'orth_resolution_joblist',

    'nextflow_config_files',

    'failed_feature_batches', 'failed_preprocessing_batches', 'failed_alignment_batches',
    'failed_orthology_batches',

    'CHAIN_FILTER_SCRIPT', 'INDEX_CHAIN_SCRIPT',
    'REF_BED_FILTER', 'CDS_TRACK_SCRIPT', 'REF_BED_TO_HDF5', 'U12_TO_HDF5_SCRIPT',
    'CONTIG_SIZE_SCRIPT', 
    'FEATURE_EXTRACTOR',
    'MODEL_TRAINER', 'FINAL_RESOLVER_SCRIPT',
    'UTR_PROJECTOR_SCRIPT', 'SCHEMA_FILE'
)

TOGA2_SLOT2ARG: Dict[str, str] = {
    'ref_2bit': 'ref_2bit',
    'query_2bit': 'query_2bit',
    'chain_file': 'chain_file',
    'ref_annotation': 'ref_annotation',
    'resume_from': 'resume_from',
    'halt_at': 'halt_at',
    'selected_feature_batches': 'selected_feature_batches', 
    'selected_preprocessing_batches': 'selected_preprocessing_batches', 
    'selected_alignment_batches': 'selected_alignment_batches',
    'skip_utr': 'no_utr_annotation',
    'isoform_file': 'isoform_file',
    'u12_file': 'u12_file',
    'spliceai_dir': 'spliceai_dir',
    'min_chain_score': 'min_chain_score',
    'min_orth_chain_score': 'min_orthologous_chain_score',
    'feature_job_num' : 'feature_jobs',
    'orthology_threshold': 'orthology_threshold',
    'se_model': 'single_exon_model',
    'me_model': 'multi_exon_model',
    'use_ld_model': 'use_long_distance_model',
    'ld_model': 'long_distance_model',
    'disable_fragment_assembly': 'disable_fragment_assembly',
    'orthologs_only': 'orthologs_only',
    'one2ones_only': 'one2ones_only',
    'enable_spanning_chains': 'enable_spanning_chains',
    'annotate_ppgenes': 'annotate_processed_pseudogenes',
    'preprocessing_job_num': 'preprocessing_jobs',
    'max_chains_per_transcript': 'max_chains_per_transcript',
    'cesar_memory_limit': 'cesar_memory_limit',
    'max_search_space_size': 'max_search_space_size',
    'extrapolation_modifier': 'extrapolation_modifier',
    'minimal_covered_fraction': 'minimal_covered_fraction',
    'exon_locus_flank': 'exon_locus_flank',
    'assembly_gap_size': 'assembly_gap_size',
    'cesar_canon_u2_acceptor': 'cesar_canon_u2_acceptor',
    'cesar_canon_u2_donor': 'cesar_canon_u2_donor',
    'cesar_non_canon_u2_acceptor': 'cesar_non_canon_u2_acceptor',
    'cesar_non_canon_u2_donor': 'cesar_non_canon_u2_donor',
    'cesar_canon_u12_acceptor': 'cesar_canon_u12_acceptor',
    'cesar_canon_u12_donor': 'cesar_canon_u12_donor',
    'cesar_non_canon_u12_acceptor': 'cesar_non_canon_u12_acceptor',
    'cesar_non_canon_u12_donor': 'cesar_canon_u12_donor',
    'cesar_first_acceptor': 'cesar_first_acceptor',
    'cesar_last_donor': 'cesar_last_donor',
    'separate_site_treat': 'separate_splice_site_treatment',
    'bigwig2wig_binary': 'bigwig2wig_binary',
    'min_splice_prob': 'min_splice_prob',
    'splice_prob_margin': 'splice_prob_margin',
    'intron_gain_check': 'intron_gain_check',
    'min_intron_gain_score': 'intron_gain_threshold',
    # 'min_intron_prob_gapped': 'min_intron_prob_gapped',
    # 'min_intron_prob_ungapped': 'min_intron_prob_ungapped',
    'min_intron_prob_trusted': 'min_intron_prob_trusted',
    'min_intron_prob_supported': 'min_intron_prob_supported',
    'min_intron_prob_unsupported': 'min_intron_prob_unsupported',
    'cesar_binary': 'cesar_binary',
    'cesar_memory_bins': 'memory_bins',
    'job_nums_per_bin': 'job_nums_per_bin',
    'allow_heavy_jobs': 'allow_heavy_jobs',
    'matrix_file': 'matrix',
    'mask_terminal_mutations': 'mask_n_terminal_mutations',
    'leave_missing_stop': 'disable_missing_stop_search',
    'consider_alt_frame': 'account_for_alternative_frame',
    'spliceai_correction_mode': 'spliceai_correction_mode',
    'accepted_loss_symbols': 'accepted_loss_symbols',
    'skip_tree_resolver': 'skip_gene_trees',
    'use_raxml': 'use_raxml',
    'max_clique_size': 'max_clique_size',
    'orth_job_num': 'orthology_jobs',
    'prank_bin': 'prank_binary',
    'tree_bin': 'tree_binary',
    'tree_cpus': 'tree_cpus',
    'utr_abs_threshold': 'utr_abs_threshold',
    'utr_rel_threshold': 'utr_rel_threshold',
    'no_utr_extrapolation': 'no_utr_boundary_extrapolation',
    'no_adjacent_utrs': 'no_adjacent_utr_extra',
    'fixed_adjacent_utrs': 'fixed_adjacent_utr_extra',
    'ref_link_file': 'link_file',
    'parallel_strategy': 'parallel_strategy',
    'nextflow_exec_script': 'nextflow_exec_script',
    'max_number_of_retries': 'max_number_of_retries',
    'nextflow_config_dir': 'nextflow_config_dir',
    'max_parallel_time': 'max_parallel_time',
    'keep_nextflow_log': 'keep_nextflow_log',
    'cluster_queue_name': 'cluster_queue_name',
    'legacy_chain_feature_extraction': 'legacy_chain_feature_extraction',
    'toga1': 'toga1_compatible',
    'toga1_plus_cesar': 'toga1_plus_corrected_cesar',
    'output': 'output',
    'keep_tmp': 'keep_temporary_files',
    'v': 'verbose',
    'email': 'email',
    'mailx_binary': 'mailx_binary'
}
