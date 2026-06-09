"""Class holding all project-wide constants."""

import os
from logging import Formatter
from typing import Dict, List, Tuple

__author__ = "Yury V. Malovichko"
__credits__ = "Bogdan M. Kirilenko"


class Constants:
    LOCATION = os.path.dirname(__file__)

    CONFIG_FORMATS: Tuple[str] = ("json", "tsv", "yaml")
    YAML: str = "yaml"

    BINARIES_TO_CHECK: Dict[str, str] = {
        "bigbedtobed_binary": "bigBedToBed",
        "bedtobigbed_binary": "bedToBigBed",
        "bigwig2wig_binary": "bigWigToWig",
        "fatotwobit_binary": "faToTwoBit",
        "twobittofa_binary": "twoBitToFa",
        "ixixx_binary": "ixIxx",
        "prank_bin": "prank",
        "tree_bin": "iqtree2",
        "mailx_binary": "mailx",
    }

    RAXML_DEFAULT_NAME: str = "raxmlHPC-PTHREADS-AVX"

    SPLICEAI_FILES: Tuple[str, ...] = (
        "spliceAiDonorPlus.bw",
        "spliceAiAcceptorPlus.bw",
        "spliceAiDonorMinus.bw",
        "spliceAiAcceptorMinus.bw",
    )

    MEM_FILE: str = "preprocessing_report"

    PROJECTION_OUTPUT: Tuple[str, ...] = (
        "query_annotation.with_discarded_exons.bed",
        "query_annotation.bed",
    )
    CESAR_REJECTION_LOG: str = "genes_rejection_reason.tsv"
    CDS_FASTA: str = "nucleotide.fa"
    CODON_ALN: str = "codon_aln.fa"
    EXON_ALN: str = "exon_aln.fa"
    EXON_META: str = "exon_meta.tsv"
    GAINED_INTRONS: str = "gained_intron_summary.tsv"
    MUTATIONS: str = "inactivating_mutations.tsv"
    QUERY_BED_CLEAN: str = "query_annotation.bed"
    QUERY_BED_RAW: str = "query_annotation.with_discarded_exons.bed"
    PROT_ALN: str = "protein_aln.fa"
    PROT_FASTA: str = "protein.fa"
    SELENO_CODONS: str = "selenocysteine_codons.tsv"
    SPLICE_SITES: str = "splice_sites.tsv"
    SPLICE_SITE_SHIFTS: str = "splice_site_shifts.tsv"
    TRANSCRIPT_META: str = "transcript_meta.tsv"
    UCSC_STUB: str = "query_annotation.for_browser.bed"
    CESAR_OUT_FILES: List[str] = [
        "query_annotation.with_discarded_exons.bed",
        "query_annotation.bed",
        "query_annotation.for_browser.bed",
        "transcript_meta.tsv",
        "exon_meta.tsv",
        "inactivating_mutations.tsv",
        "exon_aln.fa",
        "codon_aln.fa",
        "protein_aln.fa",
        "nucleotide.fa",
        "protein.fa",
        "splice_sites.tsv",
        "genes_rejection_reason.tsv",
        "gained_intron_summary.tsv",
        "splice_site_shifts.tsv",
        "selenocysteine_codons.tsv",
    ]
    CESAR_FILE_TO_DEST: Dict[str, str] = {
        CESAR_REJECTION_LOG: "alignment_rejection_log",
        CDS_FASTA: "cds_fasta_tmp",
        CODON_ALN: "codon_fasta",
        EXON_ALN: "exon_fasta",
        EXON_META: "query_exon_meta",
        GAINED_INTRONS: "gained_intron_summary",
        MUTATIONS: "mutation_report",
        PROT_ALN: "aa_fasta",
        PROT_FASTA: "prot_fasta_tmp",
        QUERY_BED_CLEAN: "query_annotation_filt",
        QUERY_BED_RAW: "query_annotation_raw",
        SELENO_CODONS: "selenocysteine_codons",
        SPLICE_SITES: "splice_sites",
        SPLICE_SITE_SHIFTS: "splice_site_shifts",
        TRANSCRIPT_META: "transcript_meta",
        UCSC_STUB: "aggr_ucsc_stub",
    }
    FINAL_UCSC_FILES: Tuple[str, ...] = (
        "{}.bb",
        "{}.ix",
        "{}.ixx",
    )  ## TODO: Devise file naming
    SCORE_CORRECTION_CMD: str = (
        "awk -F'\t' 'BEGIN{{OFS=\"\t\"}}{{$5=1000; print $0}}' {} > {}"
    )
    ALL_LOSS_SYMBOLS: Tuple[str, ...] = (
        "FI",
        "I",
        "PI",
        "UL",
        "M",
        "L",
        "PG",
        "PP",
        "N",
    )
    DEFAULT_LOSS_SYMBOLS: Tuple[str, ...] = ("FI", "I", "PI", "UL")

    FORMATTER: Formatter = Formatter(
        "[{asctime}][{toga_module}] - {levelname}: {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    RESUME_OPTIONS: List[str] = [
        "all",
        "setup",
        "feature_extraction",
        "classification",
        "preprocessing",
        "aggregate_preprocessing_res",
        "alignment",
        "aggregate_cesar_res",
        "gene_inference",
        "loss_summary",
        "orthology",
        "summarize_trees",
        "finalize",
        "ucsc_report",
    ]

    RESUME_ORDER: Dict[str, int] = {x: i for i, x in enumerate(RESUME_OPTIONS)}
    CESAR_AGGREGATION_RANK: int = RESUME_ORDER["aggregate_cesar_res"]
    OK_FILE: str = ".ok"

    CLEANUP_TARGETS: Dict[str, Tuple[str, ...]] = {
        "all": ("meta", "nextflow_dir", "tmp", "ucsc_dir"),
        "setup": ("meta", "nextflow_dir", "tmp", "ucsc_dir"),
        "feature_extraction": (
            "feature_job_dir",
            "feature_data_dir",
            "feature_res_dir",
            "feature_rejection_log",
            "feature_table",
        ),
        "classification": ("classification_dir", "pred_scores", "tr2chain_classes"),
        "preprocessing": (
            "preprocessing_job_dir",
            "preprocessing_res_dir",
            "paralog_report",
            "processed_pseudogene_report",
            "preprocessing_rejection_log",
            "spanning_chain_coords",
            "preprocessing_report",
            "fragmented_projection_list",
        ),
        "aggregate_preprocessing_res": (
            "preprocessing_rejection_log",
            "preprocessing_report",
            "spanning_chain_coords",
        ),
        "alignment": ("alignment_job_dir", "alignment_res_dir"),
        "aggregate_cesar_res": (
            "annot_dir",
            "vis_input_dir",
            "query_annotation_raw",
            "query_annotation_filt",
            "transcript_meta",
            "query_exon_meta",
            "splice_sites",
            "mutation_report",
            "aa_fasta",
            "cds_fasta",
            "codon_fasta",
            "exon_fasta",
            "exon_2bit",
            "prot_fasta",
            "cds_fasta_tmp",
            "prot_fasta_tmp",
            "codon_gzip",
            "exon_gzip",
            "prot_gzip",
            "cds_gzip",
            "splice_sites_gzip",
            "exon_meta_gzip",
            "transcript_meta_gzip",
            "alignment_rejection_log",
            "gained_intron_summary",
            "splice_site_shifts",
            "selenocysteine_codons",
        ),
        "gene_inference": (
            "all_deprecated_projs",
            "query_genes_raw",
            "query_genes_bed_raw",
            "discarded_overextended_projections",
        ),
        "loss_summary": (
            "gene_loss_summary",
            "pseudogene_annotation",
            "loss_summary_extended",
        ),
        "orthology": (
            "orthology_resolution_dir",
            "rejected_by_graph",
            "weak_ortholog_names",
            "orthology_job_dir",
            "orthology_input_dir",
            "orthology_res_dir",
            # "resolved_leaves_file",
            # "orth_resolution_report",
            # "one2zero_genes",
            "aa_hdf5",
            "orth_resolution_raw",
        ),
        "summarize_trees": (
            "resolved_leaves_file",
            "unresolved_clades_file",
            "orth_resolution_report",
            "one2zero_genes",
            "rejected_at_tree_step",
        ),
        "finalize": (
            "query_annotation_final",
            "query_annotation_with_utrs",
            "processed_pseudogene_annotation",
            "finalized_output_dir",
            "query_genes",
            "query_genes_bed",
            "summary",
            "all_discarded_projections",
            "query_genes_for_gtf",
            "query_gtf",
            "gtf_gzip",
            # 'cds_gzip', 'codon_gzip', 'exon_gzip', 'prot_gzip',
            # 'splice_sites_gzip', 'exon_meta_gzip', 'transcript_meta_gzip'
        ),
        "ucsc_report": ("ucsc_dir",),
    }

    ## Rejection-reason-pipeline-step mapping
    REJ2STEP: Dict[str, str] = {
        "ILLEGAL_NAME": "setup",
        "REJECTED_CONTIG": "setup",
        "NON_CODING": "setup",
        "OUT_OF_FRAME": "setup",
        "CHROM_UNALIGNED": "feature_extraction",
        "TRANSCRIPT_UNALIGNED": "feature_extraction",
        "NO_PROJ": "classification",
        "PPGENE_ONLY": "classification",
        "INSUFFICIENT_CHAIN_SCORE": "classification",
        "CHAIN_LIMIT_EXCEEDED": "preprocessing",
        "EXCEEDS_MEMORY": "preprocessing",
        "EXCEEDS_MEMORY+GAP": "preprocessing",
        "EXCEEDS_SPACE": "preprocessing",
        "EXCEEDS_SPACE+GAP": "preprocessing",
        "INSUFFICIENT_SEARCH_SPACE": "preprocessing",
        "INSUFFICIENT_SEARCH_SPACE+GAP": "alignment",
        "MULTIPLE_ORTHOLOGY": "preprocessing",
        "NO_EXONS_ALIGNED": "preprocessing",
        "NO_CHAINS": "preprocessing",
        "SPANNING": "preprocessing",
        "ZERO_ORTHOLOGY": "preprocessing",
        "CHIMERIC": "alignment",
        "HEAVY": "alignment",
        "REDUNDANT": "alignment",
        "SECOND_BEST": "gene_inference",
        "REDUNDANT_PARALOG": "gene_inference",
        "REDUNDANT_PPGENE": "gene_inference",
        "ALL_PARALOGS_REDUNDANT": "loss_summary",
        "ALL_ORTHS_DISCARDED": "orthology",
        "WEAK_EDGE": "orthology",
        "GENE_TREE_REJECTION": "summarize_trees",
    }

    ## Table file headers
    FILE2HEADER: Dict[str, str] = {
        "selenocysteine_codons": "SELENO_HEADER",
        "gained_intron_summary": "GAINED_INTRON_HEADER",
        "query_exon_meta": "EXON_META_HEADER",
        "preprocessing_report": "MEM_FILE_HEADER",
        "mutation_report": "MUT_FILE_HEADER",
        "transcript_meta": "TRANSCRIPT_META_HEADER",
        "splice_sites": "SPLICE_SITE_HEADER",
        "splice_site_shifts": "SPLICE_SHIFT_HEADER",
        "gene_loss_summary": "LOSS_FILE_HEADER",
        "final_rejection_log": "REJ_LOG_HEADER",
        "fragmented_projection_list": "FRAGM_PROJ_HEADER",
        "resolved_leaves_file": "RESOLVED_LEAVES_HEADER",
        "spanning_chain_coords": "SPANNING_CHAIN_HEADER",
    }

    DISCARDED_PROJECTION_FILES: Tuple[str, ...] = (
        "paralog_report",
        "processed_pseudogene_report",
        "discarded_overextended_projections",
    )

    FILES_TO_GZIP: Tuple[str, ...] = (
        "aa_fasta",
        "cds_fasta",
        "codon_fasta",
        "exon_fasta",
        "prot_fasta",
        "query_exon_meta",
        "splice_sites",
        "transcript_meta",
    )

    U12_FILE_COLS = 3
    U12_AD_FIELD = {"A", "D"}
    ISOFORMS_FILE_COLS: int = 2
    NF_DIR_NAME = "nextflow_logs"
    NEXTFLOW = "nextflow"
    CESAR_PUSH_INTERVAL: int = 30  # CESAR jobs push interval
    ITER_DURATION: int = 60  # CESAR jobs check interval
    UTF8: str = "utf-8"

    ## Nextflow configuration constants
    NEXTFLOW_SUPPORTED_EXECS: Tuple[str, ...] = (
        "awsbatch",
        "azurebatch",
        "bridge",
        "flux",
        "google-batch",
        "condor",
        "hq",
        "k8s",
        "local",
        "lsf",
        "moab",
        "nqsii",
        "oar",
        "pbs",
        "pbspro",
        "sge",
        "slurm",
    )
    ALL_PARALLEL_EXECS: Tuple[str, ...] = (*NEXTFLOW_SUPPORTED_EXECS, "para", "custom")
    NF_EXEC_SCRIPT_NAME: str = "execute_joblist.nf"
    UNIQUE_CONFIGS: Dict[str, str] = {
        "preprocessing": "preprocessing.config",
        "orthology": "orthology.config",
    }
    ALN_CONFIG: str = "alignment_{}.config"
    NEXTFLOW_STUB: str = """#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.joblist = 'NONE'  // file containing jobs

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
    if (params.joblist == "NONE") {{
        println("Usage: nextflow execute_joblist.nf  --joblist [joblist file] -c [config file]")
        System.exit(2)
    }}
    lines = Channel.fromPath(params.joblist).splitText()
    execute_jobs(lines)
}}"""

    NUM_CESAR_MEM_PRECOMP_JOBS = 500
    PARA_STRATEGIES = ("nextflow", "para", "custom")  # TODO: add snakemake

    MODULES_DIR = "modules"
    RUNNING = "RUNNING"
    CRASHED = "CRASHED"
    TEMP = "temp"

    SETUP: str = "set -eu; set -o pipefail;"

    ## projection classification features and thresholds
    SE_MODEL_FEATURES: List[str] = ["gl_exo", "flank_cov", "exon_perc", "synt_log"]
    ME_MODEL_FEATURES: List[str] = [
        "gl_exo",
        "loc_exo",
        "flank_cov",
        "synt_log",
        "intr_perc",
    ]
    LD_MODEL_FEATURES: List[str] = [
        "gl_exo",
        "flank_cov",
        "exon_perc",
        "synt_log",
        "loc_exo",
        "intr_perc",
        "score",
        "single_exon",
    ]
    PP_FEATURES: List[str] = ["clipped_exon_qlen", "clipped_intr_cover"]
    PP_CLIPPED_EXON_QLEN: float = 0.3
    PP_CLIPPED_INTRON_QLEN: float = 0.1
    DEFAULT_ORTH_THRESHOLD: float = 0.5
    SPANNING_CHAIN_PROB: float = -1.0
    PPGENE_PROB: float = -2.0

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

    DEFAULT_UCSC_PREFIX: str = "HLTOGAannot"

    ## E-mail notification
    MAILX_TEMPLATE: str = 'echo -e "{}" | {} -s "{}" {}'
    SUCCESS_EMAIL_HEADER: str = "{} - Success"
    SUCCESS_EMAIL: str = "This is an automated notification on TOGA2 project {} hosted at directory {} having finished successfully"
    PARTIAL_RUN_NOTE: str = (
        '\nThe run finished before the "{}" as requested by the user'
    )

    SANITY_CHECK_HEADER: str = ""
    SANITY_CHECK_EMAIL: str = ""

    CRASH_HEADER: str = "{} - Crashed!"
    CRASH_EMAIL: str = """
This is an automated notification on TOGA2 project {} hosted at directory {} \
having crashed with the following error:\n{}
"""
    SANITY_CHECK_PS: str = """\
If the current run's setup has known complications \
(low quality input assemblies, large evolutionary distance between references and query species, \
specific gene/genome architecture in any of the species, etc.), disregard this e-mail. \
Otherwise, you might want to halt the current run and inspect the results and/or input files.
"""
    WARNING_SUBJECT: str = "{} - Sanity check warning at {} step"
    WARNING_EMAIL: str = (
        """
This is an automated notification on TOGA2 project {} hosted at directory {} \
reporting the following warning at the {} step:\n{}
"""
        + SANITY_CHECK_PS
    )

    ## default settings for 'test' mode
    DEFAULT_CONFIG: str = os.path.join("supply", "project_args.tsv")
    DEFAULT_OUTPUT_DIR: str = "sample_output"


## rejection reason templates
class RejectionReasons:
    REJ_GENE: str = "\t".join(
        (
            "GENE",
            "{}",
            "0",
            "No (valid) transcripts found in the reference annotation",
            "ZERO_TRANSCRIPT_INPUT",
            "N",
        )
    )
    ORPHAN_TR: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "No corresponding gene found in the isoform file",
            "ZERO_GENE_INPUT",
            "N",
        )
    )
    NAME_REJ_REASON: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "Illegal character used in transcript name",
            "ILLEGAL_NAME",
            "N",
        )
    )
    CONTIG_REJ_REASON: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0Located outside of user-preferred contigs",
            "REJECTED_CONTIG",
            "N",
        )
    )
    NON_CODING_REJ_REASON: str = "\t".join(
        ("TRANSCRIPT", "{}", "0", "Does not have a coding sequence", "NON_CODING", "N")
    )
    FRAME_REJ_REASON: str = "\t".join(
        ("TRANSCRIPT", "{}", "0", "Transcript is out of frame", "OUT_OF_FRAME", "N")
    )
    UNCOV_CHROM: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "No chain corresponds to reference chromosome",
            "CHROM_UNALIGNED",
            "M",
        )
    )
    UNALIGNED_TR: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "No chain corresponds to the transcript",
            "TRANSCRIPT_UNALIGNED",
            "M",
        )
    )
    UNCLASS_REJ_REASON: str = "\t".join(
        ("TRANSCRIPT", "{}", "0", "No classifiable projections found", "NO_PROJ", "M")
    )
    PPGENE_ONLY_REASON: str = "\t".join(
        ("TRANSCRIPT", "{}", "0", "Only processed pseudogene projections found", "PPGENE_ONLY", "M")
    )
    UNDERSCORED_REJ_REASON: str = "\t".join(
        (
            "PROJECTION",
            "{}",
            "0",
            "Chain score below set threshold ({})",
            "INSUFFICIENT_CHAIN_SCORE",
            "L",
        )
    )
    LIMIT_EXCEED_REJ: str = "\t".join(
        (
            "PROJECTION",
            "{}",
            "0",
            "Number of homologous chains exceeds the set limit ({})",
            "CHAIN_LIMIT_EXCEEDED",
            "N",
        )
    )
    MULTIPLE_ORTHOLOG_REJ: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "Multiple orthologs detected",
            "MULTIPLE_ORTHOLOGY",
            "M",
        )
    )
    NO_CHAINS_REJ: str = "\t".join(
        ("TRANSCRIPT", "{}", "0", "No covering chains detected", "NO_CHAINS", "M")
    )
    NO_ALIGNED_EXON_REJ: str = "\t".join(
        ("PROJECTION", "{}", "0", "No aligned exons found", "NO_EXONS_ALIGNED", "{}")
    )
    PREPROCESSING_REJ: str = "\t".join(("PROJECTION", "{}", "{}", "{}", "{}", "{}"))
    SPANNING_CHAIN_REASON: str = "\t".join(
        ("PROJECTION", "{}", "0", "Spanning chain", "SPANNING", "{}")
    )
    ZERO_ORTHOLOGY_REJ: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "No orthologous chains detected",
            "ZERO_ORTHOLOGY",
            "M",
        )
    )
    CHIMERIC_ENTRY: str = "\t".join(
        ("PROJECTION", "{}", "0", "Potential chimeric projection", "CHIMERIC", "N")
    )
    HEAVY_ENTRY: str = "\t".join(
        ("PROJECTION", "{}", "0", "Maximum memory limit exceeded", "HEAVY", "N")
    )
    REDUNDANT_ENTRY: str = "\t".join(
        (
            "PROJECTION",
            "{}",
            "0",
            "Redundant projection to the given locus",
            "REDUNDANT",
            "N",
        )
    )
    REJ_ORTH_REASON: str = "\t".join(
        (
            "PROJECTION",
            "{}",
            "0",
            "Insufficiently covered exons in second-best projection",
            "SECOND_BEST",
            "{}",
        )
    )
    REJ_PARA_REASON: str = "\t".join(
        (
            "PROJECTION",
            "{}",
            "0",
            "Redundant paralog overlapping orthologous projections",
            "REDUNDANT_PARALOG",
            "{}",
        )
    )
    REJ_PPGENE_REASON: str = "\t".join(
        (
            "PROJECTION",
            "{}",
            "0",
            "Processed pseudogene overlapping ortholog or paralog",
            "REDUNDANT_PPGENE",
            "{}",
        )
    )
    OUTCOMPETED_PARALOG_REASON: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "All paralogous projections outcompeted by orthologous predictions of other items",
            "ALL_PARALOGS_REDUNDANT",
            "M",
        )
    )
    REMOVED_ORTH_REASON: str = "\t".join(
        (
            "TRANSCRIPT",
            "{}",
            "0",
            "No projections reached the orthology step",
            "ALL_ORTHS_DISCARDED",
            "{}",
        )
    )
    WEAK_EDGE_REASON: str = "\t".join(
        ("PROJECTION", "{}", "0", "Weak orthology graph edge", "WEAK_EDGE", "{}")
    )
    ORTH_REJ_TEMPLATE: str = "\t".join(
        (
            "PROJECTION",
            "{}",
            "0",
            "Rejected after the gene resolution step",
            "GENE_TREE_REJECTION",
            "{}",
        )
    )


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
    EXON_META_HEADER: str = (
        "\t".join(
            (
                "projection",
                "exon",
                "chain",
                "chrom",
                "start",
                "end",
                "strand",
                "exon_presence",
                "was_aligned",
                "alignment_class",
                "start_from_cesar",
                "end_from_cesar",
                "acceptor_support",
                "acceptor_prob",
                "donor_support",
                "donor_prob",
                "expected_locus",
                "found_in_expected_locus",
                "assembly_gaps",
                "spanned_by_chain",
                "nuc_id%",
                "blosum_id%",
            )
        )
        + "\n"
    )
    FRAGM_PROJ_HEADER: str = "\t".join(("transcript", "chains")) + "\n"
    GAINED_INTRON_HEADER: str = (
        "\t".join(
            (
                "projection",
                "exon",
                "big_gap_support",
                "mutation_support",
                "acceptor_prob",
                "donor_prob",
            )
        )
        + "\n"
    )
    LOSS_FILE_HEADER: str = "\t".join(("level", "entry", "status")) + "\n"
    MEM_FILE_HEADER: str = (
        "\t".join(
            (
                "transcript",
                "chain",
                "max_mem",
                "sum_mem",
                "largest_target",
                "largest_query",
                "%covered_exons",
                "chrom",
                "locus_start",
                "locus_end",
                "init_start",
                "init_end",
                "chain_start",
                "chain_end",
                "batch_path",
            )
        )
        + "\n"
    )
    MUT_FILE_HEADER: str = (
        "\t".join(
            (
                "projection",
                "exon",
                "triplet",
                "ref_codon",
                "chrom",
                "start",
                "end",
                "type",
                "description",
                "is_masked",
                "masking_reason",
                "mut_id",
            )
        )
        + "\n"
    )
    ORTHOLOGY_TABLE_HEADER: str = (
        "t_gene\tt_transcript\tq_gene\tq_transcript\torthology_class\n"
    )
    QUERY_GENE_HEADER: str = "\t".join(("query_gene", "projection")) + "\n"
    REJ_LOG_HEADER: str = (
        "\t".join(
            ("level", "item", "segment", "rejection_reason", "rej_id", "loss_status")
        )
        + "\n"
    )
    RESOLVED_LEAVES_HEADER: str = "\t".join(("reference", "query")) + "\n"
    SELENO_HEADER: str = (
        "\t".join(
            ("projection", "exon", "codon_num", "chrom", "start", "end", "query_codon")
        )
        + "\n"
    )
    SPANNING_CHAIN_HEADER: str = (
        "\t".join(("projection", "ref_chrom", "ref_start", "ref_end")) + "\n"
    )
    SPLICE_SITE_HEADER: str = (
        "\t".join(("projection", "exon", "acceptor", "donor")) + "\n"
    )
    SPLICE_SHIFT_HEADER: str = (
        "\t".join(
            (
                "projection",
                "exon",
                "exon_coords",
                "strand",
                "site",
                "shift",
                "intron_type",
                "spliceai_prob",
                "dinucleotide",
            )
        )
        + "\n"
    )
    TRANSCRIPT_META_HEADER: str = (
        "\t".join(
            (
                "projection",
                "loss_status",
                "nuc_id%",
                "blosum_id%",
                "longest_intact_fraction",
                "longest_nondeleted_fraction",
                "total_intact_fraction",
                "middle_80%_intact",
                "middle_80%_present",
            )
        )
        + "\n"
    )
    TREE_SUMMARY_HEADER: str = (
        "\t".join(("batch", "clique", "#genes", "#resolved_pairs", "model")) + "\n"
    )


class NameTemplates:
    TWOBIT: str = os.path.join("{}", "{}.2bit")
    CHAINS: str = os.path.join(
        "{}", "lastz", "vs_{}", "axtChain", "{}.{}.allfilled.chain"
    )
    CHAINS_GZ: str = os.path.join(
        "{}", "lastz", "vs_{}", "axtChain", "{}.{}.allfilled.chain.gz"
    )
    REF_ANNOT: str = os.path.join(
        "{}", "TOGA2", "currentAnnotation", "{}.toga.transcripts.bed"
    )
    REF_ISOFORMS: str = os.path.join(
        "{}", "TOGA2", "currentAnnotation", "{}.toga.isoforms.tsv"
    )
    REF_U12: str = os.path.join(
        "{}", "TOGA2", "currentAnnotation", "{}.toga.U12introns.bed"
    )
    SPLICEAI: str = os.path.join("{}", "spliceAi")
    REF_LINKS: str = os.path.join(
        "{}", "TOGA2", "currentAnnotation", "{}.toga.links.tsv"
    )
    CESAR_PROFILE_VALUES: Dict[str, str] = {
        "cesar_canon_u2_acceptor": "canon_U2_acceptor.tsv",
        "cesar_canon_u2_donor": "canon_U2_donor.tsv",
        "cesar_non_canon_u2_acceptor": "nonCanon_U2_acceptor.tsv",
        "cesar_non_canon_u2_donor": "nonCanon_U2_donor.tsv",
        "cesar_canon_u12_acceptor": "canon_U12_acceptor.tsv",
        "cesar_canon_u12_donor": "canon_U12_donor.tsv",
        "cesar_non_canon_u12_acceptor": "nonCanon_U12_acceptor.tsv",
        "cesar_non_canon_u12_donor": "nonCanon_U12_donor.tsv",
    }


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

CONTAINER_ENGINE2BIND_KEY: Dict[str, str] = {"apptainer": "--bind", "docker": "--mount"}
PRE_CLEANUP_LINE: str = "rm -rf {}/*"
IQTREE_ACCEPTED_MODELS: str = ",".join(
    (
        "JTT",
        "WAG",
        "JTTDCMut",
        "Q.LG",
        "Q.pfam",
        "Q.pfam_gb",
        "Q.mammal",  # , 'Q.bird', 'Q.insect', 'Q.plant', 'Q.yeast'
    )
)

PHYLO_NOT_FOUND: str = "{} was not found in PATH, with no defaults"

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
TOGA2_SLOTS: Tuple[str, ...] = (
    "ref_2bit",
    "query_2bit",
    "chain_file",
    "ref_annotation",
    "isoform_file",
    "no_isoform_file",
    "u12_file",
    "no_u12_file",
    "spliceai_dir",
    "no_spliceai",
    "input_dir",
    "ref_name",
    "query_name",
    "resume_from",
    "halt_at",
    "selected_feature_batches",
    "selected_preprocessing_batches",
    "selected_alignment_batches",
    "skip_utr",
    "min_chain_score",
    "min_orth_chain_score",
    "feature_job_num",
    "orthology_threshold",
    "se_model",
    "me_model",
    "use_ld_model",
    "ld_model",
    "disable_fragment_assembly",
    "orthologs_only",
    "one2ones_only",
    "annotate_paralogs",
    "paralogs_over_spanning",
    "enable_spanning_chains",
    "annotate_ppgenes",
    "preprocessing_job_num",
    "max_chains_per_transcript",
    "cesar_memory_limit",
    "max_search_space_size",
    "extrapolation_modifier",
    "minimal_covered_fraction",
    "exon_locus_flank",
    "assembly_gap_size",
    "bigwig2wig_binary",
    "bedtobigbed_binary",
    "fatotwobit_binary",
    "twobittofa_binary",
    "ixixx_binary",
    "min_splice_prob",
    "splice_prob_margin",
    "intron_gain_check",
    "max_intron_number",
    "min_intron_gain_score",
    "min_intron_prob_gapped",
    "min_intron_prob_ungapped",
    "min_intron_prob_trusted",
    "min_intron_prob_supported",
    "min_intron_prob_unsupported",
    "cesar_binary",
    "cesar_memory_bins",
    "job_nums_per_bin",
    "allow_heavy_jobs",
    "matrix_file",
    "mask_terminal_mutations",
    "leave_missing_stop",
    "consider_alt_frame",
    "spliceai_correction_mode",
    "cesar_profile_dir",
    "cesar_canon_u2_acceptor",
    "cesar_canon_u2_donor",
    "cesar_non_canon_u2_acceptor",
    "cesar_non_canon_u2_donor",
    "cesar_canon_u12_acceptor",
    "cesar_canon_u12_donor",
    "cesar_non_canon_u12_acceptor",
    "cesar_non_canon_u12_donor",
    "cesar_first_acceptor",
    "cesar_last_donor",
    "joint_site_treat",
    "accepted_loss_symbols",
    "skip_tree_resolver",
    "max_clique_size",
    "use_raxml",
    "orth_job_num",
    "prank_bin",
    "tree_bin",
    "tree_cpus",
    "utr_abs_threshold",
    "utr_rel_threshold",
    "no_utr_extrapolation",
    "no_adjacent_utrs",
    "fixed_adjacent_utrs",
    "ref_link_file",
    "ucsc_prefix",
    "parallel_strategy",
    "nextflow_exec_script",
    "max_number_of_retries",
    "nextflow_config_dir",
    "max_parallel_time",
    "project_arg_format",
    "keep_nextflow_log",
    "output",
    "keep_tmp",
    "project_name",
    "project_id",
    "timestamp",
    "cmd",
    "v",
    "debug",
    "email",
    "mailx_binary",
    "toga1",
    "toga1_plus_cesar",
    "tmp",
    "logs",
    "meta",
    "ucsc_dir",
    "nextflow_dir",
    "arg_file",
    "log_file",
    "failed_batches_file",
    "result_checker",
    "project_name",
    "local_executor",
    "logger",
    "cluster_queue_name",
    "parallel_process_names",
    "ignore_crashed_parallel_batches",
    "legacy_chain_feature_extraction",
    "container_image",
    "container_executor",
    "bindings",
    "input_data",
    "bed_file_copy",
    "ref_cds_unfilt",
    "cds_bed_file",
    "prefiltered_transcripts",
    "chain_file_copy",
    "chain_index",
    "chain_index_txt",
    "se_model",
    "me_model",
    "ld_model",
    "full_bed_hdf5",
    "cds_bed_hdf5",
    "u12_hdf5",
    "ref_contig_size_file",
    "query_contig_size_file",
    "feature_table",
    "feature_rejection_log",
    "tr2chain_classes",
    "pred_scores",
    "class_rejection_log",
    "fragmented_projection_list",
    "weak_ortholog_names",
    "discarded_overextended_projections",
    "rejected_by_graph",
    "preprocessing_report",
    "missing_transcripts",
    "resolved_leaves_file",
    "unresolved_clades_file",
    "temporary_orth_report",
    "preprocessing_rejection_log",
    "spanning_chain_coords",
    "processed_pseudogene_report",
    "paralog_report",
    "cesar_job_list_summary",
    "redundancy_rejection_log",
    "alignment_rejection_log",
    "gene_inference_rejection_log",
    "redundant_paralogs",
    "redundant_ppgenes",
    "discarded_proj_bed",
    "orth_resolution_raw",
    "transcript_meta",
    "query_annotation_raw",
    "query_annotation_filt",
    "query_gtf",
    "postoga_table",
    "final_rejection_log",
    "gene_loss_summary",
    "loss_summary_extended",
    "query_exon_meta",
    "tree_summary_table",
    "aa_fasta",
    "cds_fasta",
    "codon_fasta",
    "exon_fasta",
    "prot_fasta",
    "exon_2bit",
    "query_genes_raw",
    "query_genes_bed_raw",
    "query_genes_for_gtf",
    "query_genes",
    "query_genes_bed",
    "one2zero_genes",
    "orth_resolution_report",
    "rejected_at_tree_step",
    "mutation_report",
    "splice_sites",
    "gained_intron_summary",
    "splice_site_shifts",
    "selenocysteine_codons",
    "pseudogene_annotation",
    "aa_hdf5",
    "aa_gzip",
    "cds_gzip",
    "codon_gzip",
    "exon_gzip",
    "prot_gzip",
    "splice_sites_gzip",
    "exon_meta_gzip",
    "transcript_meta_gzip",
    "gtf_gzip",
    "query_annotation_final",
    "query_annotation_with_utrs",
    "processed_pseudogene_annotation",
    "summary",
    "decoration_track",
    "feature_job_dir",
    "feature_data_dir",
    "feature_res_dir",
    "preprocessing_job_dir",
    "preprocessing_res_dir",
    "alignment_job_dir",
    "alignment_res_dir",
    "orthology_resolution_dir",
    "orthology_job_dir",
    "orthology_input_dir",
    "orthology_res_dir",
    "orthology_results_dir",
    "rejection_dir",
    "classification_dir",
    "vis_input_dir",
    "finalized_output_dir",
    "aggr_ucsc_stub",
    "decor_stub",
    "all_deprecated_projs",
    "annot_dir",
    "cds_fasta_tmp",
    "prot_fasta_tmp",
    "postoga_tmp",
    "postoga_table_tmp",
    "all_discarded_projections",
    "feature_extraction_joblist",
    "cesar_preprocess_joblist",
    "cesar_align_joblist",
    "orth_resolution_joblist",
    "nextflow_config_files",
    "failed_feature_batches",
    "failed_preprocessing_batches",
    "failed_alignment_batches",
    "failed_orthology_batches",
    "rejection_log_cleaned",
    "CHAIN_FILTER_SCRIPT",
    "INDEX_CHAIN_SCRIPT",
    "REF_BED_FILTER",
    "CDS_TRACK_SCRIPT",
    "REF_BED_TO_HDF5",
    "U12_TO_HDF5_SCRIPT",
    "CONTIG_SIZE_SCRIPT",
    "FEATURE_EXTRACTOR",
    "MODEL_TRAINER",
    "FINAL_RESOLVER_SCRIPT",
    "FASTA_FILTER_SCRIPT",
    "UTR_PROJECTOR_SCRIPT",
    "DECORATOR_SCRIPT",
    "GTF_SCRIPT",
    "SCHEMA_FILE",
    "DECOR_SCHEMA_FILE",
)

TOGA2_SLOT2ARG: Dict[str, str] = {
    "ref_2bit": "ref_2bit",
    "query_2bit": "query_2bit",
    "chain_file": "chain_file",
    "ref_annotation": "ref_annotation",
    "isoform_file": "isoform_file",
    "no_isoform_file": "no_isoform_file",
    "u12_file": "u12_file",
    "no_u12_file": "no_u12_file",
    "spliceai_dir": "spliceai_dir",
    "no_spliceai": "no_spliceai",
    "input_dir": "input_directory",
    "ref_name": "ref_name",
    "query_name": "query_name",
    "resume_from": "resume_from",
    "halt_at": "halt_at",
    "selected_feature_batches": "selected_feature_batches",
    "selected_preprocessing_batches": "selected_preprocessing_batches",
    "selected_alignment_batches": "selected_alignment_batches",
    "skip_utr": "no_utr_annotation",
    "min_chain_score": "min_chain_score",
    "min_orth_chain_score": "min_orthologous_chain_score",
    "feature_job_num": "feature_jobs",
    "orthology_threshold": "orthology_threshold",
    "se_model": "single_exon_model",
    "me_model": "multi_exon_model",
    "use_ld_model": "use_long_distance_model",
    "ld_model": "long_distance_model",
    "disable_fragment_assembly": "disable_fragment_assembly",
    "orthologs_only": "orthologs_only",
    "one2ones_only": "one2ones_only",
    "annotate_paralogs": "annotate_paralogs",
    "paralogs_over_spanning": "paralogs_over_spanning",
    "enable_spanning_chains": "enable_spanning_chains",
    "annotate_ppgenes": "annotate_processed_pseudogenes",
    "preprocessing_job_num": "preprocessing_jobs",
    "max_chains_per_transcript": "max_chains_per_transcript",
    "cesar_memory_limit": "cesar_memory_limit",
    "max_search_space_size": "max_search_space_size",
    "extrapolation_modifier": "extrapolation_modifier",
    "minimal_covered_fraction": "minimal_covered_fraction",
    "exon_locus_flank": "exon_locus_flank",
    "assembly_gap_size": "assembly_gap_size",
    "cesar_profile_dir": "cesar_profile_dir",
    "cesar_canon_u2_acceptor": "cesar_canon_u2_acceptor",
    "cesar_canon_u2_donor": "cesar_canon_u2_donor",
    "cesar_non_canon_u2_acceptor": "cesar_non_canon_u2_acceptor",
    "cesar_non_canon_u2_donor": "cesar_non_canon_u2_donor",
    "cesar_canon_u12_acceptor": "cesar_canon_u12_acceptor",
    "cesar_canon_u12_donor": "cesar_canon_u12_donor",
    "cesar_non_canon_u12_acceptor": "cesar_non_canon_u12_acceptor",
    "cesar_non_canon_u12_donor": "cesar_non_canon_u12_donor",
    "cesar_first_acceptor": "cesar_first_acceptor",
    "cesar_last_donor": "cesar_last_donor",
    "joint_site_treat": "joint_splice_site_treatment",
    "bigwig2wig_binary": "bigwig2wig_binary",
    "min_splice_prob": "min_splice_prob",
    "splice_prob_margin": "splice_prob_margin",
    "intron_gain_check": "intron_gain_check",
    "min_intron_gain_score": "intron_gain_threshold",
    # 'min_intron_prob_gapped': 'min_intron_prob_gapped',
    # 'min_intron_prob_ungapped': 'min_intron_prob_ungapped',
    "min_intron_prob_trusted": "min_intron_prob_trusted",
    "min_intron_prob_supported": "min_intron_prob_supported",
    "min_intron_prob_unsupported": "min_intron_prob_unsupported",
    "cesar_binary": "cesar_binary",
    "cesar_memory_bins": "memory_bins",
    "job_nums_per_bin": "job_nums_per_bin",
    "allow_heavy_jobs": "allow_heavy_jobs",
    "max_intron_number": "max_intron_number",
    "matrix_file": "matrix",
    "mask_terminal_mutations": "mask_n_terminal_mutations",
    "leave_missing_stop": "disable_missing_stop_search",
    "consider_alt_frame": "account_for_alternative_frame",
    "spliceai_correction_mode": "spliceai_correction_mode",
    "accepted_loss_symbols": "accepted_loss_symbols",
    "skip_tree_resolver": "skip_gene_trees",
    "use_raxml": "use_raxml",
    "max_clique_size": "max_clique_size",
    "orth_job_num": "orthology_jobs",
    "prank_bin": "prank_binary",
    "tree_bin": "tree_binary",
    "tree_cpus": "tree_cpus",
    "utr_abs_threshold": "utr_abs_threshold",
    "utr_rel_threshold": "utr_rel_threshold",
    "no_utr_extrapolation": "no_utr_boundary_extrapolation",
    "no_adjacent_utrs": "no_adjacent_utr_extra",
    "fixed_adjacent_utrs": "fixed_adjacent_utr_extra",
    "ref_link_file": "link_file",
    "parallel_strategy": "parallel_strategy",
    "nextflow_exec_script": "nextflow_exec_script",
    "max_number_of_retries": "max_number_of_retries",
    "nextflow_config_dir": "nextflow_config_dir",
    "max_parallel_time": "max_parallel_time",
    "keep_nextflow_log": "keep_nextflow_log",
    "cluster_queue_name": "cluster_queue_name",
    "legacy_chain_feature_extraction": "legacy_chain_feature_extraction",
    "toga1": "toga1_compatible",
    "toga1_plus_cesar": "toga1_plus_corrected_cesar",
    "output": "output",
    "project_name": "project_name",
    "project_arg_format": "project_arg_format",
    "keep_tmp": "keep_temporary_files",
    "v": "verbose",
    "debug": "debug",
    "email": "email",
    "mailx_binary": "mailx_binary",
    "fatotwobit_binary": "fatotwobit_binary",
    "twobittofa_binary": "twobittofa_binary",
    "ixixx_binary": "ixixx_binary",
    "ucsc_prefix": "ucsc_prefix",
    "bedtobigbed_binary": "bedtobigbed_binary",
    "ignore_crashed_parallel_batches": "ignore_crashed_parallel_batches",
    "container_image": "container_image",
    "container_executor": "container_executor",
    "bindings": "bindings",
}

TOGA2_ARG2SLOT: Dict[str, str] = {v: k for k, v in TOGA2_SLOT2ARG.items()}
