"""
TOGA2 main class
"""

from collections import defaultdict
from contextlib import nullcontext
from .constants import Constants, TOGA2_SLOTS, TOGA2_SLOT2ARG
from .shared import (
    CommandLineManager, dir_name_by_date, get_upper_dir, hex_dir_name
)
from pathlib import Path
from shutil import copy, which
from typing import Any, Dict, List, Optional, Union
from utils.parallel_jobs_manager import (
    CustomStrategy, NextflowStrategy, ParaStrategy, ParallelJobsManager
)

import click
import logging
import subprocess
import time
import os
import sys

__author__ = 'Yury V. Malovichko'
__version__ = '2.0.2'
__year__ = '2024'
__credits__ = ('Bogdan M. Kirilenko', 'Michael Hiller')

LOCATION: str = get_upper_dir(__file__, 4)
PYTHON_DIR: str = get_upper_dir(__file__, 2)
sys.path.append(PYTHON_DIR)

class TogaMain(CommandLineManager):
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
    Mandatory arguments are:
    * REF_2BIT - a path to reference genome assembly in .2bit format;
    * QUERY_2BIT - a path to query genome assembly, also in .2bit format;
    * ALIGNMENT_CHAINS - a path to genome alignment chains, with REF_2BIT as reference and QUERY_2BIT as query. Can be compressed in .gzip format;
    * REF_ANNOTATION_BED - protein-coding gene annotation for the reference genome, in Bed12 format
    """

    __slots__ = TOGA2_SLOTS

    def __init__(
        self,
        ref_2bit: click.Path,
        query_2bit: click.Path,
        chain_file: click.Path,
        ref_annotation: click.Path,
        resume_from: Optional[str],
        halt_at: Optional[str],
        selected_feature_batches: Optional[Union[str, None]],
        selected_preprocessing_batches: Optional[Union[str, None]],
        selected_alignment_batches: Optional[Union[str, None]],
        no_utr_annotation: Optional[bool],
        isoform_file: Optional[click.Path],
        u12_file: Optional[click.Path],
        spliceai_dir: Optional[click.Path],
        min_chain_score: Optional[int],
        min_orthologous_chain_score: Optional[int],
        feature_jobs: Optional[int],
        orthology_threshold: Optional[float],
        single_exon_model: Optional[click.Path],
        multi_exon_model: Optional[click.Path],
        use_long_distance_model: Optional[bool],
        long_distance_model: Optional[click.Path],
        disable_fragment_assembly: Optional[bool],
        orthologs_only: Optional[bool],
        one2ones_only: Optional[bool],
        enable_spanning_chains: Optional[bool],
        annotate_processed_pseudogenes: Optional[bool],
        preprocessing_jobs: Optional[int],
        max_chains_per_transcript: Optional[int],
        max_search_space_size: Optional[int],
        extrapolation_modifier: Optional[int],
        minimal_covered_fraction: Optional[float],
        exon_locus_flank: Optional[int],
        assembly_gap_size: Optional[int],
        cesar_binary: Optional[Union[str, None]],
        memory_bins: Optional[str],
        job_nums_per_bin: Optional[str],
        allow_heavy_jobs: Optional[bool],
        cesar_memory_limit: Optional[int],
        cesar_canon_u2_acceptor: Optional[click.Path],
        cesar_canon_u2_donor: Optional[click.Path],
        cesar_non_canon_u2_acceptor: Optional[click.Path],
        cesar_non_canon_u2_donor: Optional[click.Path],
        cesar_canon_u12_acceptor: Optional[click.Path],
        cesar_canon_u12_donor: Optional[click.Path],
        cesar_non_canon_u12_acceptor: Optional[click.Path],
        cesar_non_canon_u12_donor: Optional[click.Path],
        cesar_first_acceptor: Optional[click.Path],
        cesar_last_donor: Optional[click.Path],
        separate_splice_site_treatment: Optional[bool],
        bigwig2wig_binary: Optional[click.Path],
        spliceai_correction_mode: Optional[int],
        min_splice_prob: Optional[float],
        splice_prob_margin: Optional[float],
        intron_gain_check: Optional[bool],
        intron_gain_threshold: Optional[float],
        min_intron_prob_trusted: Optional[bool],
        min_intron_prob_supported: Optional[bool],
        min_intron_prob_unsupported: Optional[bool],
        max_intron_number: Optional[int],
        matrix: Optional[click.Path],
        mask_n_terminal_mutations: Optional[bool],
        disable_missing_stop_search: Optional[bool],
        accepted_loss_symbols: Optional[str],
        skip_gene_trees: Optional[bool],
        max_clique_size: Optional[int],
        use_raxml: Optional[bool],
        orthology_jobs: Optional[int],
        prank_binary: Optional[Union[click.Path, None]],
        tree_binary: Optional[Union[click.Path, None]],
        tree_cpus: Optional[int],
        utr_abs_threshold: Optional[int],
        utr_rel_threshold: Optional[float],
        no_utr_boundary_extrapolation: Optional[bool],
        no_adjacent_utr_extra: Optional[bool],
        fixed_adjacent_utr_extra: Optional[bool],
        link_file: Optional[Union[click.Path, None]],
        ucsc_prefix: Optional[str],
        bedtobigbed_binary: Optional[Union[click.Path, None]],
        parallel_strategy: Optional[str],
        nextflow_exec_script: Optional[Union[click.Path, None]],
        max_number_of_retries: Optional[int],
        nextflow_config_dir: Optional[Union[click.Path, None]],
        max_parallel_time: Optional[int],
        cluster_queue_name: Optional[str],
        keep_nextflow_log: Optional[bool],
        ignore_crashed_parallel_batches: Optional[bool],
        legacy_chain_feature_extraction: Optional[bool],
        toga1_compatible: Optional[bool],
        toga1_plus_corrected_cesar: Optional[bool],
        account_for_alternative_frame: Optional[bool],
        output: Optional[click.Path],
        keep_temporary_files: Optional[bool],
        verbose: Optional[bool],
        email: Optional[Union[str, None]],
        mailx_binary: Optional[Union[str, None]],
        fatotwobit_binary: Optional[Union[click.Path, None]],
        twobittofa_binary: Optional[Union[click.Path, None]]
        # version: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.project_name: str = hex_dir_name('TOGA2')

        ## internalized attributes
        self.ref_2bit: click.Path = os.path.abspath(ref_2bit)
        self.query_2bit: click.Path = os.path.abspath(query_2bit)
        self.chain_file: click.Path = chain_file
        self.ref_annotation: click.Path = ref_annotation

        self.resume_from: str = resume_from
        self.halt_at: str = halt_at
        self.selected_feature_batches: List[str] = (
            selected_feature_batches.split(',') if selected_feature_batches is not None else None
        )
        self.selected_preprocessing_batches: List[str] = (
            selected_preprocessing_batches.split(',') if selected_preprocessing_batches is not None else None
        )
        self.selected_alignment_batches: List[str] = (
            selected_alignment_batches.split(',') if selected_alignment_batches is not None else None
        )
        self.skip_utr: bool = no_utr_annotation

        self.isoform_file: Union[click.Path, None] = isoform_file
        self.u12_file: Union[click.Path, None] = u12_file
        self.spliceai_dir: Union[click.Path, None] = spliceai_dir

        self.min_chain_score: int = min_chain_score
        self.min_orth_chain_score: int = min_orthologous_chain_score

        self.feature_job_num: int = feature_jobs
        self.orthology_threshold: float = orthology_threshold
        self.se_model: click.Path = single_exon_model
        self.me_model: click.Path = multi_exon_model
        self.use_ld_model: bool = long_distance_model
        self.ld_model: click.Path = long_distance_model

        self.disable_fragment_assembly: bool = disable_fragment_assembly
        self.annotate_ppgenes: bool = annotate_processed_pseudogenes
        self.orthologs_only: bool = orthologs_only
        self.one2ones_only: bool = one2ones_only
        self.enable_spanning_chains: bool = enable_spanning_chains

        self.preprocessing_job_num: int = preprocessing_jobs
        self.max_chains_per_transcript: int = max_chains_per_transcript
        self.cesar_memory_limit: float = cesar_memory_limit
        self.max_search_space_size: int = max_search_space_size
        self.extrapolation_modifier: float = extrapolation_modifier
        self.minimal_covered_fraction: float = (
            minimal_covered_fraction if not enable_spanning_chains else 0.0
        )
        self.exon_locus_flank: int = exon_locus_flank
        self.assembly_gap_size: int = assembly_gap_size
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

        self.bigwig2wig_binary: click.Path = bigwig2wig_binary
        self.min_splice_prob: float = min_splice_prob
        self.splice_prob_margin: float = splice_prob_margin
        self.intron_gain_check: bool = intron_gain_check
        self.min_intron_gain_score: float = intron_gain_threshold
        # self.min_intron_prob_gapped: float = min_intron_prob_gapped
        # self.min_intron_prob_ungapped: float = min_intron_prob_ungapped
        self.min_intron_prob_trusted: bool = min_intron_prob_trusted
        self.min_intron_prob_supported: bool = min_intron_prob_supported
        self.min_intron_prob_unsupported: bool = min_intron_prob_unsupported

        self.cesar_binary: Union[str, None] = cesar_binary
        self.cesar_memory_bins: str = memory_bins
        self.job_nums_per_bin: str = job_nums_per_bin
        self.allow_heavy_jobs: bool = allow_heavy_jobs
        self.matrix_file: str = matrix
        self.mask_terminal_mutations: bool = mask_n_terminal_mutations
        self.leave_missing_stop: bool = disable_missing_stop_search
        self.consider_alt_frame: bool = account_for_alternative_frame
        self.spliceai_correction_mode: int = spliceai_correction_mode

        self.accepted_loss_symbols: str = accepted_loss_symbols
        self.skip_tree_resolver: bool = skip_gene_trees
        self.max_clique_size: int = max_clique_size
        self.use_raxml: bool = use_raxml
        self.orth_job_num: int = orthology_jobs
        self.prank_bin: Union[str, None] = prank_binary
        self.tree_bin: Union[str, None] = tree_binary
        self.tree_cpus: int = tree_cpus

        self.utr_abs_threshold: int = utr_abs_threshold
        self.utr_rel_threshold: float = utr_rel_threshold
        self.no_utr_extrapolation: bool = no_utr_boundary_extrapolation
        self.no_adjacent_utrs: bool = no_adjacent_utr_extra
        self.fixed_adjacent_utrs: bool = fixed_adjacent_utr_extra

        self.ref_link_file: Union[str, None] = link_file
        self.ucsc_prefix: str = ucsc_prefix

        self.email: Union[str, None] = email
        self.mailx_binary: Union[click.Path, None] = mailx_binary

        self.bedtobigbed_binary: Union[click.Path, None] = bedtobigbed_binary
        self.fatotwobit_binary: Union[click.Path, None] = fatotwobit_binary
        self.twobittofa_binary: Union[click.Path, None] = twobittofa_binary

        self.parallel_strategy: str = parallel_strategy
        self.max_number_of_retries: int = max_number_of_retries
        self.nextflow_exec_script: str = nextflow_exec_script
        self.local_executor: str = parallel_strategy == 'local'#self.nextflow_config_dir is not None ## TODO: Add config dir content check
        self.max_parallel_time: int = max_parallel_time
        self.keep_nextflow_log: bool = keep_nextflow_log
        self.cluster_queue_name: str = cluster_queue_name
        self.ignore_crashed_parallel_batches: bool = ignore_crashed_parallel_batches
        self.legacy_chain_feature_extraction: bool = legacy_chain_feature_extraction
        self.parallel_process_names: List[str] = []

        self.output: str = Path(
            output if output else dir_name_by_date('toga2_run')
        ).absolute()
        self.keep_tmp: bool = keep_temporary_files

        ## benchmarking flags
        self.toga1: bool = toga1_compatible
        self.toga1_plus_cesar: bool = toga1_plus_corrected_cesar

        ## input genome file status flags; OBSOLETE??
        # self.ref_is_2bit: bool = False
        # self.ref_is_gzipped: bool = False
        # self.query_is_2bit: bool = False
        # self.query_is_gzipped: bool = False

        ## failed batch collections
        self.failed_feature_batches: List[str] = []
        self.failed_preprocessing_batches: List[str] = []
        self.failed_alignment_batches: List[str] = []
        self.failed_orthology_batches: List[str] = []

        ## preset Nextflow log files
        self.nextflow_config_files: Dict[str, Union[str, None]] = {}

        ## directory structure
        ## first-level
        self.tmp: str = os.path.join(self.output, 'tmp')
        self.meta: str = os.path.join(self.output, 'meta')
        # self.res: str = os.path.join(self.output, 'results')
        self.logs: str = os.path.join(self.output, 'logs')
        self.nextflow_dir: str = os.path.join(self.output, 'nextflow')
        self.nextflow_config_dir: str = (
            self.nextflow_dir if nextflow_config_dir is None else nextflow_config_dir
        )
        self.ucsc_dir: str = os.path.join(
            self.output, 'ucsc_browser_files'
        )
        self.arg_file: str = os.path.join(self.logs, f'project_args_{self.project_name}.tsv')
        self.log_file: str = os.path.join(self.logs, f'{self.project_name}.log')
        self.failed_batches_file: str = os.path.join(self.logs, f'failed_batches_{self.project_name}.tsv')

        ## temporary subdirectories
        self.input_data: str = os.path.join(self.tmp, 'input_data')
        self.feature_job_dir: str = os.path.join(
            self.tmp, 'feature_extraction_jobs'
        )
        self.feature_data_dir: str = os.path.join(
            self.tmp, 'feature_extraction_input'
        )
        self.feature_res_dir: str = os.path.join(
            self.tmp, 'feature_extraction_results'
        )
        self.classification_dir: str = os.path.join(
            self.tmp, 'projection_classification'
        )
        self.preprocessing_job_dir: str = os.path.join(
            self.tmp, 'cesar_preprocessing_jobs'
        )
        self.preprocessing_job_dir: str = os.path.join(
            self.tmp, 'cesar_preprocessing_jobs'
        )
        self.preprocessing_res_dir: str = os.path.join(
            self.tmp, 'cesar_preprocessing_res'
        )
        self.alignment_job_dir: str = os.path.join(
            self.tmp, 'cesar_alignment_jobs'
        )
        self.alignment_res_dir: str = os.path.join(
            self.tmp, 'cesar_alignment_res'
        )
        self.orthology_resolution_dir: str = os.path.join(
            self.tmp, 'orthology_resolution'
        )
        self.orthology_job_dir: str = os.path.join(
            self.tmp, 'orthology_resolution_jobs'
        )
        self.orthology_input_dir: str = os.path.join(
            self.tmp, 'orthology_resolution_input'
        )
        self.orthology_res_dir: str = os.path.join(
            self.tmp, 'orthology_resolution_res'
        )
        self.vis_input_dir: str = os.path.join(
            self.tmp, 'bed_files_for_browsers'
        )
        self.finalized_output_dir: str = os.path.join(
            self.tmp, 'finalized_output_files'
        )

        ## meta subdirectories
        self.rejection_dir: str = os.path.join(
            self.tmp, 'rejection_logs'
        )

        ## input file copies and indices
        self.bed_file_copy: str = os.path.join(
            self.input_data, 'reference_annotation.bed'
        )
        self.ref_cds_unfilt: str = os.path.join(
            self.input_data, 'reference_annotation_raw_cds.bed'
        )
        self.prefiltered_transcripts: str = os.path.join(
            self.rejection_dir, 'prefiltered_transcripts.tsv'
        )
        self.chain_file_copy: str = os.path.join(
            self.input_data, 'genome_alignment.chain'
        )
        self.chain_index: str = os.path.join(
            self.input_data, 'genome_alignment.bst'
        )
        self.chain_index_txt: str = os.path.join(
            self.input_data, 'genome_alignment.chain_ID_position'
        )
        self.full_bed_hdf5: str = os.path.join(
            self.input_data, 'reference_annotation.hdf5'
        )
        self.cds_bed_file: str = os.path.join(
            self.input_data, 'reference_annotation_cds.bed'
        )
        self.cds_bed_hdf5: str = os.path.join(
            self.input_data, 'reference_annotation_cds.hdf5'
        )
        self.u12_hdf5: str = os.path.join(
            self.input_data, 'u12_data.hdf5'
        )
        self.ref_contig_size_file: str = os.path.join(
            self.input_data, 'reference_contig_sizes.tsv'
        )
        self.query_contig_size_file: str = os.path.join(
            self.input_data, 'query_contig_sizes.tsv'
        )
        self.aa_hdf5: str = os.path.join(
            self.input_data, 'protein_sequences.hdf5'
        )
        self.annot_dir: str = os.path.join(
            self.tmp, 'annotation_raw'
        )
        self.orth_resolution_raw: str = os.path.join(
            self.annot_dir, 'orthology_classification.tsv'
        )
        self.query_annotation_filt: str = os.path.join(
            self.annot_dir, 'query_annotation.bed'
        )
        self.query_genes_raw: str = os.path.join(
            self.annot_dir, 'query_genes.tsv'
        )
        self.query_genes_bed_raw: str = os.path.join(
            self.annot_dir, 'query_genes.bed'
        )

        ## joblists
        self.feature_extraction_joblist: str = os.path.join(
            self.feature_job_dir, 'joblist'
        )
        self.cesar_preprocess_joblist: str = os.path.join(
            self.preprocessing_job_dir, 'joblist_preprocess'
        )
        self.cesar_job_list_summary: str = os.path.join(
            self.alignment_job_dir, 'joblist_description.tsv'
        )
        self.orth_resolution_joblist: str = os.path.join(
            self.orthology_job_dir, 'joblist'
        )

        ## provisional and meta files
        self.feature_table: str = os.path.join(
            self.meta, 'projection_features.tsv'
        )
        self.tr2chain_classes: str = os.path.join(
            self.meta, 'trans_to_chain_classes.tsv'
        )
        self.feature_rejection_log: str = os.path.join(
            self.rejection_dir, 'rejected_at_feature_extraction.tsv'
        )
        self.class_rejection_log: str = os.path.join(
            self.rejection_dir, 'rejected_at_classification.tsv'
        )
        self.fragmented_projection_list: str = os.path.join(
            self.meta, 'fragmented_projections.tsv'
        )
        self.preprocessing_report: str = os.path.join(
            self.meta, 'memory_requirements.tsv'
        )
        # self.spanning_chain_precedence_file: str = os.path.join(
        #     self.meta, 'spanning_chains_ref_spans.tsv'
        # )
        self.rejected_by_graph: str = os.path.join(
            self.orthology_resolution_dir, 'rejected_by_graph_reduction.tsv'
        )
        self.weak_ortholog_names: str = os.path.join(
            self.orthology_resolution_dir, 'rejected_projection_names.txt'
        )
        self.discarded_overextended_projections: str = os.path.join(
            self.meta, 'discarded_overextended_projections.txt'
        )
        self.missing_transcripts: str = os.path.join(
            self.orthology_resolution_dir, 'unclassified_transcripts.txt'
        )
        self.resolved_leaves_file: str = os.path.join(
            self.meta, 'resolved_tree_pairs.tsv'
        )
        self.unresolved_clades_file: str = os.path.join(
            self.meta, 'unresolved_tree_clades.txt'
        )
        self.temporary_orth_report: str = os.path.join(
            self.orthology_resolution_dir, 'orthology_classification.tsv'
        )
        self.one2zero_genes: str = os.path.join(
            self.orthology_resolution_dir, 'one2zero_genes.txt'
        )
        self.preprocessing_rejection_log: str = os.path.join(
            self.rejection_dir, 'rejected_at_preprocessing.tsv'
        )
        self.spanning_chain_coords: str = os.path.join(
            self.meta, 'spanning_chain_ref_coords.tsv'
        )
        self.paralog_report: str = os.path.join(
            self.meta, 'paralogous_projections_to_align.tsv'
        )
        self.processed_pseudogene_report: str = os.path.join(
            self.meta, 'processed_pseudogene_projections_to_align.tsv'
        )
        self.redundant_paralogs: str = os.path.join(
            self.meta, 'redundant_paralogs.txt'
        )
        self.redundant_ppgenes: str = os.path.join(
            self.meta, 'redundant_processed_pseudogenes.txt'
        )
        self.discarded_proj_bed: str = os.path.join(
            self.meta, 'discarded_projections.bed'
        )
        self.redundancy_rejection_log: str = os.path.join(
            self.rejection_dir, 'rejected_by_redundancy.tsv'
        )
        self.alignment_rejection_log: str = os.path.join(
            self.rejection_dir, 'rejected_at_alignment.tsv'
        )
        self.gene_inference_rejection_log: str = os.path.join(
            self.rejection_dir, 'rejected_at_gene_inference.tsv'
        )
        self.aggr_ucsc_stub: str = os.path.join(
            self.vis_input_dir, 'ucsc_bed_stub.bed23'
        )
        self.all_deprecated_projs: str = os.path.join(
            self.vis_input_dir, 'all_deprecated_projections.txt'
        )

        ## output files and directories
        self.transcript_meta: str = os.path.join(
            self.meta, 'transcript_meta.tsv'
        )
        self.splice_sites: str = os.path.join(
            self.meta, 'splice_sites.tsv'
        )
        self.gained_intron_summary: str = os.path.join(
            self.meta, 'gained_intron_summary.tsv'
        )
        self.splice_site_shifts: str = os.path.join(
            self.meta, 'splice_site_shifts.tsv'
        )
        self.selenocysteine_codons: str = os.path.join(
            self.meta, 'selenocysteine_codons.tsv'
        )
        self.query_annotation_raw: str = os.path.join(
            self.meta, 'query_annotation.with_discarded_exons.bed'
        )
        self.query_exon_meta: str = os.path.join(
            self.meta, 'exon_meta.tsv'
        )
        self.tree_summary_table: str = os.path.join(
            self.meta, 'gene_tree_job_summary.tsv'
        )
        # self.orthology_results_dir: str = os.path.join(
        #     self.tmp, 'orthology_resolution'
        # )

        ## ultimate output files
        self.pred_scores: str = os.path.join(
            self.output, 'orthology_scores.tsv'
        )
        self.query_annotation_final: str = os.path.join(
            self.output, 'query_annotation.bed'
        )
        self.query_annotation_with_utrs: str = os.path.join(
            self.output, 'query_annotation.with_utrs.bed'
        )
        self.aa_fasta: str = os.path.join(
            self.output, 'protein_aln.fa'
        )
        self.cds_fasta: str = os.path.join(
            self.output, 'nucleotide.fa'
        )
        self.codon_fasta: str = os.path.join(
            self.output, 'codon_aln.fa'
        )
        self.exon_fasta: str = os.path.join(
            self.output, 'exon_aln.fa'
        )
        self.exon_2bit: str = os.path.join(
            self.output, 'exon_seqs.2bit'
        )
        self.mutation_report: str = os.path.join(
            self.output, 'inactivating_mutations.tsv'
        )
        self.processed_pseudogene_annotation: str = os.path.join(
            self.output, 'processed_pseudogenes.bed'
        )
        self.gene_loss_summary: str = os.path.join(
            self.output, 'loss_summary.tsv'
        )
        self.final_rejection_log: str = os.path.join(
            self.output, 'rejected_items.tsv'
        )
        self.pseudogene_annotation: str = os.path.join(
            self.output, 'processed_pseudogenes.bed'
        )
        self.query_genes: str = os.path.join(
            self.output, 'query_genes.tsv'
        )
        self.query_genes_bed: str = os.path.join(
            self.output, 'query_genes.bed'
        )
        self.orth_resolution_report: str = os.path.join(
            self.output, 'orthology_classification.tsv'
        )
        self.cds_gzip: str = self.cds_fasta + '.gz'
        self.codon_gzip: str = self.codon_fasta + '.gz'
        self.exon_gzip: str = self.exon_fasta + '.gz'
        self.prot_gzip: str = self.aa_fasta + '.gz'
        self.splice_sites_gzip: str = self.splice_sites + 'gz'
        self.exon_meta_gzip: str = self.query_exon_meta + '.gz'
        self.transcript_meta_gzip: str = self.transcript_meta + '.gz'

        ## script location attributes
        ## NOTE: Most of these scripts are now obsolete since the modules are imported directly to the mothership

        self.CHAIN_FILTER_SCRIPT: str = os.path.join(
            LOCATION, 'src', 'rust', 'target', 'release', 'chain_filter'
        )
        self.INDEX_CHAIN_SCRIPT: str = os.path.join(
            PYTHON_DIR, 'modules', 'chain_bst_index.py'
        )
        self.REF_BED_FILTER: str = os.path.join(
            PYTHON_DIR, 'modules', 'filter_ref_bed.py'
        )
        self.CDS_TRACK_SCRIPT: str = os.path.join(
            LOCATION, 'src', 'rust', 'target', 'release', 'bed12ToFraction'
        ) ## TODO: Replace with the cubiculum code
        self.REF_BED_TO_HDF5: str = os.path.join(
            PYTHON_DIR, 'modules', 'bed_hdf5_index.py'
        )
        self.U12_TO_HDF5_SCRIPT: str = os.path.join(
            PYTHON_DIR, 'modules', 'intronIC_to_hdf5.py'
        )
        self.CONTIG_SIZE_SCRIPT: str = os.path.join(
            PYTHON_DIR, 'modules', 'get_contig_sizes.py'
        ) ## TODO: Alejandro must have a Rust implementation
        self.FEATURE_EXTRACTOR: str = os.path.join(
            LOCATION, 'src', 'rust', 'target', 'release', 'feature_extraction'
        )
        self.MODEL_TRAINER: str = os.path.join(
            PYTHON_DIR, 'train_model.py'
        )
        self.FINAL_RESOLVER_SCRIPT: str = os.path.join(
            LOCATION, 'modules', 'final_orthology_aggregator.py'
        )
        self.UTR_PROJECTOR_SCRIPT: str = os.path.join(
            LOCATION, 'src', 'rust', 'target', 'release', 'utr_projector'  
        )
        self.SCHEMA_FILE: str = os.path.join(
            LOCATION, 'supply', 'bb_schema_toga2.as'
        )

        self.run()

    def run(self) -> None:
        """
        Executes the TOGA 2.0 routine
        """
        self._echo('Initializing TOGA 2')
        ## Step 0: Output directory creation, input preparation, and sanity checks
        ## create the output directory
        self._echo('Creating working directory')
        self.create_wd()
        self._to_log(
            'Working directory created', 'info'
        )
        # self._echo('Working directory successfully created')
        # self._echo('Writing project metadata')
        self._to_log('Writing project metadata', 'info')
        self.write_project_meta()
        self._to_log(
            f'Project metadata dumped at {self.arg_file}', 'info'
        )
        # self._echo(f'Project metadata dumped at {self.arg_file}')

        ## check the input arguments
        # self._echo('Checking input arguments')
        self._to_log('Checking input arguments', 'info')
        self.check_arguments()
        # self._echo('All settings are correct')
        self._to_log('All settings are correct', 'info')

        ## check the availability of all the necessary binaries
        # self._echo('Checking third-part dependencies')
        self._to_log('Checking third-part dependencies', 'info')
        self.check_binaries()
        # self._echo('All dependencies were found')
        self._to_log('All dependencies were found', 'info')

        ## if SpliceAI directory was provided, check its contents
        self.check_spliceai_files()

        ## if a custom Nextflow executor script was not provided,
        ## generate one from the boilerplate, substituting the maxRetrie parameter
        if self.nextflow_exec_script is None:
            self._to_log('Generating a Nextflow execution master script')
            self._generate_nf_script()
            self._to_log('Nextflow execution script saved at %s' % self.nextflow_exec_script)

        ## if a path to custom Nextflow configuration files was provided,
        ## check the directory completeness
        self._check_nextflow_configs()

        ## check the format consistency of all the input files
        # self.check_ref_bed_file()
        # self.check_u12_file()
        # self.check_isoform_file_file()

        ## create the necessary links, indices, and HDF5 storage files
        # self._echo('Preparing input data')
        if self._execute_step('setup'):
            self._to_log('Preparing input data', 'info')
            self.prepare_data()
            self._to_log('All input data successfully checked and handled', 'info')

        ## Step 1: Extract chain features
        # self._echo('Running TOGA main procedure')
        self._to_log('Running TOGA main procedure', 'info')
        if self._execute_step('feature_extraction'):
            if not self.legacy_chain_feature_extraction:
                self._to_log('Exracting projection features for further classification')
                self.extract_chain_features()
                pass
            else:
                ## 1a: Schedule chain jobs
                if self.selected_feature_batches is None:
                    self._to_log('Scheduling projection feature extraction jobs')
                    self.schedule_feature_jobs()
                else:
                    self._to_log('Preparing selected feature extraction batches for parallel run')
                    self._schedule_selected_batches('feature_extraction')
                self._to_log('Feature extraction jobs successfully scheduled')
                ## 1b: Run the feature extraction jobs
                self._to_log('Extracting projection features (parallel step)')
                self.run_feature_extraction_jobs()
                self._to_log('All projection feature jobs finished')
                ## 1c: Aggregate chain feature data and create a classification dataset
                # self._echo('Aggregating projection features')
                self._to_log('Aggregating projection features')
                self.create_feature_dataset()
                self._to_log('Feature dataset successfully prepared')
        else:
            if self.halt_at == 'feature_extraction':
                self.notify_on_completion('feature_extraction')
                self._exit('Finish pipeline before feature extraction step as suggested')
            self._to_log('Skipping feature extraction step as suggested')

        ## Step 2: Classify chain-transcript pairs in terms of orthology
        # self._echo('Classifying projections')
        if self._execute_step('classification'):
            self._to_log('Classifying projections', 'info')
            self.classify_chains()
            # self._echo('Projections successfully classified')
            self._to_log('Projection successfully classified', 'info')
        else:
            if self.halt_at == 'classification':
                self.notify_on_completion('classification')
                self._exit('Finish pipeline before projection classification step as suggested')
            self._to_log('Skipping projection classification step as suggested')

        if self._execute_step('preprocessing'):
            ## Step 2.5: Unless explicitly disabled, stitch fragmented projections
            if not self.disable_fragment_assembly:
                # self._echo('Recovering fragmented projections')
                self._to_log('Recovering fragmented projections', 'info')
                self.recover_fragmented_projections()
                # self._echo('Fragmented projections successfully recovered')
                self._to_log('Fragmented projections successfully recovered', 'info')

            ## Step 3: Prepare and run CESAR preprocessing jobs
            ## prepare preprocessing jobs
            # self._echo('Scheduling preprocessing jobs')
            if self.selected_preprocessing_batches is None:
                self._to_log('Scheduling preprocessing jobs', 'info')
                self.schedule_preprocessing_jobs()
            else:
                self._to_log('Preparing selected preprocessing batches for parallel run')
                self._schedule_selected_batches('preprocessing')
            # self._echo('Preprocessing jobs scheduled')
            self._to_log('Preprocessing jobs scheduled', 'info')

            ## run preprocessing jobs
            self._to_log('Running preprocessing jobs (parallel step)', 'info')
            self.run_preprocessing_jobs()
            self._to_log('Preprocessing jobs finished', 'info')
        else:
            if self.halt_at == 'preprocessing':
                self.notify_on_completion('preprocessing')
                self._exit('Finish pipeline before CESAR preprocessing step as suggested')
            self._to_log('Skipping CESAR preprocessing step as suggested')

        ## Step 4: Aggregate preprocessing results
        if self._execute_step('aggregate_preprocessing_res'):
            self._to_log('Aggregating preprocessing results')
            self.aggregate_preprocessing_res()
        else:
            if self.halt_at == 'aggregate_preprocessing_res':
                self.notify_on_completion('aggregate_preprocessing_res')
                self._exit('Finish aggregatting preprocessing results step as suggested')
            self._to_log('Skipping aggregatting preprocessing results step as suggested')

        ## Step 4: Prepare and run CESAR alignment & postprocessing jobs
        if self._execute_step('alignment'):
            ## prepare alignment jobs
            # self._echo('Scheduling CESAR alignment jobs')
            if self.selected_alignment_batches is None:
                self._to_log('Scheduling CESAR alignment jobs', 'info')
                self.schedule_alignment_jobs()
            else:
                self._to_log('Preparing selected alignment batches for parallel run')
                self._schedule_selected_batches('alignment')
            # self._echo('Alignment jobs scheduler')
            self._to_log('Alignment jobs scheduled', 'info')

            ## then run them
            self._to_log('Running alignment jobs (parallel step)', 'info')
            self.run_alignment_jobs()
            # self._echo('Alignment step complete')
            self._to_log('Alignment step complete', 'info')
        else:
            if self.halt_at == 'alignment':
                self.notify_on_completion('alignment')
                self._exit('Finishing pipeline before CESAR alignment step as suggested')
            self._to_log('Skipping CESAR alignment step as suggested')

        ## Step 5: Aggregate alignment/postprocessing results
        if self._execute_step('aggregate_cesar_res'):
            self._to_log('Aggregating alignment results', 'info')
            self.aggregate_cesar_results()
            self._to_log('CESAR step results successfully aggregated', 'info')
        else:
            if self.halt_at == 'aggregate_cesar_res':
                self.notify_on_completion('aggregate_cesar_res')
                self._exit('Finishing pipeline before alignment results aggregation step as suggested')
            self._to_log('Skipping alignment results aggregation step as suggested')

        ## Step 6: Summarise ortholog presence status on projection, transcript,
        ## and (if requested) gene levels

        ## 6a: Aggregate rejection reports from all the previous steps
        # self._echo('Aggregating rejection reports')
        if self._execute_step('loss_summary'):
            self._to_log('Aggregating rejection reports', 'info')
            self.aggregate_rejection_reports()
            # self._echo(f'All rejected items are placed into {self.final_rejection_log}')
            self._to_log(
                'All rejected items are placed into %s' % self.final_rejection_log,
                'info'
            )

            ## 6b: Summarise loss data
            # self._echo('Summarizing transcript/gene conservation data')
            self._to_log('Summarizing transcript/gene conservation data', 'info')
            self.loss_status_summary()
            # self._echo('Conservation data successfully summarized')
            self._to_log('Conservation data successfully summarized', 'info')

            ## 6c: If processed pseudogenes were not considered in the previous steps,
            ## prepare a simple BED9 track
            if not self.annotate_ppgenes:
                self._to_log('Preparing a processed pseudogene BED track')
                self.prepare_pseudogene_track()
                self._to_log('Pseudogene track successfully produced')
        else:
            if self.halt_at == 'loss_summary':
                self.notify_on_completion('loss_summary')
                self._exit('Finishing pipeline before gene loss summary step as suggested')
            self._to_log('Skipping gene loss summary step as suggested')


        ## Step 7: Resolve orthology status for all reference genes
        ## collapse projections to genes
        # self._echo('Inferring query genes')
        if self._execute_step('orthology'):
            self._to_log('Inferring query genes', 'info')
            self.infer_query_genes()
            # self._echo('Query genes successfully inferred')
            self._to_log('Query genes successfully inferred', 'info')

            ## resolve orthology relationships with a graph-based method
            # self._echo('Resolving orthology relationships')
            # self._echo('Running initial step (graph-based resolution)')
            if not self.skip_tree_resolver:
                self._to_log('Converting protein annotation FASTA file into HDF5 format')
                self.convert_fasta_to_hdf5()
                self._to_log('HDF5 preparation is complete')
            self._to_log('Resolving orthology relationships', 'info')
            self._to_log('Running initial step (graph-based resolution)', 'info')
            self.resolve_orthology_initial()
            # self._echo('Graph-based step successfully finished')
            self._to_log('Graph-based step successfully finished', 'info')
            ## aggregate the results

            ## if requested, run fine orthology resolution jobs
            # self.copy_resolution_reports()
            if not self.skip_tree_resolver:
                self._to_log(
                    'Running fine orthology resolution jobs (parallel step)',
                    'info'
                )
                self.run_tree_jobs()
                self._to_log('Tree jobs finished', 'info')
                self._to_log('Updating orthology estimates with tree-inferred data', 'info')
            self._to_log('Orthology inference step complete', 'info')
        else:
            if self.halt_at == 'orthology':
                self.notify_on_completion('orthology')
                self._exit('Finishing pipeline before orthology inference step as suggested')
            self._to_log('Skipping orthology inference step as suggested')

        ## Step 8: If tree-based resolution was requested, summarize its results 
        ## and adjust orthology resolution accordingly
        if self._execute_step('summarize_trees') and not self.skip_tree_resolver:
            self._to_log('Adding gene tree step results to orthology resolution')
            self.resolve_orthology_final()
            self._to_log('Creating a summary table for the gene tree resolution step')
            self.gene_tree_summary()
            self._to_log('Tree-based orthology inference step complete')
            pass
        else:
            if self.halt_at == 'summarize_trees':
                self.notify_on_completion('summarize_trees')
                self._exit(
                    'Finishing pipeline before adding gene tree data to orthology inference as suggested'
                )
            else:
                self._to_log(
                    'Skipping adding gene tree data to orthology inference as suggested'
                )

        ## Step 9: Finalize the output
        if self._execute_step('finalize'):
            # self._to_log('Creating exon FASTA HDF5 storage for SLEASY')
            self._to_log('Aggregatig deprecated projection lists')
            self.deprecated_projection_file()
            self._to_log('Renaming & filtering final output files')
            self.filter_final_bed_files()
            if self.skip_utr:
                self._to_log('Skipping the UTR annotation step as suggested')
            else:
                self._to_log('Annotating UTR projections and adding them to the output Bed and BigBed files')
                self.annotate_utrs()
            self._to_log('Creating exon 2bit file for SLEASY')
            self.convert_exon_fasta()
            self._to_log('Compressing FASTA files from the output')
            self.gzip_fasta_files()
            self._to_log('Finalizing naming notation for query genes')
            self.rename_query_genes()
            ## before moving to the final steps, record the failed batches 
            self._to_log('Recording failed batches')
            self.write_failed_batches()
        else:
            if self.halt_at == 'finalize':
                self.notify_on_completion('finalize')
                self._exit('Finishing pipeline befiore output finalizing step as suggested')
            self._to_log('Skipping output finalizing step as suggested')

        ## Step 10: Prepare UCSC browser input files
        if self._execute_step('ucsc_report'):
            self._to_log('Preparing BigBed track for UCSC Genome Browser')
            self.prepare_bigbed_track()
            self._to_log('UCSC Genome Browser track successfully created')
        else:
            if self.halt_at == 'ucsc_report':
                self.notify_on_completion('ucsc_report')
                self._exit('Finishing pipeline before UCSC report preparation as suggested')
            self._to_log('Skipping UCSC report preparation step as suggested')

        ## Step 11: If requested, clean the temporary directory
        self._to_log('TOGA2 pipeline finished, cleaning up the temporary data')
        self.cleanup()

        ## If e-mail notification was requested, notify the user on successful pipeline run
        self.notify_on_completion('all')

        ## So long and thanks for the fish!
        self._to_log('TOGA2 pipeline successfully finished')

    def set_logging(self) -> None:
        """
        Sets up logging system for a TogaMain instance
        """
        self.logger: logging.Logger = logging.getLogger(self.project_name)
        file_handler: logging.FileHandler = logging.FileHandler(
            self.log_file, mode='a', encoding=Constants.UTF8
        )
        file_handler.setFormatter(Constants.FORMATTER)
        self.logger.addHandler(file_handler)
        if self.v:
            console_handler: logging.StreamHandler = logging.StreamHandler()
            console_handler.setFormatter(Constants.FORMATTER)
            self.logger.addHandler(console_handler)
        self.logger.propagate = False

    def _to_log(self, msg: str, level: str = 'info') -> None:
        """
        Adds the message to the TOGA2 log
        """
        getattr(self.logger, level)(msg)

    def _die(self, msg: str) -> None:
        """Error-exit with a given message"""
        if self.email is not None:
            self._email(
                Constants.CRASH_HEADER.format(self.project_name), 
                Constants.CRASH_EMAIL.format(self.project_name, self.output, msg))
        super()._die(msg)

    def _execute_step(self, step: str):
        """
        Defines whether the current step is to be executed based on 'resume' and 'halt' options
        """
        if step not in Constants.RESUME_ORDER:
            self._die('Improper step name provided')
        regular_start: bool = self.resume_from == 'all'
        regular_finish: bool = self.halt_at == 'all'
        step: int = Constants.RESUME_ORDER[step]
        resume_step: int = Constants.RESUME_ORDER[self.resume_from]
        halt_step: int = Constants.RESUME_ORDER[self.halt_at]
        exec_started: bool = resume_step <= step or regular_start
        exec_not_finished: bool = step < halt_step or regular_finish
        return exec_started and exec_not_finished
        
    def write_project_meta(self) -> None:
        """
        """
        with open(self.arg_file, 'w') as h:
            h.write(f'version\t{__version__}\n')
            for arg, option in TOGA2_SLOT2ARG.items(): ## TODO: Double-check the number
                arg_value: str = getattr(self, arg)
                value: str = (
                    arg_value.name if type(arg_value) is click.utils.LazyFile
                    else None if arg_value == ''
                    else arg_value
                )
                h.write('\t'.join(map(str, (option, value))) + '\n')

    def _rsync(self, from_: str, to_: str) -> None:
        """
        Synchonizes local object with a remote instance
        """
        cmd: str = f'rsync -a {from_} {to_}'
        _ = self._exec(
            cmd, 'Symbolic link creation failed with the following error:'
        )

    def _get_manager(self) -> ParallelJobsManager:
        """
        Creates a parallel process manager instance according to the select
        parallel process strategy
        """
        if self.parallel_strategy == 'para':
            strategy = ParaStrategy()
        elif self.parallel_strategy == 'custom':
            strategy = CustomStrategy()
        else:
            strategy = NextflowStrategy()
        jobs_manager: ParallelJobsManager = ParallelJobsManager(strategy)
        return jobs_manager

    def _run_parallel_process(
        self,
        project_name: str,
        project_path: str,
        joblist: str,
        queue: str,
        wait_for_process: bool = True,
        memory_limit: int = 6.5,
        cpu: int = 1,
        process_num: int = 1000,
        step: str = '',
        config_file: Union[str, None] = None
    ) -> None:
        """
        A general method for running parallel steps
        """
        job_manager: ParallelJobsManager = self._get_manager()
        manager_data: Dict[str, str] = {
            'project_name': project_name,
            'project_path': project_path,
            'logs_dir': project_path,
            'nextflow_dir': self.nextflow_dir,
            'NF_EXECUTE': self.nextflow_exec_script,
            'local_executor': self.local_executor,
            'keep_nf_logs': self.keep_nextflow_log,
            # 'nexflow_config_file': nextflow_config,
            'nextflow_config_dir': self.nextflow_config_dir,
            'temp_wd': self.tmp,
            'queue_name': self.cluster_queue_name,
            'logger': self.logger
        }
        if config_file is not None:
            manager_data['nexflow_config_file'] = config_file
        try:
            job_manager.execute_jobs(
                joblist,
                manager_data,
                project_name,
                wait=wait_for_process,
                memory_limit=memory_limit,
                queue_name=queue,
                clean=self.keep_tmp,
                cpu=cpu,
                process_num=process_num,
                executor=self.parallel_strategy
            )
            if not wait_for_process:
                return job_manager
            else:
                iteration: int = 1
                while job_manager.check_status() is None:
                    duration: int = iteration * Constants.ITER_DURATION
                    self._to_log(
                        'Polling %s jobs, iteration %i; already waiting %i seconds' % (
                            step, iteration, duration
                        )
                    )
                if job_manager.check_status() != 0:
                    self._to_log(
                        'Some of the monitored parallel processes at %s step died' % step,
                        'warning'
                    )
        except KeyboardInterrupt:
            # TogaUtil.terminate_parallel_processes([job_manager])
            self._to_log('Aborting the parallel step')
            self._terminate_parallel_processes([job_manager])
            self._exit()

    def _terminate_parallel_processes(self, managers: ParallelJobsManager) -> None:
        """
        Given a list of parallel job managers, 
        kills both the cluster jobs and the governing scheduler
        """
        for manager in managers:
            manager.terminate_process()

    def _monitor_jobs(
        self, 
        jobs_managers: List[ParallelJobsManager],
        step: str
    ):
        """Monitors parallel jobs if many batches run simultaneously."""
        self._to_log('Polling cluster jobs')
        iter_num: int = 0
        while True:  # Run until all jobs are done (or crashed)
            all_done: bool = True  # default val, re-define if something is not done
            for job_manager in jobs_managers:
                # check if each process is still running
                rc = job_manager.check_status()
                if rc is None:
                    all_done = False
            if all_done:
                self._to_log("All parallel jobs done ###")
                break
            else:
                if iter_num:
                    self._to_log(
                        'Polling iteration %i; already waiting %i seconds.' % (
                            iter_num, iter_num * Constants.ITER_DURATION
                        )
                    )
                time.sleep(Constants.ITER_DURATION)
                iter_num += 1

        if any(jm.check_status() != 0 for jm in jobs_managers):
            self._to_log(
                'Some of the monitored parallel processes at %s step died' % step,
                'warning'
            )

    def _schedule_selected_batches(self, step: str) -> None:
        """Prepares a short job list for selected batches for a given step"""
        jobfiles_to_execute: Dict[int, List[str]] = defaultdict(list)
        joblist2jobs: Dict[str, List[str]] = defaultdict(list)
        step2joblist: Dict[str, str] = {
            'feature_extraction': 'feature_extraction_joblist',
            'preprocessing': 'cesar_preprocess_joblist'
        }
        if step == 'feature_extraction':
            batch_ids: List[str] = self.selected_alignment_batches
            job_dir: str = self.feature_job_dir
            output_dir: str = self.feature_res_dir
            joblists: Dict[int, str] = {
                0: self.__getattribute__(step2joblist[step])
            }
        elif step == 'preprocessing':
            batch_ids: List[str] = self.selected_preprocessing_batches
            job_dir: str = self.preprocessing_job_dir
            output_dir: str = self.preprocessing_res_dir
            joblists: Dict[int, str] = {
                0: self.__getattribute__(step2joblist[step])
            }
        elif step == 'alignment':
            batch_ids: List[str] = self.selected_alignment_batches
            job_dir: str = self.alignment_job_dir
            output_dir: str = self.alignment_res_dir
            if not os.path.exists(self.cesar_job_list_summary):
                self._die(
                    (
                        'Job list summary %s for the alignment step does not exists; ' 
                        'consider re-running the alignment step altogether'
                    ) % self.cesar_job_list_summary
                )
            joblists: Dict[int, str] = {}
            with open(self.cesar_job_list_summary, 'r') as h:
                for line in h:
                    joblist, mem = line.rstrip().split('\t')
                    joblists[int(mem)] = joblist
                    with open(joblist, 'r') as h:
                        for line in h.readlines():
                            joblist2jobs[int(mem)].append(line.rstrip())
        else:
            self._die(
                'Step %s is not a valid parallel process step name' % step
            )
        for batch_id in batch_ids:
            job_file: str = os.path.join(job_dir, f'batch{batch_id}.ex')
            if not os.path.exists(job_file):
                self._die(
                    'Batch %s is missing for step %s (file %s not found)' % (batch_id, step, job_file)
                )
            if step == 'alignment':
                for joblist, jobs in joblist2jobs.items():
                    if job_file in jobs:
                        jobfiles_to_execute[joblist].append(job_file)
                        break
                else:
                    self._die(
                        'Batch %s does not belong to any of the alignment job lists from the previous runs' % job_file
                    )
            else:
                out_path: str = os.path.join(output_dir, f'batch{batch_id}')
                self._rm(out_path)
                jobfiles_to_execute[0].append(job_file)
        if step == 'alignment':
            new_summary_file: str = f'{self.cesar_job_list_summary}_partial_{self.project_name}'
            self.cesar_job_list_summary = new_summary_file
            context = open(self.cesar_job_list_summary, 'w')
        else: 
            context = nullcontext()
        with context as h:
            for i, joblist in joblists.items():
                new_joblist_name: str = f'{joblist}_partial_{self.project_name}'
                jobs: List[str] = jobfiles_to_execute[i]
                if not jobs:
                    continue
                # print(f'{i=}, {joblist=}, {jobs=}')
                if step == 'alignment':
                    h.write(f'{new_joblist_name}\t{i}\n')
                else:
                    self.__setattr__(
                        step2joblist[step], new_joblist_name
                    )
                with open(new_joblist_name, 'w') as h:
                    for job in jobs:
                        h.write(job + '\n')

    def _write_failed_batches_and_exit(self, step: str) -> None:
        """
        Produce a message informing that some of the parallel process batches died, 
        then write the batch IDs and exit
        """
        self.write_failed_batches()
        self._die(
            'Some of the monitored parallel processes at %s step died' % step 
        )

    def _create_output_stub(self, file: str) -> None:
        """
        Creates a stub file with a header for a given output file
        """
        from constants import Headers
        if file not in Constants.FILE2HEADER:
            return
        filepath: str = self.__getattribute__(file)
        header_line: str = getattr(Headers, Constants.FILE2HEADER[file])
        with open(filepath, 'w') as h:
            h.write(header_line)

    def _clean_previous_results(self) -> None:
        """Cleans the output from the previous runs potentially interfering with the new run"""
        remove_from: int = Constants.RESUME_ORDER[self.resume_from]
        # if remove_from < 2:
        #     for d in Constants.CLEANUP_TARGETS['all']:
        #         self._rmdir(d)
        #     return
        for step, step_rank in Constants.RESUME_ORDER.items():
            ## do not remove the data from the upstream steps
            if step_rank < remove_from:
                continue
            cleanup_targets: List[str] = Constants.CLEANUP_TARGETS[step]
            ## if certain batches were prompted to be re-run at parallel steps, 
            ## do not clean respective directories out; those will be cleaned up later   
            if step == 'feature_extraction' and self.selected_feature_batches is not None:
                cleanup_targets = cleanup_targets[-1:]
            if step == 'preprocessing' and self.selected_preprocessing_batches is not None:
                cleanup_targets = cleanup_targets[-2:]
            if step == 'alignment' and self.selected_alignment_batches is not None:
                continue
            for file in cleanup_targets:
                filepath: str = getattr(self, file)
                self._rm(filepath)
            # if step == 'orthology_resolution' and self.selected_orthology_batches:
            #     continue
        if remove_from > Constants.CESAR_AGGREGATION_RANK:
            self._to_log('Decompressing the annotation step results from the previous run')
            for file in Constants.FILES_TO_GZIP:
                filepath: str = f'{getattr(self, file)}.gz'
                if not os.path.exists(filepath):
                    continue
                cmd: str = f'gunzip {filepath}'
                _ = self._exec(cmd, 'Decompression failed for file %s:' % filepath)

    def _email(self, subject: str, text: str) -> None:
        """Sends an e-mail notification to the provided mail box via mailx utilitiy"""
        cmd: str = Constants.MAILX_TEMPLATE.format(
            text, self.mailx_binary, subject, self.email
        )
        subprocess.call(cmd, shell=True)
        # _ = self._exec(cmd, "E-mail notification via mailx failed: ")

    def create_wd(self) -> None:
        """Sets the output directory structure"""
        ## create the zero-level output directory
        if os.path.exists(self.output):
            if os.path.exists(self.logs):
                self.set_logging()
            else:
                self._mkdir(self.logs)
                self.set_logging()
            self._to_log(
                'Output directory does exist; cleaning up the conflicting files from the previous runs',
                'warning'
            )
            self._clean_previous_results()
        else:
            self._mkdir(self.output)
            self._mkdir(self.logs)
            self.set_logging()

        ## create the first-level directories
        ## (for results, metadata, and temporary files)
        self._mkdir(self.tmp)
        self._mkdir(self.meta)
        # self._mkdir(self.res)
        ## create subdirectories for temporary files' folder
        self._mkdir(self.input_data)
        self._mkdir(self.nextflow_dir)
        self._mkdir(self.feature_job_dir)
        self._mkdir(self.feature_res_dir)
        self._mkdir(self.preprocessing_job_dir)
        self._mkdir(self.preprocessing_res_dir)
        self._mkdir(self.alignment_job_dir)
        self._mkdir(self.alignment_res_dir)
        self._mkdir(self.orthology_resolution_dir)
        self._mkdir(self.finalized_output_dir)
        ## and the same for metadata directory
        self._mkdir(self.rejection_dir)
        self._mkdir(self.vis_input_dir)
        ## finally, create the UCSC subdirectory in the 'results' section
        self._mkdir(self.classification_dir)
        self._mkdir(self.annot_dir)
        # self._mkdir(self.orthology_results_dir)
        self._mkdir(self.ucsc_dir)

    def check_arguments(self) -> None:
        """
        Performs various argument sanity checks
        """
        ## if CESAR memory limit was set, check that no memory bin exceeds it
        if self.cesar_memory_limit:
            memory_bins_split: List[Union[int, str]] = [
                (int(x) if x != 'big' else x)
                for x in self.cesar_memory_bins.split(',') if x
            ]
            job_bins_split: List[int] = [
                int(x) for x in self.job_nums_per_bin.split(',') if x
            ]
            deprecated_bins: List[int] = []
            for i in range(len(memory_bins_split)):
                mem_bin: Union[int, str] = memory_bins_split[i]
                if memory_bins_split[i] == 'big':
                    continue
                if memory_bins_split[i] > self.cesar_memory_limit:
                    self._to_log(
                        f'Memory bin {i} requires {mem_bin} GB of memory, '
                        f'exceeding the set memory limit ({self.cesar_memory_limit}); '
                        'removing the bin',
                        'warning'
                    )
                    deprecated_bins.append(i)
            self.cesar_memory_bins = ','.join(
                map(
                    str, [
                        x for i, x in enumerate(memory_bins_split)
                        if i not in deprecated_bins
                    ]
                )
            )
            self.job_nums_per_bin = ','.join(
                map(
                    str, [
                        x for i, x in enumerate(job_bins_split)
                        if i not in deprecated_bins
                    ]
                )
            )
        ## check whether loss classes accepted for orthology resolution
        ## are consistent with the TOGA notation
        if self.accepted_loss_symbols == 'ALL':
            self.accepted_loss_symbols = ','.join(Constants.ALL_LOSS_SYMBOLS)
            
        else:
            for x in self.accepted_loss_symbols.split(','):
                if not x:
                    continue
                if x not in Constants.ALL_LOSS_SYMBOLS:
                    self._die(
                        f'ERROR: symbol {x} in the --accepted_loss_symbols '
                        'argument does not correspond to the TOGA notation. '
                        f'Accepted symbols are: {",".join(Constants.ALL_LOSS_SYMBOLS)}'
                    )

    def check_binaries(self) -> None:
        """
        Checks if all third-party binaries are accessble
        """
        ## check bigWigToWig binary
        for attr, default_name in Constants.BINARIES_TO_CHECK.items():
            if self.__getattribute__(attr) is None:
                self._to_log('Looking for %s in PATH' % default_name)
                exe_in_path: Union[str, None] = which(default_name)
                if exe_in_path is None:
                    if default_name == 'mailx':
                        if self.email is None:
                            continue
                        else:
                            self._die(
                                'E-mail notification was requested by the user but executable "mailx" '
                                'was not provided, with no defaults in $PATH'
                            )
                    self._die(
                        'Path to %s executable was not provided, with no defaults in $PATH' % default_name
                    )
                self.__setattr__(attr, exe_in_path)
        # if self.bigwig2wig_binary is None:
        #     self._to_log('Looking for bigWigToWig in PATH', 'info')
        #     bw2w_in_path: str = which('bigWigToWig')
        #     if bw2w_in_path is not None:
        #         self.bigwig2wig_binary: click.Path = os.path.abspath(bw2w_in_path)#Path(bw2w_in_path).absolute()
        #     else:
        #         if os.path.exists(HL_BW2W_PATH):
        #             self._to_log('Falling back to standard bigWigToWig', 'info')
        #             self.bigwig2wig_binary: click.Path = HL_BW2W_PATH
        #         else:
        #             self._to_log('ERROR: bigWigToWig is not accessible', 'critical')
        #             self._die('bigWigToWig binary is not accessible')
        ## TODO: add the same support for the following binaries:

        if self.cesar_binary is None:
            self._to_log('Falling back to standard CESAR2.0', 'info')
            self.cesar_binary = os.path.abspath(
                os.path.join(LOCATION, 'CESAR2.0', 'cesar')
            )
        # cds_script_in_path: str = which(self.CDS_TRACK_SCRIPT)
        # if cds_script_in_path is None:
        #     self.CDS_TRACK_SCRIPT = os.path.join(LOCATION, 'bin', 'bed12ToCDSOnly')

    def check_spliceai_files(self) -> None:
        """
        Checks the SpliceAI directory content for compliance with
        cesar_preprocess.py expected input format
        """
        if not self.spliceai_dir:
            return
        self._to_log('Checking SpliceAI directory structure', 'info')
        for file in Constants.SPLICEAI_FILES:
            file_path: str = os.path.join(self.spliceai_dir, file)
            if not os.path.exists(file_path):
                self._to_log(
                    'ERROR: file %s is missing from the %s directory' % (file, self.spliceai_dir),
                    'critical'
                )
                self._die(
                    f'ERROR: file {file} is missing '
                    f'from the {self.spliceai_dir} directory'
                )
        self._to_log('SpliceAI directory check complete', 'info')

    def _generate_nf_script(self) -> None:
        """
        Generates a Nextflow execution master script from a minimal boilerplate
        """
        if self.parallel_strategy in ('para', 'custom'):
            return
        nf_contents: str = Constants.NEXTFLOW_STUB.format(self.max_number_of_retries)
        nf_file: str = os.path.join(self.nextflow_dir, 'execute_joblist.nf')
        self.nextflow_exec_script = nf_file
        with open(nf_file, 'w') as h:
            h.write(nf_contents + '\n')

    def _check_nextflow_configs(self) -> None:
        """Checks Nextflow configuration file directory contents"""
        if self.nextflow_config_dir is None:
            return
        expected_configs: Dict[str, str] = {
            **Constants.UNIQUE_CONFIGS, 
            **{x: Constants.ALN_CONFIG.format(x) for x in self.cesar_memory_bins.split(',')}
        }
        for process, file in expected_configs.items():
            config_path: str = os.path.join(self.nextflow_config_dir, file)
            if not os.path.exists(config_path):
                self._to_log(
                    (
                        'File %s is missing from the Nextflow configuration file directory %s; '
                        'defaulting to TOGA2 Nextflow config template'
                    ) % (file, self.nextflow_config_dir),
                    'warning'
                )
                file = None
            self.nextflow_config_files[process] = file
        ## TODO: Add the file contents inspection?

    def filter_chain_file(self) -> None:
        """
        Filters the chain file by their score
        """
        cmd: str = (
            f'{self.CHAIN_FILTER_SCRIPT} -i {self.chain_file} '
            f'-s {self.min_chain_score} -o {self.chain_file_copy}'
        )
        _ = self._exec(cmd, 'Chain filtering step failed')

    def get_contig_sizes(self, ref: bool) -> None:
        """
        Infers contig sizes and saves them in a two-column tab-separated file
        """
        twobit_file: str = self.ref_2bit if ref else self.query_2bit
        dest_file: str = (
            self.ref_contig_size_file if ref else self.query_contig_size_file
        )
        twobit_cmd: str = (
            f'{Constants.SETUP} {self.twobittofa_binary} {twobit_file} stdout | '
            f'{self.CONTIG_SIZE_SCRIPT} - -o {dest_file}'
        )
        twobit_out: str = self._exec(twobit_cmd, err_msg='', die=False)
        if 'is not a twoBit file' not in twobit_out:
            return
        compress_check_cmd: str = f'{Constants.SETUP} file {twobit_file} | grep -q "compressed"'
        is_compressed: bool = bool(self._exec(compress_check_cmd), '')
        if is_compressed:
            gz_cmd: str = (
                f'{Constants.SETUP} gzip -dc {twobit_file} | '
                f'{self.CONTIG_SIZE_SCRIPT} - -o {dest_file}'
            )
            _ = self._exec(
                gz_cmd,
                'Contig size retrieval from the gz-compressed genome file failed'
            )
            return
        txt_cmd: str = (
            f'{self.CONTIG_SIZE_SCRIPT} {twobit_file} -o {dest_file}'
        )
        _ = self._exec(txt_cmd, 'Contig size retrieval from a FASTA file failed')


    def index_chains(self) -> None: ## TODO: Rustify this as soon as possible
        """
        Creates a chain ID: start byte mapping for the input chain file
        """
        cmd: str = (
            f'{self.INDEX_CHAIN_SCRIPT} {self.chain_file_copy} '
            f'{self.chain_index} -t {self.chain_index_txt}'
        )
        self._exec(cmd, 'Chain indexing failed')


    def prepare_data(self) -> None:
        """
        """
        ## filter the chain file; if filtering is disabled, simply synchronize
        ## the local copy with the original chain file
        if not self.min_chain_score:
            self._to_log('Synchronizing chain file instance', 'info')
            self._rsync(self.chain_file, self.chain_file_copy)
        else:
            self._to_log(
                'Filtering out chains with score less than %s' % self.min_chain_score,
                'info'
            )
            self.filter_chain_file()
            self._to_log('Chain filtering complete', 'info')

        ## index the chain file
        self._to_log('Indexing the chain file', 'info')
        self.index_chains()
        self._to_log('Chain indexing complete', 'info')

        ## infer contig sizes
        self._to_log('Retrieving contig sizes for the reference genome')
        self.get_contig_sizes(ref=True)
        self._to_log('Retrieving contig sizes for the query genome')
        self.get_contig_sizes(ref=False)
        self._to_log('Contig sizes retrieved')

        ## clip the UTRs from the reference Bed records 
        ## and dump the data to the HDF5 storage
        ## TODO: Repurpose all the BED-reading scripts to extract data from the HDF5 storage
        self._to_log('Formatting reference BED annotation')
        self.prepare_ref_bed_data()
        self._to_log('Reference annotation formatting complete')

        ## if --u12_file argument was provided, convert IntronIC pipeline output 
        ## into an HDF5 storage
        self.prepare_u12_data()

    def prepare_ref_bed_data(self) -> None:
        """
        Processes reference BED file for future use by downstream scripts and dumps
        the results in HDF5 format.
        By default, TOGA 2.0 creates two HDF5 storage files:
        * tmp/input_data/full_annotation.hdf5 : contains full reference transcript
          data; these are used at chain classification step;
        * tmp/input_data/cds_annotation.hdf5 : contains reference transcripts
          after removing UTR exons; these are used at CESAR alignment step
        """
        ## filter the reference BED file first ## TODO: Accommodate for the optional arguments
        from filter_ref_bed import AnnotationFilter
        args: List[str] = [
            self.ref_annotation, self.bed_file_copy, self.prefiltered_transcripts,
            '-ln', self.project_name
        ]
        self._to_log('Filtering reference BED file')
        AnnotationFilter(args, standalone_mode=False)

        ## dump the resulting data to the HDF5 storage
        from bed_hdf5_index import BedHdf5Indexer
        if self.legacy_chain_feature_extraction:
            self._to_log('Converting reference transcript BED into HDF5 format')
            BedHdf5Indexer([self.ref_annotation, self.full_bed_hdf5], standalone_mode=False)

        self._to_log('Stripping non-coding exons from the reference transcript entries')
        cmd3: str = f'{self.CDS_TRACK_SCRIPT} -i {self.bed_file_copy} -o {self.cds_bed_file} -m cds'
        _ = self._exec(cmd3, 'CDS track preparation failed')
        self._to_log('Converting reference CDS BED into HDF5 format')
        BedHdf5Indexer([self.cds_bed_file, self.cds_bed_hdf5], standalone_mode=False)
        ## strip CDS from the unfiltered BED track
        cmd4: str = f'{self.CDS_TRACK_SCRIPT} -i {self.ref_annotation} -o {self.ref_cds_unfilt} -m cds'
        _ = self._exec(cmd4, 'Unfiltered CDS track preparation failed:')


    def prepare_u12_data(self) -> None:
        """
        Converts IntronIC pipelin data into an HDF5 storage for easier further access
        """
        if self.u12_file is None:
            return
        self._to_log('Formatting U12 data')
        from intronIC_to_hdf5 import IntronIcConverter
        ## TODO: Rust implementation?
        args: List[str] = [
            self.cds_bed_hdf5, self.u12_file, self.u12_hdf5, '--hdf5_input', 
            '-ln', self.project_name
        ]
        IntronIcConverter(args, standalone_mode=False)
        self._to_log('U12 data formatting complete')

    def extract_chain_features(self) -> None:
        """
        A single-thread version of feature extraction procedure, 
        using a Rust implementation of chain_runner.py and merge_chains_output.py
        """
        cmd: str = (
            f'{self.FEATURE_EXTRACTOR} -b {self.bed_file_copy} '
            f'-c {self.chain_file_copy} -o {self.feature_table} '
        )
        if self.isoform_file is not None:
            cmd += f' -i {self.isoform_file}'
        _ = self._exec(cmd, 'Feature extraction step failed:', gather_stdout=False)

    def schedule_feature_jobs(self) -> None:
        """
        Schedule chain-transcript feature extraction jobs. 
        NOTE: This is a legacy feature borrowed from original TOGA and replaced by 
        single-core implementation in Rust. Do not use it unless you experience problems 
        with the Rust code.
        """
        from schedulers.chain_feature_scheduler import ChainFeatureScheduler
        args: List[str] = [
            self.chain_file_copy, self.bed_file_copy, 
            self.feature_job_dir, self.feature_data_dir, self.feature_res_dir, 
            '-j', f'{self.feature_job_num}',
            '-r', self.feature_rejection_log, '-ln', self.project_name
        ]
        ChainFeatureScheduler(args, standalone_mode=False)

    def run_feature_extraction_jobs(self) -> None:
        """
        Runs scheduled feature extraction jobs.
        NOTE: This is a legacy feature borrowed from original TOGA and replaced by 
        single-core implementation in Rust. Do not use it unless you experience problems 
        with the Rust code.
        """
        project_name: str = f'chain_feats_{self.project_name}'
        self.parallel_process_names.append(project_name)
        project_path: str = os.path.join(self.nextflow_dir, project_name)
        self._run_parallel_process(
            project_name, project_path, self.feature_extraction_joblist,
            queue=self.cluster_queue_name,
            step='feature extraction'
        )

    def create_feature_dataset(self) -> None:
        """
        Aggregates data over feature extracting jobs' output files and creates a
        dataset suitable for further classification
        """
        from merge_chains_output import FeatureAggregator
        for batch_id in range(0, self.feature_job_num):
            batch_name: str = f'batch{batch_id}'
            batch_job_file: str = os.path.join(self.feature_job_dir, f'{batch_name}.ex')
            if not os.path.exists(batch_job_file):
                break
            batch_path: str = os.path.join(self.feature_res_dir, batch_name)
            batch_exec_mark: str = f'{batch_name}_ok'
            batch_exec_path: str = os.path.join(self.feature_res_dir, batch_exec_mark)
            if not os.path.exists(batch_path) or not os.path.exists(batch_exec_path):
                self._to_log(
                    'Feature extraction batch %s has not finished properly' % batch_name,
                    'warning'
                )
                self.failed_feature_batches.append(batch_name)
        if self.failed_feature_batches and not self.ignore_crashed_parallel_batches:
            self._write_failed_batches_and_exit('feature extraction')
        args: List[str] = [
            self.feature_res_dir, self.bed_file_copy, self.feature_table,
            '-ln', self.project_name
        ]
        if self.isoform_file is not None:
            args.extend(('--isoforms', self.isoform_file))
        FeatureAggregator(args, standalone_mode=False)

    def classify_chains(self) -> None:
        """
        Run XGBoost classifier for projection orthologu classification
        """
        
        ## if any of the models are missing by chance, train the classifier first
        if not os.path.isfile(self.se_model) or not os.path.isfile(self.me_model):
            self._exec(self.MODEL_TRAINER)

        ## classify chains
        ## TODO: Transcripts discarded at this step are considered as missing in TOGA 1.0
        from classify_chains import ChainClassifier
        args: List[str] = [
            self.feature_table, self.classification_dir, self.se_model, self.me_model,
            '-t', f'{self.orthology_threshold}', '-ln', self.project_name,
            '-minscore', self.min_orth_chain_score
        ]
        if self.use_ld_model:
            args.extend(('--ld_model', self.ld_model))
        if self.legacy_chain_feature_extraction:
            args.append('--legacy')
        ChainClassifier(args, standalone_mode=False)

        orthology_scores_tmp: str = os.path.join(
            self.classification_dir, 'orthology_scores.tsv'
        )
        _ = self._exec(
            f'cp {orthology_scores_tmp} {self.pred_scores}',
            'Orthology score file copying failed'
        )
        tr2chain_tmp: str = os.path.join(
            self.classification_dir, 'trans_to_chain_classes.tsv'
        )
        _ = self._exec(
            f'cp {tr2chain_tmp} {self.tr2chain_classes}',
            'Projection class table copying failed'
        )
        rej_log: str = os.path.join(
            self.classification_dir, 'rejection_report.tsv'
        )
        if os.path.exists(rej_log):
            rejection_move_cmd: str = f'mv {rej_log} {self.class_rejection_log}'
            _ = self._exec(
                rejection_move_cmd,
                'Moving classification rejection log failed'
            )
        else:
            self._to_log('No items were rejected at classification step')
        ## TODO: Here goes a sanity check

    def recover_fragmented_projections(self) -> None:
        """
        Identifies fragmented projections from the orthology classifcation data,
        reports the chains used to restore these projections in a two-column
        tab-separated file
        """
        if self.disable_fragment_assembly:
            return
        from stitch_fragments import main
        args: List[str] = [
            self.chain_file_copy, self.pred_scores, self.bed_file_copy,
            '-o', self.fragmented_projection_list,
            '--orthology_threshold', f'{self.orthology_threshold}',
            '--fragmented_only'
        ]
        main(args, standalone_mode=False)

    def schedule_preprocessing_jobs(self) -> None:
        """
        Schedules CESAR preprocessing jobs for further cluster run
        """
        from schedulers.preprocess_scheduler import PreprocessingScheduler
        kwargs: Dict[str, Any] = {
            'chain_map': self.tr2chain_classes,
            'ref_annotation': self.cds_bed_file,
            'ref': self.ref_2bit,
            'query': self.query_2bit, 
            'chain_file': self.chain_index, 
            'ref_chrom_sizes': self.ref_contig_size_file,
            'query_chrom_sizes': self.query_contig_size_file, 
            'job_directory': self.preprocessing_job_dir,
            'preprocessing_directory': self.preprocessing_res_dir,
            'job_number': self.preprocessing_job_num,
            'max_chain_number': self.max_chains_per_transcript,
            'max_space_size': self.max_search_space_size,
            'extrapolation_modifier': self.extrapolation_modifier,
            'minimal_covered_fraction': self.minimal_covered_fraction,
            'exon_locus_flank': self.exon_locus_flank,
            'cesar_canon_u2_acceptor': self.cesar_canon_u2_acceptor,
            'cesar_canon_u2_donor': self.cesar_canon_u2_donor,
            'cesar_non_canon_u2_acceptor': self.cesar_non_canon_u2_acceptor,
            'cesar_non_canon_u2_donor': self.cesar_non_canon_u2_donor,
            'cesar_canon_u12_acceptor': self.cesar_canon_u12_acceptor,
            'cesar_canon_u12_donor': self.cesar_canon_u12_donor,
            'cesar_non_canon_u12_acceptor': self.cesar_non_canon_u12_acceptor,
            'cesar_non_canon_u12_donor': self.cesar_non_canon_u12_donor,
            'cesar_first_acceptor': self.cesar_first_acceptor,
            'cesar_last_donor': self.cesar_last_donor,
            'assembly_gap_size': self.assembly_gap_size,
            'paralog_report': self.paralog_report,
            'processed_pseudogene_report': self.processed_pseudogene_report,
            'rejection_report': self.preprocessing_rejection_log,
            'log_name': self.project_name,
            'verbose': True
        }
        if self.toga1 and not self.toga1_plus_cesar:
            kwargs['toga1_compatible'] = True
        if self.toga1_plus_cesar:
            kwargs['toga1_plus_cesar'] = True
        if self.cesar_memory_limit:
            kwargs['memory_limit'] = self.cesar_memory_limit
        if not self.disable_fragment_assembly and os.path.exists(self.fragmented_projection_list):
            kwargs['fragmented_projections'] = self.fragmented_projection_list
        if self.orthologs_only:
            kwargs.append['orthologs_only'] = True
        if self.one2ones_only:
            kwargs['one2one_only'] = True
        if not self.enable_spanning_chains:
            kwargs['disable_spanning_chains'] = True
        if self.u12_file is not None:
            kwargs['u12'] = self.u12_hdf5
        if self.separate_site_treat:
            kwargs['separate_splice_site_treatment'] = True
        if self.spliceai_dir is not None:
            kwargs['spliceai_dir'] = self.spliceai_dir
            kwargs['bigwig2wig_binary'] = self.bigwig2wig_binary
            kwargs['min_splice_prob'] = self.min_splice_prob
        if self.annotate_ppgenes:
            kwargs['annotate_processed_pseudogenes'] = True
        PreprocessingScheduler(**kwargs)

    def run_preprocessing_jobs(self) -> None:
        """
        Runs CESAR preprocessing jobs
        """
        project_name: str = f'cesar_preprocess_{self.project_name}'
        self.parallel_process_names.append(project_name)
        project_path: str = os.path.join(self.nextflow_dir, project_name)
        self._run_parallel_process(
            project_name, 
            project_path, 
            self.cesar_preprocess_joblist, 
            queue=self.cluster_queue_name,
            step='preprocessing',
            cpu=1,
            process_num=self.preprocessing_job_num,
            config_file=self.nextflow_config_files.get('preprocessing', None)
        )

    def aggregate_preprocessing_res(self) -> None:
        """
        Aggregates preprocessing step results
        """
        self._create_output_stub('preprocessing_report')
        # if not self.enable_spanning_chains:
        self._create_output_stub('spanning_chain_coords')
        for dir_name in os.listdir(self.preprocessing_res_dir):
            dir_path: str = os.path.join(self.preprocessing_res_dir, dir_name)
            ok_file: str = os.path.join(dir_path, Constants.OK_FILE)
            if not os.path.exists(ok_file):
                self._to_log(
                    'Preprocesssing batch %s has not finished properly' % dir_name,
                    'warning'
                )
                self.failed_preprocessing_batches.append(dir_name)
            mem_path: str = os.path.join(dir_path, 'max_memory_requirements.tsv')
            if os.path.exists(mem_path):
                aggr_cmd: str = f'{Constants.SETUP} cat {mem_path} >> {self.preprocessing_report}'
                _ = self._exec(
                    aggr_cmd, 
                    'Preprocessing aggregation data failed at batch %s' % dir_name
                )
            rej_path: str = os.path.join(dir_path, 'genes_rejection_reason.tsv')
            if os.path.exists(rej_path):
                rej_aggr_cmd: str = (
                    f'{Constants.SETUP} cat {rej_path} >> {self.preprocessing_rejection_log}'
                )
                _ = self._exec(
                    rej_aggr_cmd, 
                    'Preprocessing step rejection log aggregation failed at batch %s' % dir_name
                )
            span_path: str = os.path.join(dir_path, 'spanning_chains_ref_coords.tsv')
            if os.path.exists(span_path):
                span_aggr_cmd: str = f'{Constants.SETUP} cat {span_path} >> {self.spanning_chain_coords}'
                _ = self._exec(
                    span_aggr_cmd, 
                    'Spanning chains coordinates aggregation failed at batch %s' % dir_name
                )
        if self.failed_preprocessing_batches and not self.ignore_crashed_parallel_batches:
            self._write_failed_batches_and_exit('preprocessing')
        ## TODO: Here goes a sanity check

    def schedule_alignment_jobs(self) -> None:
        """
        Schedules CESAR preprocessing jobs for further cluster run
        """
        from schedulers.cesar_job_scheduler import CesarScheduler
        args: List[str] = [
            self.preprocessing_report, self.alignment_job_dir, self.alignment_res_dir,
            '-b', self.cesar_memory_bins, '-jb', self.job_nums_per_bin,
            '-cs', self.cesar_binary, '-m', self.matrix_file,
            '-rr', self.redundancy_rejection_log,
            '-scm', self.spliceai_correction_mode
        ]
        if self.toga1 and not self.toga1_plus_cesar:
            args.append('-t1')
        if not self.toga1 or self.toga1 and self.toga1_plus_cesar:
            args.append('-c_si') ## Ultrashort intron correction
        if os.path.exists(self.paralog_report):
            args.extend(('-pl', self.paralog_report))
        if os.path.exists(self.processed_pseudogene_report):
            args.extend(('-ppl', self.processed_pseudogene_report))
        if self.allow_heavy_jobs:
            args.append('--ahj')
        if self.mask_terminal_mutations:
            args.append('-m10m')
        if not self.leave_missing_stop:
            args.append('-rms')
        if not self.consider_alt_frame:
            args.append('-no_alt_frame')
        if self.spliceai_dir is not None:
            args.extend(('-msp', self.min_splice_prob, '-spm', self.splice_prob_margin))
        if self.intron_gain_check and self.spliceai_dir is not None:
            args.extend(
                (
                    '-ig', '-igt', self.min_intron_gain_score,
                    '-mipt', self.min_intron_prob_trusted,
                    '-mips', self.min_intron_prob_supported,
                    '-mipu', self.min_intron_prob_unsupported,
                    '-cra', self.cesar_canon_u2_acceptor,
                    '-crd', self.cesar_canon_u2_donor,
                    '-cfa', self.cesar_first_acceptor,
                    '-cld', self.cesar_last_donor
                )
            )
        CesarScheduler(args, standalone_mode=False)

    def run_alignment_jobs(self) -> None:
        """
        Runs CESAR alignment jobs
        """
        cesar_managers: List[ParallelJobsManager] = []
        try:
            with open(self.cesar_job_list_summary, 'r') as h:
                for line in h:
                    joblist, mem = line.strip().split('\t')
                    if not os.path.exists(joblist):
                        self._to_log(
                            'Alignment joblist %s does not exist' % joblist, 'warning'
                        )
                        continue
                    job_num: int = int(self.job_nums_per_bin[self.cesar_memory_bins.index(mem)])
                    project_name: str = f'cesar_align_{self.project_name}_{mem}'
                    self.parallel_process_names.append(project_name)
                    project_path: str = os.path.join(self.nextflow_dir, project_name)
                    manager: ParallelJobsManager = self._run_parallel_process(
                        project_name,
                        project_path,
                        joblist,
                        wait_for_process=False,
                        memory_limit=int(mem),
                        queue=self.cluster_queue_name,
                        process_num=job_num,
                        cpu=1,
                        config_file=self.nextflow_config_files.get(mem, None)
                    )
                    cesar_managers.append(manager)
                    time.sleep(2)
            self._monitor_jobs(cesar_managers, 'alignment')
        except KeyboardInterrupt:
            self._to_log('Keyboard interrupt - Terminating all parallel processes', 'warning')
            self._terminate_parallel_processes(cesar_managers)
            self._exit()

    def aggregate_cesar_results(self) -> None:
        """
        Aggregates the results of CESAR alignment step
        """
        batch_dirs: List[str] = os.listdir(self.alignment_res_dir)
        for dir_name in batch_dirs:
            dir_path: str = os.path.join(self.alignment_res_dir, dir_name)
            ok_file: str = os.path.join(dir_path, Constants.OK_FILE)
            if not os.path.exists(ok_file):
                self._to_log(
                    'CESAR alignment batch %s has not finished properly' % dir_name,
                    'warning'
                )
                self.failed_alignment_batches.append(dir_name)
            for out_file in Constants.CESAR_OUT_FILES:
                batch_path: str = os.path.join(dir_path, out_file)
                out_file_slot: str = Constants.CESAR_FILE_TO_DEST[out_file]
                aggr_path: str = self.__getattribute__(out_file_slot)
                if not os.path.exists(aggr_path):
                    self._create_output_stub(out_file_slot)
                if not os.path.exists(batch_path):
                    continue
                cmd: str = f'cat {batch_path} >> {aggr_path}'
                _ = self._exec(cmd, f'File aggregation failed at file {batch_path}')
        if self.failed_alignment_batches and not self.ignore_crashed_parallel_batches:
            self._write_failed_batches_and_exit('alignment')

    def aggregate_rejection_reports(self) -> None:
        """
        Aggregates rejection reports from various stages into a final report
        """
        self._create_output_stub('final_rejection_log')
        if len(os.listdir(self.rejection_dir)) == 0:
            self._to_log('No rejected items reported up to CESAR alignment step', 'warning')
            return
        cmd: str = f'cat {self.rejection_dir}/* >> {self.final_rejection_log}'
        _ = self._exec(cmd, 'Final rejection log aggregation failed')

    def loss_status_summary(self) -> None:
        """
        Summarizes sequence loss in the query on projection, 
        reference transcript, and reference gene levels
        """
        self._to_log(
            'Loss statuses considered for orthology annotation are: %s' % self.accepted_loss_symbols
        )
        from conservation_summary import main
        ## TODO: Needs a class representation for sure
        args: List[str] = [
            self.transcript_meta, '-r', self.final_rejection_log,
            '-o', self.gene_loss_summary
        ]
        if self.isoform_file is not None:
            args.extend(('-i', self.isoform_file))
        if os.path.exists(self.spanning_chain_coords):
            args.extend(('--spanning_chains_precedence_file', self.spanning_chain_coords))
        if os.path.exists(self.paralog_report):
            args.extend(('-p', self.paralog_report))
        if self.annotate_ppgenes:
            args.extend(('-pp', self.processed_pseudogene_report))
        main(args, standalone_mode=False)

    def prepare_pseudogene_track(self) -> None:
        """Prepares a BED9 track of projections classified as processed pseudogenes"""
        from prepare_pseudogene_track import PseudogeneTrackBuilder
        args: List[str] = [
            self.tr2chain_classes, self.chain_file_copy, 
            '-o', self.pseudogene_annotation, '-l', self.log_file, '-v'
        ]
        PseudogeneTrackBuilder(args, standalone_mode=False)

    def infer_query_genes(self) -> None:
        """
        Infers coding regions in the query genome, serving as proxies for query genes
        """
        from infer_query_genes import QueryGeneCollapser
        args: List[str] = [
            #self.query_annotation_filt, 
            self.query_exon_meta, self.query_genes_raw, 
            '-b', self.query_genes_bed_raw, '-ln', self.project_name,
            '-l', self.gene_loss_summary, 
            '-d', self.redundant_paralogs, '-dpp', self.redundant_ppgenes,
            '-pf', self.feature_table, '-op', self.pred_scores,
            '--insufficiently_covered_orthologs', self.discarded_overextended_projections,
            '-rl', self.gene_inference_rejection_log
        ]
        if self.isoform_file is not None:
            # args.extend(('-r', self.cds_bed_file, '-i', self.isoform_file))
            args.extend(('-r', self.ref_cds_unfilt, '-i', self.isoform_file))
        if os.path.exists(self.paralog_report):
            args.extend(
                ('-p', self.paralog_report)
            )
        if os.path.exists(self.processed_pseudogene_report):
            args.extend(
                ('-pp', self.processed_pseudogene_report)
            )
        QueryGeneCollapser(args, standalone_mode=False)
        if os.path.exists(self.gene_inference_rejection_log):
            rej_aggr_cmd: str = f'cat {self.gene_inference_rejection_log} >> {self.final_rejection_log}'
            _ = self._exec(
                rej_aggr_cmd, 
                'Adding rejection report for the gene inference step to the main rejection log failed'
            )
        else:
            self._to_log('No items were rejected at gene inference step')

    def convert_fasta_to_hdf5(self) -> None:
        """Converts TOGA2 output FASTA file into an HDF5 storage"""
        in_file: str = self.aa_fasta
        out_file: str = self.aa_hdf5
        from pairwise_fasta_to_hdf5 import FastaToHdf5Converter
        args: List[str] = [
            in_file, out_file, '-ln', self.project_name, '-v'
        ]
        FastaToHdf5Converter(args, standalone_mode=False)

    def resolve_orthology_initial(self) -> None:
        """
        Resolves orthology relationships between reference and query genes
        based on orthology graph structure; if specified by the user, schedules
        the jobs for further alignment/tree-based resolution of particularly
        entangled clades
        """

        from initial_orthology_resolver import InitialOrthologyResolver
        args: List[str] = [
            self.bed_file_copy, self.query_annotation_filt, self.gene_loss_summary,
            self.pred_scores, self.orthology_resolution_dir,
            '-qi', self.query_genes_raw, '-l', self.accepted_loss_symbols,
            '-mr', self.preprocessing_report, '-ln', self.project_name,
            # '-pf', self.feature_table
        ]
        if self.isoform_file is not None:
            args.extend(('-ri', self.isoform_file))
        if os.path.exists(self.paralog_report):
            args.extend(('-p', self.paralog_report))
        if os.path.exists(self.processed_pseudogene_report):
            args.extend(('-pp', self.processed_pseudogene_report))
        if not self.skip_tree_resolver:
            args.extend(
                (
                    '-st', '-mcs', self.max_clique_size, '-j', self.orth_job_num,
                    '-pb', self.prank_bin, '-rb', self.tree_bin,
                    '-rc', self.tree_cpus, '-f', self.aa_hdf5, '-fh',
                    '-jd', self.orthology_job_dir, '-fd', self.orthology_input_dir,
                    '-rd', self.orthology_res_dir
                )
            )
            if self.use_raxml:
                args.extend(('--use_raxml', '--tree_bootstrap', '100'))
            else:
                args.extend(('--tree_bootstrap', '5000'))
        InitialOrthologyResolver(args, standalone_mode=False)
        add_graph_rej_cmd: str = f'cat {self.rejected_by_graph} >> {self.final_rejection_log}'
        _ = self._exec(add_graph_rej_cmd, 'Adding rejected items from the orthology step failed')
        add_rej_to_discarded_cmd: str = (
            f'cat {self.weak_ortholog_names} >> {self.discarded_overextended_projections}'
        )
        _ = self._exec(add_rej_to_discarded_cmd, 'Merging files with discarded projections\' names failed')
        if self.skip_tree_resolver:
            self._to_log('Moving final orthology file')
            cmd: str = f'mv {self.temporary_orth_report} {self.orth_resolution_raw}'
            _ = self._exec(cmd, 'Moving orthology file to a finalization input directory failed')

    def run_tree_jobs(self) -> None:
        """
        Runs fine orthology resolution (PRANK + IQTree2/RAxML) jobs
        """
        project_name: str = f'tree_resolution_{self.project_name}'
        self.parallel_process_names.append(project_name)
        project_path: str = os.path.join(self.nextflow_dir, project_name)
        self._run_parallel_process(
            project_name, project_path, self.orth_resolution_joblist,
            queue=self.cluster_queue_name, 
            cpu=self.tree_cpus, 
            step='fine orthology resolution',
            process_num = self.orth_job_num,
            config_file=self.nextflow_config_files.get('orthology', None)
        )

    def resolve_orthology_final(self) -> None:
        """
        Updates the orthology resolution results, aggregating the batchwise
        results in the respective meta/ subdirectory
        """
        for batch in os.listdir(self.orthology_res_dir):
            ok_file: str = os.path.join(self.orthology_res_dir, batch, Constants.OK_FILE)
            if not os.path.exists(ok_file):
                self._to_log(
                    'Gene tree batch %s has not finished properly' % batch,
                    'warning'
                )
                self.failed_orthology_batches.append(batch)
        if self.failed_orthology_batches and not self.ignore_crashed_parallel_batches:
            self._write_failed_batches_and_exit('orthology')
        self._create_output_stub('resolved_leaves_file')
        leaf_aggr_cmd: str = (
            f'{Constants.SETUP} cat {self.orthology_res_dir}/*/resolved_pairs.tsv | cut -f1,2 >> '
            f'{self.resolved_leaves_file}'
        )
        _ = self._exec(
            leaf_aggr_cmd, 'ERROR: Resolved clades\' data aggregation failed'
        )
        unresolved_aggr_cmd: str = (
            f'cat {self.orthology_res_dir}/*/unresolved_clades.txt > '
            f'{self.unresolved_clades_file}'
        )
        _ = self._exec(
            unresolved_aggr_cmd, 'ERROR: Unresolved clades\' data aggregation failed'
        )
        ## TODO: Import the class instead!
        cmd: str = (
            f'{self.FINAL_RESOLVER_SCRIPT} {self.temporary_orth_report} '
            f'{self.resolved_leaves_file} -o {self.orth_resolution_raw} '
            f'-o2z {self.one2zero_genes} '
        )
        _ = self._exec(cmd, 'ERROR: Final orthology resolution failed')

    def gene_tree_summary(self) -> None:
        """
        Creates a summary table for the gene tree-based orthology resolution step
        """
        from gene_tree_summary import GeneTreeSummary
        args: List[str] = [
            self.orthology_input_dir, self.orthology_res_dir, 
            self.tree_summary_table, '-ln', self.project_name, '-v'
        ]
        if self.use_raxml:
            args.append('--raxml')
        GeneTreeSummary(args, standalone_mode=False)

    def deprecated_projection_file(self) -> None:
        """
        Creates a provisional list of deprecated projections 
        to be used for UCSC report production
        """
        deprecated_lists: List[str] = []
        if os.path.exists(self.discarded_overextended_projections):
            deprecated_lists.append(self.discarded_overextended_projections)
        if os.path.exists(self.redundant_paralogs):
            deprecated_lists.append(self.redundant_paralogs)
        if os.path.exists(self.redundant_ppgenes):
            deprecated_lists.append(self.redundant_ppgenes)
        if not deprecated_lists:
            return
        deprecated_aggr_cmd: str = 'cat ' + '\t'.join(deprecated_lists) + f' > {self.all_deprecated_projs}'
        _ = self._exec(deprecated_aggr_cmd, 'Aggregating deprecated projection lists failed:')

    def prepare_bigbed_track(self) -> None:
        """
        Prepares a BigBed22 track suitable for further loading to UCSC browser
        """
        ## TODO: Import the class here
        from make_ucsc_report import BigBedProducer
        args: List[str] = [
            self.aggr_ucsc_stub, self.bed_file_copy, self.feature_table,
            self.pred_scores, self.query_contig_size_file, self.SCHEMA_FILE,
            self.vis_input_dir, '-l', self.log_file, '--prefix', self.ucsc_prefix
        ]
        if self.ref_link_file:
            args.extend(['-i', self.ref_link_file])
        if os.path.exists(self.all_deprecated_projs):
            args.extend(['-d', self.all_deprecated_projs])
        if self.annotate_ppgenes:
            args.extend(['-pp', self.processed_pseudogene_report])
        if not self.skip_utr:
            args.extend(['-a', self.query_annotation_with_utrs])
        if self.v:
            args.append('-v')
        BigBedProducer(args, standalone_mode=False)

        for file in Constants.FINAL_UCSC_FILES:
            file = file.format(self.ucsc_prefix)
            from_path: str = os.path.join(self.vis_input_dir, file)
            to_path: str = os.path.join(self.ucsc_dir, file)
            copy(from_path, to_path)

    def gzip_fasta_files(self) -> None:
        """Compresses hefty output files into gzip format"""
        ## TODO: Modify ._exec() to support check_call()
        import subprocess
        for file in Constants.FILES_TO_GZIP:
            filepath: str = self.__getattribute__(file)
            if not os.path.exists(filepath):
                if os.path.exists(filepath + '.gz'):
                    self._to_log(
                        'File %s seems to be already gzipped' % filepath
                    )
                    continue
                else:
                    self._die(
                        'File %s does not exist' % filepath
                    )
            cmd: str = f'gzip -q -f -5 {filepath}'
            subprocess.check_call(cmd, shell=True)

    def annotate_utrs(self) -> None:
        """
        Projects untranslated regions (UTRs) and adds the projections to the respective transcripts 
        in the output Bed file
        """
        cmd: str = (
            f'{self.UTR_PROJECTOR_SCRIPT} -q {self.query_annotation_final} -r {self.bed_file_copy} '
            f'-c {self.chain_file_copy} -o {self.query_annotation_with_utrs} -e {self.query_exon_meta} '
            f'-a {self.utr_abs_threshold} -R {self.utr_rel_threshold}'
        )
        if self.no_utr_extrapolation:
            cmd += ' --no-extrapolation'
        if not self.no_adjacent_utrs:
            if self.fixed_adjacent_utrs:
                cmd += ' --fixed-adjacent-regions'
            else:
                cmd += ' --deduce-adjacent-regions'
        print(f'{cmd=}')
        _ = self._exec(cmd, 'UTR annotation failed: ')

    def convert_exon_fasta(self) -> None:
        """
        Creates a binary index of query exons in the gzipped exon Fasta file 
        """
        # from modules.index_gzipped_fasta import main
        from exon_fasta_to_twobit import TwoBitConverter
        if not os.path.exists(self.exon_fasta):
            if os.path.exists(self.exon_gzip):
                self._to_log(
                    'Exon FASTA file has been already gzipped; skipping conversion into 2bit',
                    'warning'
                )
                return
            self._die(
                'Exon FASTA file does not exist; dying at the exon FASTA to 2bit conversion step'
            )
        args: List[str] = [
            self.exon_fasta, self.exon_2bit, 
            '--fa2twobit', self.fatotwobit_binary,
            '-ln', self.project_name,
            '-e', self.query_exon_meta
        ]
        TwoBitConverter(args, standalone_mode=False)

    def filter_final_bed_files(self) -> None:
        """Cleans up final BED files from discarded projections and processed pseudogenes"""
        from filter_output_bed import OutputBedFilter
        args: List[str] = [
            self.query_annotation_filt, self.query_annotation_final
        ]
        if os.path.exists(self.discarded_overextended_projections):
            discarded_proj_args: List[str] = [
                '-di', self.all_deprecated_projs,
                '-do', self.discarded_proj_bed
            ]
            args.extend(discarded_proj_args)
        if self.annotate_ppgenes:
            ppgene_args: List[str] = [
                '-ppi', self.processed_pseudogene_report,
                '-ppo', self.processed_pseudogene_annotation
            ]
            args.extend(ppgene_args)
        OutputBedFilter(args, standalone_mode=False)

    def rename_query_genes(self) -> None:
        """Establishes query gene naming notation and renames gene entries in the final files"""
        from finalise_orthology_files import QueryGeneNamer
        args: List[str] = [
            self.orth_resolution_raw, self.query_genes_raw, self.finalized_output_dir,
            '-qb', self.query_genes_bed_raw, '-ln', self.project_name
        ]
        QueryGeneNamer(args, standalone_mode=False)
        for file in os.listdir(self.finalized_output_dir):
            filepath: str = os.path.join(self.finalized_output_dir, file)
            cmd: str = f'mv {filepath} {self.output}/'
            _ = self._exec(cmd, 'Moving final orthology file failed')

    def write_failed_batches(self) -> None:
        """
        Writes failed batches for feature extraction, 
        preprocessing, alignment steps to a respective metadata file
        """
        with open(self.failed_batches_file, 'w') as h:
            if self.failed_feature_batches:
                failed_feature_batches: str = ','.join(map(str, self.failed_feature_batches))
                h.write(f'FEATURE_EXTRACTION\t{failed_feature_batches}\n')
            if self.failed_preprocessing_batches:
                failed_preprocessing_batches: str = ','.join(map(str, self.failed_preprocessing_batches))
                h.write(f'PREPROCESSING\t{failed_preprocessing_batches}\n')
            if self.failed_alignment_batches:
                failed_alignment_batches: str = ','.join(map(str, self.failed_alignment_batches))
                h.write(f'ALIGNMENT\t{failed_alignment_batches}\n')
            if self.failed_orthology_batches:
                failed_orthology_batches: str = ','.join(map(str, self.failed_orthology_batches))
                h.write(f'ORTHOLOGY\t{failed_orthology_batches}\n')

    def cleanup(self) -> None:
        """Removes temporary directories, including Nextflow logs"""
        if not self.keep_tmp:
            self._to_log('Removing tmp directory %s' % self.tmp, 'info')
            self._rmdir(self.tmp)
            self._rmdir(self.vis_input_dir)
        if not self.keep_nextflow_log:
            self._to_log(
                'Removing Nextflow temporary directory %s' % self.nextflow_dir,
                'info'
            )
            self._rmdir(self.nextflow_dir)
            if self.parallel_strategy == 'para':
                self._to_log('Cleaning Para metadata', 'info')
                for project in self.parallel_process_names:
                    cleanup_cmd: str = f'para clean {project}'
                    _ = self._exec(
                        cleanup_cmd, f'ERROR: Cleanup failed for project {project}'
                    )

    def notify_on_completion(self, step: str) -> None:
        """Generates and sends an e-mail notification on successful TOGA2 run"""
        body: str = Constants.SUCCESS_EMAIL.format(self.project_name, self.output)
        if step != 'all':
            body += Constants.PARTIAL_RUN_NOTE.format(step)
        subject: str = Constants.SUCCESS_EMAIL_HEADER.format(self.project_name)
        self._email(subject, body)