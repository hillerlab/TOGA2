## v2.0.9
* `run` mode:
* `prepare-input` mode:
    * Annotation and isoform file generation from GFF3/GTF files;
* **NEW MODE**: `orthogroups` enables orthogroup modelling and generates gene family files compatible with CAFE5 input;

## v2.0.8
* `run` mode:
    * All eight CESAR2 profiles can be now provided as a single input directory with the `--cesar_profile_dir` argument.
    * `--paralogs_over_spanning` flag for swapping annotation priority
    * revised naming notation for one-to-many genes: addoitional copies for genes with more than 26 instances in the query get the binary letter suffix ('_aa', '_ab', etc.); genes with more than 300 orthologous copies lead to a hardcoded crash.
    * `--debug` mode (early access) for increased logging verbosity
    * Variable project argument formats (TSV, JSON, YAML).
* **NEW MODE**: `summary` for concise run summary generation.
* `from-config` mode:
    * Support for all accepted config formats (TSV, JSON, YAML).
* `spliceai` mode:
    * Lifting the "early access" warning (see `Minor changes`)
* `integrate` mode:
    * Paralogs overlapping orthologous projections are now retained as long as they contain enough novel sequence compared to the rest of the projections in the locus.
    * UTR-annotated input support
* `prepare-input` mode:
    * Exon sequence .2bit file for SLEASY compatibility generated along with BED and isoforms files.
    * File names are now prepended with an optional reference name prefix.
* Minor changes:
    * `run`:
        * Chimeric projections are no longer accounted for when estimated the most probable/most chain-covered items in `infer_query_genes.tsv`.
        * Transcripts with processed pseudogene projections only are now classified as *Missing* and appear in the rejection log under the `PPGENE_ONLY` label.
        * Genes with processed pseudogene projections only are also classified as *Missing*.
        * Safeguard against false loss status assignment based on the rejected items.
        * Additional loss summary updated at the `finalize` step (**NOTE**: the updated data are not reflected in `meta/loss_summary_extended.tsv`).
        * Stepwise rejection logs are now appended to `rejected_items.tsv` instead of being dumped to separate temporary files.
        * Fixing rejection level for `GENE_TREE_REJECTION` category from `TRANSCRIPT` to `PROJECTION`.
        * Default bootstrap number in `fine_ortology_resolver.py` set to 5000.
        * Resolving faulty imports from `shared.py` in scheduler scripts.
        * Temporary workarounds for conflicting paralogous/processed pseudogene projections from the rejection log in the gene loss summary module (`conservation_summary.py`).
        * Timestamps removed from the project names and moved to `projet_args.tsv` instead.
        * Modified protein sequence file (`protein.fa(.gz)`) generation
    * `spliceai`:
        * Overlapping coordinates bug fixed
        * Missing `project_name` arg fixed
    * `sequence-alignment`:
        * Error-free exit if no sequences were found across the query list for the focal transcript
        * PRANK is set as default sequence aligner
        * Fixed random seed option for PRANK
        * Query projection names added to FASTA headers if `--add_projection_names` flag if set
        * Added proper handling for exons present in one sequence only
        * Split codon phasing between the exons is disabled unless the `--phase_split_codons` flag is set
    * `prepare-input`:
        * trailing comma-insensitive parsing for BED fields 10 and 11

## v2.0.7a
* `run` mode:
    * Replacing positional arguments with keyword arguments
    * `--isoform_file`, `--u12_file`, and `--spliceai_dir` options are now "semi-mandatory"; the user is expected to provide the respective arguments unless the explicit deprecative flags are set
    * Alternative input formatting with `--input_directory`, `--ref_name`, and `--query_name` shortcuts: Format your data storage tree once and enjoy simplified command line interface
    * Postoga summary table (`toga.table.gz`) added to the output for `run` mode
    * Projections of the same reference gene/transcript overlapping by absolute coding sequence coordinates are now collapsed into a single query gene regardless of their overlap by coding exon coordinates
* **NEW MODE**: `postoga` for [Postoga](https://github.com/alejandrogzi/postoga) integration
* **NEW MODE**: `sequence-alignment` for orthologous sequence alignment across multiple same-referenced TOGA2 runs (alpha version)
* Apptainer support (see `supply/containers`):
    * Stable local execution container image
    * Batch manager-compatible image template
    * Removing `toga2.py` as a container entry point
    * Adding container support for parallel step scripts (see `supply/containers/README.md`)
* Updated local installation
    * Postoga installation
    * Conda environment support
    * Updated `bigWigToWig` version (`-bed` and `-header` options) now distributed with TOGA2
* Minor changes:
    * `run`:
        * Suppressed logging for XGBoost at `classification` step
        * Suppressed Pandas warnings at `classification` step 
        * Transcripts which do not have a single overlapping chain are now reported at `classification` step unless legacy feature extraction procedure is enabled
        * Setting default non-canonical U12 acceptor to `equiprobable_acceptor.tsv`
        * Proper non-canonical U12 donor profiel
        * Setting separate splice site treatment by default, replacing `--separate_splice_site_treatment` flag with `--joint_splice_site_treatment`
        * Fixed memory bin mem-to-jobs mapping for `alignment` step
        * Sequences in `nucleotide.fa` and `protein.fa` now contain only sequences from present (non-missing and non-deleted) exons
        * Minor gene inference speedup;
        * Overlapping projections of the same reference gene are now collapsed regardless of their overlap by chain-supported coding sequence -> reduced number of false one2many genes;
        * Replacing *N* loss status with *M* for certain rejected item categories
        * Removing the remaining 'orphan' projections after the gene tree step from the final annotation;
        * Consistent sequence selection and sorting for gene tree step input;
        * `rejected_items.tsv` is now filtered at rerunning to ensure that previous runs' results do not affect the resumed runs;
        * Proper handling for query gene rendered orphan after the gene tree step
        * Fixed random seed for PRANK
    * All modes:
        * Mandatory options introduced instead of positional arguments;
        * Grouped options for better readability

## v2.0.6
* New TOGA2 mode added: `integrate` (early access functionality)
* Query gene inference improved:
    * Gene inference step is moved from `loss_summary` to a separate workflow step now placed before gene loss summary step;
    * **All** query genes now get their names after progenitor genes in the reference, with the following prefixes:
        * intact orthologous loci do not get any prefix;
        * inactivated orthologous loci get `missing_` prefix if they were inferred from `Missing` projections alone; otherwise, the prefix is `lost_`;
        * paralogous loci get the `paralog_` prefix;
        * retrogene loci (annotated based on `Intact/Fully Intact` processed pseudogene projections) get the `retro_` prefix
* Updated loss summary files:
    * `loss_summary` step now goes after `gene_inference` to accommodate for transcript/gene status changes 
    after gene inference step;
    * top-level file `loss_summary.tsv` now contains loss statuses only for projections appearing in the final output files (`query_annotation.bed`, `query_annotation.with_utrs.bed`, UCSC browser files). It still contains loss statuses for **all** reference transcripts and genes. For loss statuses for all anotated projections, including rejected items see `meta/loss_summary_extended.tsv`
* '#paralog' postfix added for paralogous projections in the final output files (`query_annotation.bed`, `query_annotation.with_utrs.bed`, UCSC browser files)
* '#paralog' and "#retro" suffixes also added to respective projections' names in `nucleotide.fa.gz`, `protein.fa.gz`, `loss_summary.tsv`
* Entries for fragmented projections get numerical postfixes based on their order in the restored query sequence. Numbers are 1-based and are separated from the base name with dollar sign ('$')
* A GTF format copy is produced for final query annotation file (either `query_annotation.bed` or `query_annotation.with_utrs.bed`)
* Folder `nextflow_configs` contains example config files and executor script for parallel steps performed with Nextflow
* Bug fixes:
    * Third-party binaries are now sought in the `bin/` directory first.
    * Removed hardcoded instances of `bedToBigBed` and `ixIxx` in `src/python/modules/make_ucsc_report.py`
    * Error-exit if all batches for a given step failed prior to ok-file check
    * Projections discarded at gene tree filtering step now removed from the final output files
    * Added post-gene-tree orthology resolution step to the main logging channel
    * 'missing_' query gene inference;
    * Partially Intact consistently removed from accepted retrogene/trusted second-level ortholog statuses
