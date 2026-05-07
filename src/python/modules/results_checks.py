"""
Sanity check and troubleshooting manager
"""

import logging
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from .cesar_wrapper_constants import CLASS_TO_NUM
from .constants import Constants
from .shared import (
    CommandLineManager,
    get_proj2trans,
    get_upper_dir,
    parse_single_column,
    read_tab,
)

LOCATION: str = get_upper_dir(__file__, 4)
ORTH_PROB_HEADER: str = "transcript"
ORTH_PROB_FIELDS: int = 3
MEM_FILE_HEADER: str = "transcript"
ISOFORMS_FIELDS: int = 2
MEM_FILE_FIELDS: int = 15
REJ_LOG_FIELDS: int = 6
CHAIN_LIMIT_REASON: str = "CHAIN_LIMIT_EXCEEDED"
LOSS_HEADER: str = "level"
LOSS_FILE_FIELDS: int = 3
PROJECTION: str = "PROJECTION"
TRANSCRIPT: str = "TRANSCRIPT"
GENE: str = "GENE"
ORTH_CLASS_HEADER: str = "t_gene"
ORTH_CLASS_FIELDS: int = 5
IGNORED_STATUSES: Tuple[str, str, str] = ("PG", "PP", "N")
QUERY_ISOFORMS_HEADER: str = "query_gene"
REG: str = "pseudo_"
RETRO: str = "retro_"
NONE: str = "None"

MIN_ORTH_PERCENTAGE: float = 50.0
MAX_REJ_AT_PREPR_PERCENTAGE: float = (1 / 3) * 100
MIN_INTACT_PROJ_PERCENTAGE: float = 50.0
MIN_ONE2ONE_PERCENTAGE: float = 50.0
MAX_ONE2ZERO_PERCENTAGE: float = 25.0

ONE2ONE: str = "one2one"
ONE2MANY: str = "one2many"
MANY2ONE: str = "many2one"
MANY2MANY: str = "many2many"
ONE2ZERO: str = "one2zero"

EXPECTED_OUTPUT_FILE_NAMES: Dict[str, str] = {
    "orth_probs_file": "orthology_scores.tsv",
    "loss_summary": "loss_summary.tsv",
    "query_genes": "query_genes.tsv",
    "orthology_classification": "orthology_classification.tsv",
}

REF_2BIT: str = "ref_2bit"
QUERY_2BIT: str = "query_2bit"
CHAIN_FILE: str = "chain_file"
REF_ANNOTATION: str = "ref_annotation"
ORTH_THRESH: str = "orthology_threshold"
ACCEPTED_CLASSES: str = "accepted_loss_symbols"
ISOFORM_FILE: str = "isoform_file"
OUT_DIR: str = "output"
U12_FILE: str = "u12_file"
SPLICEAI_DIR: str = "spliceai_dir"
QUERY_NAME: str = "query_name"
MANDATORY_PATH_ARGS: Tuple[str, ...] = (
    REF_2BIT,
    QUERY_2BIT,
    CHAIN_FILE,
    REF_ANNOTATION,
    OUT_DIR,
)
MANDATORY_ARGS: Tuple[str, ...] = (*MANDATORY_PATH_ARGS, ORTH_THRESH, ACCEPTED_CLASSES,)
ALL_ARGS: Tuple[str, ...] = (
    *MANDATORY_ARGS,
    ISOFORM_FILE,
    U12_FILE,
    SPLICEAI_DIR,
    QUERY_NAME,
)
OUTPUT: str = "output"

BREAK_LINE: str = "#" * 100

MAIN_HEADER: str = """{br}
#### TOGA2 summary
{br}
TOGA2 annotated {num_genes} genes and {num_retro} retrogene candidates in the query genome.
In addition, TOGA2 identified {num_lost} genes classified as lost and {num_missing} genes classified as missing in the query.

Out of {ref_gene_num} reference genes, {num_with_func_orth} ({perc_with_func_orth}%) have at least one potentially functional ortholog; of these, {num_one2one} ({perc_one2one}%) are classified as 1:1 orthologs. For {num_with_func_para} ({perc_with_func_para}%) genes, TOGA2 identified only potentially functional paralogs in the query genome.

#HEADER\tQueryAssembly\tno. annotated query genes  \tno. query retrogenes\tno. lost genes in query\tno. missing genes in query\tno. ref genes with potentially functional orthologs\tno. one2one orthologs\tno. ref genes with potentially functional paralogs
#SINGLELINESUMMARY\t{query_name}\t{num_genes}\t{num_retro}\t{num_lost}\t{num_missing}\t{num_with_func_orth}\t{num_one2one}\t{num_with_func_para}

This data set was generated with TOGA2 version {version} and the following input files:
* Reference genome: {ref_2bit}
* Query genome: {query_2bit}
* Genome alignment chain file: {chains}
* Reference annotation file: {ref_annot}
* Output directory: {output}
* Isoforms file: {isoforms}
* U12/non-canonical U2 intron file: {u12_file}
* SpliceAI query directory: {spliceai_dir}{cmd}

For questions, please check whether the issue has already been addressed at https://github.com/hillerlab/TOGA2/issues. If not, please open a new issue.

If you use these data, please cite: 
Yury V. Malovichko et al. "Accurate, comprehensive gene annotation and ortholog identification across thousands of vertebrate genomes with TOGA2", in preparation

"""

CMD: str = """

TOGA2 was invoked with the following command:
{cmd}
"""

DETAILED_HEADER: str = """{br}
#### Detailed TOGA2 summary
{br}
The detailed summary lists information about the orthology classification step, statistics of projections, transcripts and genes,
and details on the orthology classification.{warnings}"""


_SUMMARY_BOILERPLATE: str = """{header}

{detailed_header}

{br}

{orthology_class_report}

{br}

{loss_summary}

{br}

{orthology_res_report}
"""


SUMMARY_BOILERPLATE: str = """
{br}
#### TOGA2 run summary
{br}

Reference genome: {ref_2bit}
Query genome: {query_2bit}
Genome alignment chain file: {chains}
Reference annotation file: {ref_annot}
Output directory: {output}

{isoforms_line}

{br}

{orthology_class_report}

{br}

{loss_summary}

{br}

{orthology_res_report}

"""

ISOFORMS_LINE: str = "Isoforms file: {}\n"
ISOFORMS_CAVEAT: str = """
\nNOTE: Reference isoform file was not provided. All gene level statistics are presented \
assuming that each reference transcript represents a separate gene
"""

CLASSIFICATION_BOILERPLATE: str = """
Orthology prediction statistics:
\tReference transcripts level::
\t\t#reference transcripts subjected to classification: {}
\t\t#reference transcripts with >=1 predicted ortholog: {} ({}%)
\t\t#reference transcripts with 1 predicted ortholog: {} ({}%)
\t\t#reference transcripts with no predicted orthologs: {} ({}%)
\t\t#reference transcripts with no classifiable projections: {} ({}%)
\t\t#predicted processed pseudogene/retrogene projections {}:
\tReference genes level:
\t\t#reference genes with >=1 predicted ortholog: {} ({}%)
\t\t#reference genes with 1 predicted ortholog: {} ({}%)
\t\t#reference genes with no predicted orthologs: {} ({}%)
\t\t#reference gene with no classifiable projections: {} ({}%)
"""
## TODO: Consider adding IGNORED_STATUSES to the summary to exclude ambiguity in the results' interpretation
LOSS_SUMMARY_BOILERPLATE: str = """
Gene loss summary statistics:
\tloss classes considered for assessing gene presence: {}
\tQuery projections level:
\t\t#Fully Intact (FI) - {} ({}%)
\t\t#Intact (I) - {} ({}%)
\t\t#Partially Intact (PI) - {} ({}%)
\t\t#Uncertain Losses (UL) - {} ({}%)
\t\t#Lost (L) - {} ({}%)
\t\t#Missing (M) - {} ({}%)
\t\t#projections considered present: {} ({}%)
\t\t#projections considered lost/missing: {} ({}%)
\tReference transcripts level:
\t\t#Fully Intact (FI) - {} ({}%)
\t\t#Intact (I) - {} ({}%)
\t\t#Partially Intact (PI) - {} ({}%)
\t\t#Uncertain Losses (UL) - {} ({}%)
\t\t#Lost (L) - {} ({}%)
\t\t#Missing (M) - {} ({}%)
\t\t#trancripts considered present: {} ({}%)
\t\t#transcripts considered lost/missing: {} ({}%)
\tReference genes level:
\t\t#Fully Intact (FI) - {} ({}%)
\t\t#Intact (I) - {} ({}%)
\t\t#Partially Intact (PI) - {} ({}%)
\t\t#Uncertain Losses (UL) - {} ({}%)
\t\t#Lost (L) - {} ({}%)
\t\t#Missing (M) - {} ({}%)
\t\t#reference genes considered present: {} ({}%)
\t\t#reference genes considered lost/missing: {} ({}%)
"""

ORTHOLOGY_BOILERPLATE: str = """
Orthology resolution:
\t#reference genes: {}
\t#query genes: {}
\t\t#with defined orthology: {} ({}%)
\t\t#lost, missing, or lacking defined orthology: {} ({}%)
\tReference gene orthology class composition:
\t\tone2one: {} ({}%)
\t\tone2many: {} ({}%)
\t\tmany2one: {} ({}%)
\t\tmany2many: {} ({}%)
\t\tone2zero: {} ({}%)
"""

CLASSIFICATION_STUB: str = """\
{num} reference {lvl} ({fraction}%) do not have any predicted orthologous projections\
"""
CLASSIFICATION_OK: str = "Sanity check passed: " + CLASSIFICATION_STUB
CLASSIFICATION_WARNING: str = (
    "Sanity check failed: "
    + CLASSIFICATION_STUB
    + ", indicating potential problems with the current TOGA2 run"
)
PREPROCESSING_STUB: str = """\
Of {} projections having reached the preprocessing step, {} projections ({}%) were discarded\
"""
PREPROCESSING_OK: str = "Sanity check passed: " + PREPROCESSING_STUB
PREPROCESSING_WARNING: str = (
    "Sanity check failed: "
    + PREPROCESSING_STUB
    + ", indicating potential problems with the current TOGA2 run"
)
LOSS_STUB: str = """\
Of {total} orthologous {lvl} having reached the gene loss summary step, \
{num_valid} {lvl} ({fraction}%) belong to loss classes \
designating gene presence ({classes})\
"""
LOSS_OK: str = "Sanity check passed: " + LOSS_STUB
LOSS_WARNING: str = (
    "Sanity check failed: "
    + LOSS_STUB
    + ", indicating potential problems with the current TOGA2 run"
)
RESOLUTION_STUB: str = """\
Of {total} reference genes, {num} ({fraction}%) genes have {type} orthology status\
"""
RESOLUTION_OK: str = "Sanity check passed: " + RESOLUTION_STUB
RESOLUTION_WARNING: str = (
    "Sanity check failed: "
    + RESOLUTION_STUB
    + ", indicating potential problems with the current TOGA2 run"
)


@dataclass
class ElementaryCheckResult:
    is_warning: bool
    log_message: str
    warning_message: str


@dataclass
class SanityCheckResult:
    step: str
    messages: List[ElementaryCheckResult]


def to_perc(numerator: int, denominator: int) -> float:
    """
    Divides numerator by denominator, multiplies by 100,
    and rounds to third digit after the dot
    """
    return 0.0 if denominator == 0 else round(numerator / denominator * 100, 3)


def arg_or_na(arg: Union[Any, None]) -> str:
    """Returns the input argument or N/A if argument was not provided

    Args:
        arg: A TOGA2 argument to check

    Returns:
        'N/A' if argument value is None or an empty string, the initial value otherwise
    """
    if isinstance(arg, (int, float, bool, complex)):
        return str(arg)
    if not arg:
        return "N/A"
    return str(arg)


class ResultChecker(CommandLineManager):
    """
    Checks output at selected TOGA2 pipeline steps;
    at the end of the pipeline, summarizes results and
    """

    __slots__ = (
        "logger",
        "orth_probs",
        "feat_rej_log",
        "class_rej_log",
        "memory_report",
        "prepr_rej_log",
        "loss_summary",
        "orth_classes",
        "query_genes_file",
        "isoforms",
        "paralogs",
        "ppgenes",
        "orth_thresh",
        "intact_classes",
        "all_initial_transcripts",
        "tr2orths",
        "gene2tr",
        "non_orthologous",
        "num_orth",
        "fr_orth",
        "num_o2o",
        "fr_o2o",
        "num_no_orth",
        "fr_no_orth",
        "num_no_chains",
        "fr_no_chains",
        "num_ppgene",
        "gene_num_orth",
        "gene_fr_orth",
        "gene_num_o2o",
        "gene_fr_o2o",
        "gene_num_no_orth",
        "gene_fr_no_orth",
        "proj2loss",
        "tr2loss",
        "gene2loss",
        "orth_status_counter",
        "gene2orth_class",
    )

    def __init__(
        self,
        logger: logging.Logger,
        output: str,
        ref_isoforms: Optional[Union[str, None]] = None,
        orthology_threshold: Optional[float] = Constants.DEFAULT_ORTH_THRESHOLD,
        intact_classes: Optional[Union[List, str]] = Constants.DEFAULT_LOSS_SYMBOLS,
    ) -> None:
        pass
        ## TODO: If pipeline was resumed from a certain step,
        ## do not forget to check the results
        self.logger: logging.Logger = logger
        self.orth_probs: str = os.path.join(output, "orthology_scores.tsv")
        self.feat_rej_log: str = os.path.join(
            output, "tmp", "rejection_logs", "rejected_at_feature_extraction.tsv"
        )
        self.class_rej_log: str = os.path.join(
            output, "tmp", "rejection_logs", "rejected_at_classification.tsv"
        )
        self.memory_report: str = os.path.join(
            output, "meta", "memory_requirements.tsv"
        )
        self.prepr_rej_log: str = os.path.join(
            output, "tmp", "rejection_logs", "rejected_at_preprocessing.tsv"
        )
        self.loss_summary: str = os.path.join(
            output, "meta", "loss_summary_extended.tsv"
        )
        self.orth_classes: str = os.path.join(
            output, "tmp", "annotation_raw", "orthology_classification.tsv"
        )
        self.query_genes_file: str = os.path.join(output, "query_genes.tsv")
        self.isoforms: Union[str, None] = ref_isoforms
        self.paralogs: str = os.path.join(
            output, "meta", "paralogous_projections_to_align.tsv"
        )
        self.ppgenes: str = os.path.join(
            output, "meta", "processed_pseudogene_projections_to_align.tsv"
        )
        self.orth_thresh: float = orthology_threshold
        if isinstance(intact_classes, str):
            if intact_classes == "ALL":
                self.intact_classes: List[str] = Constants.ALL_LOSS_SYMBOLS
            else:
                self.intact_classes: List[str] = [
                    x.strip() for x in intact_classes.split(",") if x
                ]
        else:
            self.intact_classes: List[str] = intact_classes

        self.all_initial_transcripts: Set[str] = set()
        self.tr2orths: Dict[str, List[str]] = defaultdict(list)
        self.gene2tr: Dict[str, List[str]] = defaultdict(list)

        self.non_orthologous: Set[str] = set()

        self.num_orth: int = 0
        self.fr_orth: float = 0.0
        self.num_o2o: int = 0
        self.fr_o2o: float = 0.0
        self.num_no_orth: int = 0
        self.fr_no_orth: float = 0.0
        self.num_no_chains: int = 0
        self.fr_no_chains: float = 0.0
        self.num_ppgene: int = 0
        self.gene_num_orth: int = 0
        self.gene_fr_orth: float = 0.0
        self.gene_num_o2o: int = 0
        self.gene_fr_o2o: float = 0.0
        self.gene_num_no_orth: int = 0
        self.gene_fr_no_orth: float = 0.0

        self.proj2loss: Dict[str, int] = defaultdict(int)
        self.tr2loss: Dict[str, int] = defaultdict(int)
        self.gene2loss: Dict[str, int] = defaultdict(int)

        self.orth_status_counter: Dict[str, int] = defaultdict(int)

        self.gene2orth_class: Dict[str, int] = defaultdict(int)

    def _throw_warning(
        self, text: str, to_email: Optional[bool] = False
    ) -> Union[str, None]:
        """Throw warning to the provided Logger channel"""
        pass
        if to_email:
            pass  ## return warning text

    def _parse_isoform_file(self) -> None:
        """Retrieves gene-to-transcripts mapping for reference genome annotation"""
        with open(self.isoforms, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) != ISOFORMS_FIELDS:
                    self.logger.critical(
                        (
                            "Improper formatting at reference isoform file line %i; "
                            "expected %i columns, got %i"
                        )
                        % (i, ISOFORMS_FIELDS, len(data))
                    )
                    sys.exit(1)
                gene, tr = data
                self.gene2tr[gene].append(tr)

    def _parse_non_orthologous(self) -> None:
        """Parses paralogous and retrocopy projection files"""
        if os.path.exists(self.paralogs):
            self.non_orthologous = self.non_orthologous.union(
                parse_single_column(self.paralogs)
            )
        if os.path.exists(self.ppgenes):
            self.non_orthologous = self.non_orthologous.union(
                parse_single_column(self.ppgenes)
            )

    def _get_query_gene_num(self) -> Tuple[int, int]:
        """Counts the number of query genes"""
        all_query_genes: Set[str] = set()
        if not os.path.exists(self.query_genes_file):
            self.logger.critical(
                "Reference genes file %s does not exist" % self.query_genes_file
            )
            sys.exit(1)
        with open(self.query_genes_file, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) != ISOFORMS_FIELDS:
                    self.logger.critical(
                        (
                            "Improper formatting at query isoform file line %i; "
                            "expected %i columns, got %i"
                        )
                        % (i, ISOFORMS_FIELDS, len(data))
                    )
                    sys.exit(1)
                if data[0] == QUERY_ISOFORMS_HEADER:
                    continue
                if data[0].startswith(RETRO):
                    continue
                all_query_genes.add(data[0])

        orphan: int = sum(x.startswith(REG) for x in all_query_genes)
        orth: int = len(all_query_genes) - orphan

        return (orth, orphan)

    def check_classification(self) -> SanityCheckResult:
        """
        Checks orthology prediction results. Throws a warning if:
        a) <50% of the transcripts do not have at least one ortholog;
        b) <50% reference genes do not have any isoform with at least one ortholog
        """
        if not os.path.exists(self.orth_probs):
            self.logger.critical(
                "Orthology score file %s does not exist" % self.orth_probs
            )
            sys.exit(1)
        messages: List[ElementaryCheckResult] = []
        with open(self.orth_probs, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) != ORTH_PROB_FIELDS:
                    self.logger.critical(
                        (
                            "Improper formatting at orthology probability file line %i; "
                            "expected %i columns, got %i"
                        )
                        % (i, ORTH_PROB_FIELDS, len(data)),
                    )
                    sys.exit(1)
                if data[0] == ORTH_PROB_HEADER:
                    continue
                tr: str = data[0]
                self.all_initial_transcripts.add(tr)
                chain: str = data[1]
                prob: float = float(data[2])
                proj: str = f"{tr}#{chain}"
                if prob >= self.orth_thresh:
                    self.tr2orths[tr].append(proj)
                if prob == Constants.PPGENE_PROB:
                    self.num_ppgene += 1

        discarded_at_classification: Set[str] = set()
        if os.path.exists(self.feat_rej_log):
            with open(self.feat_rej_log, "r") as h:
                for i, line in enumerate(h, start=1):
                    data: List[str] = line.strip().split("\t")
                    if not data:
                        continue
                    if len(data) != REJ_LOG_FIELDS:
                        self.logger.critical(
                            (
                                "Improper formatting at feature extraction rejection log file line %i; "
                                "expected %i columns, got %i"
                            )
                            % (i, REJ_LOG_FIELDS, len(data)),
                        )
                        sys.exit(1)
                    if data[0] != TRANSCRIPT:
                        continue
                    discarded_at_classification.add(data[1])
        if os.path.exists(self.class_rej_log):
            with open(self.class_rej_log, "r") as h:
                for i, line in enumerate(h, start=1):
                    data: List[str] = line.strip().split("\t")
                    if not data:
                        continue
                    if len(data) != REJ_LOG_FIELDS:
                        self.logger.critical(
                            (
                                "Improper formatting at classification rejection log file line %i; "
                                "expected %i columns, got %i"
                            )
                            % (i, REJ_LOG_FIELDS, len(data)),
                        )
                        sys.exit(1)
                    if data[0] != TRANSCRIPT:
                        continue
                    discarded_at_classification.add(data[1])

        self.num_no_chains = len(discarded_at_classification)
        num_tr: str = len(self.all_initial_transcripts) + self.num_no_chains
        self.num_orth = len(self.tr2orths)
        self.fr_orth = to_perc(self.num_orth, num_tr)
        self.num_o2o = sum(len(x) == 1 for x in self.tr2orths.values())
        self.fr_o2o = to_perc(self.num_o2o, num_tr)
        self.num_no_orth = sum(
            x not in self.tr2orths for x in self.all_initial_transcripts
        )
        self.fr_no_orth = to_perc(self.num_no_orth, num_tr)
        self.fr_no_chains = to_perc(self.num_no_chains, num_tr)

        ## throw a warning if <50% of transcrips do not have an ortholog
        if self.fr_no_orth > MIN_ORTH_PERCENTAGE:
            warning_message: str = CLASSIFICATION_WARNING.format(
                num=self.num_no_orth, lvl="transcripts", fraction=self.fr_no_orth
            )
            log_message = warning_message
            raise_warning = True
        else:
            warning_message: str = ""
            log_message: str = CLASSIFICATION_OK.format(
                num=self.num_no_orth, lvl="transcripts", fraction=self.fr_no_orth
            )
            raise_warning = False
        elem_trs: ElementaryCheckResult = ElementaryCheckResult(
            is_warning=raise_warning,
            warning_message=warning_message,
            log_message=log_message,
        )
        messages.append(elem_trs)

        ## parse isoforms file if present
        if self.isoforms is not None:
            self._parse_isoform_file()
            num_gene: int = sum(
                any(y in self.all_initial_transcripts for y in x) for x in self.gene2tr
            )
            self.gene_num_orth: int = sum(
                any(y in self.tr2orths for y in x) for x in self.gene2tr.values()
            )
            self.gene_fr_orth = to_perc(self.gene_num_orth, num_gene)
            self.gene_num_o2o = sum(
                all(self.tr2orths[y] == 1 for y in x) for x in self.gene2tr.values()
            )
            self.gene_fr_o2o = to_perc(self.gene_num_o2o, num_gene)
            ## throw a warning if <50% genes have an ortholog
            if self.gene_fr_no_orth > MIN_ORTH_PERCENTAGE:
                warning_message: str = CLASSIFICATION_WARNING.format(
                    num=self.gene_num_o2o, lvl="genes", fraction=self.gene_fr_o2o
                )
                log_message: str = warning_message
                raise_warning: bool = True
            else:
                warning_message: str = ""
                log_message += "\n" + CLASSIFICATION_OK.format(
                    num=self.gene_num_o2o, lvl="genes", fraction=self.gene_fr_o2o
                )
                raise_warning: bool = False
            elem_genes: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=raise_warning,
                warning_message=warning_message,
                log_message=log_message,
            )
            messages.append(elem_genes)

        result: SanityCheckResult = SanityCheckResult(
            step="classification", messages=messages
        )
        return result

    def check_preprocessing(self) -> SanityCheckResult:
        """
        Checks preprocessing rejection log.
        Throws a warning if >33% of projections having reached this point
        were discarded before alignment.
        """
        ## get the list of transcripts that reached this point
        if not os.path.exists(self.memory_report):
            self.logger.critical(
                "Memory report file %s does not exist" % self.memory_report
            )
            sys.exit(1)
        if not os.path.exists(self.prepr_rej_log):
            self.logger.warning(
                (
                    "Preprocessing step rejection log file %s does not exist; "
                    "assuming no items were discarded at the preprocessing step"
                )
                % self.prepr_rej_log
            )
            result: SanityCheckResult = SanityCheckResult(
                step="preprocessing",
                messages=[
                    ElementaryCheckResult(
                        is_warning=False,
                        warning_message="",
                        log_message="No items were discarded at preprocessing step",
                    )
                ],
            )
            return result

        accepted_num: int = 0
        with open(self.memory_report, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                if data[0] == MEM_FILE_HEADER:
                    continue
                if len(data) != MEM_FILE_FIELDS:
                    self.logger.critical(
                        (
                            "Improper formatting at memory file line %s; "
                            "expected %i fields, got %s"
                        )
                        % (i, MEM_FILE_FIELDS, len(data))
                    )
                    sys.exit(1)
                accepted_num += 1
        ## get the list of transcripts discarded at this step
        total_rej_num: int = 0
        rej_num: int = 0
        with open(self.prepr_rej_log, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) != REJ_LOG_FIELDS:
                    self.logger.critical(
                        (
                            "Improper formatting at preprocessing step rejection log file line %s; "
                            "expected %i fields, got %s"
                        )
                        % (i, REJ_LOG_FIELDS, len(data))
                    )
                total_rej_num += 1
                if data[4] == CHAIN_LIMIT_REASON:
                    continue
                rej_num += 1
        total_num: int = accepted_num + rej_num
        rej_perc: float = to_perc(rej_num, total_num)
        if rej_perc > MAX_REJ_AT_PREPR_PERCENTAGE:
            warning_message = PREPROCESSING_WARNING.format(total_num, rej_num, rej_perc)
            log_message = warning_message
            raise_warning = True
        else:
            warning_message: str = ""
            log_message = PREPROCESSING_OK.format(total_num, rej_num, rej_perc)
            raise_warning = False

        result: SanityCheckResult = SanityCheckResult(
            step="preprocessing",
            messages=[
                ElementaryCheckResult(
                    is_warning=raise_warning,
                    warning_message=warning_message,
                    log_message=log_message,
                )
            ],
        )

        return result

    def check_loss_summary(self) -> None:
        """
        Checks gene loss summary. Throws a warning if:
        a) <50% of original projections are considered present;
        b) <50% of reference transcripts are considered present;
        c) <50% of reference genes are considered present.
        """
        self._parse_non_orthologous()

        if not os.path.exists(self.loss_summary):
            self.logger.critical(
                "Gene loss summary file %s does not exist" % self.loss_summary
            )
            sys.exit(1)
        messages: List[ElementaryCheckResult] = []
        with open(self.loss_summary, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) != LOSS_FILE_FIELDS:
                    self.logger.critical(
                        (
                            "Improper formatting at loss summary file line %i; "
                            "expected %i fields, got %i"
                        )
                        % (i, LOSS_FILE_FIELDS, len(data))
                    )
                    sys.exit(1)
                if data[0] == LOSS_HEADER:
                    continue
                level, name, status = data
                if level == PROJECTION:
                    # self.proj2loss[status] += 1
                    if name in self.non_orthologous:
                        ## ignore Ns, PPgenes, and paralogs for the warning message
                        continue
                    if status == "N":
                        continue
                    self.proj2loss[status] += 1
                if level == TRANSCRIPT:
                    if status in IGNORED_STATUSES:
                        continue
                    self.tr2loss[status] += 1
                if level == GENE:
                    if status in IGNORED_STATUSES:
                        continue
                    self.gene2loss[status] += 1
        ## throw a warning for each loss summary level
        ## 1) projections
        num_proj: int = sum(self.proj2loss.values())
        num_acc_proj: int = sum(
            y for x, y in self.proj2loss.items() if x in self.intact_classes
        )
        fr_acc_proj: float = to_perc(num_acc_proj, num_proj)
        if fr_acc_proj < MIN_INTACT_PROJ_PERCENTAGE:
            warning_message: str = LOSS_WARNING.format(
                lvl="projections",
                total=num_proj,
                num_valid=num_acc_proj,
                fraction=fr_acc_proj,
                classes=",".join(self.intact_classes),
            )
            log_message = warning_message
            elem_proj: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True,
                warning_message=warning_message,
                log_message=log_message,
            )
        else:
            warning_message: str = ""
            log_message = LOSS_OK.format(
                lvl="projections",
                total=num_proj,
                num_valid=num_acc_proj,
                fraction=fr_acc_proj,
                classes=",".join(self.intact_classes),
            )
            elem_proj: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False,
                warning_message=warning_message,
                log_message=log_message,
            )
        messages.append(elem_proj)
        ## 2) reference transcripts
        num_tr: int = sum(self.tr2loss.values())
        num_acc_tr: int = sum(
            y for x, y in self.tr2loss.items() if x in self.intact_classes
        )
        fr_acc_tr: float = to_perc(num_acc_tr, num_tr)
        if fr_acc_tr < MIN_INTACT_PROJ_PERCENTAGE:
            warning_message: str = LOSS_WARNING.format(
                lvl="reference transcripts",
                total=num_tr,
                num_valid=num_acc_tr,
                fraction=fr_acc_tr,
                classes=",".join(self.intact_classes),
            )
            log_message = warning_message
            elem_tr: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True,
                warning_message=warning_message,
                log_message=log_message,
            )
        else:
            warning_message: str = ""
            log_message = LOSS_OK.format(
                lvl="reference transcripts",
                total=num_tr,
                num_valid=num_acc_tr,
                fraction=fr_acc_tr,
                classes=",".join(self.intact_classes),
            )
            elem_tr: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False,
                warning_message=warning_message,
                log_message=log_message,
            )
        messages.append(elem_tr)
        ## 3) genes (if isoform file was provided)
        if self.isoforms is not None:
            num_gene: int = sum(self.gene2loss.values())
            num_acc_gene: int = sum(
                y for x, y in self.gene2loss.items() if x in self.intact_classes
            )
            fr_acc_gene: int = to_perc(num_acc_gene, num_gene)
            if fr_acc_tr < MIN_INTACT_PROJ_PERCENTAGE:
                warning_message: str = LOSS_WARNING.format(
                    lvl="reference genes",
                    total=num_gene,
                    num_valid=num_acc_gene,
                    fraction=fr_acc_gene,
                    classes=",".join(self.intact_classes),
                )
                log_message = warning_message
                elem_gene: ElementaryCheckResult = ElementaryCheckResult(
                    is_warning=True,
                    warning_message=warning_message,
                    log_message=log_message,
                )
            else:
                warning_message: str = ""
                log_message = LOSS_OK.format(
                    lvl="reference genes",
                    total=num_gene,
                    num_valid=num_acc_gene,
                    fraction=fr_acc_gene,
                    classes=",".join(self.intact_classes),
                )
                elem_gene: ElementaryCheckResult = ElementaryCheckResult(
                    is_warning=False,
                    warning_message=warning_message,
                    log_message=log_message,
                )
            messages.append(elem_gene)

        result: SanityCheckResult = SanityCheckResult(
            step="gene loss summary", messages=messages
        )
        return result

    def check_orthology_resolution(self) -> None:
        """
        Checks orthology resolution results. Throws a warning if:
        a) <50% of genes are classified as one:one;
        b) >=20% of genes are classified as one:zero
        """
        if not os.path.exists(self.orth_classes):
            self.logger.critical(
                "Orthology summary file %s does not exist" % self.orth_classes
            )
            sys.exit(1)
        with open(self.orth_classes, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                if data[0] == ORTH_CLASS_HEADER:
                    continue
                if len(data) != ORTH_CLASS_FIELDS:
                    self.logger.critical(
                        (
                            "Improper formatting at orthology summary file line %i; "
                            "expected %i fields, got %i"
                        )
                        % (i, ORTH_CLASS_FIELDS, len(data))
                    )
                # ref_gene: str = data[0]
                # query_gene: str = data[2]
                status: str = data[4]
                self.orth_status_counter[status] += 1
                ## TODO: For sanity check, do not account for genes discarded at the initial steps
                ## (technical reasons, etc.)
                ## NOTE: those with no orthologous predictions at the initial step should still be accounted for!
        messages: List[ElementaryCheckResult] = []
        num_gene: int = sum(self.orth_status_counter.values())
        num_one2one: int = self.orth_status_counter.get("one2one", 0)
        fr_one2one: float = to_perc(num_one2one, num_gene)
        if fr_one2one < MIN_ONE2ONE_PERCENTAGE:
            warning_message: str = RESOLUTION_WARNING.format(
                total=num_gene, num=num_one2one, fraction=fr_one2one, type="1:1"
            )
            log_message: str = warning_message
            elem_one2one: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True,
                warning_message=warning_message,
                log_message=log_message,
            )
        else:
            warning_message: str = ""
            log_message: str = RESOLUTION_OK.format(
                total=num_gene, num=num_one2one, fraction=fr_one2one, type="1:1"
            )
            elem_one2one: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False,
                warning_message=warning_message,
                log_message=log_message,
            )
        messages.append(elem_one2one)
        num_one2zero: int = self.orth_status_counter.get("one2zero", 0)
        fr_one2zero: float = to_perc(num_one2zero, num_gene)
        if fr_one2zero > MAX_ONE2ZERO_PERCENTAGE:
            warning_message: str = RESOLUTION_WARNING.format(
                total=num_gene, num=num_one2zero, fraction=fr_one2zero, type="1:0"
            )
            log_message: str = warning_message
            elem_one2one: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True,
                warning_message=warning_message,
                log_message=log_message,
            )
        else:
            warning_message: str = ""
            log_message: str = RESOLUTION_OK.format(
                total=num_gene, num=num_one2zero, fraction=fr_one2zero, type="1:0"
            )
            elem_one2zero: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False,
                warning_message=warning_message,
                log_message=log_message,
            )
            messages.append(elem_one2zero)

        result: SanityCheckResult = SanityCheckResult(
            step="orthology resolution", messages=messages
        )
        return result


class LogParserForSummary(CommandLineManager):
    """
    An auxiliary project arguments file parser for summary report production
    """

    __slots__ = ("args_file", "format")

    def __init__(self, args: os.PathLike, report_format: str) -> None:
        self.v: bool = True
        self.set_logging()

        self.args_file: str = args
        self.format: str = report_format

    def extract_settings(self) -> Dict[str, Union[str, None]]:
        """Main args extraction method"""
        if self.format == "tsv":
            args: Dict[str, Union[str, None]] = self.parse_tsv()
        elif self.format == "yaml":
            args: Dict[str, Union[str, int, float, None]] = self.parse_yaml()
        elif self.format == "json":
            args: Dict[str, Union[str, None]] = self.parse_json()
        else:
            self._die("Inappropriate log file format")

        return args

    def parse_tsv(self) -> Dict[str, Union[str, None]]:
        """Extracts args from the TSV-format (legacy) project arguments file"""
        args: Dict[str, Union[str, None]] = {arg: None for arg in ALL_ARGS}
        for data in read_tab(self.args_file):
            if data[0] == REF_2BIT:
                if not os.path.exists(data[1]):
                    self._to_log(
                        (
                            "Reference 2bit file %s does not exist or "
                            "has been deleted since the run passed"
                        )
                        % data[1],
                        "warning",
                    )
                args[REF_2BIT] = data[1]
            if data[0] == QUERY_2BIT:
                if not os.path.exists(data[1]):
                    self._to_log(
                        (
                            "Query 2bit file %s does not exist or "
                            "has been deleted since the run passed"
                        )
                        % data[1],
                        "warning",
                    )
                args[QUERY_2BIT] = data[1]
            if data[0] == CHAIN_FILE:
                if not os.path.exists(data[1]):
                    self._to_log(
                        (
                            "Chain file %s does not exist or "
                            "has been deleted since the run passed"
                        )
                        % data[1],
                        "warning",
                    )
                args[CHAIN_FILE] = data[1]
            if data[0] == REF_ANNOTATION:
                if not os.path.exists(data[1]):
                    self._die(
                        (
                            "Reference annotation file %s does not exist or "
                            "has been deleted since the run passed"
                        )
                        % data[1]
                    )
                args[REF_ANNOTATION] = data[1]
            if data[0] == ISOFORM_FILE:
                if data[1] == NONE:
                    continue
                elif not os.path.exists(data[1]):
                    self._die(
                        (
                            "Reference isoforms file %s does not exist or "
                            "has been deleted since the run passed"
                        )
                        % data[1]
                    )
                args[ISOFORM_FILE] = data[1]
            if data[0] == U12_FILE:
                if data[1] == NONE:
                    continue
                if not os.path.exists(data[1]):
                    self._die(
                        (
                            "Reference U12 intron file %s does not exist or "
                            "has been deleted since the run passed"
                        ) % data[1]
                    )
                args[U12_FILE] = data[1]
            if data[0] == SPLICEAI_DIR:
                if data[1] == NONE:
                    continue
                if not os.path.exists(data[1]):
                    self._die(
                        (
                            "Query SpliceAI annotation directory %s does not exist or "
                            "has been deleted since the run passed"
                        ) % data[1]
                    )
                args[SPLICEAI_DIR] = data[1]
            if data[0] == ORTH_THRESH:
                try:
                    args[ORTH_THRESH] = float(data[1])
                except ValueError:
                    self._die(
                        "Orthology probability threshold is not a valid floating point number"
                        % data[1]
                    )
            if data[0] == ACCEPTED_CLASSES:
                accepted_loss_symbols: List[str] = [x for x in data[1].split(",") if x]
                if any(x not in Constants.ALL_LOSS_SYMBOLS for x in accepted_loss_symbols):
                    self._die(
                        '"accepted_loss_symbols" contains inappropriate loss symbols: %s'
                        % ", ".join(
                            [
                                x
                                for x in accepted_loss_symbols
                                if x not in Constants.ALL_LOSS_SYMBOLS
                            ]
                        )
                    )
                args[ACCEPTED_CLASSES] = data[1]
            if data[0] == OUT_DIR:
                if not os.path.exists(data[1]):
                    self._die(
                        (
                            "Output directory %s does not exist or "
                            "has been deleted since the run passed"
                        )
                        % data[1]
                    )
                args[OUT_DIR] = data[1]
                for attr, expected_name in EXPECTED_OUTPUT_FILE_NAMES.items():
                    expected_path: str = os.path.join(data[1], expected_name)
                    if not os.path.exists(expected_path):
                        self._die("Output file %s does not exist" % expected_name)
                    args[attr] = expected_path
            if data[0] == QUERY_NAME:
                value: Union[str, None] = None if data[1] == NONE else data[1]
                args[data[0]] = value
        missing_args: List[str] = [
            x for x, y in args.items() if y is None and x in MANDATORY_ARGS
        ]
        if missing_args:
            self._die(
                "The following arguments are missing from the project argument file: %s"
                % ", ".join(missing_args)
            )
        args["cmd"] = None
        return args

    def parse_json(self) -> Dict[str, Union[str, None]]:
        """
        Method for JSON format argument file parsing
        WARNING: This is a method stub; the method's code will be expanded
        once additional project_args formats are implemented
        """
        import json
        with open(self.args_file, "r") as h:
            all_args: Dict[str, Union[str, Dict[any]]] = json.load(h)
        args: Dict[str, Union[str, None]] = {arg: None for arg in ALL_ARGS}
        for arg in ALL_ARGS:
            value: Union[str, None] = all_args["parameters"].get(arg, None)
            if value is None:
                continue
            if arg in MANDATORY_ARGS:
                if arg in MANDATORY_PATH_ARGS:
                    if not os.path.exists(value):
                        self._die(
                            (
                                "Mandatory input file/directory %s does not exist or "
                                "has been deleted since the run passed"
                            ) % value
                        )
                elif arg == ORTH_THRESH:
                    try:
                        args[ORTH_THRESH] = float(value)
                    except ValueError:
                        self._die(
                            "Orthology probability threshold is not a valid floating point number"
                            % value
                        )
                    continue
                elif arg == ACCEPTED_CLASSES:
                    accepted_loss_symbols: List[str] = [x for x in value.split(",") if x]
                    if any(x not in Constants.ALL_LOSS_SYMBOLS for x in accepted_loss_symbols):
                        self._die(
                            '"accepted_loss_symbols" contains inappropriate loss symbols: %s'
                            % ", ".join(
                                [
                                    x
                                    for x in accepted_loss_symbols
                                    if x not in Constants.ALL_LOSS_SYMBOLS
                                ]
                            )
                        )
            args[arg] = value
        missing_args: List[str] = [
            x for x, y in args.items() if y is None and x in MANDATORY_ARGS
        ]
        if missing_args:
            self._die(
                "The following arguments are missing from the project argument file: %s"
                % ", ".join(missing_args)
            )
        output: Union[str, None] = all_args["parameters"].get(OUT_DIR)
        if output is None:
            self._die("Output directory is undefined")
        for attr, expected_name in EXPECTED_OUTPUT_FILE_NAMES.items():
            expected_path: str = os.path.join(output, expected_name)
            if not os.path.exists(expected_path):
                self._die("Output file %s does not exist")
            args[attr] = expected_path
        args["cmd"] = None
        return args

    def parse_yaml(self) -> Dict[str, Union[str, None]]:
        """
        Method for YAML format argument file parsing
        WARNING: This is a method stub; the method's code will be expanded
        once additional project_args formats are implemented
        """
        import yaml
        with open(self.args_file, "r") as h:
            all_args: Dict[str, Union[str, Dict[any]]] = yaml.safe_load(h)
        args: Dict[str, Union[str, None]] = {arg: None for arg in ALL_ARGS}
        for arg in ALL_ARGS:
            value: Union[str, None] = all_args["parameters"].get(arg, None)
            if value is None:
                continue
            if arg in MANDATORY_ARGS:
                if arg in MANDATORY_PATH_ARGS:
                    if not os.path.exists(value):
                        self._die(
                            (
                                "Mandatory input file/directory %s does not exist or "
                                "has been deleted since the run passed"
                            ) % value
                        )
                elif arg == ORTH_THRESH:
                    try:
                        args[ORTH_THRESH] = float(value)
                    except ValueError:
                        self._die(
                            "Orthology probability threshold is not a valid floating point number"
                            % value
                        )
                    continue
                elif arg == ACCEPTED_CLASSES:
                    accepted_loss_symbols: List[str] = [x for x in value.split(",") if x]
                    if any(x not in Constants.ALL_LOSS_SYMBOLS for x in accepted_loss_symbols):
                        self._die(
                            '"accepted_loss_symbols" contains inappropriate loss symbols: %s'
                            % ", ".join(
                                [
                                    x
                                    for x in accepted_loss_symbols
                                    if x not in Constants.ALL_LOSS_SYMBOLS
                                ]
                            )
                        )
            args[arg] = value
        missing_args: List[str] = [
            x for x, y in args.items() if y is None and x in MANDATORY_ARGS
        ]
        if missing_args:
            self._die(
                "The following arguments are missing from the project argument file: %s"
                % ", ".join(missing_args)
            )
        output: Union[str, None] = all_args["parameters"].get(OUT_DIR)
        if output is None:
            self._die("Output directory is undefined")
        for attr, expected_name in EXPECTED_OUTPUT_FILE_NAMES.items():
            expected_path: str = os.path.join(output, expected_name)
            if not os.path.exists(expected_path):
                self._die("Output file %s does not exist")
            args[attr] = expected_path
        args["cmd"] = None
        return args


class SummaryStat:
    """
    A small class for generating results' summary. Supersedes summary generation with ResultChecker
    """

    __slots__ = (
        "ref_2bit",
        "query_2bit",
        "chain_file",
        "ref_annotation",
        "output",
        "orthology_threshold",
        "orth_probs_file",
        "accepted_loss_symbols",
        "loss_summary",
        "query_genes",
        "orthology_classification",
        "ref_isoform_file",
        "u12_file",
        "spliceai_dir",
        "query_name",
        "version",
        "cmd",
        "commit",
    )

    def __init__(
        self,
        ref_2bit: os.PathLike,
        query_2bit: os.PathLike,
        chain_file: os.PathLike,
        ref_annotation: os.PathLike,
        output: os.PathLike,
        orthology_threshold: float,
        orth_probs_file: os.PathLike,
        accepted_loss_symbols: str,
        loss_summary: os.PathLike,
        query_genes: os.PathLike,
        orthology_classification: os.PathLike,
        isoform_file: Optional[Union[os.PathLike, None]] = None,
        u12_file: Optional[Union[os.PathLike, None]] = None,
        spliceai_dir: Optional[Union[os.PathLike, None]] = None,
        query_name: Union[str, None] = None,
        cmd: Optional[Union[str, None]] = None,
    ) -> None:
        """Entry point"""
        self.ref_2bit: os.PathLike = ref_2bit
        self.query_2bit: os.PathLike = query_2bit
        self.chain_file: os.PathLike = chain_file
        self.ref_annotation: os.PathLike = ref_annotation
        self.output: os.PathLike = output
        self.orthology_threshold: float = orthology_threshold
        self.orth_probs_file: os.PathLike = orth_probs_file
        self.accepted_loss_symbols: List[str] = [x for x in accepted_loss_symbols.split(",") if x]
        self.loss_summary: os.PathLike = loss_summary
        self.query_genes: os.PathLike = query_genes
        self.orthology_classification: os.PathLike = orthology_classification
        self.ref_isoform_file: Union[os.PathLike, None] = isoform_file
        self.u12_file: Union[os.PathLike, None] = u12_file
        self.spliceai_dir: Union[os.PathLike, None] = spliceai_dir
        self.query_name: str = (
            query_name if query_name is not None else query_2bit.split(os.sep)[-1]
        )
        self.cmd: Union[str, None] = cmd

        sys.path.append(LOCATION)
        from __version__ import __version__

        self.version: str = __version__
        sys.path.remove(LOCATION)

        # import git

        # repo = git.Repo(LOCATION, search_parent_directories=True)
        # sha = repo.head.commit.hexsha
        # self.commit: str = repo.git.rev_parse(sha, short=7)

    def summary(self) -> str:
        """Main summary method"""
        ## parse the input data
        ref_transcripts: Set[str] = set()
        ref_tr2gene: Dict[str, str] = dict()
        ref_gene2tr: Dict[str, List[str]] = defaultdict(list)
        for data in read_tab(self.ref_annotation):
            tr: str = data[3]
            ref_transcripts.add(tr)
        num_trs: int = len(ref_transcripts)
        if self.ref_isoform_file is not None:
            for data in read_tab(self.ref_isoform_file):
                tr: str = data[1]
                if tr not in ref_transcripts:  ## TODO: Track in a separate collection??
                    continue
                gene: str = data[0]
                ref_tr2gene[tr] = gene
                ref_gene2tr[gene].append(tr)
            ref_gene_num: int = len(ref_gene2tr)
        else:
            ref_gene_num: int = len(ref_transcripts)
        ## summarize the orthology classification step
        tr2chain2prob: Dict[str] = defaultdict(dict)
        for data in read_tab(self.orth_probs_file):
            tr: str = data[0]
            if tr == ORTH_PROB_HEADER:
                continue
            chain: str = data[1]
            prob: float = float(data[2])
            tr2chain2prob[tr][chain] = prob
        num_classified_trs: int = len(tr2chain2prob)
        num_orth_tr: int = sum(
            any(
                tr2chain2prob[x][y] >= self.orthology_threshold
                for y in tr2chain2prob[x]
            )
            for x in tr2chain2prob
        )
        num_zero_orth: int = len(tr2chain2prob) - num_orth_tr
        num_one2one_prob: int = sum(
            sum(
                tr2chain2prob[x][y] >= self.orthology_threshold
                for y in tr2chain2prob[x]
            )
            == 1
            for x in tr2chain2prob
        )
        num_no_proj: int = len(ref_transcripts) - len(tr2chain2prob)
        num_ppgene_pred: int = sum(
            sum(tr2chain2prob[x][y] == -2.0 for y in tr2chain2prob[x])
            for x in tr2chain2prob
        )
        if ref_tr2gene:
            gene_num_orth: int = len(
                {
                    ref_tr2gene[x]
                    for x in tr2chain2prob
                    if any(
                        y >= self.orthology_threshold for y in tr2chain2prob[x].values()
                    )
                }
            )
            gene_num_one2one: int = len(
                {
                    x
                    for x, y in ref_gene2tr.items()
                    if all(
                        sum(
                            c >= self.orthology_threshold
                            for c in tr2chain2prob[z].values()
                        )
                        < 2
                        for z in y
                    )
                    and any(
                        sum(
                            c >= self.orthology_threshold
                            for c in tr2chain2prob[z].values()
                        )
                        == 1
                        for z in y
                    )
                }
            )
            gene_num_zero: int = len(
                {
                    x
                    for x, y in ref_gene2tr.items()
                    if not any(
                        sum(
                            c >= self.orthology_threshold
                            for c in tr2chain2prob[z].values()
                        )
                        for z in y
                    )
                }
            )
            num_genes: int = len(ref_gene2tr)
            gene_num_no_proj: int = num_genes - len(
                {
                    x
                    for x, y in ref_gene2tr.items()
                    if any(len(tr2chain2prob[z]) for z in y)
                }
            )
        else:
            gene_num_orth: int = num_orth_tr
            gene_num_one2one: int = num_one2one_prob
            gene_num_zero: int = num_zero_orth
            num_genes: int = len(ref_transcripts)
            gene_num_no_proj: int = num_no_proj
        ## fetch loss statistics
        proj2loss: Dict[str, Dict[str, int]] = defaultdict(int)
        tr2loss: Dict[str, Dict[str, int]] = defaultdict(int)
        gene2loss: Dict[str, Dict[str, int]] = defaultdict(int)
        paralog2loss: Dict[str, str] = defaultdict(set)
        for data in read_tab(self.loss_summary):
            if data[0] == PROJECTION:
                proj2loss[data[2]] += 1
                if "#paralog" in data[1]:
                    tr = get_proj2trans(data[1])[0]
                    paralog2loss[tr].add(data[2])
                continue
            if data[0] == TRANSCRIPT:
                tr2loss[data[2]] += 1
                if self.ref_isoform_file is None:
                    gene2loss[data[2]] += 1
                continue
            if data[0] == GENE:
                gene2loss[data[2]] += 1
                continue
        proj_loss_classified: int = sum(proj2loss.values())
        proj_present: int = sum(proj2loss.get(x, 0) for x in self.accepted_loss_symbols)
        proj_non_present: int = proj_loss_classified - proj_present
        tr_loss_classified: int = sum(tr2loss.values())
        tr_present: int = sum(tr2loss.get(x, 0) for x in self.accepted_loss_symbols)
        tr_non_present: int = tr_loss_classified - tr_present
        genes_loss_classified: int = sum(proj2loss.values())
        genes_present: int = sum(tr2loss.get(x, 0) for x in self.accepted_loss_symbols)
        genes_non_present: int = genes_loss_classified - genes_present
        ## fetch the functional paralogs
        paralog2loss: Dict[str, str] = {
            x: max(y, key=lambda x: CLASS_TO_NUM[x]) for x, y in paralog2loss.items()
        }
        if self.ref_isoform_file:
            par_genes_found: Set[str] = {ref_tr2gene[x] for x in paralog2loss.keys()}
            num_with_func_para: int = 0
            for _gene in par_genes_found:
                best_status_per_gene: str = max(
                    [paralog2loss.get(x, "N") for x in ref_gene2tr[_gene]],
                    key=lambda x: CLASS_TO_NUM[x],
                )
                if best_status_per_gene in self.accepted_loss_symbols:
                    num_with_func_para += 1
        else:
            num_with_func_para: int = sum(
                y in self.accepted_loss_symbols for x, y in paralog2loss.items()
            )
        perc_with_func_para: float = to_perc(num_with_func_para, ref_gene_num)
        ## gene & orthology statistics
        query_genes: Set[str] = set()
        for data in read_tab(self.query_genes):
            if not data or not data[0]:
                continue
            if data[0] == QUERY_ISOFORMS_HEADER:
                continue
            query_genes.add(data[0])
        lost_gene_num: int = sum("lost_" in x for x in query_genes)
        missing_gene_num: int = sum("missing_" in x for x in query_genes)
        paralog_gene_num: int = sum("paralog_" in x for x in query_genes)
        ppgene_gene_num: int = sum("retro_" in x for x in query_genes)
        undefined_orth: int = (
            lost_gene_num + missing_gene_num + paralog_gene_num + ppgene_gene_num
        )
        defined_orth: int = len(query_genes) - undefined_orth
        gene2orth_class: Dict[str, Set[str]] = defaultdict(set)
        for data in read_tab(self.orthology_classification):
            if data[0] == ORTH_CLASS_HEADER:
                continue
            ref_gene: str = data[0]
            orth_class: str = data[4]
            gene2orth_class[orth_class].add(ref_gene)
        ## we must be all set at this point
        classification: str = CLASSIFICATION_BOILERPLATE.format(
            num_classified_trs,
            num_orth_tr,
            to_perc(num_orth_tr, num_trs),
            num_one2one_prob,
            to_perc(num_one2one_prob, num_trs),
            num_zero_orth,
            to_perc(num_zero_orth, num_trs),
            num_no_proj,
            to_perc(num_no_proj, num_trs),
            num_ppgene_pred,
            gene_num_orth,
            to_perc(gene_num_orth, num_genes),
            gene_num_one2one,
            to_perc(gene_num_one2one, num_genes),
            gene_num_zero,
            to_perc(gene_num_zero, num_genes),
            gene_num_no_proj,
            to_perc(gene_num_no_proj, num_genes),
        )
        loss_summary: str = LOSS_SUMMARY_BOILERPLATE.format(
            ", ".join(self.accepted_loss_symbols),
            proj2loss.get("FI", 0),
            to_perc(proj2loss.get("FI", 0), proj_loss_classified),
            proj2loss.get("I", 0),
            to_perc(proj2loss.get("I", 0), proj_loss_classified),
            proj2loss.get("PI", 0),
            to_perc(proj2loss.get("PI", 0), proj_loss_classified),
            proj2loss.get("UL", 0),
            to_perc(proj2loss.get("UL", 0), proj_loss_classified),
            proj2loss.get("L", 0),
            to_perc(proj2loss.get("L", 0), proj_loss_classified),
            proj2loss.get("M", 0),
            to_perc(proj2loss.get("M", 0), proj_loss_classified),
            proj_present,
            to_perc(proj_present, proj_loss_classified),
            proj_non_present,
            to_perc(proj_non_present, proj_loss_classified),
            tr2loss.get("FI", 0),
            to_perc(tr2loss.get("FI", 0), tr_loss_classified),
            tr2loss.get("I", 0),
            to_perc(tr2loss.get("I", 0), tr_loss_classified),
            tr2loss.get("PI", 0),
            to_perc(tr2loss.get("PI", 0), tr_loss_classified),
            tr2loss.get("UL", 0),
            to_perc(tr2loss.get("UL", 0), tr_loss_classified),
            tr2loss.get("L", 0),
            to_perc(tr2loss.get("L", 0), tr_loss_classified),
            tr2loss.get("M", 0),
            to_perc(tr2loss.get("M", 0), tr_loss_classified),
            tr_present,
            to_perc(tr_present, tr_loss_classified),
            tr_non_present,
            to_perc(tr_non_present, tr_loss_classified),
            gene2loss.get("FI", 0),
            to_perc(gene2loss.get("FI", 0), genes_loss_classified),
            gene2loss.get("I", 0),
            to_perc(gene2loss.get("I", 0), genes_loss_classified),
            gene2loss.get("PI", 0),
            to_perc(gene2loss.get("PI", 0), genes_loss_classified),
            gene2loss.get("UL", 0),
            to_perc(gene2loss.get("UL", 0), genes_loss_classified),
            gene2loss.get("L", 0),
            to_perc(gene2loss.get("L", 0), genes_loss_classified),
            gene2loss.get("M", 0),
            to_perc(gene2loss.get("M", 0), genes_loss_classified),
            genes_present,
            to_perc(genes_present, genes_loss_classified),
            genes_non_present,
            to_perc(genes_non_present, genes_loss_classified),
        )
        orth_summary: str = ORTHOLOGY_BOILERPLATE.format(
            ref_gene_num,
            len(query_genes),
            defined_orth,
            to_perc(defined_orth, len(query_genes)),
            undefined_orth,
            to_perc(undefined_orth, len(query_genes)),
            len(gene2orth_class.get(ONE2ONE, [])),
            to_perc(len(gene2orth_class.get(ONE2ONE, [])), ref_gene_num),
            len(gene2orth_class.get(ONE2MANY, [])),
            to_perc(len(gene2orth_class.get(ONE2MANY, [])), ref_gene_num),
            len(gene2orth_class.get(MANY2ONE, [])),
            to_perc(len(gene2orth_class.get(MANY2ONE, [])), ref_gene_num),
            len(gene2orth_class.get(MANY2MANY, [])),
            to_perc(len(gene2orth_class.get(MANY2MANY, [])), ref_gene_num),
            len(gene2orth_class.get(ONE2ZERO, [])),
            to_perc(len(gene2orth_class.get(ONE2ZERO, [])), ref_gene_num),
        )
        # if self.ref_isoform_file is not None:
        #     isoforms_line: str = ISOFORMS_LINE.format(self.ref_isoform_file)
        # else:
        #     isoforms_line: str = ISOFORMS_CAVEAT
        # summary_text: str = SUMMARY_BOILERPLATE.format(
        #     br=BREAK_LINE,
        #     ref_2bit=self.ref_2bit,
        #     query_2bit=self.query_2bit,
        #     chains=self.chain_file,
        #     ref_annot=self.ref_annotation,
        #     output=self.output,
        #     isoforms_line=isoforms_line,
        #     orthology_class_report=classification,
        #     loss_summary=loss_summary,
        #     orthology_res_report=orth_summary,
        # )

        num_with_func_orth: int = sum(
            len(y) for x, y in gene2orth_class.items() if x != ONE2ZERO
        )
        perc_with_func_orth: float = to_perc(num_with_func_orth, ref_gene_num)

        cmd_line: str = "" if self.cmd is None else CMD.format(cmd=self.cmd)

        main_header: str = MAIN_HEADER.format(
            br=BREAK_LINE,
            ref_gene_num=ref_gene_num,
            num_genes=defined_orth,
            num_retro=ppgene_gene_num,
            num_lost=lost_gene_num,
            num_missing=missing_gene_num,
            num_with_func_orth=num_with_func_orth,
            perc_with_func_orth=perc_with_func_orth,
            num_one2one=len(gene2orth_class.get(ONE2ONE, [])),
            perc_one2one=to_perc(len(gene2orth_class.get(ONE2ONE, [])), ref_gene_num),
            num_with_func_para=num_with_func_para,
            perc_with_func_para=perc_with_func_para,
            version=self.version,
            # commit=self.commit,
            ref_2bit=self.ref_2bit,
            query_2bit=self.query_2bit,
            chains=self.chain_file,
            ref_annot=self.ref_annotation,
            output=self.output,
            isoforms=arg_or_na(self.ref_isoform_file),
            u12_file=arg_or_na(self.u12_file),
            spliceai_dir=arg_or_na(self.spliceai_dir),
            query_name=self.query_name,
            cmd=cmd_line,
        )
        warnings: str = ""
        if not self.ref_isoform_file:
            warnings += "\n" + ISOFORMS_CAVEAT
        detailed_header: str = DETAILED_HEADER.format(br=BREAK_LINE, warnings=warnings)
        summary_text: str = _SUMMARY_BOILERPLATE.format(
            br=BREAK_LINE,
            header=main_header,
            detailed_header=detailed_header,
            orthology_class_report=classification,
            loss_summary=loss_summary,
            orthology_res_report=orth_summary
        )

        return summary_text
