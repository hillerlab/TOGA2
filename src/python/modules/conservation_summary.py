"""
Summarises projection classification data based on TOGA results
"""

from collections import defaultdict
from sys import stdout
from typing import Final, TextIO

import click

from .cesar_wrapper_constants import CLASS_TO_NUM
from .constants import Headers, RejectionReasons
from .shared import (
    CONTEXT_SETTINGS,
    CommandLineManager,
    base_proj_name,
    parse_single_column,
    read_tab,
)

__author__ = "Yury V. Malovichko"
__year__ = "2024"
__email__ = "yury.malovichko@senckenberg.de"
__credits__ = ("Bogdan Kirilenko", "Michael Hiller", "Virag Sharma", "David Jebb")

## define local constants
PROJECTION: Final[str] = "PROJECTION"
TRANSCRIPT: Final[str] = "TRANSCRIPT"
GENE: Final[str] = "GENE"
CESAR_STEP_REJECT: Final[str] = "No exons were aligned"
NO_EXONS: Final[str] = "NO_EXONS_FOUND"
LOW_COV: Final[str] = "INSUFFICIENT_SEQ_COVERAGE"
REDUNDANT_PARALOG: Final[str] = "REDUNDANT_PARALOG"
REDUNDANT_PPGENE: Final[str] = "REDUNDANT_PPGENE"
SECOND_BEST: Final[str] = "SECOND_BEST"
IGNORED_ITEMS: Final[tuple[str, ...]] = ("REDUNDANT_PARALOG", "REDUNDANT_PPGENE", "SECOND_BEST")
SPANNING_CLASSES: Final[tuple[str, ...]] = ("L", "M", "N")
PG: Final[str] = "PG"

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.option(
    "--transcript_meta", 
    type=click.File("r", lazy=True), 
    metavar="TRANSCRIPT_META",
    required=True,
    help="Path to the transcript meta file"
)
@click.option(
    "--rejected_projections",
    "-r",
    type=click.Path(exists=True),
    metavar="TSV",
    default=None,
    show_default=True,
    help=(
        "A six-column file containing projections discarded prior "
        "to actual CESAR alignment and respective rejection reasons"
    ),
)
@click.option(
    "--isoform_file",
    "-i",
    type=click.File("r", lazy=True),
    default=None,
    show_default=True,
    metavar="TSV",
    help="A path to two-column gene-to-isoform mapping file",
)
@click.option(
    "--paralogs",
    "-p",
    type=click.File("r", lazy=True),
    default=None,
    show_default=True,
    help=("A path to a single-column file containing paralogous projections' names"),
)
@click.option(
    "--processed_pseudogenes",
    "-pp",
    type=click.File("r", lazy=True),
    default=None,
    show_default=True,
    help=(
        "A path to a single-column file containing projection names of processed pseudogenes"
    ),
)
@click.option(
    "--spanning_chains_precedence_file",
    "-span",
    type=click.File("r", lazy=True),
    metavar="SPANNING_CHAIN_FILE",
    default=None,
    show_default=True,
    help=(
        "A file containing spanning chain gap coordinates in the reference. "
        "Nestedness of spanning gaps is used to prioritise chains for loss inference"
    ),
)
@click.option(
    "--output",
    "-o",
    type=click.File("w"),
    default=stdout,
    show_default=False,
    help="A path to store the output at [default: stdout]",
)
@click.option(
    "--log_name",
    "-ln",
    type=str,
    metavar="STR",
    default=None,
    show_default=True,
    help="Logger name to use; relevant only upon main class import",
)
@click.option(
    "--verbose", "-v", is_flag=True, default=False, help="Control logging verbosity"
)

class ConservationSummary(CommandLineManager):
    __slots__ = (
        "gene2status",
        "gene2trs",
        "isoform_file",
        "missing_paralogs",
        "output",
        "paralogs",
        "potentially_missing_paralogs", 
        "ppgenes",
        "precedence",
        "precedence_file",
        "proj2status",
        "rejected",
        "rejection_file",
        "second_bests",
        "status_by_orth",
        "tr2gene",
        "tr2projs",
        "tr2status",
        "transcript_meta",
    )

    def __init__(
        self,
        transcript_meta: click.File,
        rejected_projections: click.Path | None = None,
        isoform_file: click.File | None = None,
        paralogs: click.File | None = None,
        processed_pseudogenes: click.File | None = None,
        spanning_chains_precedence_file: click.File | None = None,
        discarded_paralogs: click.File | None = None,
        output: click.File | TextIO = stdout,
        log_name: str | None = None,
        verbose: bool = False,
    ) -> None:
        self.v: bool = verbose
        self.set_logging(name=log_name, toga_module="conservation_summary")

        self.transcript_meta: click.File = transcript_meta
        self.rejection_file: click.Path | None = rejected_projections
        self.isoform_file: click.File | None = isoform_file
        self.precedence_file: click.File | None = spanning_chains_precedence_file
        self.output: click.File | TextIO = output

        self.paralogs: set[str] = parse_single_column(paralogs)
        self.ppgenes: set[str] = parse_single_column(processed_pseudogenes)

        self.gene2status: dict[str, str] = {}
        self.gene2trs: dict[str, list[str]] = defaultdict(list)
        self.missing_paralogs: set[str] = set()
        self.potentially_missing_paralogs: set[str] = set()
        self.precedence: dict[str, str] = {}
        self.proj2status: dict[str, str] = {}
        self.rejected: set[str] = set()
        self.second_bests: set[str] = set()
        self.tr2gene: dict[str, str] = {}
        self.tr2projs: dict[str, list[str]] = defaultdict(list)
        self.tr2status: dict[str, str] = {}

        self.run()

    def run(self) -> None:
        """Entry point method"""
        ## parse isoform file
        self.parse_isoform_file()
        ## read precedence file
        self.parse_precedence_file()
        ## parse rejection log
        self.parse_rejection_file()
        ## read transcript meta and perform the high-level analysis
        self.parse_transcript_meta()
        ## perform the transcript level analysis
        self.classify_transcripts()
        ## if isoform file was provided, perform the gene loss analysis
        self.classify_genes()
        ## write the results to the file
        self.write_output()
        ## if any items were discarded, write them to the respective files
        self.write_rejected()

    def parse_isoform_file(self) -> None:
        """
        Reads reference {gene: transcript} mapping

        Arguments:
            None

        Returns:
            None

        Note:
            1. If isoform file was not provided, exits without populating class attributes
            2. Otherwise, populates ``self.tr2gene`` with {transcript -> gene} mapping, 
            and ``self.gene2trs`` with {gene -> [transcripts]} mapping
        """
        if self.isoform_file is None:
            return
        for i, data in enumerate(read_tab(self.isoform_file), start=1):
            if not data or not data[0]:
                continue
            if data[0] == "GeneID":
                continue
            if len(data) != 2:
                self._die(
                    f"Improper reference isoform file at line {i}: expected 2 fields, got {len(data)}"
                )
            gene, tr = data
            if tr in self.tr2gene:
                self._die(
                    f"Transcript {tr} appears more than once in the reference isoform file; "
                    f"first duplicate encountered at line {i}"
                )
            self.tr2gene[tr] = gene
            self.gene2trs[gene].append(tr)

    def parse_rejection_file(self) -> None:
        """
        Extracts projection and transcript annotation from rejection log file 

        Arguments:
            None

        Returns:
            None
        """
        if self.rejection_file is None:
            return
        for i, data in enumerate(read_tab(self.rejection_file), start=1):
            if not data or not data[0]:
                continue
            if data[0] == "level":
                continue
            ## workaround for compatibility with bugged older versions 
            if len(data) == 2:
                name: str = data[0]
                status: str = "N"
                self.tr2status[name] = status
                continue
            level: str = data[0]
            name: str = data[1]
            ## the entire transcript has been discarded;
            ## no reasons to expect further projections
            status: str = data[5]
            if level == TRANSCRIPT:
                self.tr2status[name] = max((status, self.tr2status.get(name, "N")), key=lambda x: CLASS_TO_NUM[x])
                continue
            if level != PROJECTION:
                raise ValueError(f"An ambiguous entry found at line {i}: {level}")
            basename: str = base_proj_name(name)
            is_paralog: bool = name in self.paralogs or basename in self.paralogs
            is_ppgene: bool = name in self.ppgenes or basename in self.ppgenes
            if is_paralog and is_ppgene:
                self._die(
                    f"Projection {name} at rejection file line {i} "
                    "is considered both a paralog and a processed pseudogene"
                )
            tr: str = "#".join(name.split("#")[:-1])
            reason: str = data[4]
            if reason == SECOND_BEST:
                self.second_bests.add(name)
            ## save the projection name
            self.rejected.add(name)
            ## save the projections loss status
            self.proj2status[name] = status
            ## add projection to the transcript -> projections mapping
            self.tr2projs[tr].append(name)

    def parse_precedence_file(self) -> None:
        """
        Parses coordinates of spanning sections of spanning chain projections.

        For each transcript covered by spanning chains, the tightest span is selected 
        for projection loss station classification. Consider the following example:
        ```
        projection    ==============[***]=====[***]================
        chain1        ----[]--------------------------------[]-----
        chain2        -----------[]-----------------------------[]-
        chain3        ---------[]--------------------[]------------
        ```
        In this example, chain 3 is selected as the best representative since its span 
        (the distance between the two aligned blocks) is the smallest of all three chains 
        considered.

        Arguments:
            None

        Returns:
            None
        """
        if self.precedence_file is None:
            return
        tr2curr_best: dict[str, tuple[str, int, int]] = {}
        for i, data in enumerate(read_tab(self.precedence_file), start=1):
            if not data or not data[0]:
                continue
            if data[0] == "projection":
                continue
            if len(data) < 4:
                self._die(
                    f"Improper precedence file formatting at line {i}: "
                    f"expected at least 4 fields, got {len(data)}"
                )
            proj: str = data[0]
            tr: str = "#".join(proj.split("#")[:-1])
            # chrom: str = data[1]
            start: int = int(data[2])
            end: int = int(data[3])
            if tr not in tr2curr_best:
                tr2curr_best[tr] = (proj, start, end)
            else:
                prev_best_start, prev_best_end = tr2curr_best[tr][1:]
                if end - start < prev_best_end - prev_best_start:
                    tr2curr_best[tr] = (proj, start, end)
        self.precedence = {k: v[0] for k, v in tr2curr_best.items()}

    def parse_transcript_meta(self) -> None:
        """
        Reads transcript metadata file and extracts loss status classification

        Arguments:
            None

        Returns:
            None
        """
        ## extract the projection loss statuses
        for i, data in enumerate(read_tab(self.transcript_meta), start=1):
            if data[0] == "projection":
                continue
            proj: str = data[0]
            basename: str = base_proj_name(proj)
            is_paralog: bool = proj in self.paralogs or basename in self.paralogs
            is_ppgene: bool = proj in self.ppgenes or basename in self.ppgenes
            if is_paralog and is_ppgene:
                self._die(
                    f"Projection {proj} at transcript meta file {i} "
                    "is considered both a paralog and a processed pseudogene"
                )
            ## save the projection's loss status
            status: str = data[1]
            self.proj2status[proj] = status
            ## add projection to the transcript -> projections mapping
            tr: str = "#".join(proj.split("#")[:-1])
            self.tr2projs[tr].append(proj)

    def classify_transcripts(self) -> None:
        """
        Classifies transcripts in terms of loss status.

        By default, transcript loss status is defined as the highest loss status among 
        its projections. Loss status precedence order is defined as followis:
        ```
        FI > I > PI > UL > L > M > N
        ```
        Note that special loss statuses PG (ParaloGous) and PP (Processed Pseudogene) do not 
        participate in the transcript status classification; moreover, orthologous projections 
        are expected never to have this loss status.

        If a transcript does not have any following projections, the loss status is defined as follows:
        * If a transcript has any spanning chains, the loss status is defined 
        based on the the most compact span (see ``parse_precedence_file``);
        * If a transcript has paralogous projections, it gets the PG status unless all the paralogous 
        projections have been discarded, otherwise it gets the Missing (M) status
        * Otherwise, the transcript gets Missing status

        Arguments:
            None

        Returns:
            None
        """
        for tr, projs in self.tr2projs.items():
            if tr in self.tr2status:
                continue
            ppgenes: list[str] = [
                x for x in projs if x in self.ppgenes or base_proj_name(x) in self.ppgenes
            ]
            paralogs: list[str] = [
                x for x in projs if x in self.paralogs or base_proj_name(x) in self.paralogs
            ]
            orthologs: list[str] = [x for x in projs if x not in ppgenes and x not in paralogs]
            ## if any ortholog is present, use them for loss status classification
            if orthologs:
                non_rejected: list[str] = [x for x in orthologs if x not in self.rejected]
                if non_rejected:
                    all_statuses: list[str] = [self.proj2status[x] for x in non_rejected]
                    max_status: str = max(all_statuses, key=lambda x: CLASS_TO_NUM[x])
                    self.tr2status[tr] = max_status
                    continue
                ## otherwise, check if there are any spanning chains recorded
                ## (they must be also present in the rejection log file)
                if tr in self.precedence:
                    self.tr2status[tr] = self.proj2status[self.precedence[tr]]
                    continue
                ## otherwise, classify based on other rejected projections
                rejected: list[str] = [x for x in orthologs if x not in non_rejected]
                non_sb: list[str] = [x for x in rejected if  x not in self.second_bests]
                if non_sb:
                    all_statuses: list[str] = [self.proj2status[x] for x in non_sb]
                    max_status: str = max(all_statuses, key=lambda x: CLASS_TO_NUM[x])
                    if max_status != "N" or not paralogs:
                        self.tr2status[tr] = max_status
                        continue
                elif rejected:
                    self.tr2status[tr] = "M"
            ## if there are any paralogs recorded, the transcript gets the "PG" status
            if paralogs:
                non_rejected: list[str] = [x for x in paralogs if x not in self.rejected]
                if non_rejected:
                    self.tr2status[tr] = "PG"
                else:
                    self.tr2status[tr] = "M"
                    self.missing_paralogs.add(tr)
                continue
            ## if there are processed pseudogenes, assign the Missing status
            if ppgenes:
                self.tr2status[tr] = "M"
                continue
            # self._to_log(f"")
            # self.tr2status[tr] = "M"
            self._to_log(f"Transcript {tr} has no projections and is missing from the rejection file")

    def classify_genes(self) -> None:
        """
        Classifies reference genes based on their transcripts' classification. 
        The process is pretty straightforward and relies on the loss status 
        precedence order only

        Arguments:
            None

        Returns:
            None
        """
        if not self.gene2trs:
            return
        for gene, trs in self.gene2trs.items():
            all_statuses: list[str] = [self.tr2status.get(x, "N") for x in trs]
            max_status: str = max(all_statuses, key=lambda x: CLASS_TO_NUM[x])
            self.gene2status[gene] = max_status

    def write_output(self) -> None:
        """
        Writes the results to the output file

        Arguments:
            None

        Returns:
            None
        """
        self.output.write(Headers.LOSS_FILE_HEADER)
        self.output.writelines(
            [f"{PROJECTION}\t{proj}\t{status}\n" for proj, status in sorted(self.proj2status.items())]
        )
        self.output.writelines(
            [f"{TRANSCRIPT}\t{tr}\t{status}\n" for tr, status in sorted(self.tr2status.items())]
        )
        self.output.writelines(
            [f"{GENE}\t{gene}\t{status}\n" for gene, status in sorted(self.gene2status.items())]
        )

    def write_rejected(self) -> None:
        """
        Appends rejection entries of transcripts with all their 
        paralogous projections missing  to ``self.rejection_file``

        Arguments:
            None

        Returns:
            None
        """
        if not self.missing_paralogs:
            return
        if self.rejection_file is not None:
            with open(self.rejection_file, "a") as h:
                h.writelines([
                    RejectionReasons.OUTCOMPETED_PARALOG_REASON.format(tr) + "\n" for tr in sorted(self.missing_paralogs)
                ])
