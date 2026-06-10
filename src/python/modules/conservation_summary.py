#!/usr/bin/env python3

"""
Summarises projection classification data based on TOGA results
"""

from collections import defaultdict
from sys import stdout
from typing import Dict, List, Optional, Set, TextIO, Tuple, Union

import click

from .cesar_wrapper_constants import CLASS_TO_NUM
from .constants import Headers, RejectionReasons
from .shared import CONTEXT_SETTINGS, base_proj_name, parse_single_column

__author__ = "Yury V. Malovichko"
__year__ = "2024"
__email__ = "yury.malovichko@senckenberg.de"
__credits__ = ("Bogdan Kirilenko", "Michael Hiller", "Virag Sharma", "David Jebb")

## define local constants
PROJECTION: str = "PROJECTION"
TRANSCRIPT: str = "TRANSCRIPT"
GENE: str = "GENE"
CESAR_STEP_REJECT: str = "No exons were aligned"
NO_EXONS: str = "NO_EXONS_FOUND"
LOW_COV: str = "INSUFFICIENT_SEQ_COVERAGE"
REDUNDANT_PARALOG: str = "REDUNDANT_PARALOG"
REDUNDANT_PPGENE: str = "REDUNDANT_PPGENE"
SECOND_BEST: str = "SECOND_BEST"
IGNORED_ITEMS: Tuple[str] = ("REDUNDANT_PARALOG", "REDUNDANT_PPGENE", "SECOND_BEST")
SPANNING_CLASSES: str = ("L", "M", "N")
PG: str = "PG"

def parse_precedence_file(file: TextIO) -> Dict[str, str]:
    """ """
    tr2curr_best: Dict[str, Tuple[str, int, int, str]] = {}
    for line in file:
        data: List[str] = line.rstrip().split("\t")
        if not data or not data[0]:
            continue
        if data[0] == "projection":
            continue
        proj: str = data[0]
        tr: str = "#".join(proj.split("#")[:-1])
        # chrom: str = data[1]
        start: int = int(data[2])
        end: int = int(data[3])
        if tr not in tr2curr_best:
            tr2curr_best[tr] = (proj, start, end)
        else:
            prev_best_start, prev_best_end = tr2curr_best[tr][1:]
            # if start >= prev_best_start and end <= prev_best_end:
            if end - start < prev_best_end - prev_best_start:
                tr2curr_best[tr] = (proj, start, end)
    tr2best: Dict[str, str] = {k: v[0] for k, v in tr2curr_best.items()}
    return tr2best


def transcript_meta_to_report(
    file: Union[str, TextIO],
    precedence: Optional[Dict[str, str]] = {},
    paralogs: Optional[Set[str]] = set(),
    ppgenes: Optional[Set[str]] = set(),
    discarded_paralogs: Optional[Set[str]] = set(),
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Parses transcript_meta.tsv output file of TOGA, inferring projection
    classification and best transcript classification
    """
    proj2status: Dict[str, str] = {}
    tr2status: Dict[str, Set[str]] = defaultdict(set)
    potentially_missing_paralogs: Set[str] = set()
    confirmed_missing_paralogs: Set[str] = set()
    if not isinstance(file, str):
        lines: List[str] = file.readlines()
    else:
        with open(file, "r") as h:
            lines: List[str] = h.readlines()
    has_ppgenes: Set[str] = set()
    ## extract the projection loss statuses
    for line in lines:
        data: List[str] = line.strip().split("\t")
        if data[0] == "projection":
            continue
        proj: str = data[0]
        basename: str = base_proj_name(proj)
        tr: str = "#".join(proj.split("#")[:-1])
        status: str = data[1]
        proj2status[proj] = status
        ## do not propagate paralog and ppgene loss status to downstream levels
        if proj in ppgenes or basename in ppgenes:
            tr2status[tr].add("N")
            has_ppgenes.add(tr)
            continue
        ## for paralogs, do not count the discarded projections towards the transcript classification
        if proj in paralogs:
            if proj in discarded_paralogs or basename in discarded_paralogs:
                potentially_missing_paralogs.add(tr)
                continue
            status = "PG"
        tr2status[tr].add(status)
    ## for paralog-only transcripts, assign the Missing status to those
    ## that have all their projections discarded
    for tr in potentially_missing_paralogs:
        if tr in tr2status:
            continue
        tr2status[tr] = "M"
        confirmed_missing_paralogs.add(tr)
    ## infer the transcript level loss statuses and add them to output
    for tr, all_classes in tr2status.items():
        max_status: str = max(all_classes, key=lambda x: CLASS_TO_NUM[x])
        ## if transcript has an established precedence order (e.g., nested spanning chains),
        ## use the top projection's status as the transcript's status estimate
        if tr in precedence and max_status in SPANNING_CLASSES:
            ## TODO: The max_status check above is a safeguard for rare quirks
            ## when both spanning and regular orthologous projections are encountered;
            ## will be likely redundant in 2.1
            preferred_proj: str = precedence[tr]
            basename: str = base_proj_name(preferred_proj)
            if preferred_proj in proj2status and (
                preferred_proj not in ppgenes and basename not in ppgenes
            ):
                preferred_loss_status: str = proj2status[preferred_proj]
                tr2status[tr] = preferred_loss_status
                continue
        ## if a transcript does not have any projections except for ppgenes, treat it as missing
        if max_status == "N" and tr in has_ppgenes:
            max_status = "M"
        tr2status[tr] = max_status

    return (proj2status, tr2status, confirmed_missing_paralogs)


# def rejection_file_to_report( ## LEGACY FORMAT
#     file: Union[str, TextIO]
# ) -> Tuple[Dict[str, str], Dict[str, str]]:
#     """
#     Parses the genes_rejection_reason.tsv file containing transcripts which were
#     rejected for one reason or another at the CESAR alignment step or upstream
#     """
#     ## TODOS:
#     ## 1) Solicit the current classification scheme with Michael
#     ##    (I suspect summing the number of exons or their lengths is better than
#     ##     summing the number of exon groups at the last classificaiton step)
#     ## 2) Double-check that all the reasons have been
#     proj2status: Dict[str, str] = {}
#     tr2status: Dict[str, str] = {}
#     if not isinstance(file, str):
#         lines: List[str] = file.readlines()
#     else:
#         with open(file, 'r') as h:
#             lines: List[str] = h.readlines()
#     for line in lines:
#         data: List[str] = line.strip().split('\t')
#         proj: str = data[0]
#         tr: str = '.'.join(proj.split('.')[:-1])
#         ## if projection was rejected prior to the CESAR step, it's marked as missing
#         if data[2] != CESAR_STEP_REJECT:
#             proj2status[proj] = 'M'
#             tr2status[proj] = 'M'
#             continue
#         reasons: List[str] = [tuple(x.split(':')) for x in data[3].split(';')]
#         ## projections for which no exons were reliably located by alignment
#         ## results are most likely lost
#         if any(map(lambda x: x[0] == NO_EXONS, reasons)):
#             proj2status[proj] = 'L'
#             tr2status[tr] = 'L'
#             continue
#         ## projections dismissed at alignment step due to insufficient fraction
#         ## of exons properly located by alignment data are (likely) missing
#         if any(map(lambda x: x[0] == LOW_COV, reasons)):
#             proj2status[proj] = 'M'
#             tr2status[tr] = 'M'
#             continue
#         ## other reasons come in two flavours: if assembly gaps were encountered
#         ## within the search space, the rejection label has the '+GAP' suffix
#         aln_gap_num: int = sum((int(x[1]) for x in reasons if '+GAP' in x[0]))
#         no_gap_num: int = sum((int(x[1]) for x in reasons if '+GAP' not in x[0]))
#         ## now, if the number of gap-containing groups ouweighs the counterpart,
#         ## mark the projection as missing unless there's already a more defined
#         if aln_gap_num > no_gap_num:
#             proj2status[proj] = 'M' if proj not in proj2status else proj2status[proj]
#             tr2status[tr] = 'M'
#         else:
#             proj2status[proj] = 'L'
#             tr2status[tr] = 'L'
#
#     return (proj2status, tr2status)


def rejection_file_to_report(
    file: Union[str, TextIO],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Parses the rejection_reasons.tsv file containing transcripts which were
    rejected for one reason or another at the CESAR alignment step or upstream
    """
    proj2status: Dict[str, str] = {}
    tr2all_statuses: Dict[str, Set[str]] = defaultdict(set)
    tr2status: Dict[str, str] = {}
    if not isinstance(file, str):
        lines: List[str] = file.readlines()
    else:
        with open(file, "r") as h:
            lines: List[str] = h.readlines()
    for line in lines:
        data: List[str] = line.rstrip().split("\t")
        if not data or not data[0]:
            continue
        if data[0] == "level":
            continue
        if len(data) == 2:  ## TEMPORARY SOLUTION TO BYPASS THE UPSTREAM BUG
            name: str = data[0]
            status: str = "N"
            tr2status[name] = status
            continue
        level: str = data[0]
        name: str = data[1]
        status: str = data[5]
        if level == TRANSCRIPT:
            tr2status[name] = status
            continue
        if level != PROJECTION:
            raise ValueError(f"An ambiguous entry found: {line}")
        tr: str = "#".join(name.split("#")[:-1])
        reason: str = data[4]
        if reason in IGNORED_ITEMS:
            continue
        proj2status[name] = status
        tr2all_statuses[tr].add(status)
    for tr in tr2all_statuses:
        if tr in tr2status:
            continue
        all_classes: Set[str] = tr2all_statuses[tr]
        max_status: str = max(all_classes, key=lambda x: CLASS_TO_NUM[x])
        tr2status[tr] = max_status

    return (proj2status, tr2status)


def add_rejection_data(
    orig_proj2status: Dict[str, str],
    orig_tr2status: Dict[str, str],
    reject_proj2status: Dict[str, str],
    reject_tr2status: Dict[str, str],
    precedence: Optional[Dict[str, str]] = {},
    paralogs: Optional[Set[str]] = set(),
    ppgenes: Optional[Set[str]] = set(),
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Updates projection- and transcript-level loss status dictionaries with
    the data for the beforehand-rejected projections. Adds the new instance
    if no instance of this projection or transcript was subjected to alignment,
    leaves the best loss status otherwise
    """
    tr_from_rej_log: Set[str] = set()
    for proj, rej_proj_status in reject_proj2status.items():
        tr: str = "#".join(proj.split("#")[:-1])
        if proj not in orig_proj2status:
            orig_proj2status[proj] = rej_proj_status
        else:
            orig_proj_status: str = orig_proj2status[proj]
            orig_proj2status[proj] = max(
                (orig_proj_status, rej_proj_status), key=lambda x: CLASS_TO_NUM[x]
            )
        ## do not propagate rejected paralog statuses
        if proj in paralogs or proj in ppgenes:
            continue
        rej_tr_status: str = reject_tr2status[tr]
        if tr not in orig_tr2status:
            tr_from_rej_log.add(tr)
            orig_tr2status[tr] = rej_tr_status
        ## do not account for items non-spanning items
        elif rej_proj_status not in SPANNING_CLASSES:
            continue
        else:
            orig_tr_status: str = orig_tr2status[tr]
            orig_tr2status[tr] = max(
                (orig_tr_status, rej_tr_status), key=lambda x: CLASS_TO_NUM[x]
            )
        if tr in precedence and orig_tr2status in SPANNING_CLASSES:
            preferred_proj: str = precedence[tr]
            if preferred_proj in orig_proj2status:
                preferred_loss_status: str = orig_proj2status[preferred_proj]
                orig_tr2status[tr] = preferred_loss_status
            elif preferred_proj in reject_proj2status:
                preferred_loss_status: str = reject_proj2status[preferred_proj]
                orig_tr2status[tr] = preferred_loss_status
            else:
                raise KeyError(
                    f"Top-precedent projection {preferred_proj} is absent "
                    "from both transcript meta and rejection reports"
                )
            continue
    for rej_tr, rej_tr_status in reject_tr2status.items():
        if rej_tr not in orig_tr2status:
            orig_tr2status[rej_tr] = rej_tr_status
        else:
            orig_tr_status: str = orig_tr2status[tr]
            ## TODO: This check is a safeguard for rare quirks
            ## when both spanning and regular orthologous projections are encountered;
            ## will be likely redundant in 2.1
            if orig_tr_status == PG:
                continue
            orig_tr2status[tr] = max(
                (orig_tr_status, rej_tr_status), key=lambda x: CLASS_TO_NUM[x]
            )
    return (orig_proj2status, orig_tr2status)


def gene_loss_report(
    file: Union[str, TextIO], tr_report: Dict[str, str]
) -> Dict[str, str]:
    """
    Parses the TOGA-formatted isoform file to collapse transcript loss report
    to gene level
    """
    gene2status: Dict[str, Set[str]] = defaultdict(set)
    if not isinstance(file, str):
        lines: List[str] = file.readlines()
    else:
        with open(file, "r") as h:
            lines: List[str] = h.readlines()
    for line in lines:
        if line.startswith("GeneID"):
            continue
        gene, isoform = line.strip().split("\t")
        gene2status[gene].add(tr_report.get(isoform, "N"))
    for gene in gene2status:
        all_classes: Set[str] = gene2status[gene]
        max_class: str = max(all_classes, key=lambda x: CLASS_TO_NUM[x])
        gene2status[gene] = max_class

    return gene2status


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    "transcript_meta", type=click.File("r", lazy=True), metavar="TRANSCRIPT_META"
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
    "--discarded_paralogs",
    "-rp",
    type=click.File("r", lazy=True),
    metavar="DISCARDED_PARALOG_FILE",
    default=None,
    show_default=True,
    help=(
        "A single-column file containing names of discarded paralogous projections. "
        "For transcripts which are annotated with paralogous projections alone, "
        "only those which retain at least one non-discarded projections are kept for "
        "transcript- and gene-level loss estimation"
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
def main(
    transcript_meta: click.File,
    rejected_projections: Optional[click.File],
    isoform_file: Optional[Union[click.File, None]],
    paralogs: Optional[Union[click.File, None]],
    processed_pseudogenes: Optional[Union[click.File, None]],
    spanning_chains_precedence_file: Optional[Union[click.File, None]],
    discarded_paralogs: Optional[Union[click.File, None]],
    output: Optional[click.File],
) -> None:
    """
    Summarises projection classification data based on TOGA results. Arguments are:\n
    \tTRANSCRIPT_META is a projection metadata output TSV file of the CESAR wrapper
    (transcript_meta.tsv for CESAR wrapper or combined_projection_meta.tsv for TOGA)\n
    """
    if spanning_chains_precedence_file:
        spanning_status: Dict[str, str] = parse_precedence_file(
            spanning_chains_precedence_file
        )
    else:
        spanning_status: Dict[str, str] = {}
    paralogs: Set[str] = parse_single_column(paralogs)
    ppgenes: Set[str] = parse_single_column(processed_pseudogenes)
    discarded_paralogs: Set[str] = parse_single_column(discarded_paralogs)
    proj2status, tr2status, confirmed_rejected_paralogs = transcript_meta_to_report(
        transcript_meta,
        spanning_status,
        paralogs=paralogs,
        ppgenes=ppgenes,
        discarded_paralogs=discarded_paralogs,
    )
    if rejected_projections is not None:
        with open(rejected_projections, "r") as h:
            r_proj2status, r_tr2status = rejection_file_to_report(h)
            proj2status, tr2status = add_rejection_data(
                proj2status, 
                tr2status, 
                r_proj2status,
                r_tr2status, 
                spanning_status,
                paralogs=paralogs,
                ppgenes=ppgenes,
            )
    output.write(Headers.LOSS_FILE_HEADER)
    for proj, status in proj2status.items():
        output.write("\t".join((PROJECTION, proj, status)) + "\n")
    for tr, status in tr2status.items():
        output.write("\t".join((TRANSCRIPT, tr, status)) + "\n")
    ## TODO: Make transcripts appearing in the isoform_file but not in the TOGA results be listed in output as N/M
    if isoform_file is not None:
        gene2status: Dict[str, str] = gene_loss_report(isoform_file, tr2status)
        for gene, status in gene2status.items():
            output.write("\t".join((GENE, gene, status)) + "\n")
    if confirmed_rejected_paralogs and rejected_projections is not None:
        with open(rejected_projections, "a") as h:
            for transcript in confirmed_rejected_paralogs:
                h.write(RejectionReasons.OUTCOMPETED_PARALOG_REASON.format(transcript) + "\n")


if __name__ == "__main__":
    main()
