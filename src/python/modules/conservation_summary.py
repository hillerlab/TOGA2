#!/usr/bin/env python3

"""
Summarises projection classification data based on TOGA results
"""
from collections import defaultdict
from constants import Headers
from shared import CONTEXT_SETTINGS, parse_single_column
from cesar_wrapper_constants import CLASS_TO_NUM
from sys import stdout
from typing import Dict, List, Optional, Set, TextIO, Tuple, Union

import click

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__version__ = '2.0'
__email__ = 'yury.malovichko@senckenberg.de'
__credits__ = ('Bogdan Kirilenko', 'Michael Hiller', 'Virag Sharma', 'David Jebb')

## define local constants
##
# NUM_TO_CLASS: Dict[int, str] = {
#     -1: 'N', 0: 'PP', 1: 'PG', 2: 'PM', 3: 'M', 4: 'L', 5: 'UL', 6: 'PI', 7: 'I', 8: 'FI'
# }
PROJECTION: str = 'PROJECTION'
TRANSCRIPT: str = 'TRANSCRIPT'
GENE: str = 'GENE'
# CLASS_TO_NUM: Dict[str, int] = {v: k for k, v in NUM_TO_CLASS.items()}
# ## set color code for the HTML report
# CLASS_TO_COL: Dict[str, str] = {
#     N_: BLACK,
#     PG: BROWN,
#     PM: GREY,
#     L: LIGHT_RED,
#     M: GREY,
#     UL: SALMON,
#     PI: LIGHT_BLUE,
#     I: BLUE,
# }
CESAR_STEP_REJECT: str = 'No exons were aligned'
NO_EXONS: str = 'NO_EXONS_FOUND'
LOW_COV: str = 'INSUFFICIENT_SEQ_COVERAGE'

def parse_precedence_file(file: TextIO) -> Dict[str, str]:
    """
    """
    tr2curr_best: Dict[str, Tuple[str, int, int, str]] = {}
    for line in file:
        data: List[str] = line.rstrip().split('\t')
        if not data or not data[0]:
            continue
        if data[0] == 'projection':
            continue
        proj: str = data[0]
        tr: str = '#'.join(proj.split('#')[:-1])
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
    tr2best: Dict[str, str] = {k: v[0] for k,v in tr2curr_best.items()}
    # print(f'{tr2best["ENST00000409539.INMT"]=}')
    # print(f'{tr2best["ENST00000013222.INMT"]=}')
    return tr2best

def transcript_meta_to_report(
    file: Union[str, TextIO],
    precedence: Optional[Dict[str, str]] = {},
    paralogs: Optional[Set[str]] = set(),
    ppgenes: Optional[Set[str]] = set()
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Parses transcript_meta.tsv output file of TOGA, inferring projection
    classification and best transcript classification
    """
    proj2status: Dict[str, str] = {}
    tr2status: Dict[str, Set[str]] = defaultdict(set)
    if not isinstance(file, str):
        lines: List[str] = file.readlines()
    else:
        with open(file, 'r') as h:
            lines: List[str] = h.readlines()
    for line in lines:
        data: List[str] = line.strip().split('\t')
        if data[0] == 'projection':
            continue
        proj: str = data[0]
        tr: str = '#'.join(proj.split('#')[:-1])
        status: str = data[1]
        proj2status[proj] = status
        ## do not propagate paralog and ppgene loss status to downstream levels
        if proj in ppgenes:
            continue
        if proj in paralogs:
            status = 'PG'
        tr2status[tr].add(status)
    for tr in tr2status:
        if tr in precedence:
            preferred_proj: str = precedence[tr]
            if preferred_proj in proj2status and preferred_proj not in ppgenes:
                preferred_loss_status: str = proj2status[preferred_proj]
                tr2status[tr] = preferred_loss_status
                continue
        all_classes: Set[str] = tr2status[tr]
        max_status: str = max(all_classes, key=lambda x: CLASS_TO_NUM[x])
        tr2status[tr] = max_status

    return (proj2status, tr2status)


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
    file: Union[str, TextIO]
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Parses the genes_rejection_reason.tsv file containing transcripts which were
    rejected for one reason or another at the CESAR alignment step or upstream
    """
    proj2status: Dict[str, str] = {}
    tr2all_statuses: Dict[str, Set[str]] = defaultdict(set)
    tr2status: Dict[str, str] = {}
    if not isinstance(file, str):
        lines: List[str] = file.readlines()
    else:
        with open(file, 'r') as h:
            lines: List[str] = h.readlines()
    for line in lines:
        data: List[str] = line.rstrip().split('\t')
        if not data or not data[0]:
            continue
        if data[0] == 'level':
            continue
        if len(data) == 2: ## TEMPORARY SOLUTION TO BYPASS THE UPSTREAM BUG
            name: str = data[0]
            status: str = 'N'
            tr2status[name] = status
            continue
        level: str = data[0]
        name: str = data[1]
        status: str = data[5]
        if level == TRANSCRIPT:
            tr2status[name] = status
            continue
        if level != PROJECTION:
            raise ValueError(f'An ambiguous entry found: {line}')
        tr: str = '#'.join(name.split('#')[:-1])
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
    precedence: Optional[Dict[str, str]] = {}
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Updates projection- and transcript-level loss status dictionaries with
    the data for the beforehand-rejected projections. Adds the new instance
    if no instance of this projection or transcript was subjected to alignment,
    leaves the best loss status otherwise
    """
    for proj, rej_proj_status in reject_proj2status.items():
        tr: str = '#'.join(proj.split('#')[:-1])
        if proj not in orig_proj2status:
            orig_proj2status[proj] = rej_proj_status
        else:
            orig_proj_status: str = orig_proj2status[proj]
            orig_proj2status[proj] = max(
                (orig_proj_status, rej_proj_status), key=lambda x: CLASS_TO_NUM[x]
            )
        if tr in precedence:
            preferred_proj: str = precedence[tr]
            if preferred_proj in orig_proj2status:
                preferred_loss_status: str = orig_proj2status[preferred_proj]
                orig_tr2status[tr] = preferred_loss_status
            elif preferred_proj in reject_proj2status:
                preferred_loss_status: str = reject_proj2status[preferred_proj]
                orig_tr2status[tr] = preferred_loss_status
            else:
                raise KeyError(
                    f'Top-precedent projection {preferred_proj} is absent '
                    'from both transcript meta and rejection reports'
                )
            continue
        rej_tr_status: str = reject_tr2status[tr]
        if tr not in orig_tr2status:
            orig_tr2status[tr] = rej_tr_status
        else:
            orig_tr_status: str = orig_tr2status[tr]
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
        with open(file, 'r') as h:
            lines: List[str] = h.readlines()
    for line in lines:
        if line.startswith('GeneID'):
            continue
        gene, isoform = line.strip().split('\t')
        gene2status[gene].add(tr_report.get(isoform, 'N'))
    for gene in gene2status:
        all_classes: Set[str] = gene2status[gene]
        max_class: str = max(all_classes, key=lambda x: CLASS_TO_NUM[x])
        gene2status[gene] = max_class

    return gene2status


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'transcript_meta',
    type=click.File('r', lazy=True),
    metavar='TRANSCRIPT_META'
)
@click.option(
    '--rejected_projections',
    '-r',
    type=click.File('r', lazy=True),
    metavar='TSV',
    default=None,
    show_default=True,
    help=(
        'A six-column file containing projections discarded prior '
        'to actual CESAR alignment and respective rejection reasons'
    )
)
@click.option(
    '--isoform_file',
    '-i',
    type=click.File('r', lazy=True),
    default=None,
    show_default=True,
    metavar='TSV',
    help='A path to two-column gene-to-isoform mapping file'
)
@click.option(
    '--paralogs',
    '-p',
    type=click.File('r', lazy=True),
    default=None,
    show_default=True,
    help=(
        'A path to a single-column file containing paralogous projections\' names'
    )
)
@click.option(
    '--processed_pseudogenes',
    '-pp',
    type=click.File('r', lazy=True),
    default=None,
    show_default=True,
    help=(
        'A path to a single-column file containing projection names of processed pseudogenes'
    )
)
@click.option(
    '--spanning_chains_precedence_file',
    '-span',
    type=click.File('r', lazy=True),
    metavar='SPANNING_CHAIN_FILE',
    default=None,
    show_default=True,
    help=(
        'A file containing spanning chain gap coordinates in the reference. '
        'Nestedness of spanning gaps is used to prioritise chains for loss inference'
    )
)
@click.option(
    '--output',
    '-o',
    type=click.File('w'),
    default=stdout,
    show_default=False,
    help='A path to store the output at [default: stdout]'
)

def main(
    transcript_meta: click.File,
    rejected_projections: Optional[click.File],
    isoform_file: Optional[Union[click.File, None]],
    paralogs: Optional[Union[click.File, None]],
    processed_pseudogenes: Optional[Union[click.File, None]],
    spanning_chains_precedence_file: Optional[Union[click.File, None]],
    output: Optional[click.File]
) -> None:
    """
    Summarises projection classification data based on TOGA results. Arguments are:\n
    \tTRANSCRIPT_META is a projection metadata output TSV file of the CESAR wrapper
    (transcript_meta.tsv for CESAR wrapper or combined_projection_meta.tsv for TOGA)\n
    """
    if spanning_chains_precedence_file:
        spanning_status: Dict[str, str] = parse_precedence_file(spanning_chains_precedence_file)
    else:
        spanning_status: Dict[str, str] = {}
    paralogs: Set[str] = parse_single_column(paralogs)
    ppgenes: Set[str] = parse_single_column(processed_pseudogenes)
    proj2status, tr2status = transcript_meta_to_report(
        transcript_meta, 
        spanning_status, 
        paralogs=paralogs, 
        ppgenes=ppgenes
    )
    if rejected_projections:
        r_proj2status, r_tr2status = rejection_file_to_report(rejected_projections)
        proj2status, tr2status = add_rejection_data(
            proj2status, tr2status, r_proj2status, r_tr2status, spanning_status
        )
    output.write(Headers.LOSS_FILE_HEADER)
    for proj, status in proj2status.items():
        output.write('\t'.join((PROJECTION, proj, status)) + '\n')
    for tr, status in tr2status.items():
        output.write('\t'.join((TRANSCRIPT, tr, status)) + '\n')
    ## TODO: Make transcripts appearing in the isoform_file but not in the TOGA results be listed in output as N/M
    if isoform_file is not None:
        gene2status: Dict[str, str] = gene_loss_report(isoform_file, tr2status)
        for gene, status in gene2status.items():
            output.write('\t'.join((GENE, gene, status)) + '\n')

if __name__ == '__main__':
    main()
