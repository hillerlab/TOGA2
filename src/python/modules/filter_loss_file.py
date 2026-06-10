#!/usr/bin/env python3

""" """

from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple, Union

import click

from .cesar_wrapper_constants import CLASS_TO_NUM
from .constants import Headers
from .shared import (
    CONTEXT_SETTINGS,
    CommandLineManager,
    base_proj_name,
    get_proj2trans,
    parse_single_column,
)

HEADER: str = "level"
PROJECTION: str = "PROJECTION"
TRANSCRIPT: str = "TRANSCRIPT"
GENE: str = "GENE"
MISSING_STATS: Tuple[str, ...] = ("N", "M", "L")
MISSING: str = "M"

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    "init_loss_summary", type=click.File("r", lazy=True), metavar="INIT_LOSS_SUMMARY"
)
@click.argument(
    "query_bed_file", type=click.File("r", lazy=True), metavar="FINAL_BED_FILE"
)
@click.argument(
    "final_loss_summary", type=click.File("w", lazy=True), metavar="FINAL_LOSS_SUMMARY"
)
@click.option(
    "--paralogs",
    "-p",
    type=click.File("r", lazy=True),
    metavar="PARALOG_FILE",
    default=None,
    show_default=True,
    help="A single-column file containing paralogous projections' identifiers",
)
@click.option(
    "--processed_pseudogenes",
    "-pp",
    type=click.File("r", lazy=True),
    metavar="PROCESSED_PSEUDOGENE_FILE",
    default=None,
    show_default=True,
    help="A single-column file containing processed pseudogene projections' identifiers",
)
@click.option(
    "--isoform_file",
    "-i",
    type=click.Path(exists=True),
    metavar="ISOFORM_FILE",
    default=None,
    show_default=True,
    help="A two-column file containing reference gene-to-transcript mapping"
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
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls execution verbosity",
)
class LossFileFilter(CommandLineManager):
    """
    A straightforward loss summary filtering class. Given a full loss summary file
    and a final annotation BED file, picks loss records only for projections present
    in the annotation. Transcript- and gene-level records are not filtered and reported
    as in the original file.

    \b
    Arguments are:
    \t* INIT_LOSS_SUMMARY is a full loss summary file produced by TOGA2 (`meta/loss_summary_extended.tsv` by default);
    \t* FINAL_BED_FILE is a final query annotation file in BED12 format (`query_annotation.tsv` or `query_annotation.with_utrs.bed`);
    \t* FINAL_LOSS_SUMMARY is a path to the final, filtered loss summary file (`loss_summary.tsv` in TOGA2 by default)
    """

    __slots__ = ("loss_in", "bed_file", "loss_out", "paralogs", "ppgenes", "isoform_file")

    def __init__(
        self,
        init_loss_summary: click.File,
        query_bed_file: click.File,
        final_loss_summary: click.File,
        paralogs: Optional[click.File],
        processed_pseudogenes: Optional[click.File],
        isoform_file: Optional[click.Path],
        log_name: Optional[str],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)
        self.loss_in: click.File = init_loss_summary
        self.bed_file: click.File = query_bed_file
        self.loss_out: click.File = final_loss_summary
        self.paralogs: Set[str] = parse_single_column(paralogs)
        self.ppgenes: Set[str] = parse_single_column(processed_pseudogenes)
        self.isoform_file: Union[str, None] = isoform_file

        self.run()

    def run(self) -> None:
        """Entry point"""
        projections: Set[str] = set()
        for i, line in enumerate(self.bed_file, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 12:
                self._die(
                    (
                        "Improper final Bed file formatting at line %i; "
                        "expected 12 fields, got %i"
                    )
                    % (i, len(data))
                )
            name: str = base_proj_name(data[3])
            projections.add(name)
        gene2trs: Dict[str, List[str]] = defaultdict(list)
        if self.isoform_file is not None:
            with open(self.isoform_file, "r") as h:
                for line in h:
                    data: List[str] = line.strip().split("\t")
                    if not data or not data[0]:
                        continue
                    gene: str = data[0]
                    tr: str = data[1]
                    gene2trs[gene].append(tr)
        tr2loss: Dict[str, str] = {}
        gene2loss: Dict[str, str] = {}
        tr2best_loss: Dict[str, str] = {}
        self.loss_out.write(Headers.LOSS_FILE_HEADER)
        for i, line in enumerate(self.loss_in, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 3:
                self._die(
                    (
                        "Improper loss file formatting at line %i; "
                        "expected 3 fields, got %i"
                    )
                    % (i, len(data))
                )
            if data[0] == HEADER:
                continue
            # if data[0] != PROJECTION:
            #     self.loss_out.write(line)
            #     continue
            ## filter the projection lines, leaving only the ones 
            ## corresponding to the final output list
            if data[0] == PROJECTION:
                if (
                    data[1] not in projections
                    and base_proj_name(data[1]) not in projections
                ):
                    continue
                tr: str = get_proj2trans(data[1])[0]
                if data[1] in self.paralogs:
                    data[1] += "#paralog"
                    line = "\t".join(data) + "\n"
                    tr2best_loss[tr] = "PG"
                elif data[1] in self.ppgenes:
                    data[1] += "#retro"
                    line = "\t".join(data) + "\n"
                    tr2best_loss[tr] = max(
                        (tr2best_loss.get(tr, "N"), "M"), key=lambda x: CLASS_TO_NUM[x]
                    )
                else:
                    tr2best_loss[tr] = max(
                        (tr2best_loss.get(tr, "N"), data[2]), key=lambda x: CLASS_TO_NUM[x]
                    )
                self.loss_out.write(line)
            ## for the other two levels, save their reported loss statuses
            elif data[0] == TRANSCRIPT:
                tr2loss[data[1]] = data[2]
            elif data[0] == GENE:
                gene2loss[data[1]] = data[2]
        ## for transcripts, check if they still have any projections
        self._to_log("Double-checking the transcript level loss classification")
        for tr, status in tr2loss.items():
            ## the transcript has any orhtologous projections reported: proceed
            if tr in tr2best_loss:
                status: str = tr2best_loss[tr]
                tr2loss[tr] = status
                self.loss_out.write(f"{TRANSCRIPT}\t{tr}\t{status}\n")
                continue
            if status in MISSING_STATS:
                self.loss_out.write(f"{TRANSCRIPT}\t{tr}\t{status}\n")
                continue
            status = MISSING
            tr2loss[tr] = status
            self.loss_out.write(f"{TRANSCRIPT}\t{tr}\t{status}\n")
        self._to_log("Double-checking the gene level loss classification")
        for gene, status in gene2loss.items():
            if any(x in tr2best_loss for x in gene2trs.get(gene, [])):
                best_losses: List[Union[str, None]] = [
                    tr2best_loss.get(x) for x in gene2trs.get(gene, []) if x in tr2best_loss
                ]
                if best_losses:
                    status: str = max(best_losses,  key=lambda x: CLASS_TO_NUM[x])
                self.loss_out.write(f"{GENE}\t{gene}\t{status}\n")
                continue
            if status in MISSING_STATS:
                self.loss_out.write(f"{GENE}\t{gene}\t{status}\n")
                continue
            all_trs: List[str] = gene2trs.get(gene, [])
            if not all_trs:
                self._die("No transcripts found for gene %s" % gene)
            all_statuses: List[str] = []
            for transcript in all_trs:
                tr_stat: Union[str, None] = tr2loss.get(transcript)
                if tr_stat is None:
                    self._to_log(
                        "No loss status recorded for transcript %s" % transcript, "warning"
                    )
                    continue
                all_statuses.append(tr_stat)
            max_tr_stat: str = max(all_statuses, key=lambda x: CLASS_TO_NUM[x])
            self.loss_out.write(f"{GENE}\t{gene}\t{max_tr_stat}\n")



if __name__ == "__main__":
    LossFileFilter()
