"""Filters the final loss files from discarded items"""

from collections import defaultdict
from typing import Final

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

HEADER: Final[str] = "level"
PROJECTION: Final[str] = "PROJECTION"
TRANSCRIPT: Final[str] = "TRANSCRIPT"
GENE: Final[str] = "GENE"
MISSING_STATS: Final[tuple[str, ...]] = ("N", "M", "L")
MISSING: Final[str] = "M"

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

    __slots__ = ("bed_file", "isoform_file", "loss_in", "loss_out", "paralogs", "ppgenes")

    def __init__(
        self,
        init_loss_summary: click.File,
        query_bed_file: click.File,
        final_loss_summary: click.File,
        paralogs: click.File | None = None,
        processed_pseudogenes: click.File | None = None,
        isoform_file: click.Path | None = None,
        log_name: str | None = None,
        verbose: bool = False,
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)
        self.loss_in: click.File = init_loss_summary
        self.bed_file: click.File = query_bed_file
        self.loss_out: click.File = final_loss_summary
        self.paralogs: set[str] = parse_single_column(paralogs)
        self.ppgenes: set[str] = parse_single_column(processed_pseudogenes)
        self.isoform_file: str | None = isoform_file

        self.run()

    def run(self) -> None:
        """Entry point"""
        projections: set[str] = set()
        for i, line in enumerate(self.bed_file, start=1):
            data: list[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 12:
                self._die(
                    f"Improper final Bed file formatting at line {i}; "
                    f"expected 12 fields, got {len(data)}"
                )
            name: str = base_proj_name(data[3])
            projections.add(name)
        gene2trs: dict[str, list[str]] = defaultdict(list)
        if self.isoform_file is not None:
            with open(self.isoform_file, "r") as h:
                for line in h:
                    data: list[str] = line.strip().split("\t")
                    if not data or not data[0]:
                        continue
                    gene: str = data[0]
                    tr: str = data[1]
                    gene2trs[gene].append(tr)
        tr2loss: dict[str, str] = {}
        gene2loss: dict[str, str] = {}
        tr2best_loss: dict[str, str] = {}
        self.loss_out.write(Headers.LOSS_FILE_HEADER)
        for i, line in enumerate(self.loss_in, start=1):
            data: list[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 3:
                self._die(
                    f"Improper loss file formatting at line {i}; "
                    f"expected 3 fields, got {len(data)}"
                )
            if data[0] == HEADER:
                continue
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
                    tr2best_loss[tr] = max(
                        (tr2best_loss.get(tr, "N"), "PG"), key=lambda x: CLASS_TO_NUM[x]
                    )
                    # tr2best_loss[tr] = "PG"
                elif data[1] in self.ppgenes:
                    data[1] += "#retro"
                    line = "\t".join(data) + "\n"
                    tr2best_loss[tr] = max(
                        (tr2best_loss.get(tr, "N"), "PP"), key=lambda x: CLASS_TO_NUM[x]
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
                if status == "PP":
                    status = "M"
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
                best_losses: list[str, None] = [
                    tr2best_loss.get(x) for x in gene2trs.get(gene, []) if x in tr2best_loss
                ]
                if best_losses:
                    status: str = max(best_losses,  key=lambda x: CLASS_TO_NUM[x])
                self.loss_out.write(f"{GENE}\t{gene}\t{status}\n")
                continue
            if status in MISSING_STATS:
                self.loss_out.write(f"{GENE}\t{gene}\t{status}\n")
                continue
            all_trs: list[str] = gene2trs.get(gene, [])
            if not all_trs:
                self._die(f"No transcripts found for gene {gene}")
            all_statuses: list[str] = []
            for transcript in all_trs:
                tr_stat: str | None = tr2loss.get(transcript)
                if tr_stat is None:
                    self._to_log(
                        f"No loss status recorded for transcript {transcript}", "warning"
                    )
                    continue
                all_statuses.append(tr_stat)
            max_tr_stat: str = max(all_statuses, key=lambda x: CLASS_TO_NUM[x])
            self.loss_out.write(f"{GENE}\t{gene}\t{max_tr_stat}\n")



if __name__ == "__main__":
    LossFileFilter()
