#!/usr/bin/env python3

"""
Creates a provisional isoform file for BED-to-GTF conversion
"""

from collections import defaultdict
from os import PathLike
from typing import Dict, List, Optional, Union

import click

from .shared import CONTEXT_SETTINGS, CommandLineManager, segment_base

HEADER: str = "query_gene"


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument("isoform_file", type=click.Path(exists=True), metavar="IN_ISOFORM_FILE")
@click.argument("bed_file", type=click.Path(exists=True), metavar="FINAL_BED_FILE")
@click.argument("output", type=click.Path(exists=False), metavar="OUT_ISOFORM_FILE")
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
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls execution verbosity",
)
class ProvisionalIsoformMapper(CommandLineManager):
    """
    An auxiliary script used for GTF file preparation by creating a modified temporary version
    of the query isoform file ('query_genes.tsv' by default). The major change is adding all the fragments of fragmented projections
    to the isoform mapping.

    \b
    Arguments are:
    * IN_ISOFORM_FILE is a final isoform file produced by TOGA2 ('query_genes.tsv' by default);
    * FINAL_BED_FILE is a BED12 file produced by TOGA2 for which the GTF conversion is requested (either 'query_annotation.bed' or 'query_annotation.with_utrs.bed');
    * OUT_ISOFORM_FILE is a path to modified isoform file
    """

    def __init__(
        self,
        isoform_file: Union[str, PathLike],
        bed_file: Union[str, PathLike],
        output: Union[str, PathLike],
        log_name: Optional[Union[str, None]] = None,
        verbose: Optional[bool] = False,
    ) -> None:
        self.v: bool = verbose
        self.set_logging(name=log_name, toga_module="isoforms_for_gtf")
        self.run(isoform_file, bed_file, output)

    def run(
        self,
        isoform_file: Union[str, PathLike],
        bed_file: Union[str, PathLike],
        output: Union[str, PathLike],
    ) -> None:
        """Entry point function"""

        ## record the fragmented projections' names
        fragmented_names: Dict[str, int] = defaultdict(int)
        with open(bed_file, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) != 12:
                    self._die(
                        (
                            "Improper formatting at BED file line %i; "
                            "expected 12 fields, got %i"
                        )
                        % (i, len(data))
                    )
                name: str = data[3]
                if "$" not in name:
                    continue
                basename: str = segment_base(name)
                fragmented_names[basename] += 1

        ## read the isoform file, check each line in it;
        ## fragmented projection mapping is to be expanded to accommodate for each fragment,
        ## the rest of the lines are written as-is
        with open(isoform_file, "r") as ih, open(output, "w") as oh:
            for i, line in enumerate(ih, start=1):
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) != 2:
                    self._die(
                        (
                            "Improper formatting at isoform file line %i; "
                            "expected 2 fields, got %i"
                        )
                        % (i, len(data))
                    )
                gene, proj = data
                if gene == HEADER:
                    continue
                if proj not in fragmented_names:
                    oh.write(line)
                    continue
                for num in range(1, fragmented_names[proj] + 1):
                    fragm_name: str = f"{proj}${num}"
                    oh.write(gene + "\t" + fragm_name + "\n")


if __name__ == "__main__":
    ProvisionalIsoformMapper()
