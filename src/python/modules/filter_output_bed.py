#!/usr/bin/env python3

"""
Filters TOGA2 output BED files and adds the necessary postfixes
"""

from typing import List, Optional, Set, TextIO, Union

import click

from .shared import (
    CONTEXT_SETTINGS,
    CommandLineManager,
    base_proj_name,
    parse_single_column,
)

FI_COL: str = "0,0,100"
I_COL: str = "0,0,200"


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument("input_bed", type=click.File("r", lazy=True), metavar="INPUT_BED")
@click.argument("output_bed", type=click.File("w", lazy=True), metavar="OUTPUT_BED")
@click.option(
    "--deprecated_projection_list",
    "-di",
    type=click.File("r", lazy=True),
    metavar="PROJECTION_LIST_FILE",
    default=None,
    show_default=True,
    help=(
        "A single-column file containing projections to filter out from the input file"
    ),
)
@click.option(
    "--deprecated_projection_bed",
    "-do",
    type=click.File("w", lazy=False),
    metavar="DISCARDED_PROJECTION_BED",
    default=None,
    show_default=True,
    help=("A path to write BED entries for deprecated projections to"),
)
@click.option(
    "--paralog_list",
    "-pi",
    type=click.File("r", lazy=True),
    metavar="PARALOG_LIST_FILE",
    default=None,
    show_default=True,
    help=(
        "A single-column list of paralogous projections "
        "fo filter out from the input file"
    ),
)
@click.option(
    "--processed_pseudogene_list",
    "-ppi",
    type=click.File("r", lazy=True),
    metavar="PSEUDOGENE_LIST_FILE",
    default=None,
    show_default=True,
    help=(
        "A single-column list of processed pseudogene projections "
        "fo filter out from the input file"
    ),
)
@click.option(
    "--processed_pseudogene_bed",
    "-ppo",
    type=click.File("w", lazy=True),
    metavar="PSEUDOGENE_PROJECTION_BED",
    default=None,
    show_default=True,
    help=("A path to write BED entries for processed pseudogene projections to"),
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
class OutputBedFilter(CommandLineManager):
    __slots__ = (
        "discarded_projs",
        "discarded_bed_lines",
        "discarded_file",
        "paralog_projs",
        "ppgene_projs",
        "ppgene_bed_lines",
        "ppgene_file",
    )

    def __init__(
        self,
        input_bed: click.File,
        output_bed: click.File,
        deprecated_projection_list: Optional[Union[click.File, None]],
        deprecated_projection_bed: Optional[Union[click.File, None]],
        paralog_list: Optional[Union[click.File, None]],
        processed_pseudogene_list: Optional[Union[click.File, None]],
        processed_pseudogene_bed: Optional[Union[click.File, None]],
        log_name: Optional[str],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.set_logging(name=log_name, toga_module="query_fasta_filter")

        self.discarded_projs: Set[str] = parse_single_column(deprecated_projection_list)
        self.discarded_bed_lines: List[str] = []
        self.discarded_file: Union[click.File, None] = deprecated_projection_bed

        self.paralog_projs: Set[str] = parse_single_column(paralog_list)
        self.ppgene_projs: Set[str] = parse_single_column(processed_pseudogene_list)
        self.ppgene_bed_lines: List[str] = []
        self.ppgene_file: Union[click.File, None] = processed_pseudogene_bed

        self.parse_and_filter_bed(input_bed, output_bed)
        self.write_discarded_projections()
        self.write_pseudogenes()

    def parse_and_filter_bed(self, input_: TextIO, output: TextIO) -> None:
        """
        Reads input BED file, filters out deprecated projections,
        and writes the filtered output to another file
        """
        for line in input_:
            data: List[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            name: str = data[3]
            basename: str = base_proj_name(name)
            if name in self.discarded_projs or basename in self.discarded_projs:
                self.discarded_bed_lines.append(line)
                continue
            if name in self.ppgene_projs or basename in self.ppgene_projs:
                self.ppgene_bed_lines.append(line)
                if data[8] not in (FI_COL, I_COL):
                    continue
                if "#retro" not in data[3]:
                    data[3] = f"{data[3]}#retro"
                line = "\t".join(data) + "\n"
            if name in self.paralog_projs or basename in self.paralog_projs:
                if "#paralog" not in data[3]:
                    data[3] = f"{data[3]}#paralog"
                line = "\t".join(data) + "\n"
            output.write(line)

    def write_discarded_projections(self) -> None:
        """Writes BED lines for projections deemed as deprecated"""
        if not self.discarded_bed_lines:
            return
        if self.discarded_file is None:
            return
        for line in self.discarded_bed_lines:
            self.discarded_file.write(line)

    def write_pseudogenes(self) -> None:
        """Writes BED lines for processed pseudogene projections"""
        if not self.ppgene_bed_lines:
            return
        if self.ppgene_file is None:
            return
        for line in self.ppgene_bed_lines:
            self.ppgene_file.write(line)


if __name__ == "__main__":
    OutputBedFilter()
