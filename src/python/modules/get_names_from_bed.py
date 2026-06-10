#!/usr/bin/env python3
"""Just extract names from toga output bed file.

Works like xenoRefGenelx.pl"""

from typing import Set
import click
import sys
from shared import CommandLineManager, CONTEXT_SETTINGS, read_tab

# from version import __version__

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.option(
    "--input",
    "-i",
    type=click.File("r", lazy=True),
    metavar="INPUT",
    default=sys.stdin,
    show_default=False,
    help="Path to input BED file [default: stdin]"
)
@click.option(
    "--output",
    "-o",
    type=click.File("w", lazy=True),
    metavar="OUTPUT",
    default=sys.stdout,
    show_default=False,
    help="Path to output index file [default: stdout]"
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

class BedNameRetriever(CommandLineManager):
    """
    Extracts names from the BED3+ file for UCSC index creation. Imitates xenoRefGenelx.pl
    """
    __slots__ = ("v",)

    def __init__(
        self, 
        input: click.File, 
        output: click.File,
        log_name: str,
        verbose: bool,
    ) -> None:
        self.v: bool = verbose
        self.set_logging(name=log_name, toga_module="bed_name_extractor")
        self._to_log("Starting BedNameRetriever")

        out_set: Set[str] = set()
        for i, data in enumerate(read_tab(input), start=1):
        # with open(input, "r") as h:
        # for i, line in enumerate(input, start=1):
            # data: List[str] = line.strip().split("\t")
            # if not data or not data[0]:
            #     continue
            if len(data) < 4:
                self._die(
                    (
                        "Improper number of fields in the BED file at line %i; "
                        "expected at least 4 fields, got %i"
                    ) % (i, len(data))
                )
            id_field: str = data[3]
            no_chain_id = "#".join(id_field.split("#")[:-1])
            dot_split = no_chain_id.split("#")
            if len(dot_split) > 1:
                to_out = [id_field, no_chain_id] + dot_split
            else:
                to_out = [id_field, no_chain_id]
            line = "\t".join(to_out)
            # output.write(line)
            out_set.add(line)
        for line in sorted(out_set):
            output.write(line)


if __name__ == "__main__":
    BedNameRetriever()


# if len(sys.argv) != 2:
#     to_read = None
#     sys.exit(f"Usage: {sys.argv[0]} [query_annotation.bed] | sort -u > ix.txt")
# else:
#     to_read = sys.argv[1]

# f = open(to_read, "r")
# for line in f:
#     id_field = line.rstrip().split("\t")[3]
#     no_chain_id = "#".join(id_field.split("#")[:-1])
#     dot_split = no_chain_id.split("#")
#     if len(dot_split) > 1:
#         to_out = [id_field, no_chain_id] + dot_split
#     else:
#         to_out = [id_field, no_chain_id]
#     line = "\t".join(to_out)
#     print(line)
# f.close()
