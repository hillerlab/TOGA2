#!/usr/bin/env python3

"""
Shared functionality across the scripts
"""

import ctypes
import fcntl
import logging
import os
import subprocess
import sys
import time
from datetime import datetime
from logging import Formatter
from shutil import copy2, move, rmtree
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Set,
    TextIO,
    Tuple,
    Union,
)

import click
import click_option_group
import gzip
import networkx as nx
from click_option_group import OptionGroup

## Constants
CONTEXT_SETTINGS: Dict[str, Any] = {
    "help_option_names": ["-h", "-help", "--help"],
    "ignore_unknown_options": True,
    "allow_extra_args": True,
    "max_content_width": 150,
}
## borrowed from splitFile
SPLIT_JOB_HEADER: Tuple[str, str, str] = ("#!/bin/bash", "set -eu", "set -o pipefail")
SLIB_NAME = "chain_bst_lib.so"
UTF8: str = "utf-8"
FORMATTER: Formatter = Formatter(
    # "[{asctime}][{filename}] - {levelname}: {message}",
    "[{asctime}][{toga_module}] - {levelname}: {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M:%S",
)

## Sequence handling & score calculation data
COMPLEMENT: Dict[str, str] = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "n": "n",
    "-": "-",
}

## Types
Numeric = Union[float, int]


class PrettyGroup(OptionGroup):
    def get_help_record(self, ctx: click.Context) -> Optional[Tuple[str, str]]:
        init_help: Union[Tuple[str, str], None] = super().get_help_record(ctx)
        if init_help is None:
            return None
        return "\n" + init_help[0], init_help[1]


class DependentOption(click_option_group.GroupedOption):
    """
    An expansion of the original GroupedOption class with mutex and dependency functionlities.
    Mutex functionality is based on the solution from https://github.com/pallets/click/issues/257
    """

    def __init__(self, *args, **kwargs) -> None:
        ## mutex functionality: deprecate the joint use of mutually exclusive options
        self.competes_with: Union[List[str], None] = kwargs.pop("competes_with", None)
        self.required_mutex: bool = kwargs.pop("required_mutex", False)
        ## dependence functionality: enforce the use of the option in tandem with its dependencies
        self.requires: Union[List[str], None] = kwargs.pop("requires", None)
        ## 'soft dependence' functionality: enforce the use of the option unless alternative is provided
        self.not_required_if: Union[List[str], None] = kwargs.pop(
            "not_required_if", None
        )
        # if kwargs["help"] and self.competes_with is not None:
        if self.competes_with or self.not_required_if or self.requires:
            kwargs["help"] = "" if not kwargs["help"] else kwargs["help"] + ". "
            if self.competes_with is not None:
                kwargs["help"] += (
                    "Mutually exclusive with the following options: "
                    + ", ".join(self.competes_with)
                    + ". "
                )
            if self.requires is not None:
                kwargs["help"] += (
                    "Requires the following options: " + ", ".join(self.requires) + ". "
                )
            if self.not_required_if is not None:
                kwargs["help"] += (
                    "Not required if the following options are provided: "
                    + ", ".join(self.not_required_if)
                    + ". "
                )
        super(DependentOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args) -> None:
        current_opt: bool = self.name in opts
        alternative_found: bool = False
        if self.competes_with:
            for mutex_opt in self.competes_with:
                if mutex_opt in opts:
                    if current_opt:
                        raise click.UsageError(
                            "Options %s and %s are mutually exclusive"
                            % (self.name, mutex_opt)
                        )
                    else:
                        self.prompt = None
                    alternative_found = True
        if self.not_required_if:
            for alternative in self.not_required_if:
                if alternative in opts:
                    alternative_found = True
                    self.required = False
                    self.prompt = None
                    break
        if self.requires and current_opt:
            for requirement in self.requires:
                if requirement not in opts:
                    raise click.UsageError(
                        "Option %s requires option %s, which is missing"
                        % (self.name, requirement)
                    )
        if not current_opt and not alternative_found and self.required_mutex:
            raise click.UsageError(
                "One of the following options is required for TOGA2 execution: %s"
                % (self.name + ", " + ",".join(self.competes_with))
            )
        return super(DependentOption, self).handle_parse_result(ctx, opts, args)


# class MutexOption(click_option_group.GroupedOption):
#     """
#     Mutually exclusive Click option class.
#     Based on the solution from https://github.com/pallets/click/issues/257
#     """

#     def __init__(self, *args, **kwargs) -> None:
#         if "competes_with" not in kwargs:
#             raise click.UsageError(
#                 "Argument \"competes_with=\" mandatory for MutexOption class "
#                 "is missing for option %s"
#                 % self.name
#             )
#         self.competes_with: List[str] = kwargs.pop("competes_with")
#         if self.competes_with is None or not self.competes_with:
#             raise click.UsageError(
#                 "No competing options for an MutexOption instance %s, "
#                 "with no defaults. Please provide a list of competing options "
#                 'with "competes_with" argument'
#                 % self.name
#             )
#         self.required_mutex: bool = (
#             kwargs.pop("required_mutex") if "required_mutex" in kwargs else False
#         )
#         self.not_required_if: Union[List[str], None] = (
#             kwargs.pop("not_required_if") if "not_required_if" in kwargs else None
#         )
#         kwargs["help"] = (
#             "" if kwargs.get("help") is None else (kwargs.get("help", "") + ". ") +
#             "Mutually exclusive with the following options: " +
#             ", ".join(self.competes_with) +
#             "."
#         )
#         super(MutexOption, self).__init__(*args, **kwargs)

#     def handle_parse_result(self, ctx, opts, args) -> None:
#         current_opt: bool = self.name in opts
#         alternative_found: bool = False
#         for mutex_opt in self.competes_with:
#             if mutex_opt in opts:
#                 if current_opt:
#                     raise click.UsageError(
#                         "Options %s and %s are mutually exclusive"
#                         % (self.name, mutex_opt)
#                     )
#                 else:
#                     self.prompt = None
#                 alternative_found = True
#         if self.not_required_if is not None:
#             for alternative in self.not_required_if:
#                 if alternative in opts:
#                     self.prompt = None
#                     alternative_found = True
#                     break
#         if not current_opt and not alternative_found and self.required_mutex:
#             raise click.UsageError(
#                 "One of the following options is required for TOGA2 execution: %s"
#                 % (self.name + ", " + ",".join(self.competes_with))
#             )
#         return super(MutexOption, self).handle_parse_result(ctx, opts, args)


# class DependentOption(click_option_group.GroupedOption):
#     """
#     A grouped option that requires one or more other options for the client to start.
#     A rework of the MutexClass above.
#     """
#     def __init__(self, *args, **kwargs) -> None:
#         if "requires" not in kwargs:
#             raise click.UsageError(
#                 "Argument \"requires=\" mandatory for DependentOption class "
#                 "is missing for option %s"
#                 % kwargs["name"]
#             )
#         self.requires: List[str] = kwargs.pop("requires")
#         if self.requires is None or not self.requires:
#             raise click.UsageError(
#                 "No required options for an DependentOption instance %s, "
#                 "with no defaults. Please provide a list of required options "
#                 'with "requires" argument'
#                 % kwargs["name"]
#             )
#         kwargs["help"] = (
#             "" if kwargs.get("help") is None else (kwargs.get("help", "") + ". ") +
#             "Requires the following options: " +
#             ", ".join(self.requires) +
#             "."
#         )
#         super(DependentOption, self).__init__(*args, **kwargs)

#     def handle_parse_result(self, ctx, opts, args) -> None:
#         if self.name in opts:
#             for requirement in self.requires:
#                 if requirement not in opts:
#                     raise click.UsageError(
#                         "Option %s requires option %s, which is missing"
#                         % (self.name, requirement)
#                     )
#         return super(MutexOption, self).handle_parse_result(ctx, opts, args)


class TogaDirConfig:
    ## NOTE: Slots are currently deprecated due to Postoga using vars() when logging starting args
    ## TO be changed after the respective change is pushed to Postoga
    # __slots__ = (
    #     "togadir",
    #     "orthology_class",
    #     "loss_status",
    #     "orthology_score",
    #     "to",
    #     "bed_type",
    #     "min_paralog_score",
    #     "with_isoforms",
    #     "extract",
    #     "outdir",
    #     "only_table",
    #     "only_convert",
    #     "log_level",
    #     "depure",
    # )

    def __init__(
        self,
        togadir: os.PathLike,
        by_loss_status: Optional[Union[str, None]] = None,
        by_orthology_class: Optional[Union[str, None]] = None,
        by_orthology_score: Optional[Union[float, None]] = None,
        to: Optional[str] = "bed",
        target: Optional[str] = "utr",
        by_paralog_score: Optional[Union[float, None]] = None,
        with_isoforms: Optional[Union[os.PathLike, None]] = None,
        extract: Optional[Union[bool, str]] = False,
        outdir: Optional[Union[os.PathLike, None]] = None,
        only_table: Optional[bool] = True,
        only_convert: Optional[bool] = False,
        level: Optional[str] = "warning",
        depure: Optional[bool] = True,
    ) -> None:
        self.togadir: os.PathLike = togadir
        self.loss_status: Union[str, None] = by_loss_status
        self.orthology_class: Union[str, None] = by_orthology_class
        self.orthology_score: Union[float, None] = by_orthology_score
        self.to: str = to
        self.bed_type: str = target
        self.min_paralog_score: Union[float, None] = by_paralog_score
        self.with_isoforms: Union[os.PathLike, None] = with_isoforms
        self.extract: Union[bool, str] = extract
        self.outdir: Union[os.PathLike, None] = outdir
        self.only_table: bool = only_table
        self.only_convert: bool = only_convert
        self.log_level: str = level
        self.depure: bool = depure


## Executables
def timestamp() -> str:
    """Returns the current date- and timestamp"""
    return datetime.now().strftime("%H:%M_%d.%m.%y")


def dir_name_by_date(prefix: str) -> str:
    """Returns current date and time preceded by a given prefix"""
    return f"{prefix}_{timestamp()}"


def hex_code() -> str:
    """
    Generates a random five-digit decimal number
    and converts it into hexadecimal
    """
    return os.urandom(5).hex()


def hex_dir_name(prefix: str) -> str:
    """A combination of the two functions above"""
    # return f"{dir_name_by_date(prefix)}_{hex_code()}"
    return f"{prefix}_{hex_code()}"


def die(message: str) -> int:
    click.echo(message)
    sys.exit(1)


def get_upper_dir(file: str, lvl: int):
    """
    Get absolute path to a paternal directory of 'lvl' levels above the file.
    Setting lvl to zero will return the current directory for the file.
    """
    curr: str = os.path.abspath(file)
    for i in range(lvl):
        curr = os.path.dirname(curr)
    return curr


def aggregate(input_dir: str, output_file: str) -> None:
    """
    Aggregate the contents of all files in the specificed directory into one file
    """
    subprocess.run(f"cat {input_dir}/* > {output_file}", shell=True)


def parts(lst: Iterable[Any], n: int = 3):
    """Split an iterable into parts with size n."""
    return [lst[i : i + n] for i in iter(range(0, len(lst), n))]


def base_proj_name(projection: str) -> str:
    """Removes the metadata suffixes from the projection name"""
    return projection.split("$")[0].replace("#paralog", "").replace("#retro", "")


def segment_base(projection: str) -> str:
    """Returns the fragmented projection's name, ignoring the fragment number"""
    return projection.split("$")[0]


def get_proj2trans(projection: str) -> Tuple[str, str]:
    """Safely extract transcript name from the chain projection name"""
    data: List[str] = projection.split("$")[0].split("#")
    if data[-1] == "retro" or data[-1] == "paralog":
        data = data[:-1]
    return "#".join(data[:-1]), data[-1]


def safe_div(dividend: Union[float, int], divisor: Union[float, int]) -> float:
    """Divides dividend by divisor, return zero if divisor equals zero"""
    if not divisor:
        return 0.0
    return dividend / divisor


def nn(num: int) -> int:
    """Return the maximum between the provided number and zero"""
    return max(0, num)


def flatten(lst):
    """Flat list out of list of lists."""
    return [item for sublist in lst for item in sublist]


def safe_make_dir(dir: str) -> None:
    """
    Create directory if it does not exist,
    raise error if the name is already reserved for a regular file,
    do nothing otherwise
    """
    if os.path.isfile(dir) and not os.path.isdir(dir):
        die(
            f"{dir} was passed as a directory name "
            "but this name is already occupied by a regular file"
        )
    os.makedirs(dir) if not os.path.isdir(dir) else None


def intersection(start1: int, end1: int, start2: int, end2: int) -> int:
    """
    Returns intersection length between two segments with coordinates
    (start1, end1) and (start2, end2); values smaller than one indicate no intersection
    """
    return min(end1, end2) - max(start1, start2)


def chain_extract_id(index_file, chain_id, chain_file=None):
    """Extract chain text using index file."""
    # within TOGA should be fine:
    chain_file = chain_file if chain_file else index_file.replace(".bst", ".chain")
    if not os.path.isfile(chain_file):
        # need this check anyways
        sys.exit(f"chain_extract_id error: cannot find {chain_file} file")
    # connect shared library
    # .so must be there: in the modules/ dir
    script_location = os.path.dirname(__file__)
    slib_location = os.path.join(script_location, SLIB_NAME)
    sh_lib = ctypes.CDLL(slib_location)
    sh_lib.get_s_byte.argtypes = [
        ctypes.c_char_p,
        ctypes.c_uint64,
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
    ]
    sh_lib.get_s_byte.restype = ctypes.c_uint64

    # call library: find chain start byte and offset
    c_index_path = ctypes.c_char_p(index_file.encode())
    c_chain_id = ctypes.c_uint64(int(chain_id))
    c_sb = ctypes.c_uint64(0)  # write results in c_sb and c_of
    c_of = ctypes.c_uint64(0)  # provide them byref -> like pointers

    _ = sh_lib.get_s_byte(
        c_index_path, c_chain_id, ctypes.byref(c_sb), ctypes.byref(c_of)
    )

    if c_sb.value == c_of.value == 0:
        # if they are 0: nothing found then, raise Error
        sys.stderr.write(f"Error, chain {chain_id} ")
        sys.stderr.write("not found\n")
        sys.exit(1)

    # we got start byte and offset, extract chain from the file
    f = open(chain_file, "rb")
    f.seek(c_sb.value)  # jump to start_byte_position
    chain = f.read(c_of.value).decode("utf-8")  # read OFFSET bytes
    f.close()
    return chain


def make_cds_track(line):
    """Trim UTRs from a bed track."""
    line_data = line.rstrip().split("\t")
    if len(line_data) != 12:
        sys.exit(f"Error! Bed line:\n{line}\nis a not bed-12 formatted line!")
    # parse bed12 line according to the specification
    chrom = line_data[0]
    chrom_start = int(line_data[1])
    # chromEnd = int(line_data[2])
    name = line_data[3]  # gene_name usually
    name += "_CDS"  # mark that UTRs are trimmed
    bed_score = int(line_data[4])  # never used
    strand = line_data[5]
    thick_start = int(line_data[6])
    thick_end = int(line_data[7])
    item_rgb = line_data[8]  # never used
    block_count = int(line_data[9])
    # chrom start and end define the entire transcript location
    # this includes both UTRs and CDS
    # thick start and end limit the CDS only
    block_sizes = [int(x) for x in line_data[10].split(",") if x != ""]
    block_starts = [int(x) for x in line_data[11].split(",") if x != ""]
    block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
    # block starts are given in the relative coordinates -> need to convert them
    # into absolute coordinates using chrom start
    block_abs_starts = [block_starts[i] + chrom_start for i in range(block_count)]
    block_abs_ends = [block_ends[i] + chrom_start for i in range(block_count)]
    # arrays for blocks with trimmed UTRs
    block_new_starts, block_new_ends = [], []

    for block_num in range(block_count):
        # go block-by-block
        blockStart = block_abs_starts[block_num]
        blockEnd = block_abs_ends[block_num]

        # skip the block if it is entirely UTR
        if blockEnd <= thick_start:
            continue
        elif blockStart >= thick_end:
            continue

        # if we are here: this is not an entirely UTR exon
        # it might intersect the CDS border or to be in the CDS entirely
        # remove UTRs: block start must be >= CDS_start (thick_start)
        # block end must be <= CDS_end (thick_end)
        block_new_start = blockStart if blockStart >= thick_start else thick_start
        block_new_end = blockEnd if blockEnd <= thick_end else thick_end
        # save blocks with updated coordinates
        # also convert them back to relative coordinates with - thick_start
        # after the update thick_start/End are equal to chrom_start/End
        block_new_starts.append(block_new_start - thick_start)
        block_new_ends.append(block_new_end - thick_start)

    # block_count could change due to entirely UTR exons
    block_new_count = len(block_new_starts)
    # this is also a subject to change
    blockNewSizes = [
        block_new_ends[i] - block_new_starts[i] for i in range(block_new_count)
    ]

    # save the updated bed line with trimmed UTRs
    new_track = [
        chrom,
        thick_start,
        thick_end,
        name,
        bed_score,
        strand,
        thick_start,
        thick_end,
        item_rgb,
        block_new_count,
        ",".join([str(x) for x in blockNewSizes]) + ",",
        ",".join([str(x) for x in block_new_starts]) + ",",
    ]
    new_line = "\t".join([str(x) for x in new_track])
    return new_line


def get_connected_components(graph: nx.Graph) -> List[nx.Graph]:
    """ """
    ## check the NetworkX version
    nx_v: str = nx.__version__
    v_split: List[str] = [x for x in nx_v.split(".") if x.isnumeric()]
    if len(v_split) > 1:
        NX_VERSION: float = float(f"{v_split[0]}.{v_split[1]}")
    else:
        NX_VERSION: float = float(v_split[0])
    if NX_VERSION < 2.4:
        raw_components = list(nx.connected_component_subgraphs(graph))
    else:
        raw_components = [graph.subgraph(c) for c in nx.connected_components(graph)]
    return raw_components


def parse_single_column(file: TextIO) -> Set[str]:
    """
    Parse single-column file into a string list. Simple as.
    """
    output: Set[str] = set()
    if file is None:
        return output
    for line in file:
        line = line.rstrip()
        if not line:
            continue
        output.add(line)
    return output


def parse_fasta(file: str) -> Dict[str, str]:
    """
    A simple FASTA parser.
    Input is the file content.
    Input is a dictionary with sequence names as keys and aligned sequences as values.
    """
    output: Dict[str, str] = dict()
    key, seq, is_entry = "", "", False
    for line in file.split("\n"):
        if not line:
            continue
        if line[0] == ">":
            if is_entry:
                output[key] = seq
                seq = ""
            else:
                is_entry = True
            key = line.strip().split()[0][1:]
        else:
            seq += line.strip()
    if seq:
        output[key] = seq
    return output


def parse_score_file(filter_file: TextIO, score_column: int) -> Dict[str, Numeric]:
    out_dict: Dict[str, Numeric] = {}
    for line in filter_file.readlines():
        line: str = line.strip()
        line_data: List[str] = line.split("\t")
        if len(line_data) < 2:
            continue
        exon_name: str = line_data[0]
        sc: int = score_column - 1 if score_column <= len(line_data) else 1
        try:
            exon_filt_score: int = int(float(line_data[sc]))
            out_dict[exon_name] = exon_filt_score
        except ValueError:
            out_dict[exon_name] = 0
    return out_dict


def parse_one_column(file: Union[str, TextIO]) -> List[str]:
    """
    Parse single-column file into a string list. Simple as.
    """
    if isinstance(file, str):
        return list(map(lambda x: x.strip("\n\r\t"), file.readlines()))
    else:
        with open(file, "r") as h:
            return list(map(lambda x: x.strip("\n\r\t"), h.readlines()))


def read_tab(file: str) -> Iterable[str]:
    """Read a TSV file line by line, yield a generator of field-split lines"""
    if type(file) is str:
        gzipped: bool = file.endswith(".gz") or file.endswith("gzip")
        opener: Callable[str] = gzip.open if gzipped else open
        with opener(file, "r") as h:
            for line in h:
                if gzipped:
                    line = line.decode("utf8")
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                yield data
    else:
        for line in file:
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            yield data


def reverse_complement(seq: str) -> str:
    """Returns a reverse complement of a standard alphabet nucleotide sequence"""
    out_seq: str = "".join([COMPLEMENT[x] for x in seq[::-1]])
    return out_seq


def is_locked(file: str) -> bool:
    if not os.path.isfile(file):
        return False
    try:
        with open(file, "a"):
            pass
    except IOError:
        return True

    dummy: str = file + ".dummy"
    try:
        copy2(file, dummy)
        os.remove(dummy)
        return False
    except WindowsError:
        return True


class Lock:
    """
    A homebrew lock manager to coordinate multiple scripts sharing resources
    """

    def __init__(
        self, lockfile: str, retry: float = 0.25, timeout: float = 20.0
    ) -> None:
        self.lockfile: str = lockfile
        self.retry: float = retry
        self.timeout: float = timeout
        self.filehandle: Union[TextIO, None] = None

    def __enter__(self) -> TextIO:
        print(self.lockfile, self.retry, self.timeout)
        # is_available: bool = not is_locked(self.lockfile) ## check if lockfile exists
        # print(f'{is_available = }')
        fh: TextIO = os.open(self.lockfile, os.O_RDWR | os.O_CREAT)
        is_available: bool = False
        curr_time: float = time.time()
        start_time: float = curr_time
        while True:
            if curr_time > start_time + self.timeout:
                break
            try:
                # is_available = not is_locked(self.lockfile)
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
                is_available = True
                break
            except (IOError, OSError):
                print("Lock not acquired; waiting")
                time.sleep(self.retry)  ## sleep the expected amount of time
                curr_time = time.time()

        if is_available:  ## create the lockfile
            print("File lock successfully acquired")
            self.filehandle = fh
        else:
            print("Failed to obtain the block")
            os.close(fh)

        return self.filehandle

    def __exit__(self, exc_type: Any, exc_value: Any, exc_tb: Any) -> None:
        if self.filehandle:
            fcntl.flock(self.filehandle, fcntl.LOCK_UN)
            os.close(self.filehandle)


class CommandLineManager:
    __slots__ = ("v", "debug", "logger")
    """
    A minimal functionality class to be decorated with Click functionality and
    further extended to suit particular scripts' needs
    """

    def set_logging(self, name: str = __name__, toga_module: str = None) -> None:
        """
        Sets up logging system for a TogaMain instance
        """
        logging.root.handlers = []
        if name is None:
            name = __name__
        if toga_module is None:
            toga_module = "unknown_module"
        self.logger: logging.Logger = logging.getLogger(name)
        if self.logger.handlers:
            self.logger = logging.LoggerAdapter(
                self.logger, {"toga_module": toga_module}
            )
            return
        if hasattr(self, "debug") and self.debug:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)
        if hasattr(self, "log_file") and self.log_file:
            file_handler: logging.FileHandler = logging.FileHandler(
                self.log_file, mode="a", encoding=UTF8
            )
            file_handler.setFormatter(FORMATTER)
            self.logger.addHandler(file_handler)
        if hasattr(self, "v") and self.v:
            console_handler: logging.StreamHandler = logging.StreamHandler()
            console_handler.setFormatter(FORMATTER)
            self.logger.addHandler(console_handler)
        self.logger = logging.LoggerAdapter(self.logger, {"toga_module": toga_module})

    def _to_log(self, msg: str, level: Optional[str] = "info") -> None:
        """Logs a message at a given level

        Args:
            msg: A message to report in the log channel/file
            level: Logging level to report the message at

        Raises:
            Returns without logging anything if the current class has no 'logger' attribute
        """
        if not hasattr(self, "logger"):
            return
        getattr(self.logger, level)(msg)

    def _debug(self, msg: str) -> None:
        """Logs message at DEBUG level

        Args:
            msg: A message to report in the log channel/file

        Raises:
            Returns without logging anything if `debug` attribute does not exist
        or is set to False or if `logger` attribute does not exist
        """
        if not hasattr(self, "debug") or not self.debug:
            return
        self._to_log(msg=msg, level="debug")

    def _echo(self, msg: str) -> None:
        """Report a line to standard output if verbosity is enabled"""
        click.echo(msg) if self.v else None

    def _stderr(self, msg: str) -> None:
        """Report a line to standard error stream regardless of verbosity settings"""
        sys.stderr.write(msg + "\n")

    def _exit(self, msg: str = None) -> None:
        """Safe exit witn zero return code and a given message"""
        self._to_log(msg) if msg is not None else None
        sys.exit(0)

    def _die(self, msg: str = None) -> None:
        """Error-exit with a given message"""
        if msg is not None:
            self._to_log(msg, "critical")
        sys.exit(1)

    def _cp(self, from_: os.PathLike, to_: os.PathLike) -> None:
        """Copies file from from_ to to_"""
        try:
            copy2(from_, to_)
        except Exception as e:
            self._die(
                ("Unexpected behavior when trying to copy from %s to %s: \n")
                % (from_, to_)
                + e
            )

    def _mkdir(self, d: os.PathLike) -> None:
        """Safe directory creation method"""
        try:
            os.makedirs(d)
        except FileExistsError:
            pass

    def _rm(self, f: os.PathLike) -> None:
        """File deletion method"""
        try:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                self._rmdir(f)
        except FileNotFoundError:
            pass
        except Exception:
            self._die("Unexpected behavior observed while trying to remove %s" % f)

    def _rmdir(self, d: os.PathLike) -> None:
        """Recursive directory deletion method"""
        try:
            if os.path.islink(d):
                path: os.PathLike = os.readlink(d)
                rmtree(path, ignore_errors=True)
                os.remove(d)
            else:
                rmtree(d, ignore_errors=True)
        except FileNotFoundError:
            pass
        except Exception:
            self._die("Unexpected behavior observed while trying to remove %s" % d)

    def _mv(self, file: os.PathLike, dest: os.PathLike) -> None:
        """Moves a file or directory to a new location"""
        move(file, dest)

    def _abspath(self, path: os.PathLike) -> str:
        """Checks whether a path is absolute, prepends root prefix if not"""
        if os.path.isabs(path):
            return path
        return os.path.abspath(path)

    def _exec(
        self,
        cmd: Union[str, List[str]],
        err_msg: str,
        input_: bytes = None,
        die: bool = True,
        shun_verbosity: bool = True,
        gather_stdout: bool = True,
    ) -> str:
        """Runs subprocesses, handles the exceptions"""
        ## TODO: Consider slowly moving to shell=False
        # if isinstance(cmd, str):
        #     cmd: List[str] = [x for x in cmd.split(' ') if x]
        pr: subprocess.Pipe = subprocess.Popen(
            cmd,
            shell=True,
            executable="bash",
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE if gather_stdout else None,
            stderr=subprocess.PIPE,
        )
        if self.v and not shun_verbosity and not gather_stdout:
            for line in pr.stdout:
                if not line:
                    continue
                self._echo(line)  ##what was this supposed to mean?
        stdout, stderr = pr.communicate(input=input_)
        rc: int = pr.returncode
        if rc != 0:
            err: str = stderr.decode("utf8")
            if not die:
                return err
            msg: str = f"{err_msg}:\n{err}"
            self._die(msg)
        if gather_stdout:
            return stdout.decode("utf8")

    def _symlink(self, from_: os.PathLike, to_: os.PathLike) -> None:
        cmd: str = f"ln -s {from_} {to_}"
        self._exec(cmd, "Symlink creation failed:")
