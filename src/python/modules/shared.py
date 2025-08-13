#!/usr/bin/env python3

"""
Shared functionality across the scripts
"""

from click_option_group import OptionGroup
from .constants import Constants
from datetime import datetime
from shutil import copy2, rmtree
from typing import (
    Any, Dict, Iterable, List, Optional, TextIO, Set, Tuple, Union
)

import click
import fcntl
import logging
import networkx as nx
import os
import subprocess
import sys
import time

## Constants
CONTEXT_SETTINGS: Dict[str, Any] = {
    'help_option_names': [
        '-h', '-help', '--help'
    ],
    'ignore_unknown_options': True,
    'allow_extra_args': True,
    'max_content_width': 150
}
## borrowed from splitFile
SPLIT_JOB_HEADER: Tuple[str] = (
    "#!/bin/bash",
    "set -eu",
    "set -o pipefail"
)
TMP_DIR_SUBS: Dict[str, List[str]] =  {
    'joblists': [],
    'job_partitions': [],
    'coords': [],
    'lastz': [
        'search_space_fasta',
        'exon_fasta',
        'lastz',
        'missing_exon_fasta',
        'missing_lastz',
        'bed'
    ],
    'hmm': [],
    'cesar': ['input', 'output']
}

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
    "-": "-"
}

## Types
Numeric = Union[float, int]

class PrettyGroup(OptionGroup):
    def get_help_record(self, ctx: click.Context) -> Optional[Tuple[str, str]]:
        init_help: Union[Tuple[str, str], None] = super().get_help_record(ctx)
        if init_help is None:
            return None
        return '\n' + init_help[0], init_help[1]


## Executables
def dir_name_by_date(prefix: str) -> str:
    """Returns current date and time preceded by a given prefix"""
    return f"{prefix}_{datetime.now().strftime('%H:%M_%d.%m.%y')}"


def hex_code() -> str:
    """
    Generates a random five-digit decimal number 
    and converts it into hexadecimal
    """
    return os.urandom(5).hex()


def hex_dir_name(prefix: str) -> str:
    """A combination of the two functions above"""
    return f'{dir_name_by_date(prefix)}_{hex_code()}'


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
    subprocess.run(f'cat {input_dir}/* > {output_file}', shell=True)


def parts(lst: Iterable[Any], n: int =3):
    """Split an iterable into parts with size n."""
    return [lst[i : i + n] for i in iter(range(0, len(lst), n))]


def get_proj2trans(projection: str) -> Tuple[str, str]:
    """Safely extract transcript name from the chain projection name"""
    data: List[str] = projection.split('#')
    return '#'.join(data[:-1]), data[-1]


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
            f'{dir} was passed as a directory name '
            'but this name is already occupied by a regular file'
        )
    os.makedirs(dir) if not os.path.isdir(dir) else None


def intersection(
    start1: int,
    end1: int,
    start2: int,
    end2: int
) -> int:
    """
    Returns intersection length between two segments with coordinates 
    (start1, end1) and (start2, end2); values smaller than one indicate no intersection
    """
    return (min(end1, end2) - max(start1, start2))


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
    """
    """
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
        raw_components = [
            graph.subgraph(c) for c in nx.connected_components(graph)
        ]
    return raw_components


def parse_single_column(file: TextIO) -> Set[str]:
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
    key, seq, is_entry = '', '', False
    for line in file.split('\n'):
        if not line:
            continue
        if line[0] == '>':
            if is_entry:
                output[key] = seq
                seq = ''
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
        line_data: List[str] = line.split('\t')
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
        return list(map(lambda x: x.strip('\n\r\t'), file.readlines()))
    else:
        with open(file, 'r') as h:
            return list(map(lambda x: x.strip('\n\r\t'), h.readlines()))


def reverse_complement(seq: str) -> str:
    """Returns a reverse complement of a standard alphabet nucleotide sequence"""
    out_seq: str = ''.join([COMPLEMENT[x] for x in seq[::-1]])
    return out_seq


def is_locked(file: str) -> bool:
    if not os.path.isfile(file):
        return False
    try:
        with open(file, 'a'):
            pass
    except IOError:
        return True

    dummy: str = file + '.dummy'
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
                print('Lock not acquired; waiting')
                time.sleep(self.retry) ## sleep the expected amount of time
                curr_time = time.time()

        if is_available: ## create the lockfile
            print('File lock successfully acquired')
            self.filehandle = fh
        else:
            print('Failed to obtain the block')
            os.close(fh)

        return self.filehandle


    def __exit__(self, exc_type: Any, exc_value: Any, exc_tb: Any) -> None:
        if self.filehandle:
            fcntl.flock(self.filehandle, fcntl.LOCK_UN)
            os.close(self.filehandle)


class CommandLineManager:
    __slots__ = ('v', 'logger')
    """
    A minimal functionality class to be decorated with Click functionality and
    further extended to suit particular scripts' needs
    """
    def set_logging(self, name: str = __name__) -> None:
        """
        Sets up logging system for a TogaMain instance
        """
        if name is None:
            name = __name__
        self.logger: logging.Logger = logging.getLogger(name)
        if self.logger.handlers:
            return
        self.logger.setLevel(logging.DEBUG)
        if hasattr(self, 'log_file') and self.log_file:
            file_handler: logging.FileHandler = logging.FileHandler(
                self.log_file, mode='a', encoding=Constants.UTF8
            )
            file_handler.setFormatter(Constants.FORMATTER)
            self.logger.addHandler(file_handler)
        if hasattr(self, 'v') and self.v:
            console_handler: logging.StreamHandler = logging.StreamHandler()
            console_handler.setFormatter(Constants.FORMATTER)
            self.logger.addHandler(console_handler)

    def _to_log(self, msg: str, level: str = 'info'):
        """
        Adds the message to the log
        """
        if not hasattr(self, 'logger'):
            return
        getattr(self.logger, level)(msg)

    def _echo(self, msg: str) -> None:
        """Report a line to standard output if verbosity is enabled"""
        click.echo(msg) if self.v else None

    def _stderr(self, msg: str) -> None:
        """Report a line to standard error stream regardless of verbosity settings"""
        sys.stderr.write(msg + '\n')

    def _exit(self, msg: str = None) -> None:
        """Safe exit witn zero return code and a given message"""
        self._to_log(msg) if msg is not None else None
        sys.exit(0)

    def _die(self, msg: str = None) -> None:
        """Error-exit with a given message"""
        if msg is not None:
            self._to_log(msg, 'critical')
        sys.exit(1)

    def _cp(self, from_: str, to_: str) -> None:
        """Copies file from from_ to to_"""
        try:
            copy2(from_, to_)
        except Exception as e:
            self._die(
                (
                    'Unexpected behavior when trying to copy from %s to %s: \n'
                ) % (from_, to_) + e
            )

    def _mkdir(self, d: str) -> None:
        """Safe directory creation method"""
        try:
            os.makedirs(d)
        except FileExistsError:
            pass

    def _rm(self, f: str) -> None:
        """File deletion method"""
        try:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                rmtree(f)
        except FileNotFoundError:
            pass
        except Exception:
            self._die('Unexpected behaviour observed while trying to remove %s' % f)

    def _rmdir(self, d: str) -> None:
        """Recursive directory deletion method"""
        try:
            rmtree(d)
        except FileNotFoundError:
            pass

    def _exec(
        self,
        cmd: Union[str, List[str]],
        err_msg: str,
        input_: bytes = None,
        die: bool = True,
        shun_verbosity: bool = True,
        gather_stdout: bool = True
    ) -> str:
        """Runs subprocesses, handles the exceptions"""
        ## TODO: Consider slowly moving to shell=False
        # if isinstance(cmd, str):
        #     cmd: List[str] = [x for x in cmd.split(' ') if x]
        pr: subprocess.Pipe = subprocess.Popen(
            cmd, 
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE if gather_stdout else None,
            stderr=subprocess.PIPE
        )
        if self.v and not shun_verbosity:
            for line in pr.stdout:
                if not line:
                    continue
                # self._echo(line)
        stdout, stderr = pr.communicate(input=input_)
        rc: int = pr.returncode
        if rc != 0:
            err: str = stderr.decode('utf8')
            if not die:
                return err
            msg: str = f'{err_msg}:\n{err}'
            self._die(msg)
        if gather_stdout:
            return stdout.decode('utf8')
