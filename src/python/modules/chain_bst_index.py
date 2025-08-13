#!/usr/bin/env python3
"""Save chain index file.

Chain index is just a BST saved into a separate file.
Using this index file for any chain ID we can extract:
1) Start byte of this chain in the chain file
2) Length of this chain in the file.
And then simply extract it.
"""
import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from modules.common import to_log
from modules.shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Optional, Union
# from version import __version__

import click
import ctypes

__author__ = "Bogdan M. Kirilenko"
__credits__ = ('Yury V. Malovichko')

SLIB_NAME = "chain_bst_lib.so"


def chain_bst_index(chain_file, index_file, txt_index=None):
    """Create index file for chain."""
    # assume that shared lib is in the same dir
    script_location = os.path.dirname(__file__)
    slib_location = os.path.join(script_location, SLIB_NAME)
    if not os.path.isfile(slib_location):
        sys.exit("chain_bst_lib.so lib not found, call ./configure to compile")
    # connect shared lib
    sh_lib = ctypes.CDLL(slib_location)
    sh_lib.make_index.argtypes = [
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.c_uint64,
        ctypes.c_char_p,
    ]
    sh_lib.make_index.restype = ctypes.c_int

    # read chain file, get necessary data: start bytes and offsets
    chain_ids = [0]
    start_bytes = [0]
    offsets = [0]

    byte_num = 0
    offset = 0

    f = open(chain_file, "rb")
    for line in f:
        if not line.startswith(b"chain"):
            # just count these bytes
            byte_num += len(line)
            offset += len(line)
            continue
        # if we're here -> this is a chain header
        offsets.append(offset)
        header = line.decode("utf-8").rstrip()
        chain_id = header.split()[-1]
        chain_ids.append(int(chain_id))
        start_bytes.append(int(byte_num))
        # byte num is absolute and offset (chain size) is relative
        byte_num += len(line)  # si just continue incrementing byte num
        offset = len(line)  # and reset (to line length) offset
    f.close()

    offsets.append(offset)
    arr_size = len(chain_ids)
    del offsets[0]

    if arr_size == 0:
        to_log(f"chain_bst_index: ERROR: No chains found. Abort")
        sys.exit(1)
    to_log(f"chain_bst_index: indexing {arr_size} chains")

    if txt_index is not None:
        # save text (non-binary) dict for chain ids and positions in the file
        # in some cases this is more efficient that extracting data from BST
        # via shared library
        f = open(txt_index, "w")
        for x in zip(chain_ids, start_bytes, offsets):
            f.write(f"{x[0]}\t{x[1]}\t{x[2]}\n")
        f.close()

    # call shared lib
    c_chain_ids = (ctypes.c_uint64 * (arr_size + 1))()
    c_chain_ids[:-1] = chain_ids
    c_s_bytes = (ctypes.c_uint64 * (arr_size + 1))()
    c_s_bytes[:-1] = start_bytes
    c_offsets = (ctypes.c_uint64 * (arr_size + 1))()
    c_offsets[:-1] = offsets

    c_arr_size = ctypes.c_uint64(arr_size)
    c_table_path = ctypes.c_char_p(str(index_file).encode())
    _ = sh_lib.make_index(c_chain_ids, c_s_bytes, c_offsets, c_arr_size, c_table_path)
    to_log(f"chain_bst_index: Saved chain {chain_file} index to {index_file}")

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'input',
    type=click.Path(exists=True),
    metavar='CHAIN_FILE'
)
@click.argument(
    'bst_output',
    type=click.Path(exists=False),
    metavar='BST_INDEX_OUTPUT'
)
@click.option(
    '--text_output',
    '-t',
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help='If provided, doubles the index file in plain text mode to the provided path'
)
@click.option(
    '--log_name',
    '-ln',
    type=str,
    metavar='STR',
    default=None,
    show_default=True,
    help='Logger name to use; relevant only upon main class import'
)
@click.option(
    '--verbose',
    '-v',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)

class ChainIndexer(CommandLineManager):
    __slots__ = ('input', 'bst_output', 'text_output', 'log_name', 'v')

    def __init__(
        self,
        input: click.Path,
        bst_output: click.Path,
        text_output: Optional[click.Path],
        log_name: Optional[int],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)

        chain_bst_index(input, bst_output, txt_index=text_output)

if __name__ == "__main__":
    ChainIndexer()
