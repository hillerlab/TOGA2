#!/usr/bin/env python3
"""Convert text bed file to hdf5.

This allows TOGA to extract bed track
for a particular transcript ID immediately.
"""
import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from modules.common import to_log
from modules.shared import CommandLineManager, CONTEXT_SETTINGS
from numpy import bytes_
from typing import List, Optional

import click
import h5py
# from version import __version__

__author__ = 'Yury V. Malovichko'
__credits__ = ('Bogdan M. Kirilenko')


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'in_bed',
    type=click.File('r', lazy=True),
    metavar='BED'
)
@click.argument(
    'out_db',
    type=click.Path(exists=False),
    metavar='HDF5'
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

class BedHdf5Indexer(CommandLineManager):
    """
    Converts BED file into a string HDF5 collection with sequence IDs as keys and full BED records as values.
    Arguments are:\n
    * BED is an input BED file; by default TOGA expects an uncompressed file with at least twelve columns;\n
    * HDF5 is a path to output HDF5 file
    """
    __slots__ = ('log_file')

    def __init__(
        self, 
        in_bed: click.File, 
        out_db: click.Path,
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)

        i: int = 0
        with h5py.File(out_db, 'w') as f:
            for i, line in enumerate(in_bed, start=1):
                # line = line.rstrip()
                data: List[str] = line.split('\t')
                if not data or not data[0]:
                    continue
                if len(data) < 12:
                    self._to_log(
                        'Input BED file contains less than 12 fields at line %i' % i,
                        'warning'
                    )
                name: str = data[3]
                f.create_dataset(name, data=bytes_(line))
        if not i:
            self._die('ERROR - Empty BED file provided for HDF5 indexing')


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could fail
    BedHdf5Indexer()
