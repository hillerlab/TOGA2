#!/usr/bin/env python3

"""
Turns the IntronIC results into an HDF5 storage containing data
on reference genome's non-canonical introns
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from collections import defaultdict
from numpy import array, str_
from shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Dict, Iterable, List, Optional, Tuple, TypeVar

import h5py
import click

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__include__ = (None)

ACCEPTOR: str = 'A'
DONOR: str = 'D'

CoordKey: TypeVar = TypeVar('CoordKey', bound='Tuple[str, int, int, str]')
Record: TypeVar = TypeVar('Record', bound='List[str]')

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'bed12',
    type=click.Path(exists=True),
    metavar='BED_OR_HDF5'
)
@click.argument(
    'intronic',
    type=click.Path(exists=True),
    metavar='INTRON_IC_OUT'
)
@click.argument(
    'output',
    type=click.Path(exists=False),
    metavar='OUTPUT'
)
@click.option(
    '--hdf5_input',
    '-hdf5',
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, BED_OR_HDF5 is treated as a transcript:BED_line HDF5 storage'
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

class IntronIcConverter(CommandLineManager):
    __slots__ = (
        'bed12', 'intronic', 'output',
        'hdf5_input', 'recorded_sites', 'sites_to_write'
    )

    def __init__(
        self,
        bed12: click.Path,
        intronic: click.Path,
        output: click.Path,
        hdf5_input: Optional[bool],
        log_name: Optional[int],
        verbose: Optional[bool]
    ) -> None:
        self.v  = verbose
        self.set_logging(log_name)

        self.bed12: click.Path = bed12
        self.intronic: click.Path = intronic
        self.output: click.Path = output
        self.hdf5_input: bool = hdf5_input

        self.recorded_sites: Dict[CoordKey, Record] = {}
        self.sites_to_write: Dict[str, List[Record]] = defaultdict(list)

        self.run()

    def run(self) -> None:
        """
        """
        self.parse_intronic_output()
        if self.hdf5_input:
            with h5py.File(self.bed12, 'r') as f:
                for tr in f:
                    data: Iterable[str] = f[tr][()].decode('utf-8')
                    self.process_transcript_instance(data)
        else:
            with open(self.bed12, 'r') as h:
                for line in h:
                    self.process_transcript_instance(line)
        self.write_output()

    def parse_intronic_output(self) -> None:
        """
        """
        with open(self.intronic, 'r') as h:
            for line in h:
                data: List[str] = line.strip().split('\t')
                if not line:
                    continue
                if len(data) < 6:
                    self._die(
                        'ERROR: IntronIC output contains less than six columns. '
                        'Make sure that a proper IntronIC pipeline output in the '
                        'BED6 format was provided'
                    )
                chrom: str = data[0]
                start: int = int(data[1])
                end: int = int(data[2])
                strand: str = data[5]
                key: CoordKey = (chrom, start, end, strand)
                name: str = data[3]
                _, intron_class, dinuc = name.split('_') ## TODO: Must always correspond to the IntronIC pipeline's notation
                donor_dinuc, acc_dinuc = dinuc.split('-')
                self.recorded_sites[key] = (intron_class, donor_dinuc, acc_dinuc)

    def process_transcript_instance(self, line: str) -> None:
        """
        """
        data: List[str] = line.strip().split('\t')
        if not data:
            return
        if len(data) != 12:
            self._die(
                'ERROR: Reference annotation passed is improperly formatted. '
                'Make sure that you provided a valid BED12 file or a HDF5 '
                'storage containing BED12 data'
            )
        chrom: str = data[0]
        name: str = data[3]
        strand: str = data[5]
        is_pos: bool = strand == '+'
        cds_start: int = int(data[6])
        cds_end: int = int(data[7])
        exon_num: int = int(data[9])
        exon_sizes: int = list(map(int, data[10].split(',')[:-1]))
        exon_starts: int = list(map(int, data[11].split(',')[:-1]))
        for i in range(exon_num - 1):
            ex_num: int = i + 1 if is_pos else exon_num - i
            next_ex_num: int = ex_num + (1 if is_pos else -1)
            exon_start: int = cds_start + exon_starts[i]
            exon_end: int = exon_start + exon_sizes[i]
            next_exon_start: int = cds_start + exon_starts[i+1]
            intron_key: CoordKey = (chrom, exon_end, next_exon_start, strand)
            if intron_key in self.recorded_sites:
                intron: str = ex_num if is_pos else next_ex_num
                intron_data: Tuple[str] = self.recorded_sites[intron_key]
                intron_value: Record = (str(intron), *intron_data)
                self.sites_to_write[name].append(intron_value)

    def write_output(self) -> None:
        """
        """
        with h5py.File(self.output, 'w') as f:
            for name, record in self.sites_to_write.items():
                dataset: array = array(record).astype(str_)
                try:
                    ds = f.create_dataset(
                        name,
                        shape=dataset.shape,
                        dtype=h5py.string_dtype(encoding='utf-8'),
                        chunks=True
                    )
                    ds[:] = dataset
                except ValueError:
                    del f[name]
                    ds = f.create_dataset(
                        name,
                        shape=(1, 4),
                        dtype=h5py.string_dtype(encoding='utf-8'),
                        chunks=True
                    )
                    ds[:] = array(list(record[1:])).astype(str_)


if __name__ == '__main__':
    IntronIcConverter()
