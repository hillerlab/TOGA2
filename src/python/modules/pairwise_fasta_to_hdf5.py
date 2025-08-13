#!/usr/bin/env python3

"""
Converts pairwise FASTA file into an HDF5 storage 
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from collections import defaultdict
from modules.shared import CommandLineManager, CONTEXT_SETTINGS
from numpy import bytes_
from typing import Dict, List, Optional, TextIO

import click
import h5py

HEADER_START: str = '>'
REFERENCE: str = 'REFERENCE'
REF_SOURCE: str = '_ref'
QUERY_SOURCE: str = '_query'
REF_EXON: str = 'reference_exon'

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'input',
    type=click.File('r', lazy=True),
    metavar='INPUT_FASTA'
)
@click.argument(
    'output',
    type=click.Path(exists=False),
    metavar='OUTPUT_HDF5'
)
@click.option(
    '--exon_fasta',
    '-e',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, expects TOGA2 exon alignment FASTA file as input. Exon number will be added '
        'to HDF5 sequences identifiers'
    )
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
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)

class FastaToHdf5Converter(CommandLineManager):
    """
    Converts pairwise FASTA into HDF5 storage. Sequences are stored under the key of 
    {projection_name}_{source}, where {source} is "_ref" for reference and 
    "_query" for query sequences, respectively.\n
    This script is used by TOGA2 to convert pairwise protein file into HDF5 used for gene tree step 
    input construction; it is not intended to be used outside of the TOGA2 pipeline.\n
    Arguments are:\n
    * INPUT is input pairwise FASTA file;\n
    * OUTPUT is a path to the output HDF5 storage
    """
    __slots__ = ()

    def __init__(
        self,
        input: click.File,
        output: click.Path,
        exon_fasta: Optional[bool],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)

        if exon_fasta:
            self.write_exons_for_sleasy(input, output)
        else:
            self.write_proteins_for_phylo(input, output)

    def write_proteins_for_phylo(self, input: TextIO, output: str) -> None:
        """
        Writes pairwise protein Fasta records to HDF5, one dataset for each entry
        """
        with h5py.File(output, 'w') as f:
            # header_encountered: bool = False
            header: str = ''
            hdf5_id: str = ''
            seq: str = ''
            for line in input:
                line  = line.rstrip()
                if not line:
                    continue
                if line[0] == HEADER_START:
                    if header:
                        # print(hdf5_id)
                        try:
                            f.create_dataset(
                                hdf5_id, 
                                data=bytes_(f'{header}\n{seq}')
                                # dtype=string_,
                                # shape=(len(header) + len(seq) + 1,),
                                # chunks=True
                            )
                            header = ''
                            hdf5_id = ''
                            seq = ''
                        except ValueError as e:
                            print(e)
                            self._die(
                                'Sequence %s occurs twice in the input FASTA file' % hdf5_id
                            )
                    header = line
                    header_split: List[str] = header.split(' | ')
                    proj: str = header_split[0].lstrip(HEADER_START)
                    source: str = REF_SOURCE if header_split[2] == REFERENCE else QUERY_SOURCE
                    hdf5_id = f'{proj}{source}'
                else:
                    seq += line
            if header:
                f.create_dataset(
                    hdf5_id, 
                    data=bytes_(f'{header}\n{seq}')
                    # dtype=string_,
                    # shape=(len(header) + len(seq) + 1,),#len(header) + len(seq) + 1,
                    # chunks=True
                )

    def write_exons_for_sleasy(self, input: TextIO, output: str) -> None:
        """
        Writes query exons to the HDF5 file. Datasets are organized projection-wise, 
        sequences are stored as ordered lists of variable-length strings
        """
        exon_seq_dict: Dict[str, List[str]] = defaultdict(list)
        # exon_seqs: List[str] = []
        # exon_names: List[str] = []
        proj: str = ''
        exon_num: int = 0
        seq: str = ''
        for line in input:
                line  = line.rstrip()
                if not line:
                    continue
                if line[0] == HEADER_START:
                    if proj:
                        seq = seq.replace('-','')
                        if not seq:
                            seq = 'DELETED'
                        exon_seq_dict[proj].append((exon_num, seq))
                        # exon_seqs.append(seq)
                        # exon_names.append(f'{proj}_{exon_num}')
                    seq = ''
                    header_split: List[str] = line.split(' | ')
                    if header_split[3] == REF_EXON:
                        proj = ''
                        exon_num = 0
                        continue
                    proj: str = header_split[0].lstrip(HEADER_START)
                    exon_num: int = int(header_split[1])
                    continue
                if not proj:
                    continue
                seq += line.rstrip().replace('-','')
        if proj:
            exon_seq_dict[proj].append((exon_num, seq))
            # exon_seqs.append(seq)
            # exon_names.append(f'{proj}_{exon_num}')

        with h5py.File(output, 'w') as f:
            # f.create_dataset(
            #     'sequence_names', 
            #     data=exon_names, 
            #     dtype=h5py.string_dtype(encoding='utf-8')
            # )
            # f.create_dataset(
            #     'sequences', 
            #     data=exon_seqs, 
            #     dtype=h5py.string_dtype(encoding='utf-8'), 
            #     compression='gzip', 
            #     compression_opts=9
            # )

        
            for projection, exons in exon_seq_dict.items():
                sorted_exons: List[str] = [
                    bytes_(x[1]) for x in sorted(exons, key=lambda y: y[0])
                ]
                shape: int = len(sorted_exons)
                f.create_dataset(
                    projection, 
                    data=sorted_exons,
                    shape=shape,
                    compression='gzip',
                    compression_opts=9
                    # dtype=bytes_,
                    # shape=(len(header) + len(seq) + 1,),#len(header) + len(seq) + 1,
                    # chunks=True
                )


if __name__ == '__main__':
    FastaToHdf5Converter()