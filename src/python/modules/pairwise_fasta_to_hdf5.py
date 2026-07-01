#!/usr/bin/env python3

"""
Converts pairwise FASTA file into an HDF5 storage
"""

from typing import Dict, List, Tuple, Optional, TextIO

import click
import h5py
from numpy import array, int32
from .shared import CONTEXT_SETTINGS, CommandLineManager, get_proj2trans

HEADER_START: str = ">"
REFERENCE: str = "REFERENCE"
REF_SOURCE: str = "_ref"
QUERY_SOURCE: str = "_query"
REF_EXON: str = "reference_exon"


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument("input", type=click.File("r", lazy=True), metavar="INPUT_FASTA")
@click.argument("output", type=click.Path(exists=False), metavar="OUTPUT_HDF5")
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
class FastaToHdf5Converter(CommandLineManager):
    """
    Converts pairwise FASTA into HDF5 storage. Sequences are stored under the key of
    {projection_name}_{source}, where {source} is "_ref" for reference and
    "_query" for query sequences, respectively.\n
    \n
    The datasets in the output HDF5 file are:\n
    \t* "name": Sequence name (projection ID + "_ref"/"_query" suffix); guaranteed to be sorted;\n
    \t* "seq_num": Number of the the sequence in the third dataset ("seq") 
    corresponding to the name in the first dataset ("name")\n
    \t* "seq": non-redundant protein sequences\n
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
        log_name: Optional[str],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.set_logging(name=log_name, toga_module="fasta2hdf5")

        self.write_proteins_for_phylo(input, output)

    def write_proteins_for_phylo(self, input: TextIO, output: str) -> None:
        """
        Writes pairwise protein Fasta records to HDF5 as three parallel variable-length
        string datasets:
            * "name": Sequence name (projection ID + "_ref"/"_query" suffix)
            * "seq_num": Number of the the sequence in the third dataset ("seq") 
            corresponding to the name in the first dataset ("name")
            * "seq": non-redundant protein sequences

        Args:
            input: file handle for the input FASTA file
            output: path to the output HDF5 storage

        Returns:
            None
        """
        # records: Dict[str, str] = {}
        name2id: Dict[str, int] = {}
        seq2id: Dict[str, int] = {}
        header: str = ""
        hdf5_id: str = ""
        curr_id: int = 0
        seq_parts: List[str] = []

        for line in input:
            line = line.rstrip()
            if not line:
                continue
            if line[0] == HEADER_START:
                if header:
                    if hdf5_id in name2id:#records:
                        self._die(
                            "Sequence %s occurs twice in the input FASTA file" % hdf5_id
                        )
                    seq: str = "".join(seq_parts).replace("-", "")
                    if seq in seq2id:
                        seq_num: int = seq2id[seq]
                        name2id[hdf5_id] = seq_num
                    else:
                        seq2id[seq] = curr_id
                        name2id[hdf5_id] = curr_id
                        curr_id += 1
                    # records[hdf5_id] = f"{header}\n{seq}"
                    seq_parts = []
                header = line
                header_split: List[str] = header.split(" | ")
                proj: str = header_split[0].lstrip(HEADER_START)
                source: str = (
                    REF_SOURCE if header_split[2] == REFERENCE else QUERY_SOURCE
                )
                if source == REF_SOURCE:
                    proj = get_proj2trans(proj)[0]
                hdf5_id = f"{proj}{source}"
            else:
                seq_parts.append(line)

        if header:
            seq: str = "".join(seq_parts).replace("-", "")
            if seq in seq2id:
                seq_num: int = seq2id[seq]
                name2id[hdf5_id] = seq_num
            else:
                seq2id[seq] = curr_id
                name2id[hdf5_id] = curr_id
            # records[hdf5_id] = f"{header}\n{''.join(seq_parts)}"

        dt = h5py.string_dtype()
        # sorted_items: List[Tuple[str, str]] = sorted(records.items())
        sorted_items: List[Tuple[str, int]] = sorted(name2id.items())
        sorted_seqs: List[Tuple[str, int]] = sorted(seq2id.items(), key=lambda x: x[1])
        with h5py.File(output, "w") as f:
            # f.create_dataset(
            #     "keys",
            #     data=array([k for k, _ in sorted_items], dtype=object),
            #     dtype=dt,
            # )
            # f.create_dataset(
            #     "values",
            #     data=array([v for _, v in sorted_items], dtype=object),
            #     dtype=dt,
            # )
            f.create_dataset(
                "name",
                data=array([k for k, _ in sorted_items], dtype=dt),
                dtype=dt,
            )
            f.create_dataset(
                "seq_num",
                data=array([v for _, v in sorted_items], dtype=int32),
                dtype=int32,
            )
            f.create_dataset(
                "seq",
                data=array([k for k, _ in sorted_seqs], dtype=dt),
                dtype=dt,
            )
            f.attrs["sorted"] = True


if __name__ == "__main__":
    FastaToHdf5Converter()
