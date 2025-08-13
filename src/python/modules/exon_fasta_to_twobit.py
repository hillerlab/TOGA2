#!/usr/bin/env python3

"""

"""

from shared import CommandLineManager, CONTEXT_SETTINGS, hex_code
from shutil import which
from typing import Dict, List, Optional, TextIO, Tuple, Union

import click
import os

__author__ = 'Yury V. Malovichko'
__year__ = '2025'

FA2TWOBIT: str = 'faToTwoBit'
REF_EXON: str = 'reference_exon'
SEP_DUMMY: str = 'n'
ERR: str = 'FASTA to 2bit conversion failed'
EXON_HEADER_SEP: str = ' | '

def prepare_seq(seq: str) -> str:
    """Removes gaps and gained introns from the query exon sequence"""
    return ''.join([x for x in seq if x.isalpha() and x.isupper()])


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'input',
    type=click.File('r',  lazy=True),
    metavar='INPUT'
)
@click.argument(
    'output',
    type=click.Path(exists=False),
    metavar='OUTPUT'
)
@click.option(
    '--tmp_dir',
    '-t',
    type=click.Path(exists=True),
    default=os.getcwd(),
    show_default=False,
    help=(
        'A directory to store the temporary FASTA file to. Ideally has to be something '
        'subjected to automatic cleanup since temporary files from the failed runs '
        'will not be deleted'
    )
)
@click.option(
    '--exon_meta',
    '-e',
    type=click.File('r', lazy=True),
    metavar='EXON_META_FILE',
    default=None,
    show_default=True,
    help=(
        'TEMPORARY: An exon metadata file to retrieve the exon presence status from. In future the data '
        'should be retrieved from the exon FASTA header'
    )
)
@click.option(
    '--fa2twobit',
    type=click.Path(exists=True),
    metavar='FA2TWOBIT_EXEC',
    default=None,
    show_default=True,
    help=(
        'A path to UCSC faToTwoBit executable. If not set, '
        'the executable will be sought for in PATH'
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
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)

class TwoBitConverter(CommandLineManager):
    """
    """
    __slots__ = ('fa2twobit', 'exon2status')

    @staticmethod
    def parse_exon_header(header: str) -> Tuple[str, int, str]:
        """
        Retrieves transcript name, exon number, and referenve/query attribution 
        """
        split_header: List[str] = header.rstrip().lstrip('>').split(EXON_HEADER_SEP)
        proj: str = split_header[0].replace(',', '.')
        exon_num: int = int(split_header[1])
        source: str = split_header[-1]
        return (proj, exon_num, source)

    def __init__(
        self,
        input: click.Path,
        output: click.Path,
        tmp_dir: Optional[click.Path],
        exon_meta: Optional[Union[click.Path, None]],
        fa2twobit: Optional[Union[click.Path, None]],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)
        self.fa2twobit: Union[str, None] = None
        self.set_twobit_path(fa2twobit)

        self.exon2status: Dict[Tuple[str, int], str] = {}
        self.parse_exon_meta(exon_meta)

        tmp_file: str = 'exon_fasta_for_compression' + hex_code()
        tmp_path: str = os.path.join(tmp_dir, tmp_file)

        with open(tmp_path, 'w') as h:
            transcript: str = ''
            prev_tr: str = ''
            exon_num: int = 0
            prev_exon_num: int = 0
            sequences: str = ''
            # aggr_results: str = ''
            for line in input:
                line = line.rstrip()
                if line[0] == '>':
                    transcript, exon_num, source = self.parse_exon_header(line)
                    # print(f'{transcript=}, {exon_num=}, {source=}, {prev_tr=}')
                    if prev_tr:
                        if prev_tr != transcript:
                            # print('Proceeding to another transcript')
                            out_header: str = '>' + prev_tr
                            h.write(out_header + '\n' + sequences + '\n')
                            prev_exon_num = 0
                            sequences = ''
                    prev_tr = transcript
                    # transcript, exon_num, source = self.parse_exon_header()
                    if source == REF_EXON:
                        transcript = ''
                        exon_num = 0
                        prev_exon_num = 0
                        continue
                    if prev_exon_num and prev_exon_num != exon_num - 1:
                        self._die(
                            'Exon FASTA files was improperly sorted'
                        )
                    ## new exon started -> add a separator
                    if exon_num > 1:
                        sequences += 'n'
                    prev_exon_num = exon_num
                    continue
                if not transcript:
                    continue
                if self.exon2status:
                    if self.exon2status.get((transcript, exon_num), 'D') != 'I':
                        continue
                sequences += prepare_seq(line)#line.replace('-','')
            if transcript:
                out_header: str = '>' + transcript
                # aggr_results += out_header + '\n' + sequences + '\n'
                h.write(out_header + '\n' + sequences + '\n')
        # print(aggr_results)
        # aggr_results = aggr_results.encode('utf8')
        cmd = f'faToTwoBit {tmp_path} {output}'
        # print(cmd)
        _ =  self._exec(cmd, ERR)
        self._rm(tmp_path)

    def set_twobit_path(self, twobit_exe: str) -> None:
        """
        Sets UCSC faToTwoBit executable; if not provided, will search for executable in PATH
        """
        if twobit_exe is not None:
            self.fa2twobit: str = twobit_exe
            return
        exe_in_path: str = which(FA2TWOBIT)
        if exe_in_path is not None:
            self._to_log('Found faToTwoBit in PATH at %s' % exe_in_path)
            self.fa2twobit: str = exe_in_path
            return
        self._die(
            'faToTwoBit executable was not found in PATH, with no default path provided'
        )

    def parse_exon_meta(self, file: TextIO) -> None:
        """Retrieves exon loss status from the exon_meta.tsv file"""
        if file is None:
            return
        for line in file:
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if data[0] == 'projection':
                continue
            proj: str = data[0].replace(',', '.')
            exon_num: int = int(data[1])
            status: str = data[7]
            self.exon2status[(proj, exon_num)] = status


if __name__ == '__main__':
    TwoBitConverter()