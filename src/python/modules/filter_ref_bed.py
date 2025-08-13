#!/usr/bin/env python3

"""
Filters reference BED file for further use within the TOGA 2.0 pipeline
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Dict, List, Optional, Set, Tuple, Union

import click

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = ('Bogdan M. Kirilenko')
__all__ = (None)

ALLOWED_CHARSET: Tuple[int] = tuple(
    [
        35, 45, 46, *range(48, 58), *range(65, 91), 95, *range(97, 123), 124
    ]
)
NAME_REJ_REASON: str = (
    'TRANSCRIPT\t{}\t0\tIllegal character used in transcript name\tILLEGAL_NAME\tN'
)
CONTIG_REJ_REASON: str = (
    'TRANSCRIPT\t{}\t0\tLocated outside of user-preferred contigs\tREJECTED_CONTIG\tN'
)
NON_CODING_REJ_REASON: str = (
    'TRANSCRIPT\t{}\t0\tDoes not have a coding sequence\tNON_CODING\tN'
)
FRAME_REJ_REASON: str = (
    'TRANSCRIPT\t{}\t0\tTranscript is out of frame\tOUT_OF_FRAME\tN'
)

def consistent_name(name: str) -> bool:
    """
    Checks whether only the allowed symbols are used in the transcript name
    """
    return all(ord(x) in ALLOWED_CHARSET for x in name)

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'ref_bed',
    type=click.File('r', lazy=True),
    metavar='REF_BED'
)
@click.argument(
    'output',
    type=click.File('w', lazy=True),
    metavar='OUTPUT'
)
@click.argument(
    'rejection_log',
    type=click.File('w', lazy=True),
    metavar='REJECTION_LOG'
)
@click.option(
    '--contigs',
    '-c',
    type=str,
    metavar='CONTIG_LIST',
    default=None,
    show_default=True,
    help='A comma-separated list of contigs to restrict the output transcripts to'
)
@click.option(
    '--excluded_contigs',
    '-e',
    type=str,
    metavar='CONTIG_LIST',
    default=None,
    show_default=True,
    help='A comma-separated list of contigs to exclude from the output'
)
@click.option(
    '--disable_frame_filter',
    '-f',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, does not filter the out-of-frame transcripts out of the the BED file'
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

class AnnotationFilter(CommandLineManager):
    """
    Filter the provided BED12 file of the likely artifacts, namely:\n
    * Non-coding transcripts;\n
    * Transcripts containing symbols other than latin letters, numbers, dots (.),
    underscores (_), and hyphens (-);\n
    * Transcripts with incomplete reading frames (unless the respective flag is set).\n\n
    Arguments are:\n
    * REF_BED is an input BED12 file. The script will ignore blank lines but crash
    if an inconsistent number of fields is encountered;\n
    * OUTPUT is a path to the filtered BED12 file;\n
    * REJECTION_LOG is a path to a two-column tab-separated file containing
    discarded transcripts identifiers and reasons for their rejection.
    """

    __slots__ = [
        'ref_bed', 'output', 'rejection_log', 'approved_contigs', 'excluded_contigs',
        'no_frame_filter', 'v'
    ]

    def __init__(
        self,
        ref_bed: click.File,
        output: click.File,
        rejection_log: click.File,
        contigs: str,
        excluded_contigs: str,
        disable_frame_filter: bool,
        log_name: Optional[str],
        verbose: bool
    ) -> None:

        self.v: bool = verbose
        self.set_logging(log_name)

        self.ref_bed: click.File = ref_bed
        self.output: click.File = output
        self.rejection_log: click.File = rejection_log
        self.approved_contigs: Union[List[str], None] = (
            [x for x in contigs.split(',') if x] if contigs else None
        )
        self.excluded_contigs: Union[List[str], None] = (
            [x for x in excluded_contigs.split(',') if x]
            if excluded_contigs else None
        )
        self.no_frame_filter: bool = disable_frame_filter
        self.run()

    def run(self) -> None:
        """
        """
        transcript_names: Set[str] = set()
        for i, line in enumerate(self.ref_bed, start=1):
            data: List[str] = line.strip().split('\t')
            ## ignore blank lines
            if not data:
                continue
            ## die if inconsistently formatted line is encountered
            if len(data) != 12:
                self._die(
                    f'Input file contains inconsistent number of columns at line {i}. '
                    'Please make sure you are using a tab-separated BED12 file'
                )
            name: str = data[3]
            ## reject the transcripts
            if not consistent_name(name):
                self.rejection_log.write(
                    NAME_REJ_REASON.format(name) + '\n'
                )
                continue
            if name in transcript_names:
                self._die(
                    f'Transcript name {name} appears at least twice in the file. '
                    'Please make sure that the file contains no duplicate entries'
                )
            chrom: str = data[0]
            ## if entries were restricted to specific contigs,
            ## apply the respective filters
            if self.approved_contigs is not None and chrom not in self.approved_contigs:
                self._to_log('Transcript %s is located on a deprecated contig' % chrom, 'warning')
                self.rejection_log.write(
                    CONTIG_REJ_REASON.format(name) + '\n'
                )
                continue
            if self.excluded_contigs is not None and chrom in self.excluded_contigs:
                self._to_log('Transcript %s is located on a deprecated contig' % chrom, 'warning')
                self.rejection_log.write(
                    CONTIG_REJ_REASON.format(name) + '\n'
                )
                continue
            thin_start: int = int(data[1])
            thin_end: int = int(data[2])
            cds_start: int = int(data[6])
            cds_end: int = int(data[7])
            if cds_end < cds_start:
                self._die(
                    'Input file contains an entry with confused start-stop '
                    f'coordinates at line {i}'
                )
            if cds_start == cds_end:
                self.rejection_log.write(
                    NON_CODING_REJ_REASON.format(name) + '\n'
                )
                continue
            ## iterate over exon entries to infer the CDS length
            frame_length: int = 0
            sizes: List[int] = list(map(int, data[10].split(',')[:-1]))
            starts: List[int] = list(map(int, data[11].split(',')[:-1]))
            for start, size in zip(starts, sizes):
                start += thin_start
                end: int = start + size
                # print(f'{start=}, {size=}')
                if start < cds_start:
                    if end > cds_start:
                        # print(f'Partially coding 5-prime block; adding {end - cds_start} to the frame counter')
                        size -= (cds_start - start)
                    else:
                        continue
                if end > cds_end:
                    if start < cds_end:
                        # print(f'Partially coding 5-prime block; adding {cds_end - start} to the frame counter')
                        # frame_length += (cds_end - start)
                        size -= (end - cds_end)
                    else:
                        continue
                frame_length += size
            # print(f'{frame_length=}')
            if frame_length % 3 and not self.no_frame_filter:
                self.rejection_log.write(
                    FRAME_REJ_REASON.format(name) + '\n'
                )
                continue
            self.output.write(line)
            transcript_names.add(name)


if __name__ == '__main__':
    AnnotationFilter()
