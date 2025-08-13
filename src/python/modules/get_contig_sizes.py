#!/usr/bin/env python3

"""
Infers contig size from the FASTA file
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from modules.shared import CONTEXT_SETTINGS

import click

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__all__ = (None)

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'input_',
    type=click.File('r', lazy=True),
    metavar='INPUT_FASTA',
)
@click.option(
    '--output',
    '-o',
    type=click.File('w', lazy=True),
    metavar='OUTPUT_TSV',
    default=sys.stdout,
    show_default=None,
    help='A path to write the output to [default: stdout]'
)

def main(input_: click.File, output: click.File) -> None:
    """
    Retrieves contig sizes from the FASTA file, writing the results
    to a two-column tab-separated file
    """
    if input_ == '-':
        input_ = sys.stdin
    contig_name: str = ''
    counter: int = 0
    for line in input_:
        line = line.rstrip()
        if not line:
            continue
        if line[0] == '>':
            if contig_name and counter:
                output.write(f'{contig_name}\t{counter}\n')
            contig_name = line[1:]
            counter = 0
            continue
        counter += len(line)
    if contig_name and counter:
        output.write(f'{contig_name}\t{counter}\n')

if __name__ == '__main__':
    main()
