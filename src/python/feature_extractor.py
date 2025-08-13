#!/usr/bin/env python3

"""
Extracts chain-transcript pair features for orthology classification
"""

from dataclasses import dataclass
from modules.shared import CONTEXT_SETTINGS
from typing import Dict, List

import click
import os

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = [
    'Bogdan M. Kirilenko', 'Michael Hiller', 'Virag Sharma', 'David Jebb'
]
__all__ = [None]

@dataclass
class TranscriptFeatures:
    __slots__ = []


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'chain_id',
    type=str,
    metavar='CHAIN_ID'
)
@click.argument(
    'transcripts',
    type=str,
    metavar='TRANSCRIPT_ID(S)'
)
@click.argument(
    'chain_file',
    type=click.Path(exists=True),
    metavar='CHAIN_FILE'
)
@click.argument(
    'bed_file',
    type=click.Path(exists=True),
    metavar='BED_FILE'
)


class FeatureExtractor:
    __slots__ = ['chain_id', 'transcripts', 'chain_file', 'bed_file']

    def __init__(
        chain_id: str,
        transcripts: str,
        chain_file: click.Path,
        bed_file: click.Path
    ) -> None:
        """
        """
        self.chain_id: str = chain_id
        self.transcripts: List[str] = [x for x in transcripts.split(',') if x]
        self.chain_file: click.Path = chain_file
        self.bed_file: click.Path = bed_file

        self.run()

    def run(self) -> None:
        """
        """
        ## extract chain string, get global chain statistics from the header line
        self.get_chain_string()
        self.extract_chain_features()

        ## extract transcript data from the BED file (to be replace with HDF5 in the future)
        self.parse_annotation_bed()

        ## organize reference transcript data

        ##
