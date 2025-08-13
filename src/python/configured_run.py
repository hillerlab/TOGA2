#!/usr/bin/env python3

"""
Runs TOGA2 with project arguments listed in the config file
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from modules.constants import TOGA2_SLOT2ARG
from modules.shared import CommandLineManager, CONTEXT_SETTINGS
from toga2 import __version__, TogaMain
from typing import Tuple, Dict, List

import click

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__all__ = (None)

TOGA2: str = os.path.join(LOCATION, 'toga2.py')
COL_NUM: int = 2
REF: str = 'ref_2bit'
QUERY: str = 'query_2bit'
CHAIN: str = 'chain_file'
ANNOT: str = 'ref_annotation'
MANDATORY_ARGS: Tuple[str] = (REF, QUERY, CHAIN, ANNOT)
TRUE: str = 'True'
FALSE: str = 'False'
NONE: str = 'None'
PROJ_ARG_FILE: str = 'project_args.json'
VERSION: str = 'version'

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'config_file',
    type=click.File('r', lazy=True)
)

class Toga2ConfiguredLauncher(CommandLineManager):
    """
    Runs local TOGA2 version with settings listed in a configuration file provided.\n
    The sole argument is a two-column configuration file with TOGA2 arguments/options
    in the first column and their values in the second column.\n
    The following arguments are mandatory and must be present
    in the configuration file: {REF}, {QUERY}, {CHAIN}, {ANNOT}.\n
    For flag options, provide '{TRUE}' to set the flag, '{FALSE}' otherwise.\n
    For disabled non-flag options and options expected to be used with default values,
    provide '{NONE}'.\n
    By default, TOGA2 outputs project setup file '{PROJ_ARG_FILE}' which can be used
    to reproduce previous runs. Note that argument number and names may change
    between TOGA2 releases
    """

    __slots__ = ('v')

    def __init__(self, config_file: click.File) -> None:
        self.v: bool = True
        args: Dict[str, str] = {}
        # options: str = ''
        options: List[str] = []
        for line in config_file:
            data: List[str] = line.rstrip().split('\t')
            if len(data) != COL_NUM:
                continue
            arg, value = data
            if arg == VERSION:
                self._echo(f'Specified TOGA2 version is {arg}')
                if arg != __version__:
                    self._echo(
                        f'WARNING: Local TOGA2 version is {__version__} but '
                        f'the config file was prepared for using with version {arg}; '
                        'certain settings might be incompatible between the versions!'
                    )
                continue
            if arg not in TOGA2_SLOT2ARG.values():
                self._echo(f'WARNING: Argument/option {arg} is not recognized; skipping')
                continue
            if arg in MANDATORY_ARGS:
                args[arg] = value
                continue
            if value == FALSE:
                continue
            if value == NONE:
                continue
            if value == TRUE:
                # options += f' --{arg}'
                options.append(f'--{arg}')
                continue
            # options += f' --{arg} {value}'
            options.extend([f'--{arg}', value])
        missing_args: List[str] = [x for x in MANDATORY_ARGS if x not in args]
        if missing_args:
            missing_arg_line: str = ','.join(missing_args)
            self._die(
                'ERROR: the following mandatory arguments are missing from '
                f'the configuration file: {missing_arg_line}'
            )
        cmd: str = (
            f'{TOGA2} {args[REF]} {args[QUERY]} {args[CHAIN]} {args[ANNOT]} {options}'
        )
        cmd_args: List[str] = [args[REF], args[QUERY], args[CHAIN], args[ANNOT]] + options

        self._echo(f'Running the following TOGA2 command:\n{cmd}')
        TogaMain(cmd_args)
        # self._exec(
        #     cmd,
        #     'ERROR: TOGA2 run failed:',
        #     shun_verbosity=False
        # )

if __name__ == '__main__':
    Toga2ConfiguredLauncher()
