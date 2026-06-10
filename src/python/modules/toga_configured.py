#!/usr/bin/env python3

"""
Wrapper over TogaMain for file-configured runs
"""

import sys
from typing import Any, Dict, List, Optional, TextIO, Tuple, Union

from .constants import TOGA2_ARG2SLOT, TOGA2_SLOT2ARG
from .shared import CommandLineManager, get_upper_dir

__author__ = "Yury V. Malovichko"
__year__ = "2024"
__all__ = None

LOCATION: str = get_upper_dir(__file__, 4)
COL_NUM: int = 2
REF: str = "ref_2bit"
QUERY: str = "query_2bit"
CHAIN: str = "chain_file"
ANNOT: str = "ref_annotation"
MANDATORY_ARGS: Tuple[str, ...] = (REF, QUERY, CHAIN, ANNOT)
TRUE: str = "True"
FALSE: str = "False"
NONE: str = "None"
PROJ_ARG_FILE: str = "project_args.json"
VERSION: str = "version"
START_TIME: str = "start_time"
PARAMETERS: str = "parameters"


class Toga2ConfiguredLauncher(CommandLineManager):
    __slots__ = ("v", "config_file", "config_format", "override", "version")

    def __init__(
        self,
        config_file: TextIO,
        config_format: Optional[str],
        override: Optional[Union[str, None]],
    ) -> None:
        self.v: bool = True
        self.set_logging()
        self.logger.propagate = False
        self.config_file: TextIO = config_file
        self.config_format: str = config_format
        self.override: Union[str, None] = override

        sys.path.append(LOCATION)
        from __version__ import __version__
        self.version: str = __version__
        sys.path.remove(LOCATION)

    def run(self) -> Dict[str, str]:
        if self.config_format == "tsv":
            args: Dict[str, Any] = self.from_tsv()
        elif self.config_format == "yaml":
            args: Dict[str, Any] = self.from_yaml()
        elif self.config_format == "json":
            args: Dict[str, Any] = self.from_json()
        else:
            self._die("Illegal configuration file format: %s" % self.config_format)
        if self.override is not None:
            args = self.parse_override(args)
        return args

    def from_tsv(self) -> Dict[str, Any]:
        cmd_args: Dict[str, str] = {}
        for i, line in enumerate(self.config_file, start=1):
            data: List[str] = line.strip().split("\t")
            if not len(data):
                continue
            if len(data) != 2:
                self._die(
                    "Improper formtting at configuration file line %i; expected 2 columns, got %i"
                    % (i, len(data))
                )
            arg, value = data
            if arg == VERSION:
                self._echo("Specified TOGA2 version is %s" % value)
                if value != self.version:
                    self._to_log(
                        (
                            "Current TOGA2 version is %s but the configuration file was written by/for "
                            "version %s. Certain features listed in the provided file might differ or be "
                            "inaccessible in the current version"
                        )
                        % (self.version, value),
                        "warning",
                    )
                continue
            if arg == START_TIME:
                continue
            if arg not in TOGA2_SLOT2ARG.values():
                self._to_log(
                    "WARNING: Argument/option %s is not recognized; skipping" % arg,
                    "warning",
                )
                continue
            # if arg in MANDATORY_ARGS:
            #     cmd_args[arg] = value
            #     continue
            if value == FALSE:
                value = False
            elif value == TRUE:
                value = True
            elif value == NONE:
                value = None
            elif value.isdigit():
                value = int(value)
            elif value.replace(".", "").isdigit():
                value = float(value)
            # if value == TRUE:
            #     cmd_args[arg] = True
            cmd_args[arg] = value
        return cmd_args

    def from_json(self) -> Dict[str, Any]:
        """Parses the JSON project argument file"""
        import json
        args: Dict[str, Union[str, Dict[any]]] = json.load(self.config_file)
        if PARAMETERS not in args:
            self._die("\"parameters\" is not in the argument file %s" % self.config_file.name)
        return args[PARAMETERS]

    def from_yaml(self) -> Dict[str, Any]:
        """Parses the YAML project argument file"""
        import yaml
        try:
            args: Dict[str, Union[str, Dict[any]]] = yaml.safe_load(self.config_file)
        except yaml.YAMLError as exc:
            self._die("YAML file reading failed with the following error: \n%s" % exc)
        if PARAMETERS not in args:
            self._die("\"parameters\" is not in the argument file %s" % self.config_file.name)
        return args[PARAMETERS]

    def parse_override(self, args: Dict[str, Any]) -> Dict[str, Any]:
        override = self.override.lstrip("\"'").rstrip("\"'")
        overriden_args: List[str] = override.split(" ")
        over_arg_num: int = len(overriden_args)
        curr: int = 0
        while curr < len(overriden_args):
            curr_arg: str = overriden_args[curr]
            curr_arg_name: str = curr_arg.lstrip("-")
            if curr_arg_name not in TOGA2_ARG2SLOT:
                self._echo(
                    "Option %s is not supported by TOGA2 %s"
                    % (curr_arg_name, self.version)
                )
                recognized_arg: bool = False
            else:
                recognized_arg: bool = True
                # slot_name: str = TOGA2_ARG2SLOT[curr_arg_name]
                slot_name: str = curr_arg_name
            ## must be a flag option unless something went wrong
            if curr == over_arg_num - 1:
                if recognized_arg:
                    args[slot_name] = True
                break
            next_arg: str = overriden_args[curr + 1]
            if next_arg.startswith("-"):
                ## next item in the string is also an option, therefore this one must be a flag
                if recognized_arg:
                    args[slot_name] = True
                curr += 1
                continue
            ## otherwise this is an option with an argument
            if recognized_arg:
                if next_arg.isdigit():
                    next_arg = int(next_arg)
                elif next_arg.replace(".", "").isdigit():
                    next_arg = float(next_arg)
                args[slot_name] = next_arg
            curr += 2
            continue
        return args
