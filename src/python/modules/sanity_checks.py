#!/usr/bin/env python3

"""
Sanity check and troubleshooting manager
"""

from shared import CommandLineManager
from typing import Union

class SanityChecker(CommandLineManager):
    __slots__ = []

    def __init__(self) -> None:
        pass

    def check_classification(self) -> None:
        pass

    def check_preprocessing(self) -> None:
        pass

    def check_loss_summary(self) -> None:
        pass

    def check_orthology_resolution(self) -> Union[str, None]: ## TODO: Think how to organize it properly
        pass

    def _email(self) -> None:
        pass


CLASSIFICATION_BOILERPLATE: str = ""
PREPROCESSING_BOILERPLATE: str = ""
LOSS_SUMMARY_BOILERPLATE: str = ""
ORTHOLOGY_BOILERPLATE: str = ""