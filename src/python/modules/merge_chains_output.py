#!/usr/bin/env python3
"""Parse raw chain runner output.

Chain features extraction steps results in numerous files.
This script merges these files and then
builds s table containing chain features.
"""

import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, TextIO, Union

# from version import __version__
import click

from .shared import CONTEXT_SETTINGS, CommandLineManager

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

__author__ = ""
__credits__ = "Bogdan M. Kirilenko"

GENES: str = "genes"
CHAIN: str = "chain"
TIMESTAMP: str = "#estimated"
HEADER: str = "\t".join(
    (
        "transcript",
        "gene_overs",
        "chain",
        "synt",
        "gl_score",
        "gl_exo",
        "chain_len",
        "exon_qlen",
        "loc_exo",
        "exon_cover",
        "intr_cover",
        "gene_len",
        "ex_num",
        "ex_fract",
        "intr_fract",
        "flank_cov",
        "clipped_exon_qlen",
        "clipped_intr_cover",
    )
)
OK: str = "ok"


@dataclass
class TranscriptFeatures:
    __slots__ = ("gene_len", "exon_fraction", "intron_fraction", "exons_num")
    gene_len: int
    exon_fraction: float
    intron_fraction: float
    exons_num: int


@dataclass
class ChainFeatures:
    __slots__ = (
        "synteny",
        "global_score",
        "global_exons",
        "chain_q_len",
        "chain_len",
        "local_exons",
        "coverages",
        "introns",
        "flank_cov",
        "clipped_chain_qlen",
        "clipped_exons",
        "clipped_introns",
    )
    synteny: str
    global_score: str
    global_exons: str
    chain_q_len: str
    chain_len: str
    local_exons: Dict[str, str]
    coverages: Dict[str, str]
    introns: Dict[str, str]
    flank_cov: Dict[str, str]
    clipped_chain_qlen: str
    clipped_exons: Dict[str, str]
    clipped_introns: Dict[str, str]


def parse_pairs(pairs: str) -> Dict[str, str]:
    """Parse lines like X=5,Y=56, and returns a dict."""
    return {
        k.strip(): v.strip() for k, v in (x.split("=") for x in pairs.split(",") if x)
    }


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    "results_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    metavar="RESULTS_DIR",
)
@click.argument("bed_file", type=click.File("r", lazy=True), metavar="BED_FILE")
@click.argument("output", type=click.File("w", lazy=True))
@click.option(
    "--isoforms",
    "-i",
    type=click.File("r", lazy=True),
    metavar="FILE",
    default=None,
    show_default=True,
    help=("A file containing reference gene-to-isoform mapping"),
)
@click.option(
    "--exon_cov_chains",
    "-e",
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        "If set, only the projections with at least one reference base covered "
        "are considered for further analysis"
    ),
)
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
    metavar="FLAG",
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls the execution verbosity",
)
class FeatureAggregator(CommandLineManager):
    """
    Merges output files from the projection feature extraction step. Arguments are:\n
    * RESULTS_DIR is a directory containing output files from the feature extraction step;\n
    * BED_FILE is a reference annotation in BED format;\n
    * OUTPUT is a path to the output file
    """

    __slots__ = (
        "ref_bed",
        "isoform2gene",
        "gene2chains",
        "chain_feature_data",
        "only_covering",
    )

    def __init__(
        self,
        results_dir: click.Path,
        bed_file: click.File,
        output: click.File,
        isoforms: Optional[Union[click.Path, None]],
        exon_cov_chains: Optional[bool],
        log_name: Optional[str],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.set_logging(name=log_name, toga_module="merge_chain_features")

        self.only_covering: bool = exon_cov_chains
        self.ref_bed: Dict[str, TranscriptFeatures] = {}
        self._to_log("Parsing reference BED file")
        self.parse_ref_bed(bed_file)
        self.isoform2gene: Dict[str, List[str]] = {}
        self.parse_ref_isoforms(isoforms)
        self.gene2chains: Dict[str, str] = defaultdict(list)
        self.chain_feature_data: Dict[str, List[str]] = {}
        self.load_features(results_dir)
        self._to_log("Feature data successfully aggregated")
        self._to_log("Writing output data")
        self.write_combined_output(output)

    def parse_ref_bed(self, file: TextIO) -> None:
        """Retrieves the necessary features from the reference BED file"""
        """Get the necessary data from the bed file."""
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 12:
                self._die(
                    "Reference BED file contains less than twelve fiels at line %i" % i
                )
            if data[1] == data[2] or data[6] == data[7]:
                self._to_log(
                    "Reference BED file contains an entry with no coding sequence at line %i"
                    % i,
                    "warning",
                )
            thin_start: int = int(data[1])
            thin_end: int = int(data[2])
            name: str = data[3]
            cds_start: int = int(data[6])
            cds_end: int = int(data[7])
            exon_num: int = 0
            exon_sum: int = 0
            cds_sum: int = 0
            block_sizes: List[int] = [int(x) for x in data[10].split(",") if x]
            block_starts: List[int] = [int(x) for x in data[11].split(",") if x]
            if len(block_sizes) != len(block_starts):
                self._die(
                    "Number of block starts does not equal the number of block sizes at line %i"
                    % i
                )
            for i in range(len(block_sizes)):
                block_start: int = thin_start + block_starts[i]
                block_end: int = block_start + block_sizes[i]
                exon_sum += block_sizes[i]
                exon_num += 1
                if block_end <= cds_start or block_start >= cds_end:
                    continue
                block_start = max(block_start, cds_start)
                block_end = min(block_end, cds_end)
                cds_sum += block_end - block_start
            gene_len: int = thin_end - thin_start
            intron_sum: int = gene_len - exon_sum
            self.ref_bed[name] = TranscriptFeatures(
                gene_len, cds_sum, intron_sum, exon_num
            )
        self._to_log(
            "Extracted classification data for %i reference transcripts"
            % len(self.ref_bed)
        )

    def parse_ref_isoforms(self, file: Union[TextIO, None]) -> None:
        """Parses reference isoform file"""
        if file is None:
            self._to_log("No isoform mapping file was provided")
            return
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 2:
                self._die(
                    "Isoforms file contained less than two columns at line %i" % i
                )
            gene: str = data[0]
            tr: str = data[1]
            self.isoform2gene[tr] = gene

    def load_features(self, results_dir: click.Path) -> None:
        """Aggregates the results over multiple chain_runner.py runs stored in a single directory"""
        for file in os.listdir(results_dir):
            ## ignore the successful execution stamps
            if file.split("_")[-1] == OK:
                continue
            filepath: str = os.path.join(results_dir, file)
            with open(filepath, "r") as h:
                for i, line in enumerate(h, start=1):
                    data: List[str] = line.rstrip().split("\t")
                    if not data or not data[0]:
                        continue
                    if data[0] == GENES:
                        split_gene_data: List[str] = [
                            x.split("=") for x in data[1:] if "=" in x
                        ]
                        for gene, chain in split_gene_data:
                            self.gene2chains[gene].append(chain)
                    elif data[0] == CHAIN:
                        chain_id: str = data[1]
                        synteny: str = data[2]
                        global_score: str = data[3]
                        global_exon: str = data[4]
                        chain_q_len: str = data[5]
                        chain_len: str = data[10]
                        local_exons: Dict[str, str] = parse_pairs(data[6])
                        coverages: Dict[str, str] = parse_pairs(data[7])
                        introns: Dict[str, str] = parse_pairs(data[8])
                        flank_cov: Dict[str, str] = parse_pairs(data[9])
                        clipped_chain_qlen: str = data[11]
                        clipped_exons: Dict[str, str] = parse_pairs(data[12])
                        clipped_introns: Dict[str, str] = parse_pairs(data[13])
                        self.chain_feature_data[chain_id] = ChainFeatures(
                            synteny,
                            global_score,
                            global_exon,
                            chain_q_len,
                            chain_len,
                            local_exons,
                            coverages,
                            introns,
                            flank_cov,
                            clipped_chain_qlen,
                            clipped_exons,
                            clipped_introns,
                        )
                    elif TIMESTAMP in line:
                        continue
                    else:
                        self._die(
                            "Erroneous formatting at line %i in file %s" % (i, filepath)
                        )

    def write_combined_output(self, output: TextIO) -> None:
        """Combines the loaded data into output suitable for TOGA2 use"""
        output.write(HEADER + "\n")
        for chain, chain_features in self.chain_feature_data.items():
            trs: List[str] = chain_features.local_exons.keys()
            synteny: str = (
                chain_features.synteny
                if not self.isoform2gene
                else str(self._get_synteny(trs))
            )
            for tr in trs:
                tr_features: Union[TranscriptFeatures, None] = self.ref_bed.get(
                    tr, None
                )
                if tr_features is None:
                    self._to_log(
                        "Transcript %s is missing from the reference annotation data"
                        % tr
                    )
                    continue
                local_exon_score: str = chain_features.local_exons[tr]
                exon_coverage: str = str(chain_features.coverages[tr])
                intron_coverage: str = str(chain_features.introns[tr])
                flank_cov: str = str(chain_features.flank_cov[tr])
                # get a number of chains that covert this gene
                # else: if you need the chains that overlap EXONS
                gene_overs = (
                    len(self.gene2chains[tr])
                    if not self.only_covering
                    else len([x for x in self.gene2chains[tr] if x != "None"])
                )
                ## get additional features for processed pseudogene classification
                clipped_introns: str = chain_features.clipped_introns[tr]
                line_data: List[str] = [
                    tr,
                    gene_overs,
                    chain,
                    synteny,
                    chain_features.global_score,
                    chain_features.global_exons,
                    chain_features.chain_len,
                    chain_features.chain_q_len,
                    local_exon_score,
                    exon_coverage,
                    intron_coverage,
                    tr_features.gene_len,
                    tr_features.exons_num,
                    tr_features.exon_fraction,
                    tr_features.intron_fraction,
                    flank_cov,
                    chain_features.clipped_chain_qlen,
                    clipped_introns,
                ]
                output.write("\t".join(map(str, line_data)) + "\n")

    def _get_synteny(self, isoforms: List[str]) -> int:
        """Returns the number of unique isoforms for a given set of genes"""
        return len({self.isoform2gene[i] for i in isoforms if i in self.isoform2gene})


if __name__ == "__main__":
    # main()
    FeatureAggregator()
