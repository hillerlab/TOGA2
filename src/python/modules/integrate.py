#!/usr/bin/env python3

"""
Integrates TOGA2 results obtained with multiple references for a single query
to produce a combined, multi-reference TOGA2 annotation
"""

import gzip
import json
import os
from collections import defaultdict
from dataclasses import dataclass
from shutil import which
from typing import Dict, Iterable, List, Optional, Set, TextIO, Tuple, Union

import logging
import networkx as nx

from .cesar_wrapper_constants import CLASS_TO_COL, CLASS_TO_NUM
from .constants import Constants
from .results_checks import BREAK_LINE, to_perc
from .shared import (
    CommandLineManager,
    base_proj_name,
    get_proj2trans,
    get_upper_dir,
    make_cds_track,
    intersection,
    read_tab,
)

TOGA2_ROOT: str = get_upper_dir(__file__, 4)
BIN: str = os.path.join(TOGA2_ROOT, "bin")
DEFAULT_BIGBED2BED: str = os.path.join(BIN, "bigBedToBed")
LOCATION: str = os.path.dirname(os.path.abspath(__file__))
MAKE_IX_SCRIPT: str = os.path.join(LOCATION, "get_names_from_bed.py")
TOGA2_DIR: str = get_upper_dir(__file__, 4)
GTF_SCRIPT: str = os.path.join(TOGA2_DIR, "bed2gtf", "target", "release", "bed2gtf")

ALL: str = "ALL"
BED_FIELD_NUM: int = 12
COL2CLASS: Dict[str, str] = {y: x for x, y in CLASS_TO_COL.items()}
EXON_HEADER: str = "projection"
INTACT_EXON: str = "I"
INTACT_PROJECTIONS: Tuple[str, str] = ("FI", "I")
SUPPORTED: str = "CHAIN_SUPPORTED"
PLUS: str = "+"

SUMMARY: str = """{br}
#### TOGA2 summary
{br}
The combined TOGA2 annotation contains {num_genes} genes and {num_retro} retrogene candidates.
In addition, the annotation contains {num_lost} genes classified as lost and {num_missing} genes classified as missing in the query, and {num_para} additional genes were annotated using potentially functional paralogous transcript annotations.

#HEADER\tQueryAssembly\tno. annotated query genes\tno. query retrogenes\tno. lost genes in query\tno. missing genes in query\tno. query paralogs
#SINGLELINESUMMARY\t{query}\t{num_genes}\t{num_retro}\t{num_lost}\t{num_missing}\t{num_para}

The following reference species were used for the integrated annotation preparation:
{refs}

For questions, please check whether the issue has already been addressed at https://github.com/hillerlab/TOGA2/issues. If not, please open a new issue.

If you use these data, please cite: 
Yury V. Malovichko et al. "Accurate, comprehensive gene annotation and ortholog identification across thousands of vertebrate genomes with TOGA2", in preparation
"""

REF_SUPPORT_LINE: str = "* {ref}: appears in {ref_genes}/{all_genes} ({gene_perc}%) orthologous query loci\n"
GENE_TSV_HEADER: str = "query_gene\tprojection\n"

## check the NetworkX version
nx_v: str = nx.__version__
v_split: List[str] = [x for x in nx_v.split(".") if x.isnumeric()]
if len(v_split) > 1:
    NX_VERSION: float = float(f"{v_split[0]}.{v_split[1]}")
else:
    NX_VERSION: float = float(v_split[0])

logging.root.handlers = []


def safe_div(numerator: int, denominator: int) -> float:
    """
    Divides numerator by denominator, obviating the division by zero scenario.

    Args:
        numerator: a number to divide
        denominator: a number to divide by

    Returns:
        Zero if denominator is zero, numerator divided by denominator otherwise
    """
    if denominator == 0:
        return 0.0
    return numerator / denominator


class ReferenceBundle:  ## initialized from a JSON object
    """
    Data class for reference- and query annotation data
    bound to a certain TOGA2 reference
    """

    __slots__ = (
        "query_bed",
        "exon_meta",
        "ref_isoforms",
        "ref_bed",
        "ucsc_bigbed",
        "decorator_bigbed",
        "protein_file",
        "nucleotide_file",
        "priority",
    )

    def __init__(
        self,
        kwargs: Dict[str, Union[str, os.PathLike, int]],
        priority: Optional[int] = None,
    ) -> None:
        self.query_bed: str = kwargs["query_bed"]
        self.exon_meta: str = kwargs["exon_meta"]
        self.ref_isoforms: Union[str, None] = kwargs.get("reference_isoforms", None)
        self.ref_bed: Union[str, None] = kwargs.get("reference_bed", None)
        self.ucsc_bigbed: Union[str, None] = kwargs.get("ucsc_bigbed", None)
        self.decorator_bigbed: Union[str, None] = kwargs.get("decorator_bigbed", None)
        self.protein_file: Union[str, None] = kwargs.get("protein_file", None)
        self.nucleotide_file: Union[str, None] = kwargs.get("nucleotide_file", None)
        self.priority: int = kwargs.get("priority", priority)


@dataclass
class BedRecord:
    """
    Data class containing data for a single BED file record
    """

    __slots__ = (
        "name",
        "ref",
        "chrom",
        "start",
        "end",
        "strand",
        "loss_status",
        "lines",
        "cds_lines",
        "exons",
    )

    def __init__(
        self,
        name: str,
        ref: str,
        chrom: str,
        start: int,
        end: int,
        strand: bool,
        loss_status: str,
        # lines: List[str],
        lines: Dict[str, str],
    ) -> None:
        self.name: str = name
        self.ref: str = ref
        self.chrom: str = chrom
        self.start: int = start
        self.end: int = end
        self.strand: bool = strand
        self.loss_status: str = loss_status
        self.lines: Dict[str, str] = lines
        self.cds_lines: Dict[str, str] = {
            x: make_cds_track(y) for x, y in lines.items()
        }
        self.exons: List[str] = []

    def return_bed_line(self, prefix: Union[str] = "") -> Iterable[str]:
        """Returns the initial BED line for the projection"""
        for num, line in self.lines.items():
            name = self.name
            if num != 0:
                name += f"${num}"
            yield line.format(
                (prefix + "." + name) if prefix else name,
                CLASS_TO_COL[self.loss_status],
            )

    def coords(self) -> Tuple[int, int]:
        return (self.start, self.end)


@dataclass
class ExonRecord:
    """Exon coordinate record"""

    ## TODO: Salvage a similar class from preprocessing code
    __slots__ = ("projection", "num", "chrom", "start", "end", "strand")

    projection: str
    num: int
    chrom: str
    start: int
    end: int
    strand: bool

    def length(self) -> int:
        return self.end - self.start

    def coords(self) -> Tuple[int, int]:
        return (self.start, self.end)


class AnnotationIntegrator(CommandLineManager):
    __slots__ = (
        "ref_data",
        "query_projections",
        "query_proj2ref",
        "query_annotation",
        "query_name",
        "ref_proj2gene",
        "intersecting_ref_genes",
        "paralog_pool",
        "ppgene_pool",
        "graph",
        "discarded_items",
        "final_projections",
        "accepted_statuses",
        "paralog_rel_novelty_threshold",
        "paralog_abs_novelty_threshold",
        "lost_rel_novelty_threshold",
        "lost_abs_novelty_threshold",
        "output",
        "gene_tsv",
        "gene_bed",
        "projection_bed",
        "protein_file",
        "nucleotide_file",
        "ucsc_dir",
        "bigbed_stub",
        "decorator_stub",
        "bigbed",
        "decorator",
        "ix",
        "ixx",
        "prefix",
        "skip_ucsc",
        "has_ucsc_data",
        "bigbedtobed_binary",
        "bedtobigbed_binary",
        "ixixx_binary",
        "schema",
        "decorator_schema",
        "chrom_sizes",
        "bed_index",
        "query_gtf",
        "query_genes_for_gtf",
        "summary",
        "num_genes",
        "num_retro",
        "num_lost",
        "num_missing",
        "num_para",
        "gene2species",
        "support_matrix_file",
        "species2num_item",
        "species2num_genes",
        "v",
    )

    def __init__(
        self,
        ref_data: Union[str, os.PathLike],
        output: Union[str, os.PathLike],
        query_name: Union[str],
        accepted_statuses: str,
        paralog_rel_novelty_threshold: float,
        paralog_abs_novelty_threshold: int,
        lost_rel_novelty_threshold: float,
        lost_abs_novelty_threshold: int,
        prefix: str,
        skip_ucsc: bool,
        chrom_sizes: Union[str, os.PathLike, None],
        bigbedtobed_binary: Union[str, os.PathLike],
        bedtobigbed_binary: Union[str, os.PathLike],
        ixixx_binary: Union[str, os.PathLike],
        verbose: Optional[bool] = False,
    ) -> None:
        self.v: bool = verbose
        self.set_logging()

        ## parse the input dictionary
        self.ref_data: Dict[str, ReferenceBundle] = {}
        with open(ref_data, "r") as h:
            bundles: Dict[str, Dict[str, str]] = json.loads(h.read())
            for i, (species, data) in enumerate(bundles.items()):
                self.ref_data[species] = ReferenceBundle(data, priority=i)
        if accepted_statuses == ALL:
            self.accepted_statuses: List[str] = Constants.ALL_LOSS_SYMBOLS
        else:
            self.accepted_statuses: List[str] = [
                x for x in accepted_statuses.split(",") if x
            ]  ## TODO: Add sanity checks
        self.paralog_rel_novelty_threshold: float = paralog_rel_novelty_threshold
        self.paralog_abs_novelty_threshold: int = paralog_abs_novelty_threshold
        self.lost_rel_novelty_threshold: float = lost_rel_novelty_threshold
        self.lost_abs_novelty_threshold: int = lost_abs_novelty_threshold
        self.query_name: str = query_name
        self.query_projections: Dict[str, BedRecord] = {}
        self.query_proj2ref: Dict[str, str] = {}
        self.query_annotation: Dict[str, List[str]] = defaultdict(list)
        self.ref_proj2gene: Dict[str, str] = {}
        self.intersecting_ref_genes: Dict[str, Set[str]] = defaultdict(set)
        self.paralog_pool: Set[str] = set()
        self.ppgene_pool: Set[str] = set()
        self.graph: nx.Graph = nx.Graph()
        self.discarded_items: Set[str] = set()
        self.final_projections: Set[str] = set()
        self.gene2species: Dict[str, Set[str]] = {}

        self.skip_ucsc: bool = skip_ucsc
        self.bigbedtobed_binary: Union[str, os.PathLike, None] = bigbedtobed_binary
        self.bedtobigbed_binary: Union[str, os.PathLike, None] = bedtobigbed_binary
        self.ixixx_binary: Union[str, os.PathLike, None] = ixixx_binary
        self.chrom_sizes: Union[str, os.PathLike, None] = chrom_sizes
        if not skip_ucsc:
            self._check_binaries()
        if self.chrom_sizes is None and not (
            self.skip_ucsc or all(x.ucsc_bigbed is None for x in self.ref_data.values())
        ):
            self._die("Query chromosome size file was not provided")
        self.schema: str = os.path.join(TOGA2_ROOT, "supply", "bb_schema_toga2.as")
        self.decorator_schema: str = os.path.join(TOGA2_ROOT, "supply", "decoration.as")
        self.has_ucsc_data: List[str] = []

        self.output: os.PathLike = output
        self.gene_tsv: os.PathLike = os.path.join(output, "query_genes.tsv")
        self.gene_bed: os.PathLike = os.path.join(output, "query_genes.bed")
        self.projection_bed: os.PathLike = os.path.join(output, "query_annotation.bed")
        self.query_gtf: os.PathLike = os.path.join(output, "query_annotation.gtf")
        self.query_genes_for_gtf: os.PathLike = os.path.join(output, "query_genes.for_gtf.tsv")
        self.summary: os.PathLike = os.path.join(output, "summary.txt")
        self.support_matrix_file: os.PathLike = os.path.join(output, "query_genes_by_reference.tsv")
        self.ucsc_dir: os.PathLike = os.path.join(output, "ucsc_browser_files")

        self.protein_file: str = os.path.join(output, "protein.fa.gz")
        self.nucleotide_file: str = os.path.join(output, "nucleotide.fa.gz")

        self.bigbed_stub: str = os.path.join(self.ucsc_dir, f"{prefix}.bed")
        self.decorator_stub: str = os.path.join(
            self.ucsc_dir, f"{prefix}.decorator.bed"
        )
        self.bigbed: str = os.path.join(self.ucsc_dir, f"{prefix}.bb")
        self.decorator: str = os.path.join(self.ucsc_dir, f"{prefix}.decorator.bb")
        self.ix: str = os.path.join(self.ucsc_dir, f"{prefix}.ix")
        self.ixx: str = os.path.join(self.ucsc_dir, f"{prefix}.ixx")
        self.bed_index: str = os.path.join(self.ucsc_dir, f"{prefix}.ix.txt")

        self.num_genes: int = 0
        self.num_retro: int = 0
        self.num_lost: int = 0
        self.num_missing: int = 0
        self.num_para: int = 0
        self.species2num_item: Dict[str, int] = defaultdict(int)
        self.species2num_genes: Dict[str, int] = defaultdict(int)

    def run(self) -> None:
        """Main execution method"""
        self._mkdir(self.output)
        self.process()
        self.infer_genes()
        self.pick_best_isoforms()
        self.prepare_gtf()
        self.generate_support_matrix()
        self.generate_summary()
        if not self.skip_ucsc:
            self._mkdir(self.ucsc_dir)
        self.prepare_ucsc_file()
        self.prepare_sequence_files()

    def process(self) -> None:
        """Extract data for every reference directory"""
        for species in self.ref_data:
            self.read_annotation(species)
            self.read_exon_meta(species)
            self.read_ref_isoforms(species)
            self.get_overlapping_genes(species)

    def _check_binaries(self) -> None:
        """
        Checks the provided UCSC binaries,
        looks for an alternative if needed
        """
        for attr, default_name in Constants.BINARIES_TO_CHECK.items():
            if attr not in self.__slots__:
                continue
            binary: Union[str, None] = getattr(self, attr)
            if binary is not None:
                self._to_log("Testing %s binary at %s" % (default_name, binary))
                if os.access(binary, os.X_OK):
                    self._to_log(
                        "The provided binary is executable; using the stated %s instance"
                        % default_name
                    )
                    continue
                else:
                    self._to_log(
                        (
                            "%s binary at %s does not seem to be executable; "
                            "looking for alternatives"
                        )
                        % (default_name, binary),
                        "warning",
                    )
            else:
                self._to_log(
                    (
                        "No %s executable was not provided; "
                        "looking for alternatives in $PATH"
                    )
                    % default_name
                )
            self._to_log(
                "Checking the default %s executable at TOGA2/bin" % default_name
            )
            default_path: str = os.path.join(BIN, default_name)
            if os.path.exists(default_path):
                self._to_log(
                    "Found default %s executable at %s" % (default_name, default_path)
                )
                if os.access(default_path, os.X_OK):
                    self._to_log(
                        "The default UCSC %s binary is executable; using the default version"
                        % default_name
                    )
                    self.__setattr__(attr, default_path)
                    continue
                else:
                    self._to_log(
                        (
                            "Default UCSC %s binary at %s does not seem to be executable; "
                            "looking for alternatives in $PATH"
                        )
                        % (default_name, default_path),
                        "warning",
                    )
            binary_in_path: Union[str, None] = which(default_name)
            if binary_in_path is not None:
                self._to_log(
                    (
                        "Found bigBedTobed instance at %s; "
                        "checking the execution permissions"
                    )
                    % binary_in_path
                )
                if os.access(binary_in_path, os.X_OK):
                    self._to_log(
                        "The found binary is executable; using the found %s instance"
                        % default_name
                    )
                    self.__setattr__(attr, binary_in_path)
                    continue
                self._die(
                    (
                        "The %s binary found in $PATH at %s is not executable; "
                        "check your $PATH or provide a valid bigBedToBed instance"
                    )
                    % (default_name, binary_in_path)
                )
            self._die(
                (
                    "No %s binary found in $PATH; "
                    "check your $PATH or provide a valid bigBedToBed instance"
                )
                % default_name
            )

    def read_annotation(self, species: str) -> None:
        """
        Parses the reference BED file, filtering the records by their loss status
        """
        self._to_log("Reading annotation file for reference %s" % species)
        file: str = self.ref_data[species].query_bed
        if not os.path.exists(file):
            self._die("File %s does not exist" % file)
        with open(file, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) != BED_FIELD_NUM:
                    self._die(
                        (
                            "Improper formatting at BED file %s line %i; "
                            "expected %i fields, got %i"
                        )
                        % (file, i, BED_FIELD_NUM, len(data))
                    )
                color: str = data[8]
                status: str = COL2CLASS.get(color)
                if status is None:
                    self._die(
                        ("Invalid loss status color at BED file %s line %i: %s")
                        % (file, i, color)
                    )
                # if status not in self.accepted_statuses:
                #     continue
                # name: str = data[3]
                name: str = f"{species}.{data[3]}"
                if "," in name:
                    # name = name.split("$")[0]
                    name, segment = name.split("$")
                    segment = int(segment)
                else:
                    segment = 0
                chrom: str = data[0]
                start: int = int(data[6])
                end: int = int(data[7])
                strand: bool = data[5] == "+"
                line_template: str = "\t".join(
                    [*data[:3], "{}", *data[4:8], "{}", *data[9:]]
                )
                if name in self.query_projections:
                    if "," in name:
                        self.query_projections[name].lines[segment] = line_template
                    else:
                        self._die(
                            "Duplicated non-fragmented entry for reference %s: %s"
                            % (species, name)
                        )
                else:
                    record: BedRecord = BedRecord(
                        name,
                        species,
                        chrom,
                        start,
                        end,
                        strand,
                        status,
                        # [line_template],
                        {segment: line_template},
                    )
                    self.query_projections[name] = record
                self.query_proj2ref[name] = species
                self.query_annotation[chrom].append(name)
                if "#paralog" in name:
                    self.paralog_pool.add(base_proj_name(name))
                    # self.paralog_pool.add(name)
                if "#retro" in name:
                    self.ppgene_pool.add(base_proj_name(name))
                    # self.ppgene_pool.add(name)

    def read_exon_meta(self, species) -> None:
        """
        Extracts exon coordinates from the exon metadata file for a given species,
        removing missing/deleted and extrapolated exons
        """
        self._to_log("Reading exon metadata file for reference %s" % species)
        file: str = self.ref_data[species].exon_meta
        gzipped: bool = file.endswith(".gz")
        if not os.path.exists(file):
            self._die("File %s does not exist" % file)
        with gzip.open(file, "rb") if gzipped else open(file, "r") as h:
            for line in h:
                if gzipped:
                    line: str = line.decode("utf8")
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                ## TODO: Array length check
                if data[0] == EXON_HEADER:
                    continue
                # proj: str = data[0]
                proj: str = f"{species}.{data[0]}"
                if proj in self.paralog_pool:
                    proj += "#paralog"
                elif proj in self.ppgene_pool:
                    proj += "#retro"
                ## do not include projections from transcripts outside of the final annotation files
                if proj not in self.query_projections:
                    continue
                status: str = data[7]
                ## do not account for missing and deleted exons
                if status != INTACT_EXON:
                    continue
                support: str = data[19]
                ## do not account for extrapolated exons
                if support != SUPPORTED:
                    continue
                chrom: str = data[3]
                start: int = int(data[4])
                end: int = int(data[5])
                strand: bool = data[6] == PLUS
                exon_num: int = int(data[1])
                record: ExonRecord = ExonRecord(
                    proj, exon_num, chrom, start, end, strand
                )
                self.query_projections[proj].exons.append(record)

    def read_ref_isoforms(self, species: str) -> None:
        """Reads the transcript-to-gene mapping for a given reference"""
        file: Union[str, None] = self.ref_data[species].ref_isoforms
        if file is None:
            self._to_log(
                "No reference isoform file provided for reference %s" % species
            )
            return
        if not os.path.exists(file):
            self._to_log("File %s does not exist; skipping" % file)
            return
        with open(file, "r") as h:
            for line in h:
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                # gene, tr = data
                gene, tr = map(lambda x: f"{species}.{x}", data)
                self.ref_proj2gene[tr] = gene

    def get_overlapping_genes(self, species: str) -> None:
        """Extracts data on non-synonymous genes overlapping in the reference"""
        file: Union[str, None] = self.ref_data[species].ref_bed
        if file is None:
            self._to_log(
                "No reference annotation BED file provided for reference %s" % species
            )
            return
        isoforms: Union[str, None] = self.ref_data[species].ref_isoforms
        if isoforms is None:
            self._to_log(
                (
                    "No reference isoform file provided for reference %s; "
                    "ignoring the reference annotation BED file"
                )
                % species
            )
            return
        if not os.path.exists(file):
            self._to_log("File %s does not exist; skipping" % file)
            return
        chrom2tr_coords: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        for data in read_tab(file):
            # tr: str = data[3]
            tr: str = f"{species}.{data[3]}"
            chrom: str = data[0]
            start: int = int(data[6])
            end: int = int(data[7])
            chrom2tr_coords[chrom].append((tr, start, end))
        for chrom, trs in chrom2tr_coords.items():
            trs = sorted(trs, key=lambda x: (x[1], x[2]))
            for i, tr1 in enumerate(trs):
                gene1: str = self.ref_proj2gene.get(tr1[0], "")
                if not gene1:
                    continue
                for tr2 in trs[i + 1 :]:
                    if tr2[1] >= tr1[2]:
                        break
                    gene2: str = self.ref_proj2gene.get(tr2[0], "")
                    if not gene2:
                        continue
                    if gene1 == gene2:
                        continue
                    inter: int = intersection(
                        *tr1[1:], *tr2[1:]
                    )  ## TODO: Redundant unless we impose overlap threshold

                    if inter <= 0:
                        continue
                    self.intersecting_ref_genes[gene1].add(gene2)
                    self.intersecting_ref_genes[gene2].add(gene1)

    def read_paralogs(
        self, species: str
    ) -> None:  ## TODO: Can be replaced with postfixes
        """DEPRECATED: Extracts paralogous projections' names"""
        file: Union[str, None] = self.ref_data[species].paralogs
        if file is None:
            self._to_log(
                "No paralogous projection list provided for reference %s" % species
            )
            return
        if not os.path.exists(file):
            self._to_log("File %s does not exist; skipping" % file, "warning")
            return
        self._to_log("Reading paralogous projection list for reference %s" % species)
        with open(file, "r") as h:
            for line in h:
                line = line.strip()
                self.paralog_pool.add(line)

    def read_ppgenes(
        self, species: str
    ) -> None:  ## TODO: Can be replaced with postfixes
        """DEPRECATED: Extracts processed pseudogene/retrogene projections' names"""
        file: Union[str, None] = self.ref_data[species].ppgenes
        if file is None:
            self._to_log(
                "No retrogene projection list provided for reference %s" % species
            )
            return
        if not os.path.exists(file):
            self._to_log("File %s does not exist; skipping" % file, "warning")
            return
        self._to_log("Reading retrogene projection list for reference %s" % species)
        with open(file, "r") as h:
            for line in h:
                line = line.strip()
                self.ppgene_pool.add(line)

    def infer_genes(self) -> None:
        """
        Infers query genes from combined reference annotations.
        The logic is effectively borrowed from infer_query_genes.py :
        any two projections overlapping by at least one coding base
        on the same strand are collapsed into a single gene,
        with the following exceptions:
        1) paralogs overlapping orthologs are not included
        in the gene composition and are further discarded;
        2) retrogenes overlapping either orthologs or paralogs
        are not included in the gene composition and are further discarded
        """
        self._to_log(
            "Inferring query genes from the combined annotation; might take a while"
        )
        ## create the graph structure
        ## iterate over chromosomes
        for chrom, projections in self.query_annotation.items():
            ## sort the projections
            projections.sort(
                key=lambda x: (
                    self.query_projections[x].start,
                    self.query_projections[x].end,
                )
            )
            ## start iterating over projections
            for i, name_out in enumerate(projections):
                proj_out: BedRecord = self.query_projections[name_out]
                out_start: int = proj_out.start
                out_end: int = proj_out.end
                out_tr: str = get_proj2trans(name_out)[0]
                # out_paralog: bool = name_out in self.paralog_pool
                # out_ppgene: bool = name_out in self.ppgene_pool
                out_paralog: bool = base_proj_name(name_out) in self.paralog_pool
                out_valid_paralog: bool = (
                    out_paralog and proj_out.loss_status in INTACT_PROJECTIONS
                )
                out_ppgene: bool = base_proj_name(name_out) in self.ppgene_pool
                out_ortholog: bool = not (out_paralog or out_ppgene)
                discarded: bool = False
                edges: List[str] = []
                for j, name_in in enumerate(projections[i + 1 :], start=i):
                    proj_in: BedRecord = self.query_projections[name_in]
                    ## if an item has been already discarded, skip it
                    if name_in in self.discarded_items:
                        continue
                    in_tr: str = get_proj2trans(name_in)[0]
                    in_start: int = proj_in.start
                    in_end: int = proj_in.end
                    ## as long as transcripts are sorted properly,
                    ## it's safe to break the inner loop here
                    if in_start >= out_end:
                        break
                    ## if proj_in lies upstream, it must have been already processed
                    if in_end <= out_start:
                        continue
                    ## ignore projections encoded by different strands
                    if proj_in.strand != proj_out.strand:
                        continue
                    ## ignore projections encoded by different strands
                    # in_paralog: bool = name_in in self.paralog_pool
                    # in_ppgene: bool = name_in in self.ppgene_pool
                    in_paralog: bool = base_proj_name(name_in) in self.paralog_pool
                    in_valid_paralog: bool = (
                        in_paralog and proj_in.loss_status in INTACT_PROJECTIONS
                    )
                    in_ppgene: bool = base_proj_name(name_in) in self.ppgene_pool
                    in_ortholog: bool = not (in_paralog or in_ppgene)
                    has_intersection: bool = False
                    for ex_out in proj_out.exons:
                        for ex_in in proj_in.exons:
                            if (
                                intersection(
                                    ex_out.start, ex_out.end, ex_in.start, ex_in.end
                                )
                                > 0
                            ):
                                has_intersection = True
                                break
                        if has_intersection:
                            # if out_paralog or in_paralog:
                            break
                    if has_intersection:
                        ## check for intersections in the reference
                        out_gene: str = self.ref_proj2gene.get(out_tr)
                        in_gene: str = self.ref_proj2gene.get(in_tr)
                        if out_gene and in_gene:
                            if out_gene in self.intersecting_ref_genes.get(
                                in_gene, []
                            ) or in_gene in self.intersecting_ref_genes.get(
                                out_gene, []
                            ):
                                continue
                        ## ortholog + paralog/pp: discard the non-orthologous prediction
                        if out_ortholog and not in_ortholog and not in_valid_paralog:
                            self.discarded_items.add(name_in)
                            continue
                        if (
                            not out_ortholog and not out_valid_paralog and in_ortholog
                        ):  # not out_ortholog and in_ortholog:
                            self.discarded_items.add(name_out)
                            discarded = True
                            break
                        ## ppgene + paralog: discard the ppgene
                        if out_paralog and in_ppgene:
                            self.discarded_items.add(name_in)
                            continue
                        if out_ppgene and in_paralog:
                            self.discarded_items.add(name_out)
                            discarded = True
                            break
                        ## same orthology class: add the edge
                        edges.append(name_in)
                if discarded:
                    continue
                ## add the edges to the graph
                if name_out not in self.graph.nodes:
                    self.graph.add_node(name_out)
                for name_in in edges:
                    self.graph.add_edge(name_out, name_in)

    def pick_best_isoforms(self) -> None:
        """
        For each query gene (=overlap graph clique), pick the representatives:
        *   If at least one isoform corresponds to the user-defined loss statuses,
            ignore all other isoforms;
            *   If two or more isoforms are completely identical in terms of coordinates,
                pick the one with the better loss status;
            *   If two or more isoforms have identical layout and loss status, pick the one
                corresponding to a higher priority reference/run;
            *   Otherwise, pick the first encountered representative for this layout-status-priority combination
        *   If no projections with the requested status have been encountered, pick the one
        """
        # _all_names: Set[str] = set()
        gene_name2counter: Dict[str, int] = defaultdict(int)
        if NX_VERSION < 2.4:
            components = list(nx.connected_component_subgraphs(self.graph))
        else:
            components = [
                self.graph.subgraph(c) for c in nx.connected_components(self.graph)
            ]
        with (
            open(self.gene_tsv, "w") as gt,
            open(self.gene_bed, "w") as gb,
            open(self.projection_bed, "w") as qb,
        ):
            gt.write(GENE_TSV_HEADER)
            for component in components:
                # v = any("XM_047428684.1#NBPF6#165" in x for x in component)
                ## initialize a temporary storage for initial candidates
                selected: Dict[str, str] = {}
                name2lines_selected: Dict[str, List[str]] = defaultdict(list)
                ## define the best class
                best_status: int = 0
                ## set a semaphore for whether the user-defined loss classes
                ## have been encountered in the clique
                allowed_class_found: bool = False
                ## keep track of the paralogs present in this clique
                paralogs: Set[str] = set()
                lost: Set[str] = set()
                all_paralogs: bool = all(
                    base_proj_name(x) in self.paralog_pool for x in component
                )
                for name in component:
                    is_paralog: bool = base_proj_name(name) in self.paralog_pool
                    ## paralogs are to be handled later
                    if is_paralog and not all_paralogs:
                        paralogs.add(name)
                        continue
                    proj: BedRecord = self.query_projections[name]
                    status: int = CLASS_TO_NUM[proj.loss_status]
                    if not is_paralog and status == CLASS_TO_NUM["L"]:
                        lost.add(name)
                    allowed_status: bool = proj.loss_status in self.accepted_statuses
                    prev_is_better: bool = False
                    if allowed_status:  ## let this projection in
                        if not allowed_class_found:
                            ## this is the first representative of the allowed loss classes;
                            ## wipe out the previous instances
                            selected.clear()
                            name2lines_selected.clear()
                        allowed_class_found = True
                    else:
                        ## do not let the items worse than the current best
                        if status < best_status:
                            prev_is_better = True
                            ## skip this item and proceed further
                            continue
                        elif status > best_status:
                            ## not the allowed class but already better
                            ## than what has been already encountered; clear the selected lists
                            selected.clear()
                            name2lines_selected.clear()
                            # continue
                    ## at this point, this is a likely candidate
                    ## however, chances are an item in exactly the same coordinates
                    ## has been already found
                    if any(x in selected for x in proj.cds_lines.values()):
                        for line in proj.cds_lines.values():
                            if line not in selected:
                                continue
                            prev_name: str = selected[line]
                            prev_proj: BedRecord = self.query_projections[prev_name]
                            prev_status: int = CLASS_TO_NUM[prev_proj.loss_status]
                            ## pick the one with the best loss status
                            if prev_status > status:
                                prev_is_better = True
                                break
                            ## if it is a tie, go for the more preferred reference
                            elif prev_status == status:
                                species: str = self.query_proj2ref[name]
                                priority: int = self.ref_data[species].priority
                                prev_species: str = self.query_proj2ref[prev_name]
                                prev_priority: int = self.ref_data[
                                    prev_species
                                ].priority
                                ## if the previous prediction's priority is higher (lower) or equal, keep it
                                if prev_priority < priority:
                                    prev_is_better = True
                                    break
                                ## literally no reason to pick over another, so go by the alphabet order
                                elif prev_priority == priority:
                                    if prev_name < name:
                                        prev_is_better = True
                                        break
                                # else:
                            for prev_line in name2lines_selected[prev_name]:
                                del selected[prev_line]
                        if prev_is_better:
                            continue
                    ## the existing item is better; proceed further
                    if prev_is_better:
                        continue
                    ## otherwise, the new prediction is the winner
                    best_status = max(best_status, status)
                    for line in proj.cds_lines.values():
                        selected[line] = name
                    name2lines_selected[name] = list(proj.cds_lines.values())
                ## now, process the paralogs
                ## first, retrieve all the exon records
                valid_exons: Dict[str, List[ExonRecord]] = defaultdict(list)
                for proj_name in selected.values():
                    proj: BedRecord = self.query_projections[proj_name]
                    for exon in proj.exons:
                        valid_exons[exon.chrom].append(exon)
                ## sort them chromwise
                for chrom in valid_exons:
                    valid_exons[chrom].sort(key=lambda x: (x.start, x.end))
                paralogs = sorted(
                    paralogs, key=lambda x: self.query_projections[x].coords()
                )
                added_paralogs: Set[str] = set()
                for paralog in paralogs:
                    proj: BedRecord = self.query_projections[paralog]
                    status: int = CLASS_TO_NUM[proj.loss_status]
                    species: str = self.query_proj2ref[paralog]
                    priority: int = self.ref_data[species].priority
                    ## ignore if exactly the same item has been already encountered
                    ## (by default, paralogs are expected to have one 'line' alone)
                    prev_is_better: bool = False
                    for line in proj.cds_lines.values():
                        if line not in selected:
                            continue
                        prev_name: str = selected[line]
                        ## if the previous candidate is not a paralog, drop the current paralog
                        if base_proj_name(prev_name) not in self.paralog_pool:
                            prev_is_better = True
                            break
                        ## loss status filter - same as for orthologs
                        prev_proj: BedRecord = self.query_projections[prev_name]
                        prev_status: int = CLASS_TO_NUM[prev_proj.loss_status]
                        if prev_status > status:
                            prev_is_better = True
                            break
                        elif prev_status == status:
                            ## species filter - same as for orthologs
                            prev_species: str = self.query_proj2ref[prev_name]
                            prev_priority: int = self.ref_data[prev_species].priority
                            if prev_priority < priority:
                                prev_is_better = True
                                break
                            elif prev_priority == priority:
                                ## alphabet priority - same as for orthologs
                                if prev_name < paralog:
                                    prev_is_better = True
                                    break
                    if prev_is_better:
                        continue
                    already_present: bool = all(
                        x in selected for x in proj.cds_lines.values()
                    )
                    if already_present:
                        for line in proj.cds_lines.values():
                            selected[line] = paralog
                        name2lines_selected[paralog] = list(proj.cds_lines.values())
                        continue
                    if not selected:
                        for line in proj.cds_lines.values():
                            selected[line] = paralog
                        name2lines_selected[paralog] = list(proj.cds_lines.values())
                        continue
                    ## a projection with all segments already encountered,
                    ## it will definitely not add any new exons
                    to_add: bool = False
                    for paralog_exon in proj.exons:
                        ## record the minimal overlap
                        min_abs: int = paralog_exon.length()
                        for valid_exon in valid_exons[proj.chrom]:
                            if paralog_exon.end < valid_exon.start:
                                break
                            if paralog_exon.start > valid_exon.end:
                                continue
                            ## find the intersection size
                            inter_size: int = intersection(
                                paralog_exon.start,
                                paralog_exon.end,
                                valid_exon.start,
                                valid_exon.end,
                            )
                            inter_size = max(inter_size, 0)
                            min_abs = min(min_abs, paralog_exon.length() - inter_size)
                        ## if at least one exon meets the requirements, add it to the output
                        min_rel: float = safe_div(min_abs, paralog_exon.length())
                        if (
                            min_abs >= self.paralog_abs_novelty_threshold
                            and min_rel >= self.paralog_rel_novelty_threshold
                        ):
                            to_add = True
                            break
                    if to_add:
                        for line in proj.cds_lines.values():
                            selected[line] = paralog
                        name2lines_selected[paralog] = list(proj.cds_lines.values())
                ## last round: process the losses
                if allowed_class_found:
                    ## first, add the exons coming from the newly included paralogs
                    for added_paralog in added_paralogs:
                        proj: BedRecord = self.query_projections[added_paralog]
                        for exon in proj.exons:
                            valid_exons[exon.chrom].append(exon)
                    for chrom in valid_exons:
                        valid_exons[chrom].sort(key=lambda x: (x.start, x.end))
                    ## now, the autism round starts
                    for lost_proj in lost:
                        proj: BedRecord = self.query_projections[lost_proj]
                        status: int = CLASS_TO_NUM[proj.loss_status]
                        species: str = self.query_proj2ref[lost_proj]
                        priority: int = self.ref_data[species].priority
                        if all(x in selected for x in proj.cds_lines.values()):
                            continue
                        prev_is_better: bool = False
                        for line in proj.cds_lines.values():
                            if line not in selected:
                                continue
                            prev_name: str = selected[line]
                            prev_proj: BedRecord = self.query_projections[prev_name]
                            prev_status: int = CLASS_TO_NUM[prev_proj.loss_status]
                            if prev_status > status:
                                prev_is_better = True
                                break
                            elif prev_status == status:
                                prev_species: str = self.query_proj2ref[prev_name]
                                prev_priority: int = self.ref_data[
                                    prev_species
                                ].priority
                                if prev_priority < priority:
                                    prev_is_better = True
                                    break
                                elif prev_priority == priority:
                                    if prev_name < lost_proj:
                                        prev_is_better = True
                                        break
                        if prev_is_better:
                            continue
                        to_add: bool = False
                        for lost_exon in proj.exons:
                            ## record the minimal overlap
                            min_abs: int = lost_exon.length()
                            if min_abs == 0:
                                continue
                            for valid_exon in valid_exons[proj.chrom]:
                                if lost_exon.end < valid_exon.start:
                                    break
                                if lost_exon.start > valid_exon.end:
                                    continue
                                ## find the intersection size
                                inter_size: int = intersection(
                                    lost_exon.start,
                                    lost_exon.end,
                                    valid_exon.start,
                                    valid_exon.end,
                                )
                                inter_size = max(inter_size, 0)
                                min_abs = min(min_abs, lost_exon.length() - inter_size)
                            ## if at least one exon meets the requirements, add it to the output
                            min_rel: float = safe_div(min_abs, lost_exon.length())
                            if (
                                min_abs >= self.lost_abs_novelty_threshold
                                and min_rel >= self.lost_rel_novelty_threshold
                            ):
                                to_add = True
                                break
                        if to_add:
                            for line in proj.cds_lines.values():
                                selected[line] = lost_proj
                            name2lines_selected[lost_proj] = list(
                                proj.cds_lines.values()
                            )
                ## all the projections have been processed; name the gene and define its coordinates
                filtered_component: nx.Graph = component.copy()
                nodes_to_remove: Set[str] = {
                    x for x in filtered_component.nodes() if x not in selected.values()
                }
                filtered_component.remove_nodes_from(nodes_to_remove)
                if NX_VERSION < 2.4:
                    minor_components = list(
                        nx.connected_component_subgraphs(filtered_component)
                    )
                else:
                    minor_components = [
                        filtered_component.subgraph(c)
                        for c in nx.connected_components(filtered_component)
                    ]
                for minor_component in minor_components:
                    # all_projs: List[BedRecord] = [
                    #     self.query_projections[x] for x in selected.values()
                    # ]
                    # self.num_genes += 1
                    all_projs: List[BedRecord] = [
                        self.query_projections[x] for x in minor_component.nodes()
                    ]
                    ## define the coordinates; that's the easy part
                    chrom: str = all_projs[0].chrom
                    strand: str = "+" if all_projs[0].strand else "-"
                    start: int = min([x.start for x in all_projs])
                    end: int = max([x.end for x in all_projs])
                    ## define the name; a slightly trickier part
                    ## first, define how many gene participate in the locus annotation
                    genes: Set[str] = {
                        self.ref_proj2gene.get(get_proj2trans(x.name)[0], x.name)
                        for x in all_projs
                    }
                    ## define the main name; for simplicity, do not bother with chain numbers
                    genes = sorted(genes)
                    if len(genes) == 1:  ## single gene; assign its name to query locus
                        gene: str = genes.pop()
                    elif (
                        len(genes) < 4
                    ):  ## up to three genes; combine their names separated by comma
                        gene: str = ",".join(genes)
                    else:  ## more than three genes; gene+
                        gene: str = genes.pop() + "+"
                    ## define prefix
                    if all(
                        base_proj_name(x.name) in self.paralog_pool for x in all_projs
                    ):
                        # if any(x.name in self.paralog_pool for x in all_projs):
                        prefix: str = "paralog_"
                        best_loss: str = max(
                            [x.loss_status for x in all_projs], key=lambda x: CLASS_TO_NUM[x]
                        )
                        if best_loss in self.accepted_statuses:
                            self.num_para += 1
                    elif all(
                        base_proj_name(x.name) in self.ppgene_pool for x in all_projs
                    ):
                        # elif any(x.name in self.ppgene_pool for x in all_projs):
                        prefix: str = "retro_"
                        self.num_retro += 1
                    elif not allowed_class_found:
                        if best_status > CLASS_TO_NUM["M"]:
                            prefix: str = "lost_"
                            self.num_lost += 1
                        else:
                            prefix: str = "missing_"
                            self.num_missing += 1
                    else:
                        prefix: str = ""
                        self.num_genes += 1
                        # print(f"Gene: {gene}")
                        # if gene in _all_names:
                        #     print(f"DUPLICATE: {gene=}, {minor_component=}, {minor_components=}")
                        # _all_names.add(gene)
                    gene_name2counter[gene] += 1
                    if gene_name2counter[gene] > 1:
                        gene += f"_{gene_name2counter[gene]}"
                    name: str = prefix + gene
                    all_species: Set[str] = {self.query_proj2ref[x] for x in minor_component.nodes()}
                    self.gene2species[name] = all_species
                    for _spec in all_species:
                        self.species2num_item[_spec] += 1
                        if not prefix:
                            self.species2num_genes[_spec] += 1
                    ## done! now, write the output
                    ## for gene BED, it's a single line
                    gene_bed: str = "\t".join(
                        map(str, [chrom, start, end, name, 0, strand])
                    )
                    gb.write(gene_bed + "\n")
                    ## then, iterate over projection
                    _num_proj: int = 0
                    for proj in all_projs:
                        basename: str = base_proj_name(proj.name)
                        self.final_projections.add(basename)
                        ## get the reference name and add it as a prefix
                        species: str = self.query_proj2ref[proj.name]
                        ## gene isoform file; write the gene and projections names
                        gt.write(name + "\t" + proj.name + "\n")
                        ## projection BED file; write the original BED line
                        # for proj_bed in proj.return_bed_line(prefix=species):
                        for proj_bed in proj.return_bed_line():
                            qb.write(proj_bed + "\n")
                        _num_proj += 1
                    if _num_proj == 0:
                        self._to_log("No projections for gene %s" % name, "warning")

    def prepare_gtf(self) -> None:
        """
        Converts the BED annotation file into the GTF format (output/query_annotation.gtf).
        The results are compressed into gzip format by default.

        Args:
            None

        Returns:
            None
        """
        from .isoforms_for_gtf import ProvisionalIsoformMapper

        ProvisionalIsoformMapper(
            [
                self.gene_tsv,
                self.projection_bed,
                self.query_genes_for_gtf,
                "-ln",
                "integrate",
            ],
            standalone_mode=False,
        )
        gtf_cmd: str = (
            f"{GTF_SCRIPT} -b {self.projection_bed} -i {self.query_genes_for_gtf} "
            f"-o {self.query_gtf} && gzip -5 -f {self.query_gtf}"
        )
        _ = self._exec(
            gtf_cmd, "Converting integrated annotation into GTF format failed:"
        )
        self._rm(self.query_genes_for_gtf)

    def generate_support_matrix(self) -> None:
        """
        Generates a reference support matrix for all the inferred genes
        
        Args:
            None

        Returns:
            None
        """
        sorted_refs: List[str] = sorted(self.ref_data.keys(), key=lambda x: self.ref_data[x].priority)
        header: str = "gene\t" + "\t".join(sorted_refs)
        with open(self.support_matrix_file, "w") as h:
            h.write(header + "\n")
            for gene, species in self.gene2species.items():
                line_stub: str = f"{gene}\t"
                ref_presence: List[str] = []
                for ref in sorted_refs:
                    ref_presence.append("+" if ref in species else "-")
                line_stub += "\t".join(ref_presence) + "\n"
                h.write(line_stub)

    def generate_summary(self) -> None:
        """
        Generates a short summary file (output/summary.txt)

        Args:
            None

        Returns:
            None
        """
        ref_list: List[str] = sorted(
            self.ref_data.keys(), key=lambda x: self.ref_data[x].priority
        )
        # refs: str = "* " + "\n* ".join(ref_list)
        refs: str = ""
        for ref in ref_list:
            ref_genes: int = self.species2num_genes[ref]
            gene_perc: int = to_perc(ref_genes, self.num_genes)
            refs += REF_SUPPORT_LINE.format(
                ref=ref, 
                ref_genes=ref_genes, 
                all_genes=self.num_genes, 
                gene_perc=gene_perc,
            )
        with open(self.summary, "w") as h:
            summary_txt: str = SUMMARY.format(
                br=BREAK_LINE,
                query=self.query_name,
                num_genes=self.num_genes,
                num_retro=self.num_retro,
                num_lost=self.num_lost,
                num_missing=self.num_missing,
                num_para=self.num_para,
                refs=refs,
            )
            h.write(summary_txt)

    def prepare_ucsc_file(self) -> None:
        """
        Parses the UCSC BigBed file, extracting the lines for integrated annotation BigBed
        as well as sequences for final sequence files
        """
        if self.skip_ucsc:
            self._to_log("Skipping the joint UCSC BigBed file creation as suggested")
            return

        from .get_names_from_bed import BedNameRetriever

        nuc_seqs: List[str] = []
        prot_seqs: List[str] = []
        longest_word: int = 0
        bb_input_found: bool = False
        with open(self.bigbed_stub, "w") as h:
            for ref, data in self.ref_data.items():
                seqs_found: Set[str] = set()
                file: str = data.ucsc_bigbed
                if file is None:
                    self._to_log(
                        "No UCSC BigBed file provided for reference %s; skipping" % ref,
                        "warning",
                    )
                    continue
                if not os.path.exists(file):
                    self._die("UCSC BigBed file %s does not exist" % file)
                bb_input_found = True
                self.has_ucsc_data.append(ref)
                read_cmd: str = f"{self.bigbedtobed_binary} {file} stdout"
                output: TextIO = self._exec(
                    read_cmd, "Attempt to read BigBed file %s failed: " % file
                )
                for line in output.split("\n"):
                    # line = line.decode('utf8')
                    data: List[str] = line.strip().split("\t")
                    if not data or not data[0]:
                        continue
                    # name: str = data[3]
                    name: str = f"{ref}.{data[3]}"
                    basename: str = base_proj_name(name)
                    if basename not in self.final_projections:
                        continue
                    # name already has ref prefix from line above
                    data[3] = name
                    line = "\t".join(data)
                    h.write(line + "\n")
                    longest_word = max(longest_word, len(name))
                    if basename in seqs_found:
                        continue
                    nuc_seq: str = data[34]
                    nuc_seqs.append(f">{name}\n{nuc_seq}")
                    prot_seq: str = data[35]
                    prot_seqs.append(f">{name}\n{prot_seq}")
                    seqs_found.add(basename)
        if not bb_input_found:
            self._to_log(
                (
                    "No BigBed files were provided for any of the annotation; "
                    "skipping the genome browser track creation step"
                ),
                "warning",
            )
            return
        ## sort the resulting file, convert it into BigBed format
        bb_cmd: str = (
            f"sort -k1,1 -k2,2n -o {self.bigbed_stub} {self.bigbed_stub} && "
            f"{self.bedtobigbed_binary} -type=bed12+26 {self.bigbed_stub} "
            f"{self.chrom_sizes} {self.bigbed} "
            f"-tab -extraIndex=name -as={self.schema}"
        )
        _ = self._exec(bb_cmd, "BigBed generation failed: ")
        self._to_log("BigBed file successfully created")
        # bed_ix_cmd: str = (
        #     f"{MAKE_IX_SCRIPT} {self.bigbed_stub} | sort -u > {self.bed_index}"
        # )
        # _ = self._exec(bed_ix_cmd, "BED file indexing failed")
        BedNameRetriever(
            (
                "--input",
                self.bigbed_stub,
                "--output",
                self.bed_index,
                "--log_name",
                "integrate",
            ),
            standalone_mode=False,
        )
        bigbed_ix_cmd: str = (
            f"{self.ixixx_binary} {self.bed_index} {self.ix} {self.ixx} "
            f"-maxWordLength={longest_word}"
        )
        _ = self._exec(bigbed_ix_cmd, "Index creation for UCSC BigBed file failed: ")
        self._to_log("BigBed index successfully created; removing the temporary files")
        self._rm(self.bigbed_stub)
        self._rm(self.bed_index)

        self._to_log("Preparing the decoration track")
        self.prepare_decorator_track()
        self._to_log(
            "Decorator track successfully created; removing the temporary files"
        )

        if prot_seqs:
            with gzip.open(self.protein_file, "wb") as h:
                for prot in prot_seqs:
                    h.write((prot + "\n").encode("utf8"))
        if nuc_seqs:
            with gzip.open(self.nucleotide_file, "wb") as h:
                for nuc in nuc_seqs:
                    h.write((nuc + "\n").encode("utf8"))

    def prepare_decorator_track(self) -> None:
        """Prepares a companion decorator track for genome browser"""
        longest_word: int = 0
        bb_input_found: bool = False
        with open(self.decorator_stub, "w") as h:
            for ref, data in self.ref_data.items():
                file: str = data.decorator_bigbed
                if file is None:
                    self._to_log(
                        "No decorator BigBed file provided for reference %s; skipping"
                        % ref,
                        "warning",
                    )
                    continue
                if not os.path.exists(file):
                    self._die("Decorator BigBed file %s does not exist" % file)
                bb_input_found = True
                read_cmd: str = f"{self.bigbedtobed_binary} {file} stdout"
                output: TextIO = self._exec(
                    read_cmd, "Attempt to read BigBed file %s failed: " % file
                )
                for data in read_tab(output.split("\n")):
                    name_field: str = data[12]
                    name: str = ":".join(name_field.split(":")[2:])
                    basename: str = base_proj_name(name)
                    if basename not in self.final_projections:
                        continue
                    name = f"{ref}.{name}"
                    # data[3] = name
                    data[12] = ":".join([*name_field.split(":")[:2], name])
                    line = "\t".join(data)
                    h.write(line + "\n")
                    longest_word = max(longest_word, len(name))
        if not bb_input_found:
            self._to_log(
                (
                    "No decotator BigBed files were provided for any of the annotation; "
                    "skipping the decorator track creation step"
                ),
                "warning",
            )
            return
        decor_cmd: str = (
            f"sort -o {self.decorator_stub} -k1,1 -k2,2n {self.decorator_stub} && "
            f"{self.bedtobigbed_binary} -type=bed12+ -as={self.decorator_schema} "  ## ADD!!
            f"{self.decorator_stub} {self.chrom_sizes} {self.decorator}"
        )
        _ = self._exec(decor_cmd, "BigBed generation failed: ")
        self._rm(self.decorator_stub)
        # _ = self._exec(decor_cmd, "Decoration track production failed:")
        # decor_bb_cmd: str = (
        #     f"{self.bedtobigbed_binary} -type=bed12+ -as={self.DECOR_SCHEMA_FILE} "
        #     f"{self.decor_stub} {self.query_contig_size_file} {self.decoration_track}"
        # )

    def prepare_sequence_files(self) -> None:
        """p"""
        seq_files_recorded: bool = bool(self.has_ucsc_data)
        for ref, data in self.ref_data.items():
            if ref in self.has_ucsc_data:
                continue
            prot_file: str = data.protein_file
            if prot_file is None or not os.path.exists(prot_file):
                self._die(
                    (
                        "Reference %s has neither valid UCSC BigBed file "
                        "nor valid protein sequence file"
                    )
                    % ref
                )
            nuc_file: str = data.nucleotide_file
            if nuc_file is None or not os.path.exists(nuc_file):
                self._die(
                    (
                        "Reference %s has neither valid UCSC BigBed file "
                        "nor valid nucleotide sequence file"
                    )
                    % ref
                )
            prot_gzipped: bool = prot_file.endswith(".gz")
            with (
                gzip.open(prot_file, "rb") if prot_gzipped else open(prot_file, "r")
            ) as h:
                prot_seqs: List[str] = self._read_fasta(h, ref)
            nuc_gzipped: bool = nuc_file.endswith(".gz")
            with gzip.open(nuc_file, "rb") if nuc_gzipped else open(nuc_file, "r") as h:
                nuc_seqs: List[str] = self._read_fasta(h, ref)
            write_mode: str = "a" if seq_files_recorded else "w"
            with open(self.protein_file, write_mode) as h:
                for entry in prot_seqs:
                    h.write(entry + "\n")
            with open(self.nucleotide_file, write_mode) as h:
                for entry in nuc_seqs:
                    h.write(entry + "\n")
            seq_files_recorded = True

    def _read_fasta(self, handle: TextIO, ref: str) -> List[str]:
        """
        Reads the sequence file, return the list of FASTA entries corresponding to integrated annotation
        """
        header: str = ""
        seq: str = ""
        sequences: List[str] = []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entry: str = f">{ref}.{header}\n{seq}"
                    sequences.append(entry)
                header = ""
                seq = ""
                proj: str = base_proj_name(line.split()[0]).lstrip(">")
                if proj not in self.final_projections:
                    continue
                header = proj
            elif header:
                seq += line
        if seq:
            entry: str = f">{ref}.{header}\n{seq}"
            sequences.append(entry)
        return sequences

    def set_logging(self) -> None:
        """Sets logging and disables logging propagation"""
        super().set_logging(name="integrate", toga_module="integrate")
        self.logger.propagate = False
