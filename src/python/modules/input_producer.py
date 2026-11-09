"""
Reference input annotation preparation functionality
"""

__author__ = "Yury V. Malovichko"
__credits__ = ("Michael Hiller", "Bernhard Bein")
__year__ = "2025"

import logging
import os
from collections import Counter, defaultdict
from collections.abc import Callable, Sequence
from contextlib import nullcontext
from dataclasses import dataclass
from enum import Enum, IntEnum
from functools import partial
from shutil import which
from types import ModuleType
from typing import Final

import click

from .constants import Headers, RejectionReasons
from .gxf import Gxf
from .shared import (
    CommandLineManager,
    cluster_suffix,
    consistent_name,
    dir_name_by_date,
    get_upper_dir,
    hex_dir_name,
    intersection,
    read_tab,
    replace_illegal_characters,
)

logging.root.handlers = []

DEFAULT_PREFIX: Final[str] = "TOGA2_ref_annotation"
LEVEL: Final[str] = "level"
GENE: Final[str] = "GENE"
TRANSCRIPT: Final[str] = "TRANSCRIPT"
ILLEGAL_NAME: Final[str] = "ILLEGAL_NAME"
REJECTED_CONTIG: Final[str] = "REJECTED_CONTIG"
NON_CODING: Final[str] = "NON_CODING"
OUT_OF_FRAME: Final[str] = "OUT_OF_FRAME"
NUMERIC_FIELDS: Final[tuple[int]] = (1, 2, 6, 7, 9, 10, 11)
U12: Final[str] = "U12"
U2: Final[str] = "U2"
CANON_SITES: Final[tuple[str, ...]] = ("GT-AG", "GC-AG")
U12_CANON_SITES: Final[str] = "GT-AG"
MIN_INTRON_LENGTH_FOR_CLASSIFICATION: Final[int] = 30
MIN_INTRON_LENGTH_FOR_PROFILES: Final[int] = 70
ENTRY_START: Final[str] = ">"
ACC: Final[str] = "acceptor"
DONOR: Final[str] = "donor"
ACC_PROFILE_LEN: Final[int] = 22
DONOR_PROFILE_LEN: Final[int] = 6
ACC_START_POS: Final[int] = -22
DONOR_START_POS: Final[int] = 1
NUCS: Final[tuple[str, ...]] = ("A", "T", "C", "G")
N: Final[str] = "N"
CANON: Final[str] = "canon"
NONCANON: Final[str] = "nonCanon"
SEP_DUMMY: Final[str] = "n"
NAME: Final[str] = "name"
TYPE_ID: Final[str] = "type_id"
DEFAULT_CLASS_COL: Final[int] = 12

GENE_BED_LINE: Final[str] = "{chrom}\t{start}\t{end}\t{name}\t1000\t{strand}\n"
BED_LINE: Final[str] = (
    "{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t{cds_start}\t{cds_end}\t0,0,0,\t{exon_num}\t{sizes}\t{starts}\n"
)
ISOFORM_LINE: Final[str] = "{gene}\t{tr}\n"

TOGA2_ROOT: Final[str] = get_upper_dir(__file__, 4)
DEFAULT_TWOBITTOFA: Final[str] = os.path.join(TOGA2_ROOT, "bin", "twoBitToFa")
DEFAULT_BED2FRACTION: Final[str] = os.path.join(
    TOGA2_ROOT, "src", "rust", "target", "release", "bed12ToFraction"
)
PROFILE_DIR: Final[str] = "CESAR2.0_profiles"
EQUI_ACC: Final[str] = "equiprobable_acceptor.tsv"
EQUI_DONOR: Final[str] = "equiprobable_donor.tsv"
DEFAULT_MEMORY_LIMIT: Final[int] = 24
EXTRACTION_ERR_MSG: Final[str] = "ERROR: twoBitToFa call failed"
CONVERSION_ERR_MSG: Final[str] = "ERROR: FASTA to 2bit conversion failed"
BED12TO6_ERR: Final[str] = "BED12 to BED6 conversion failed:"

TRANSCRIPTS: Final[str] = "toga.transcripts.bed"
ISOFORMS: Final[str] = "toga.isoforms.tsv"
U12_FILE: Final[str] = "toga.U12introns.bed"
GENE_BED: Final[str] = "toga.genes.bed"
IMMUNE_BED: Final[str] = "toga.immunoglobulins.bed"
IMMUNE_ISOFORMS: Final[str] = "toga.immunoglobulin_isoforms.tsv"
IMMUNE_GENE_BED: Final[str] = "toga.immunoglobulin_genes.bed"
SLEASY: Final[str] = "sleasy.exons.2bit"
REJ_LOG: Final[str] = "rejected_items.tsv"
RENAME_MAP: Final[str] = "rename_map.tsv"
FISSION_MAP: Final[str] = "fission_map.tsv"
FUSION_MAP: Final[str] = "fusion_map.tsv"
EXON_BED6: Final[str] = "all_exons.bed6"
FA_FOR_SLEASY: Final[str] = "all_exons.fasta"
PROV_BED: Final[str] = "transcripts.from_gxf.raw.bed"
PROV_ISOFORM: Final[str] = "isoforms.from_gxf.raw.bed"

ATTR2BIN: Final[dict[str, str]] = {
    "twobittofa_binary": "twoBitToFa",
    "fatotwobit_binary": "faToTwoBit",
    "intronic_binary": "intronIC",
}
DEFAULT_TWOBIT2FA: Final[str] = os.path.join(TOGA2_ROOT, "bin", "twoBitToFa")
DEFAULT_FA2TWOBIT: Final[str] = os.path.join(TOGA2_ROOT, "bin", "faToTwoBit")
DEFAULT_INTRONIC: Final[str] = os.path.join(
    TOGA2_ROOT, "bin", "intronIC", "intronIC", "intronIC.py"
)
BIN2DEFAULT: Final[dict[str, str]] = {
    "twoBitToFa": DEFAULT_TWOBIT2FA,
    "faToTwoBit": DEFAULT_FA2TWOBIT,
    "intronIC": DEFAULT_INTRONIC,
}

NEWTAB: Final[str] = "\n\t"


## Enums
class BridgeType(Enum):
    SELF  = "Self-isoform bridge"
    FUS_ISO = "Fusion isoform"
    FUSION = "Fusion gene"
    COMBINED = "Fusion + fission"


class IsNMD(IntEnum):
    """Levels of nonsense-mediated decay (NMD) signatures detected"""
    NO = 0
    WEAK = 1
    STRONG = 2


class RejEnum(IntEnum):
    """Rejection reason enumeration"""
    NON_CODING = 0
    FRAME = 1
    CONTIG = 2
    NAME = 3
    FISSION_GENE = 4
    FISSION_TR = 5
    FUSION_GENE = 6
    FUSION_TR = 7
    ORPHAN = 8
    SHORT_EXON = 9
    SHORT_INTRON = 10
    NMD = 11
    LENGTH_FILTER = 12
    GENE = 13


## Classes
@dataclass(frozen=True)
class FissionReport:
    gene: str
    units: list[tuple[str, int]]
    bridges: list[str]


@dataclass(frozen=True)
class FusionReport:
    node: tuple[str, int]
    bridge_type: BridgeType
    fusion_map: dict[str, list[list[str]]]
    split_genes: list[str]


@dataclass(frozen=True)
class FaultyModelReport:
    fissions: list[FissionReport]
    fusions: list[FusionReport]


class NmdConstants:
    STOP_JUNCTION_DIST: int = 55
    WEAK_STOP_JUNCTION_DIST: int = 80
    START_STOP_DIST: int = 100
    BIG_EXON_JUNCTION_DIST: int = 400


@dataclass(frozen=True)
class Transcript:
    chrom: str
    start: int
    end: int
    strand: bool
    exons: list[tuple[int, int]]
    line: str

    def cds_len(self) -> int:
        """
        Reports coding sequence length as the sum of all coding exons' lengths

        Arguments:
            None

        Returns:
            int, coding sequence length
        """
        total: int = 0
        for e in self.exons:
            e_start, e_end = e
            if e_end <= self.start:
                continue
            if e_start >= self.end:
                break
            e_start = max(self.start, e_start)
            e_end = min(self.end, e_end)
            total += (e_end - e_start)
        return total

    def cds_span(self) -> int:
        """
        Returns coding sequence span as difference between CDS start 
        and end coordinates

        Arguments:
            None

        Returns:
            int, coding sequence span
        """
        return self.end - self.start


## Functions
def add_prefix(template: str, prefix: str | None = "") -> str:
    """Adds a prefix to the output file name

    Arguments:
        prefix: A prefix to add to the file name. Can be empty.
        template: A file name template.

    Returns:
        A filename in {prefix}.{template} if prefix is non-empty, {template} otherwise.
    """
    if not prefix:
        return template
    return f"{prefix}.{template}"


def cds_exon_overlap(a: Transcript, b: Transcript, min_inter: int = 1) -> bool:
    """
    Computes coding sequence overlap. 

    Arguments:
        a, b: Transcript, transcript to assess the overlap between
        min_inter: int, minimal overlap size between two exons, in bp

    Returns:
        bool, indication whether the sufficient overlap has been found
    """
    len_a: int = len(a.exons)
    len_b: int = len(b.exons)
    i = j = 0
    while i < len_a and j < len_b:
        s1, e1 = a.exons[i]
        s2, e2 = b.exons[j]
        if s1 > e2:
            j += 1
            continue
        if e1 < s2:
            i += 1
            continue
        if intersection(s1, e1, s2, e2) >= min_inter:
            return True
        if e1 < e2:
            i += 1
        else:
            j += 1
    return False


class InputProducer(CommandLineManager):
    """
    Main class for input annotation preparation
    """

    __slots__ = (
        "all_intron_bed", 
        'annot', 
        'annot_format', 
        'bed2fraction_binary', 
        'bed6_exons', 
        'contigs',
        'disable_cesar_profiles', 
        'disable_intron_classification',
        "disable_length_filter",
        "disable_models",
        "do_not_correct_names",
        'excluded_contigs', 
        'fasta_for_sleasy', 
        'fatotwobit_binary', 
        'filtered_annotation', 
        'filtered_isoforms',
        "fission_map_file",
        "fusion_map_file",
        "fusions",
        "gene2fissions",
        "gene2trs",
        'gene_bed',
        'ic_cores', 
        'immune', 
        'immune_annotation', 
        'immune_gene_bed', 
        'immune_isoforms', 
        'intron2class', 
        'intron2coords', 
        'intron_file', 
        'intronic', 
        'intronic_binary', 
        'intronic_cores', 
        'isoforms', 
        'keep_tmp', 
        'min_int_exon_length', 
        'min_intron_length', 
        'min_intron_length_cesar', 
        'min_intron_length_intronic', 
        'nmd_level',
        "no_gene_suffix",
        'output', 
        'profile_dir', 
        'profiles', 
        'provisional_bed', 
        'provisional_isoforms', 
        'rejection_log',
        "rejection_reasons",
        "rel_thresh",
        "remove_fissions",
        "remove_fusions",
        "rename_map",
        "rename_map_file",
        'sleasy_2bit', 
        'strip_versions', 
        'tmp_dir', 
        "tr2annot", 
        "tr2gene",
        'twobit', 
        'twobittofa_binary', 
        'v',
    )

    @staticmethod
    def nmd_filter(
        thin_start: int,
        thin_end: int,
        cds_start: int,
        cds_end: int,
        sizes: list[int],
        starts: list[int],
        strand: bool,
    ) -> IsNMD:
        """
        Filter nonsense-mediated decay (NMD) targets from the final annotation.

        The logic goes as follows:
            1. The code tracks the distance between the end of the coding sequence
            and the last (first if the strand is negative) UTR-UTR exon junction.
            In practical terms this means all exons past the CDS end up to N-1 exon
            included (or all exons before the CDS start starting from the second
            in case of the negative strand).

            2. Cumulative UTR length in this interval is recorded and compared to
            the threshold value (55 bases by default, according to the "55 bp rule");

            3. If 3'-UTR region contains at least one exon junction (=intron), and
            the distance between the end of the coding sequence and the last UTR
            exon junction is more than 55 bp, the exon is considered a potential NMD target;

            4. The following conditions might hinder NMD, so fulfilling either of them
            downgrades the evidence for NMD to the "weak evidence" category:
                * CDS-UTR junction distance is  <=80 nucleotides long (relaxed "55 bp rule");
                * Stop-containing exon consists of both CDS and UTR, and the UTR
                sequence is at least 400 bases long;
                * Coding sequence is <=100 bases long.

        Arguments:
            * thin_start: int - sequence start in the genome, including UTRs;
            * thin_end: int - sequence end in the genome, including UTRs;
            * cds_start: int - coding sequence start in the genome;
            * cds_end: int - coding sequence end in the genome;
            * sizes: List[int] - list of exon sizes (BED12 field 11);
            * starts: List[int] - list of exon relative start coordinates (BED12 field 12);
            * strand: bool - True if strand is positive.

        Returns:
            * enum IsNMD:
                * 0 (NO) if transcript is not an NMD target;
                * 1 (WEAK) if evidence for NMD is considered tentative;
                * 2 (STRONG) if there is enough evidence for NMD targeting for the transcript.
        """
        junction_counter: int = 0
        cds_length: int = 0
        utr_length: int = 0
        distance_to_utr_junction: int = 0
        exon_num: int = len(sizes)
        for (
            i,
            start,
        ) in enumerate(starts):
            exon_start: int = thin_start + start
            exon_end: int = exon_start + sizes[i]
            ## add the coding portion's length to the total CDS length
            exon_cds_start: int = max(cds_start, exon_start)
            exon_cds_end: int = min(cds_end, exon_end)
            length: int = exon_cds_end - exon_cds_start
            cds_length += max(length, 0)
            ## ignore coding/g'-UTR exons
            if (strand and exon_end < cds_end) or (
                not strand and exon_start > cds_start
            ):
                continue
            if not distance_to_utr_junction:
                if strand and exon_start < cds_end < exon_end:
                    distance_to_utr_junction = exon_end - exon_cds_end
                elif not strand and exon_start < cds_start < exon_end:
                    distance_to_utr_junction = exon_cds_start - exon_start
            if not ((strand and i == exon_num - 1) or (not strand and i == 0)):
                if strand:
                    utr_length += exon_end - max(exon_start, cds_end)
                else:
                    utr_length += min(exon_end, cds_start) - exon_start
                junction_counter += 1
        status: int = IsNMD.NO
        if junction_counter > 0 and utr_length > NmdConstants.STOP_JUNCTION_DIST:
            ## "big exon condition": UTR part of the CDS:UTR exon is bigger than the threshold value
            big_exon: bool = (
                distance_to_utr_junction >= NmdConstants.BIG_EXON_JUNCTION_DIST
            )
            ## start codon proximity: stop codon is within first 100 coding bases
            short_cds: bool = cds_length <= NmdConstants.START_STOP_DIST
            ## distance to the last UTR junction is <= 80 bp
            short_distance: int = utr_length <= NmdConstants.WEAK_STOP_JUNCTION_DIST
            if big_exon or short_cds or short_distance:
                status = IsNMD.WEAK
            else:
                status = IsNMD.STRONG
        return status

    @staticmethod
    def strip_v(name: str) -> str:
        """
        Removes the version identifier from the item (chromosome/contig/scaffold/transcript/gene) name.
        If a name contains no dot symbols, returns the input string. If a name starts with leading dots,
        strips them beforehand.

        Arguments:
            name: item name

        Returns:
            str item name without the version identifier

        Example:
            base: str = InputProducer.strip_v("XM_035385144.1")
            assert(base == "XM_035385144")
        """
        comps: list[str] = name.lstrip(".").split(".")
        if len(comps) == 1:
            return name
        if not comps[-1].isdigit():
            return name
        return ".".join(comps[:-1])

    @staticmethod
    def _sides(node: tuple[str, int], adj: Sequence[str]) -> list[set[tuple[str, int]]]:
        """
        Returns subgraphs of the focal node's connected component after the node removal. 

        Arguments:
            node: tuple, name of the focal node (e.g., articulation node)
            adj: dictionary or NetworkX Graph object

        Returns:
            list of connected components disjointed by removal of the focal node, 
            represented by sets of node labels

        Note:
            Mimics the following NetworkX code:
            ```
            comp = nx.node_connected_component(adj, node)
            return list(nx.connected_components(adj.subgraph(comp - {node})))
            ```
        """
        seen, sides = {node}, []
        for n in adj[node]:
            if n in seen: 
                continue
            seen.add(n)
            comp, stack = [n], [n]
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y not in seen: 
                        seen.add(y)
                        comp.append(y)
                        stack.append(y)
            sides.append(comp)
        return sides

    @staticmethod
    def _classify_bridge(
        node: tuple[str, int], 
        sides: list[set[tuple[str, int]]], 
        members: dict[tuple[str, int], list[str]],
        cross_unit: dict[str, set[tuple[str, int]]],
    ) -> FusionReport:
        """
        Assesses whether the focal bridge transcript represents a fusion event

        Bridging transcript is a transcript overlapping two or more transcripts that do not 
        mutually overlap, e.g.:

        ====tr1[***]====[***]============[***]====[***]tr1===
        ====tr2[***]====[***]tr2======tr3[***]====[***]tr3===

        In the overlap graph, bridging transcripts (or multi-transcript units) form articulation nodes; 
        their removal from the graph creates a set of disjointed connected components.


        Arguments:
            node: str or tuple, an articulation node of an overlap graph (string or unit-like tuple)
            side: list, two or more connected components disjointed by the ``node`` removal
            members: dict, {unit -> [unit_members]} mapping between the gene unit 
                and transcripts contained within it
            cross_unit: {transcript -> {units}} mapping between transcripts and units 
                of other genes they overlap

        Returns:
            FusionReport instance with the following fields:
                * ``node`` - focal node name
                * ``bridge_type`` - one of four BridgeType values
                * ``fusion_map`` - dictionary mapping individual fusion isoforms 
                    to units they bridge
                * ``split`` - list of genes split between the bridged loci 
                    (potential gene model fission)
        """
        ## expand "sides" into {transcript: subgraph} mapping
        tr2side: dict[str, int] = {u: i for i, s in enumerate(sides) for u in s}
        ## detect transcript-level fusions: 
        ## transcripts overlapping other gene units in more than one subgraph
        fusion_trs: list[str] = sorted(
            t for t in members[node] if len({tr2side[u] for u in cross_unit.get(t, ())}) > 1
        )
        ## get the units these transcripts bridge
        fusion_map: dict[str, list[list[str]]] = {
            x: [members[u] for u in sorted(cross_unit[x])] for x in fusion_trs
        }
        ## whole-node bridge condition: all isoforms overlap at least one other gene
        has_independent: bool = any(t not in cross_unit for t in members[node])
        per_side: dict[int, int] = Counter(g for s in sides for g in {u[0] for u in s})
        ## get genes occurring in more than one subgraph -> fissions coupled with fusion
        split: list[str] = sorted(g for g, n in per_side.items() if n > 1)
        if not has_independent:
            bridge_type: BridgeType = BridgeType.COMBINED if split else BridgeType.FUSION
        else:
            bridge_type: BridgeType = BridgeType.FUS_ISO if fusion_trs else BridgeType.SELF
        return FusionReport(node, bridge_type, fusion_map, split)

    @classmethod
    def _detect_faulty_models(
        cls, 
        members: dict[tuple[str, int], list[str]],
        cross_unit: dict[str, set[tuple[str, int]]],
        articulation_nodes: set[tuple[str, int]],
        sides_fn: Callable[[tuple[str, int]], list[set[tuple[str, int]]]],
    ) -> FaultyModelReport:
        """
        Arguments:
            //

        Returns:
            FaultyModeReport with the following fields:
                * ``fissions``: list of FissionReport instances
                * ``fusions``: list of FusionReport instances
        """
        fusions: list[FusionReport] = [
            cls._classify_bridge(c, sides_fn(c), members, cross_unit) for c in articulation_nodes
        ]
        ## expand {unit -> transcript} adjacency dict to {gene -> transcript}
        units_of: dict[str, list[str]] = defaultdict(list)
        for u in members: 
            units_of[u[0]].append(u)
        bridged: dict[str, tuple[str, int]] = defaultdict(list)
        for f in fusions:
            for g in f.split_genes: 
                bridged[g].append(f.node)
        ##
        fissions: list[FissionReport] = [
            FissionReport(g, sorted(us), sorted(bridged[g]))
            for g, us in units_of.items() if len(us) > 1
        ]
        return FaultyModelReport(
            fissions=fissions, 
            fusions=sorted(fusions, key=lambda f: f.node),
        )

    def __init__(
        self,
        ref_2bit: click.Path,
        ref_annot: click.Path,
        annot_format: str = Gxf.BED,
        ref_isoforms: click.Path = None,
        output: click.Path | None = None,
        prefix: str = "",
        do_not_add_gene_suffix: bool = False,
        disable_transcript_filtering: bool = False,
        contigs: str | None = None,
        excluded_contigs: str | None = None,
        do_not_correct_names: bool = False,
        min_internal_exon_length: int = 0,
        min_intron_length: int = 0,
        nmd_filter_level: bool = 0,
        disable_model_analysis: bool = False,
        remove_fusions: bool = False,
        remove_fissions: bool = False,
        disable_length_filter: bool = False, ## TODO: Add flag to toga2.py
        relative_length_threshold: bool = 0.8, ## TODO: Add flag to toga2.py
        disable_intron_classification: bool = False,
        disable_cesar_profiles: bool = False,
        intronic_binary: click.Path | None = None,
        intronic_cores: int = 1,
        min_intron_length_intronic: int = MIN_INTRON_LENGTH_FOR_CLASSIFICATION,
        twobittofa_binary: click.Path | None = None,
        fatotwobit_binary: click.Path | None = None,
        min_intron_length_cesar: int = MIN_INTRON_LENGTH_FOR_PROFILES,
        keep_temporary: bool = False,
        strip_versions: bool = False,
    ) -> None:
        self.v: bool = True
        self.set_logging()

        if annot_format != Gxf.BED and ref_isoforms is not None:
            self._die(
                f"Ref annotation format was set to {annot_format}, which is incompatible "
                "with the external isoform mapping. Please remove the isoform file "
                "from the settings or provide a BED annotation file instead. "
            )

        self.twobit: click.Path = ref_2bit
        self.annot: click.Path = ref_annot
        self.isoforms: click.Path | None = ref_isoforms
        self.contigs: str | None = contigs
        self.excluded_contigs: str | None = excluded_contigs
        self.do_not_correct_names: bool = do_not_correct_names
        self.min_int_exon_length: int = min_internal_exon_length
        self.min_intron_length: int = min_intron_length
        self.nmd_level: bool = nmd_filter_level
        self.disable_length_filter: bool = disable_length_filter
        self.rel_thresh: float = relative_length_threshold
        self.disable_models: bool = disable_model_analysis
        self.remove_fusions: bool = remove_fusions
        self.remove_fissions: bool = remove_fissions
        self.annot_format: str = annot_format
        self.disable_intron_classification: bool = disable_intron_classification
        self.disable_cesar_profiles: bool = disable_cesar_profiles
        self.min_intron_length_intronic: int = min_intron_length_intronic
        self.min_intron_length_cesar: int = min_intron_length_cesar
        self.output: str = (
            output if output is not None else hex_dir_name(DEFAULT_PREFIX)
        )
        self.strip_versions: bool = strip_versions
        self.no_gene_suffix: bool = do_not_add_gene_suffix
        self.tmp_dir: str = os.path.join(self.output, dir_name_by_date("tmp"))
        self.keep_tmp: bool = keep_temporary

        self.filtered_annotation: os.PathLike = os.path.join(
            self.output, add_prefix(TRANSCRIPTS, prefix)
        )
        self.filtered_isoforms: os.PathLike = os.path.join(
            self.output, add_prefix(ISOFORMS, prefix)
        )
        self.gene_bed: os.PathLike = os.path.join(
            self.output, add_prefix(GENE_BED, prefix)
        )
        self.sleasy_2bit: os.PathLike = os.path.join(
            self.output, add_prefix(SLEASY, prefix)
        )
        self.immune_annotation: str = os.path.join(
            self.output, add_prefix(IMMUNE_BED, prefix)
        )
        self.immune_isoforms: str = os.path.join(
            self.output, add_prefix(IMMUNE_ISOFORMS, prefix)
        )
        self.immune_gene_bed: str = os.path.join(
            self.output, add_prefix(IMMUNE_GENE_BED, prefix)
        )
        self.rejection_log: str = os.path.join(self.output, REJ_LOG)
        self.rename_map_file: str = os.path.join(self.output, RENAME_MAP)
        self.fission_map_file: str = os.path.join(self.output, FISSION_MAP)
        self.fusion_map_file: str = os.path.join(self.output, FUSION_MAP)

        self.bed6_exons: os.PathLike = os.path.join(self.tmp_dir, EXON_BED6)
        self.fasta_for_sleasy: os.PathLike = os.path.join(self.tmp_dir, FA_FOR_SLEASY)

        self.intron_file: os.PathLike = os.path.join(
            self.output, add_prefix(U12_FILE, prefix)
        )
        self.all_intron_bed: os.PathLike = os.path.join(self.tmp_dir, "all_introns.bed")
        self.provisional_bed: str = os.path.join(self.tmp_dir, PROV_BED)
        self.provisional_isoforms: str = os.path.join(self.tmp_dir, PROV_ISOFORM)

        self.twobittofa_binary: os.PathLike | None = twobittofa_binary
        self.fatotwobit_binary: os.PathLike | None = fatotwobit_binary
        self.intronic_binary: os.PathLike | None = intronic_binary
        self.bed2fraction_binary: str = DEFAULT_BED2FRACTION
        self.intronic_cores: int = intronic_cores

        self.tr2annot: dict[str, Transcript] = {}
        self.gene2trs: dict[str, list[str]] = {}
        self.tr2gene: dict[str, str] = {}
        self.gene2fissions: dict[str, tuple[str]] = {}
        self.fusions: dict[str, list[list[str]]] = {}
        self.rejection_reasons: dict[str, dict[str, list[RejEnum]]] = {
            GENE: defaultdict(list), TRANSCRIPT: defaultdict(list)
        }
        self.immune: set[str] = set()
        self.rename_map: dict[str, str] = {}

        self.run()

    def run(self) -> None:
        """Entry point"""
        ## checking all the necessary binaries
        self._to_log("Checking the necessary binaries")
        ## create output directory
        self._to_log("Checking binaries")
        self.check_binaries()
        self._to_log("Creating output directory")
        self._mkdir(self.output)
        self._mkdir(self.tmp_dir)
        ## step 0: if output is GTF/GFF3, extract the data from the file
        if self.annot_format != Gxf.BED:
            self._to_log("Extracting the annotation data from the GTF/GFF3 file")
            self.from_gxf()
        ## step 1: annotation file check
        self._to_log("Refining the reference annotation BED file")
        self.check_annotation()
        ## step 2, optional: isoform file check, potential further annotation filtering
        if self.isoforms is not None:
            self._to_log("Refining the input isoform file")
            self.read_isoforms()
            if not self.disable_models:
                self._to_log("Checking the input annotation for gene model fission/fusion")
                self.report_corrupted_models()
            if not self.disable_length_filter:
                self.isoform_length_filter()
        ## write the results for steps 1 and 2
        self._to_log("Writing the annotation results")
        self.write_annotation()
        self.write_isoforms()
        self.write_rejection_log()
        self.write_rename_map()
        self.write_fission_map()
        self.write_fusion_map()
        ## step 3: exon 2bit file
        self._to_log("Creating 2bit exon sequence storage for SLEASY")
        self.create_sleasy_input()
        ## step 4: intron classification, U12 input file preparation
        if not self.disable_intron_classification:
            if not self.disable_cesar_profiles:
                self.intron2class: dict[str, tuple[str, str]] = {}
                self.intron2coords: dict[str, tuple[str, int, int, bool]] = {}
            self._to_log("Classifying reference introns")
            self.intron_classifier()

        ## step 5: CESAR2 profile generation; will not shoot if step 3 is disabled
        if not (self.disable_intron_classification or self.disable_cesar_profiles):
            self._to_log("Generating CESAR profiles")
            self.profiles: dict[tuple[str, bool, str], dict[int, dict[str, int]]] = (
                defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
            )
            self.profile_dir: str = os.path.join(self.output, PROFILE_DIR)
            self._mkdir(self.profile_dir)
            self.generate_cesar_profiles()

    def check_binaries(self) -> None:
        """Checks binary availability and execution permissions for all third-party programs.

        Since the input is provided from the click Python CLI, the existence
        of non-empty option value is basically guaranteed at this point;
        nevertheless, the code further checks whether the provided binary
        instance is executable.

        If no value is provided, the method searches for binary availability in $PATH.
        Once it is found, the execution permissions are further ensured.

        The method throws error if no (executable)  instance is found. The only exception is
        intronIC absence if `--disable_intron_classification` flag was set
        """
        for attr, exp_name in ATTR2BIN.items():
            binary: str | None = getattr(self, attr)
            if binary is not None:
                self._to_log(f"Testing {exp_name} binary at {binary}")
                if os.access(binary, os.X_OK):
                    self._to_log(
                        f"The provided binary is executable; using the stated {exp_name} instance"
                    )
                    setattr(self, attr, binary)
                    continue
                else:
                    self._to_log(
                        (
                            f"{exp_name} binary at {binary} does not seem executable; "
                            "looking for alternatives"
                        ),
                        "warning",
                    )
            else:
                self._to_log(
                    f"No {exp_name} executable was provided; looking for alternatives"
                )
            ## check for the default version in bin/
            default_binary: str = BIN2DEFAULT[exp_name]
            if os.path.exists(default_binary):
                self._to_log(
                    f"Found {exp_name} instance at {default_binary}; checking the execution permissions"
                )
                if os.access(default_binary, os.X_OK):
                    self._to_log(
                        f"The found binary is executable; using the TOGA2-supplied {default_binary} instance"
                    )
                    setattr(self, attr, default_binary)
                    continue
                self._to_log(
                    f"TOGA2-supplied {exp_name} at {default_binary} is not executable; "
                    "looking for alternatives in $PATH"
                )
            else:
                self._to_log(
                    f"{exp_name} is missing at {default_binary}; looking for alternatives in $PATH"
                )
            binary_in_path: str | None = which(exp_name)
            if binary_in_path is not None:
                self._to_log(
                    f"Found {exp_name} instance at {binary_in_path}; checking the execution permissions"
                )
                if os.access(binary_in_path, os.X_OK):
                    self._to_log(
                        f"The found binary is executable; using the found {binary_in_path} instance"
                    )
                    setattr(self, attr, binary_in_path)
                    continue
                self._die(
                    f"The {exp_name} binary found in $PATH at {binary_in_path} is not executable; "
                    f"check your $PATH or provide a valid {exp_name} instance"
                )
            if exp_name == "intronIC" and self.disable_intron_classification:
                self._to_log("No available intronIC instance was found; skipping")
                continue
            self._die(
                f"No {exp_name} binary found in $PATH; "
                f"check your $PATH or provide a valid {exp_name} instance"
            )

    def from_gxf(self) -> None:
        """
        Prepares annotation from a GTF/GFF3 file

        Arguments:
            None

        Returns:
            None
        """
        from .gxf import GxfTranscript, parse_gxf_attrs

        tr_objects: dict[str, GxfTranscript] = {}
        for i, data in enumerate(read_tab(self.annot), start=1):
            if data[0].startswith("#"):
                continue
            if len(data) < 9:
                self._to_log(
                    f"Less than nine fields encountered in the GXF file at line {i}",
                    "warning",
                )
                continue
            chrom: str = data[0]
            level: str = data[2]
            strand: bool = data[6] == "+"
            attrs: dict[str, str] = parse_gxf_attrs(
                data[8], gff3=(self.annot_format == Gxf.GFF3)
            )
            start: int = int(data[3]) - 1  ## GXF is 1-based, BED is 0-based
            end: int = int(data[4])
            if level in Gxf.CHILD_TYPES:
                parent: str = Gxf.get_exon_parent(attrs)
                ## an orphan exon encountered; continue
                if not parent or parent not in tr_objects:
                    continue
                if level == Gxf.CDS:
                    ## add the coding exon to the transcript 
                    ## and update coding sequence boundaries
                    tr_objects[parent].add_cds(start, end)
                elif level == Gxf.STOP_CODON:
                    ## add stop codon to the sequence
                    tr_objects[parent].add_stop(start, end)
                else:
                    ## add the exon to the sequence
                    tr_objects[parent].add_exon(start, end)
                continue
            if level not in Gxf.CODING_TAGS and level not in Gxf.VDJ_TAGS:
                continue
            biotype: str | None = Gxf.get_biotype(attrs)
            tr_id: str = (
                attrs.get(Gxf.TR_ID, "")
                if self.annot_format == Gxf.GTF
                else attrs.get(Gxf.ID, "")
            )
            tr_id = Gxf.remove_tr_prefices(tr_id)
            if not tr_id:
                continue
            ## check if the item is a regular gene or an immune gene fragment
            is_coding = biotype not in Gxf.VDJ_TAGS and level not in Gxf.VDJ_TAGS
            gene: str = Gxf.get_transcript_parent(attrs)
            tr_obj: GxfTranscript = GxfTranscript(
                chrom,
                start,
                end,
                strand,
                tr_id,
                gene,
                is_coding,
            )
            tr_objects[tr_id] = tr_obj
        has_isoforms: bool = False
        with (
            open(self.provisional_bed, "w") as bed,
            open(self.provisional_isoforms, "w") as iso,
        ):
            for tr_obj in tr_objects.values():
                if not tr_obj.exons:
                    continue
                thin_start: int = tr_obj.exons[0][0]
                thin_end: int = tr_obj.exons[-1][1]
                name: str = tr_obj.name
                if not tr_obj.is_coding:
                    cds_start: int = (
                        thin_start if tr_obj.cds_start is None else tr_obj.cds_start
                    )
                    cds_end: int = (
                        thin_start if tr_obj.cds_end is None else tr_obj.cds_end
                    )
                    self.immune.add(name)
                else:
                    if not tr_obj.has_cds:
                        continue
                    cds_start: int = tr_obj.cds_start
                    cds_end: int = tr_obj.cds_end
                sizes: list[str] = []
                starts: list[str] = []
                count: int = 0
                for e_start, e_end in tr_obj.exons:
                    start: str = str(e_start - thin_start)
                    size: str = str(e_end - e_start)
                    starts.append(start)
                    sizes.append(size)
                    count += 1
                sizes: str = ",".join(sizes) + ","
                starts: str = ",".join(starts) + ","
                strand: str = "+" if tr_obj.strand else "-"
                if self.strip_versions:
                    name = self.strip_v(tr_obj.name)
                    # tr_obj.parent = self.strip_v(tr_obj.parent)
                bed_line: str = BED_LINE.format(
                    chrom=tr_obj.chrom,
                    start=thin_start,
                    end=thin_end,
                    name=name,
                    strand=strand,
                    cds_start=cds_start,
                    cds_end=cds_end,
                    exon_num=count,
                    sizes=sizes,
                    starts=starts,
                )
                bed.write(bed_line)
                if tr_obj.parent is not None:
                    iso_line: str = ISOFORM_LINE.format(gene=tr_obj.parent, tr=name)
                    iso.write(iso_line)
                    has_isoforms = True
            if not has_isoforms:
                self._to_log("No transcript parents found in the input GxF file", "warning")
            else:
                self.isoforms = self.provisional_isoforms

    def check_annotation(self) -> None:
        """
        Checks the transcripts in the BED file, filtering out offending items.

        Filtering criteria applied are divided into 'hard' and 'soft' filters:

        1. Violating 'hard' filtering criteria renders transcripts mechanically 
        incompatible with TOGA2 annotation procedure; such transcrtipts are removed 
        from the annotation with no option to salvage them. 'Hard' criteria include: 
        * absence of coding sequence;
        * shifted frame in the coding sequence;
        * location on a deprecated chromosome (contig, scaffold, etc.);
        * use of deprecated symbols in transcript name if ``--do_not_correct_names`` 
        flag was set
        2. In contrast, transcripts violating `soft` filtering criteria can be used 
        as TOGA2 source material; however, their biological plausability is questionable. 
        All `soft` criteria can be modified or inactivated with command line parameters. 
        Such criteria include:
        * transcripts with short internal (non-terminal) coding exons; 
        these are often introduced by gene annotation software to compensate for 
        base-level assembly error or faulty gene models;
        * short coding sequence introns; their origin is similar to that of 
        short internal exons;
        * potentional nonsense-mediated decay (NMD) targets
        In addition, soft-filtered transcripts are retained in the final annotation if 
        their progenitor gene has no intact transcripts.

        Unless the ``--do_not_correct_names`` flag is set, transcripts with illegal characters 
        in their names get renamed via invoking ``replace_illegal_characters()``. Original names 
        are used internally up until the point of writing the results to the output file(s).

        Arguments:
            None

        Returns:
            None

        Raises:
            Error-exits if any non-empty line:
            * contains less than 12 fields 
            * or has non-numeric values in coordinate fields (field 2, 3, 7, 8, 10, 11, 12)
            * or has start coordinate greater than end coordinate
        """
        illegal_name: set[str] = set()
        rejected_contigs: set[str] = set()
        non_coding: set[str] = set()
        short_exons: set[str] = set()
        short_introns: set[str] = set()
        out_of_frame: set[str] = set()
        nmd_targets: set[str] = set()
        annot_input: str = (
            self.provisional_bed if self.annot_format != Gxf.BED else self.annot
        )
        for i, data in enumerate(read_tab(annot_input), start=1):
            if len(data) != 12:
                self._die(
                    f"Improper formatting at reference annotation file line {i}; "
                    f"expected 12 fields, got {len(data)}"
                )
            if any(x.strip() == "" for x in data):
                self._die(
                    f"Improper formatting at reference annotation file line {i}; "
                    "empty fields encountered"
                )
            for field in NUMERIC_FIELDS:
                if not data[field].replace(",", "").isdigit():
                    self._die(
                        f"Improper formatting at reference annotation file line {i}; "
                        f"field {field} contains non-numeric data"
                    )
            name: str = data[3]
            if self.strip_versions and self.annot_format == Gxf.BED:
                name = self.strip_v(name)
            is_coding: bool = name not in self.immune
            ## remove the transcripts with improperly formatted names
            if not consistent_name(name):
                illegal_name.add(name)
                if self.do_not_correct_names:
                    continue
                else:
                    self.rename_map[name] = replace_illegal_characters(name)
            chrom: str = data[0]
            if self.strip_versions:
                chrom = self.strip_v(chrom)
                data[0] = chrom
            line = "\t".join(data)
            ## if entries were restricted to specific contigs,
            ## apply the respective filters
            if self.contigs and chrom not in self.contigs:
                rejected_contigs.add(name)
                continue
            if self.excluded_contigs and chrom in self.excluded_contigs:
                rejected_contigs.add(name)
                continue
            ## check coding sequence presence and frame intactness
            thin_start: int = int(data[1])
            thin_end: int = int(data[2])
            strand: bool = data[5] == "+"
            cds_start: int = int(data[6])
            cds_end: int = int(data[7])
            if cds_end < cds_start:
                self._die(
                    f"Improper formatting at reference annotation file line {i}; "
                    "coding sequence start coordinate greater than the start coordinate"
                )
            if cds_start == cds_end:
                non_coding.add(name)
                continue
            ## iterate over exon entries to infer the CDS length
            frame_length: int = 0
            sizes: list[int] = [int(x) for x in data[10].split(",") if x]
            starts: list[int] = [int(x) for x in data[11].split(",") if x]
            hard_filtered: bool = False
            # soft_filtered: bool = False
            if self.nmd_level != IsNMD.NO and is_coding:
                nmd_status: int = self.nmd_filter(
                    thin_start,
                    thin_end,
                    cds_start,
                    cds_end,
                    sizes,
                    starts,
                    strand,
                )
                if nmd_status == IsNMD.STRONG or (
                    nmd_status == IsNMD.WEAK and self.nmd_level == IsNMD.STRONG
                ):
                    nmd_targets.add(name)
                    # soft_filtered = True
            has_short_exons: bool = False
            prev_end: int | None = None
            has_short_introns: bool = False
            exons: list[tuple[int, int]] = []
            for ex, (start, size) in enumerate(zip(starts, sizes)):
                start += thin_start
                end: int = start + size
                if start < cds_start:
                    if end > cds_start:
                        size -= cds_start - start
                    else:
                        continue
                if end > cds_end:
                    if start < cds_end:
                        size -= end - cds_end
                    else:
                        continue
                if size < self.min_int_exon_length and not (
                    start <= cds_start or end >= cds_end
                ):
                    has_short_exons = True
                frame_length += size
                exons.append((max(start, cds_start), min(end, cds_end)))
                if prev_end is not None:
                    intron_size: int = start - prev_end
                    if intron_size < self.min_intron_length:
                        has_short_introns = True
                prev_end = end
            if frame_length % 3:
                out_of_frame.add(name)
                hard_filtered = True
            if has_short_exons and is_coding:
                short_exons.add(name)
                # soft_filtered = True
            if has_short_introns and is_coding:
                short_introns.add(name)
                # soft_filtered = True
            if hard_filtered:
                continue
            self.tr2annot[name] = Transcript(
                chrom, cds_start, cds_end, strand, exons, line 
            )

        ## report the discarded items
        if non_coding:
            self._log_item_list(
                    "The following transcripts were filtered out "
                    "due to absence of coding sequence:",
                non_coding,
            )
            for i in non_coding:
                self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.NON_CODING)
        if out_of_frame:
            self._log_item_list(
                "The following transcripts were filtered out "
                "due to shifted reading frame:",
                out_of_frame,
            )
            for i in out_of_frame:
                self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.FRAME)
        if rejected_contigs:
            self._log_item_list(
                "The following transcripts were filtered out "
                "due to their location in deprecated contigs:",
                rejected_contigs
            )
            for i in rejected_contigs:
                self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.CONTIG)
        if illegal_name:
            if self.do_not_correct_names:
                self._log_item_list(
                    "The following transcripts were filtered out "
                    "due to illegal symbols used in their names:",
                    illegal_name
                )
                for i in illegal_name:
                    self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.NAME)
            else:
                self._log_item_list(
                    "Names for the following transcripts were corrected "
                    "due to illegal symbols found:",
                    illegal_name,
                )
        if short_exons:
            self._log_item_list(
                f"The following transcripts contain short (<{self.min_int_exon_length} bp) "
                "internal exons and might be filtered out:",
                short_exons,
            )
            for i in short_exons:
                self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.SHORT_EXON)
        if short_introns:
            self._log_item_list(
                f"The following transcripts contain short (<{self.min_intron_length} bp) "
                "CDS introns and might be filtered out:",
                short_introns,
            )
            for i in short_introns:
                self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.SHORT_INTRON)
        if nmd_targets:
            self._log_item_list(
                "The following transcripts are potential nonsense-mediated "
                "decay targets and might be filtered out:",
                nmd_targets,
            )
            for i in nmd_targets:
                self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.NMD)
        ## proceed further

    def read_isoforms(self) -> None:
        """
        Reads isoform file, populating mutual mapping between genes and transcripts

        Arguments:
            None

        Returns:
            None

        Raises:
            Error-exists if:
            * any of the non-empty lines in the isoform file contains less or more than 2 columns;
            * any isoform appears more than once in the isoform file

        Note:
            Populates ``self.gene2trs`` ({gene -> [transcripts]}) and ``self.tr2gene`` 
            ({transcript -> gene}) mapping.
            Launches ``self.check_isoforms()`` at the end of the method
        """
        iso_file: str = (
            self.isoforms if self.annot_format == Gxf.BED else self.provisional_isoforms
        )
        for i, data in enumerate(read_tab(iso_file), start=1):
            if len(data) != 2:
                self._die(
                    f"Improper formatting at isoforms file line {i}; "
                    f"expecting 2 columns, got {len(data)}"
                )
            gene, tr = data
            if self.strip_versions and self.annot_format == Gxf.BED:
                ## stripping genes from the dot identifier is discouraged
                ## since some RefSeq genes (sadly!) contain dots in gene symbols
                # gene, tr = map(self.strip_v, data)
                tr = self.strip_v(tr)
            if gene not in self.gene2trs:
                self.gene2trs[gene] = []
            if tr in self.tr2gene:
                self._die(
                    f"Isoform {tr} appears at least twice in the isoform file; "
                    f"first duplicate found at line {i}"
                )
            self.gene2trs[gene].append(tr)
            self.tr2gene[tr] = gene
        self.check_isoforms()

    def check_isoforms(self) -> None:
        """
        Filters the reference isoform (gene-to-transcript) mapping file
        by removing genes for which transcripts (isoforms) were discarded
        at any other step.

        If a gene ends up with no valid transcripts, the method tries to salvage 
        "soft-filtered" isoforms (isoforms with short internal coding exons, 
        isoforms with short internal coding introns, and nonsense-mediated decay targets) ,
        removing them from 

        At the same time, all the transcripts recorded at the annotation filter step
        and missing gene mapping in the isoforms file are further excluded
        from the final annotation.

        Arguments:
            None

        Returns:
            None
        """
        rejected_genes: list[str] = []
        illegal_name: list[str] = []
        for gene, trs in list(self.gene2trs.items()):
            if not consistent_name(gene) and not (gene in self.rename_map):
                if self.do_not_correct_names:
                    illegal_name.append(gene)
                    continue
                self.rename_map[gene] = replace_illegal_characters(gene)
                ## internally, refer to the gene by its old name!
                ## the corrected name will be written to the output files
                # gene = self.rename_map[gene]
            ## filter out the rejected isoforms 
            ## unless removing the "soft"-filtered items 
            ## leaves the gene without any isoforms 
            _trs = [x for x in trs if x not in self.rejection_reasons[TRANSCRIPT]]
            if not _trs:
                ## if there are none, resuscitate the soft-filtered items
                _trs = [
                    x for x in trs if min(self.rejection_reasons[TRANSCRIPT][x]) > RejEnum.ORPHAN
                ]
                ## if all items were hard-filtered, discard the gene
                if not _trs:
                    rejected_genes.append(gene)
                    self.gene2trs.pop(gene)
                    continue
                ## otherwise, retain the transcripts
                for x in _trs:
                    self.rejection_reasons[TRANSCRIPT].pop(x) 
            self.gene2trs[gene] = _trs

        trs_found: set[str] = {x for y in self.gene2trs.values() for x in y}

        ## report genes which ended up having no transcripts
        if rejected_genes:
            log_msg: str = "\n\t".join(map(str, rejected_genes))
            self._to_log(
                (
                    "The following genes were removed from the "
                    "isoforms file because all the respective "
                    f"transcripts were removed from the annotation:\n\t{log_msg}"
                ),
                "warning",
            )
            for i in rejected_genes:
                self.rejection_reasons[GENE][i].append(RejEnum.GENE)
        ## report genes subject to/discarded du to illegal symbols in their names
        if illegal_name:
            log_msg: str = "\n\t".join(illegal_name)
            if self.do_not_correct_names:
                self._to_log(
                    (
                        "The following genes were filtered out "
                        f"due to illegal symbols used in their names:\n\t{log_msg}"
                    ),
                    "warning",
                )
                for i in illegal_name:
                    self.rejection_reasons[GENE][i].append(RejEnum.NAME)
            else:
                self._to_log(
                    (
                        "Names for the following genes were corrected "
                        f"due to illegal symbols found:\n\t{log_msg}"
                    ),
                    "warning",
                )
        ## report transcripts for which genes were not found in the isoform file
        rejected_transcripts: list[str] = [
            x for x in self.tr2annot if x not in trs_found and 
            x not in self.rejection_reasons[TRANSCRIPT]
        ]
        if rejected_transcripts:
            log_msg: str = "\n\t".join(map(str, rejected_transcripts))
            self._to_log(
                (
                    "The following transcripts were removed from the "
                    "annotation because they were not mapped to any gene "
                    f"in the isoform file:\n\t{log_msg}"
                ),
                "warning",
            )
            for i in rejected_transcripts:
                self.rejection_reasons[TRANSCRIPT][i].append(RejEnum.ORPHAN)
                self.tr2gene.pop(i)
                self.tr2annot.pop(i)

    def report_corrupted_models(self) -> None:
        """
        Checks the input annotation for fusions and fissions

        1. Fusions are defined as isoforms or genes overlapping two or more 
        separate genes that do not overlap each other:
        ```
        gene1     ====g1[***]====[***]============[***]====[***]g1===
        genes2,3  ====g2[***]====[***]g2========g3[***]====[***]g2===
        ```
        In strict terms, TOGA2 defines fusions as articulation nodes in the 
        overlap graph between gene "units" of same-gene transcripts.
        If `remove_fusion` is set to True, potential fusion genes and isoforms 
        are removed from the final results. Otherwise, they are retained in the 
        annotation files and listed in the output file "fusion_map.tsv".

        2. Fissions are defined as genes (gene models) containing of two or more 
        clusters of non-overlapping isoforms:
        ```
        isoA     [***]=====[***]=====[***]
        isoB                                    [***]===[***]=======[***]
        ```
        Each subunit resulting from the gene model fission is assigned its own name 
        derived from the progenitor gene's name with suffix added ('_A', ..., '_Z', '_AA', ...), 
        If `remove_fissions` option set to True, split gene models are removed together with all 
        their isoforms. Otherwise,  affected isoforms are then reassigned to their 
        respective subunits in the output files; mapping between original name and 
        subunit names is listed in "fission_map.tsv".

        Arguments:
            None

        Returns:
            None
        """
        import networkx as nx

        gene_coords: dict[tuple[str, bool], list[tuple[...]]] = defaultdict(list)
        ## sort transcripts by chromosome and strand
        for name, tr in self.tr2annot.items():
            if name not in self.rejection_reasons[TRANSCRIPT]:
                gene_coords[(tr.chrom, tr.strand)].append((name, tr))
        fusion_genes: set[str] = set()
        fusion_isoforms: set[str] = set()
        split_genes: list[str | tuple[str, list[str]]] = []
        for transcripts in gene_coords.values():
            transcripts.sort(key=lambda x: (x[1].start, x[1].end, x[0]))
            ## get all fissions and fusions
            report, members = self._infer_fusions_and_fissions(transcripts, nx)
            ## process and report fusions first
            flagged, genes = self._process_fusions(report.fusions, members)
            fusion_isoforms |= flagged
            fusion_genes |= genes
            ## if fusions are requested to be removed, 
            ## run the analysis for the second time to recover fission events
            ## covered by fusions
            if flagged and self.remove_fusions:
                transcripts = [x for x in transcripts if x[0] not in flagged]
                report, members = self._infer_fusions_and_fissions(transcripts, nx)
            split_genes.extend(self._process_fissions(report.fissions, members))

        ## log all finding
        if split_genes:
            if self.remove_fissions:
                self._log_item_list(
                    "The following fragmented gene models will be deleted, "
                    "alongside with all their isoforms:",
                    sorted(split_genes),
                )
            else:
                self._log_item_list(
                    "The following fragmented gene models were split into new individual genes:",
                    sorted([x + ": " + ", ".join(self.gene2fissions[x]) for x in split_genes]),
                )
        if fusion_genes:
            if self.remove_fusions:
                self._log_item_list(
                    "The following genes were removed from the annotation "
                    "as potential fusions, together with all their isoforms:",
                    sorted(fusion_genes),
                )
            else:
                self._log_item_list(
                    "The following genes were classified as potential fusions:",
                    sorted(fusion_genes),
                )
        if fusion_isoforms:
            if self.remove_fusions:
                self._log_item_list(
                    "The following isoforms were removed from the annotation "
                    "as potential fusions:",
                    sorted(fusion_isoforms)
                )
            else:
                self._log_item_list(
                    "The following isoforms were marked as potential fusions:",
                    sorted(fusion_isoforms),
                )
        ## cleanup: remove faulty genes and their isoforms from the annotation
        self.check_isoforms()

    def isoform_length_filter(self) -> None:
        """
        For each gene, applies relative isoform length fitler
        Arguments:
            None

        Returns:
            None
        """
        discarded: list[str] = []
        for gene, trs in self.gene2trs.items():
            tr_lengths: dict[str, int] = {tr: self.tr2annot[tr].cds_len() for tr in trs}
            greatest_len: int = max(tr_lengths.values())
            len_thresh: int = int(greatest_len * self.rel_thresh)
            passed_trs: list[str] = [x for x, y in tr_lengths.items() if y >= len_thresh]
            discarded_trs: list[str] = [x for x in trs if x not in passed_trs]
            for tr in discarded_trs:
                self.rejection_reasons[TRANSCRIPT][tr].append(RejEnum.LENGTH_FILTER)
            discarded.extend(discarded_trs)
            self.gene2trs[gene] = passed_trs
        self._log_item_list(
            f"The following transcripts are shorted than {self.rel_thresh} times the length "
            "of their genes' longest isoforms and are therefore discarded", 
            discarded,
        )


    def write_annotation(self) -> None:
        """
        Writes the filtered reference annotation to the file

        Arguments:
            None

        Returns:
            None
        """
        if not self.tr2annot:
            self._die("All transcripts were filtered out for various reasons")
        self._to_log(f"Writing the filtered annotation to {self.filtered_annotation}")
        with (
            open(self.filtered_annotation, "w") as h,
            open(self.immune_annotation, "w") as ih,
        ):
            for name, tr in self.tr2annot.items():
                if name in self.rejection_reasons[TRANSCRIPT]:
                    continue
                corr_name: str = self.rename_map.get(name, name)
                if not self.no_gene_suffix:
                    corr_name = f"{corr_name}#{self.tr2gene[name]}"
                line: str = tr.line
                line = line.replace(name, corr_name)
                if name in self.immune:
                    ih.write(line + "\n")
                else:
                    h.write(line + "\n")

    def write_isoforms(self) -> None:
        """
        Writes the filtered isoform mapping and coordinates to the respective files

        Arguments:
            None

        Returns:
            None

        Note:
            Creates and populates four different files:
            * ``self.filtered_isoforms``: main isoform mapping file (TSV)
            * ``self.gene_bed``: main gene coordinates file (BED6)
            * ``self.immune_isoforms``: immunoglobulin isoform mapping file (TSV)
            * ``self.immune_gene_bed``: immunoglobulin gene coordinates file (BED6)
        """
        with (
            open(self.filtered_isoforms, "w") as iso,
            open(self.gene_bed, "w") as bed,
            open(self.immune_isoforms, "w") as i_iso,
            open(self.immune_gene_bed, "w") as i_bed,
        ):
            for gene, trs in self.gene2trs.items():
                gene_name: str = self.rename_map.get(gene, gene)
                chrom: str | None = None
                start: int | None = None
                end: int | None = None
                strand: bool | None = None
                immune: bool = False
                for tr in trs:
                    tr_name: str = self.rename_map.get(tr, tr)
                    if not self.no_gene_suffix:
                        tr_name = f"{tr_name}#{gene}"
                    if tr in self.immune:
                        i_iso.write(f"{gene_name}\t{tr_name}\n")
                        immune = True
                    else:
                        iso.write(f"{gene_name}\t{tr_name}\n")
                    _chrom: str = self.tr2annot[tr].chrom
                    if chrom is None:
                        chrom = _chrom
                    elif chrom != _chrom:
                        self._die(
                            f"Conflicting chromosomes encountered for gene {gene}: {chrom} and {_chrom}"
                        )
                    _start: int = self.tr2annot[tr].start
                    start = _start if start is None else min(start, _start)
                    _end: int = self.tr2annot[tr].end
                    end = _end if end is None else max(end, _end)
                    _strand: bool = self.tr2annot[tr].strand
                    if strand is None:
                        strand = _strand
                    elif strand != _strand:
                            self._die(f"Conflicting strand encountered for gene {gene}")
                gene_bed_line: str = GENE_BED_LINE.format(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=gene_name,
                    strand=("+" if strand else "-"),
                )
                if immune:
                    i_bed.write(gene_bed_line)
                else:
                    bed.write(gene_bed_line)

    def write_rejection_log(self) -> None:
        """
        Writes the rejected items to the file

        Arguments:
            None

        Returns:
            None
        """
        if not any(x for x in self.rejection_reasons.values()):
            return
        self._to_log(f"Writing the rejected items to {self.rejection_log}")
        with open(self.rejection_log, "w") as h:
            h.write(Headers.REJ_LOG_HEADER)
            h.writelines([
                self._build_rejection_line(name, level) 
                for level, items in self.rejection_reasons.items()
                for name in items
            ])

    def write_rename_map(self) -> None:
        """
        Writes the original_name -> corrected_name mapping to the file

        Arguments:
            None

        Returns:
            None
        """
        if not self.rename_map:
            return
        with open(self.rename_map_file, "w") as h:
            h.writelines([f"{k}\t{v}\n" for k,v in self.rename_map.items()])

    def write_fission_map(self) -> None:
        """
        Writes mapping between the original split gene model 
        and the new inferred genes to the file.

        The resulting file is a "wide-data" tab-separated table 
        with original gene names in the first column and the new split 
        in the second column

        Arguments:
            None

        Returns:
            None
        """
        if not self.gene2fissions:
            return
        with open(self.fission_map_file, "w") as h:
            h.writelines([
                f"{self.rename_map.get(k, k)}\t{x}\n" for k,v in self.gene2fissions.items() 
                for x in v if x not in self.rejection_reasons[GENE]
            ])

    def write_fusion_map(self) -> None:
        """
        Writes all fusion isoforms and genes they bridge to the file

        Arguments:
            None

        Returns:
            None

        Note:
            * In the output table, the first column contains fusion isoforms; 
            the column number corresponds to the number of bridged genes and 
            might therefore differ between the lines.
            * For each gene, all of its isoforms are listed, regardless of whether 
            they overlap the bridging isoform
        """
        if not self.fusions:
            return
        with open(self.fusion_map_file, "w") as h:
            for fusion in sorted(self.fusions):
                out_line: list[str] = [fusion]
                for fused_locus in sorted(self.fusions[fusion]):
                    fused_locus = [x for x in fused_locus if x not in self.rejection_reasons[TRANSCRIPT]]
                    if not fused_locus:
                        continue
                    fused_gene: str = self.tr2gene[fused_locus[0]]
                    if fused_gene in self.rejection_reasons[GENE]:
                        continue
                    fused_locus: str = self.rename_map.get(fused_gene, fused_gene) + ":" + ",".join([x for x in fused_locus])
                    out_line.append(fused_locus)
                if len(out_line) < 2:
                    continue
                h.write("\t".join(out_line) + "\n")
            # h.writelines(
            #     [x + "\t" +"\t".join([",".join(z) for z in y]) + "\n" for x, y in self.fusions.items()] 
            # )

    def create_sleasy_input(self) -> None:
        """
        Creates a 2bit exon storage for SLEASY compatibility

        Arguments:
            None

        Returns:
            None

        Raises:
            Error-exits if any of the launched subprocesses 
            (bed12ToFraction, twoBitToFa, faToTwoBit) fails or if provisional BED file 
            was improperly sorted
        """

        self._to_log("Writing temporary BED6 file")
        bed6_cmd: str = (
            f"{self.bed2fraction_binary} -i {self.filtered_annotation} -m cds -b | "
            "sort -k4,4 -k5,5n | "
            f"awk -F'\t' '{{print $1,$2,$3,$4\"_exon\"$5,$5,$6}}' > {self.bed6_exons}"
        )
        _ = self._exec(bed6_cmd, BED12TO6_ERR)
        self._to_log("Writing temporary BED6 file complete")

        self._to_log("Extracting sequences from the 2bit file; might take time")
        twobit2fa_cmd: str = (
            f"{self.twobittofa_binary} -bed={self.bed6_exons} {self.twobit} stdout"
        )
        fasta_lines: str = self._exec(
            twobit2fa_cmd, EXTRACTION_ERR_MSG, gather_stdout=True
        )
        self._to_log("Sequence extraction complete")

        self._to_log("Writing temporary FASTA file")
        with open(self.fasta_for_sleasy, "w") as h:
            name: str = ""
            exon: int = 0
            prev_name: str = ""
            prev_exon: int = 0
            curr_seq: str = ""
            new_exon: bool = True
            for line in fasta_lines.split("\n"):
                if not line:
                    continue
                if line[0] == ">":
                    name, exon = line.lstrip(">").split("_exon")
                    if not prev_name:
                        prev_name = name
                    exon: int = int(exon)
                    if name != prev_name and prev_name:
                        h.write(">" + prev_name + "\n" + curr_seq + "\n")
                        curr_seq = ""
                        prev_name = name
                        prev_exon = 0
                    if prev_exon >= exon:
                        self._die(
                            "Temporary BED file was improperly sorted; "
                            "check that twoBitToFa returned sequences as presented in the BED file. "
                            f"Troublemaker: {name}, exons {prev_exon}-{exon}"
                        )
                    prev_exon = exon
                    new_exon: bool = exon > 1
                    continue
                if new_exon:
                    curr_seq += SEP_DUMMY
                    new_exon = False
                curr_seq += line.upper()
            if curr_seq:
                h.write(">" + prev_name + "\n" + curr_seq + "\n")
        self._to_log("Writing temporary FASTA file complete")

        self._to_log("Converting reference exons into 2bit format")
        fa2twobit_cmd: str = (
            f"{self.fatotwobit_binary} {self.fasta_for_sleasy} {self.sleasy_2bit}"
        )
        _ = self._exec(fa2twobit_cmd, CONVERSION_ERR_MSG)
        self._to_log("Execution successfully completed; cleaning up and exiting")

    def intron_classifier(self) -> None:
        """
        Classifies the introns in the filtered annotation, separating them
        into U2 and U12 classes according to the intronIC predictions and adding
        terminal dinucleotide data for further classification into canonical (GT-AG)
        and non-canonical (any other dinucleotide combination) subclasses.

        The resulting file is a Bed6 file of class- and dinucleotide-annotated reference introns.
        This file can be further used with TOGA2 runs provided with the --ut12_file/-u2 option.
        If CESAR2 profile annotation is requested, the same predictions are further used
        to generate reference-specific profiles.

        Arguments:
            None

        Returns:
            None

        Raises:
            Error-exists if any of subprocesses 
            (twoBitToFa, bed12ToFraction, BED sorting, intronIC) fails
        """
        ## temporarily decompress the reference .2bit genome file
        self._to_log("Converting .2bit genome into Fasta format")
        genome_fasta: str = os.path.join(self.tmp_dir, "genome.fa")
        decompr_cmd: str = f"{self.twobittofa_binary} {self.twobit} {genome_fasta}"
        _ = self._exec(decompr_cmd, "twoBitToFa conversion failed:")
        ## get the intron Bed6 file
        self._to_log("Extracting the unique coding introns")
        raw_intron_file: str = os.path.join(self.tmp_dir, "raw_introns.bed")
        intron_bed_cmd: str = (
            f"{self.bed2fraction_binary} -i {self.filtered_annotation} -o {raw_intron_file} "
            "-m cds -n -b"
        )
        _ = self._exec(intron_bed_cmd, "Intron Bed6 extraction failed:")
        ## select unique introns and anonimize them
        ## TODO: Must be merged with the previous command after writing to stdout is fixed in bed12ToFraction
        intron_bed: str = os.path.join(self.tmp_dir, "intronic_input.bed")
        uniq_cmd: str = (
            f"cut -f1-3,5,6 {raw_intron_file} | sort -u | "
            'awk \'BEGIN{OFS="\\t"}{print $1,$2,$3,"intron"NR,0,$5}\' > '
            f"{intron_bed}"
        )
        _ = self._exec(uniq_cmd, "Intron deduplication failed:")
        ## run intronIC
        self._to_log("Running intronIC")
        ic_output: str = os.path.join(self.tmp_dir, "output")
        intronic_cmd: str = (
            f"{self.intronic_binary} -g {genome_fasta} -b {intron_bed} -n {ic_output} "
            f"-p {self.intronic_cores}  --min_intron_len {self.min_intron_length_intronic} "
            "--no_nc_ss_adjustment --no_abbreviate"
        )
        _ = self._exec(intronic_cmd, "intronIC run failed:")
        ## now, prepare the final file
        ## parse the output bed file
        ## score field contains probability, in per cent value
        ## for the final Bed file, those values are multiplied by 10 (to confine within 0<=x<=1000)
        ## and rounded to the closest integer
        out_bed_file: str = f"{ic_output}.bed.iic"
        intron2coords: dict[str, str] = {}
        with open(out_bed_file, "r") as h:
            for i, line in enumerate(h, start=1):
                data: list[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) < 6:
                    self._die(
                        f"Improper intronIC bed file formatting at line {i}; "
                        f"expected 6 fields, got {len(data)}"
                    )
                name: str = data[3].split(";")[0]
                ## convert the probability into a Bed-compatible score
                score: int = (
                    0
                    if (data[4] == "." or data[4] == "NA")
                    else int(float(data[4]) * 10)
                )
                intron2coords[name] = (
                    data[0],
                    data[1],
                    data[2],
                    "{}",
                    str(score),
                    data[5],
                )
        ## all is left is intron class and terminal dinucleotides
        ## extract those from the .meta.iic file, then write the results to file
        out_meta_file: str = f"{ic_output}.meta.iic"
        with (
            open(out_meta_file, "r") as ih,
            open(self.intron_file, "w") as oh,
            (
                nullcontext
                if self.disable_cesar_profiles
                else open(self.all_intron_bed, "w")
            ) as ah,
        ):
            for i, line in enumerate(ih, start=1):
                data: list[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) < 14:
                    self._die(
                        f"Improper intronIC file formatting at line {i}; "
                        f"expected at least 14 fields, got {len(data)}"
                    )
                if data[0] == NAME:
                    try:
                        class_col: int = data.index(TYPE_ID)
                    except ValueError:
                        class_col: int = DEFAULT_CLASS_COL
                    continue
                name: str = data[0].split(";")[0]
                dinuc: str = data[2]
                intron_class: str = data[class_col].upper()
                upd_name: str = f"{name}_{intron_class}_{dinuc}"
                out_line: str = "\t".join(intron2coords[name]).format(upd_name)
                ## for TOGA2 input, only U12 and non-canonical U2 introns are required
                is_canon: bool = (intron_class == U2 and dinuc in CANON_SITES) or (
                    intron_class == U12 and dinuc == U12_CANON_SITES
                )
                if intron_class == U12 or not is_canon:
                    oh.write(out_line + "\n")
                ## if CESAR2 profile generation was requested,
                ## save all the introns to the provisional file
                ## and store their class data
                if not self.disable_cesar_profiles:
                    ah.write(out_line + "\n")
                    chrom, start, end, _, _, strand = intron2coords[name]
                    start = int(start)
                    end = int(end)
                    strand: bool = strand == "+"
                    if end - start < self.min_intron_length_cesar:
                        continue
                    intron_key: tuple[str, int, int, bool] = (chrom, start, end, strand)
                    num: str = name.replace("intron", "")
                    self.intron2class[num] = (intron_class, is_canon)
                    self.intron2coords[num] = intron_key

        ## if profile generation is disabled and/or not requested to leave,
        ## remove the temporary directory
        if not self.keep_tmp or self.disable_cesar_profiles:
            self._rmdir(self.tmp_dir)

    def generate_cesar_profiles(self) -> None:
        """
        Given the intronIC classification for filtered annotation reference introns, 
        creates reference-specific CESAR2 profiles

        Arguments:
            None

        Returns:
            None

        Raises:
            Error-exits if the launched twoBitToFa subprocess fails
        """
        ## step 1: prepare profile sequence coordinates in Bed format
        bed_string: str = ""
        ## extract coordinates
        for num, coords in self.intron2coords.items():
            chrom, start, end, strand = coords
            ## define donor and acceptor sequence coordinates
            if strand:
                donor_start: int = start
                donor_end: int = min(donor_start + DONOR_PROFILE_LEN, end)
                acc_start: int = max(end - ACC_PROFILE_LEN, start)
                acc_end: int = end
            else:
                donor_start: int = max(end - DONOR_PROFILE_LEN, start)
                donor_end: int = end
                acc_start: int = start
                acc_end: int = min(start + ACC_PROFILE_LEN, end)
            bed_strand: str = "+" if strand else "-"
            donor_name: str = f"{num}_{DONOR}"
            acc_name: str = f"{num}_{ACC}"
            donor_bed_line: str = "\t".join(
                map(str, [chrom, donor_start, donor_end, donor_name, 0, bed_strand])
            )
            acc_bed_line: str = "\t".join(
                map(str, [chrom, acc_start, acc_end, acc_name, 0, bed_strand])
            )
            bed_string += donor_bed_line + "\n" + acc_bed_line + "\n"
        ## then extract the sequences in Fasta format
        cmd: str = f"{self.twobittofa_binary} -bed=/dev/stdin {self.twobit} stdout"
        res: str = self._exec(
            cmd, err_msg=EXTRACTION_ERR_MSG, input_=bed_string.encode("utf8")
        )

        ## step 2: parse the resulting Fasta and update the resulting profiles
        self._to_log("Inferring positional letter probabilities")
        header: str = ""
        seq: str = ""
        for line in res.split("\n"):
            line = line.rstrip()
            if not line:
                continue
            if line[0] == ENTRY_START:
                if header:
                    num, site_type = header.split("_")
                    seq = seq.upper().strip("\n")
                    intron_class, canon = self.intron2class[num]
                    ## update the respective CESAR2 profile
                    self._update_profile(seq, intron_class, canon, site_type)
                    seq = ""
                header = line[1:]
            else:
                seq += line
        if seq:
            num, site_type = header.split("_")
            seq = seq.upper()
            intron_class, canon = self.intron2class[num]
            ## again, update the respective profile for the final item
            self._update_profile(seq, intron_class, canon, site_type)

        ## step 3: compute the letter probabilities per position from the recorded frequencies
        ## and write the resulting profiles
        self._to_log(
            f"Computing CESAR2 profiles and writing them to {self.profile_dir}" 
        )
        header: str = "\t".join(NUCS)
        for key in self.profiles:
            intron_class, canon, site_type = key
            canon_line: str = CANON if canon else NONCANON
            filename: str = f"{canon_line}_{intron_class}_{site_type}.tsv"
            with open(os.path.join(self.profile_dir, filename), "w") as h:
                h.write(header + "\n")
                for pos in sorted(self.profiles[key].keys()):
                    pos_sum: float = float(sum(self.profiles[key][pos].values()))
                    frequencies: list[float] = []
                    for nuc in NUCS:
                        freq: float = round(
                            self.profiles[key][pos].get(nuc, 0.0) / pos_sum, 3
                        )
                        frequencies.append(freq)
                    h.write("\t".join(map(str, frequencies)) + "\n")

        ## generate equiprobable acceptor and donor profiles
        ## equiprobable acceptor is recommended for non-canonical U12 in mammals
        ## equiprobable donors have not been tested for any purpose yet
        ## but it's nice to have both options generated automatically
        equi_line: str = "\t".join(map(str, [0.25] * 4))# + "\n"
        equi_acc: str = os.path.join(self.profile_dir, EQUI_ACC)
        with open(equi_acc, "w") as h:
            h.write(header + "\n")
            h.writelines([equi_line + "\n"] * ACC_PROFILE_LEN)
        equi_donor: str = os.path.join(self.profile_dir, EQUI_DONOR)
        with open(equi_donor, "w") as h:
            h.write(header + "\n")
            h.writelines([equi_line + "\n"] * DONOR_PROFILE_LEN)

        ## remove temporary directory unless requested to leave
        if not self.keep_tmp:
            self._rmdir(self.tmp_dir)

    def _update_profile(
        self, seq: str, intron_class: str, canon: bool, site: str,
    ) -> None:
        """
        Updates the CESAR2 profile with respect to spliceosomal class, 'canonicity', 
        and splice site

        Arguments:
            seq: str, nucleotide sequence
            intron_class: str, either U2 or U12 class
            canon: bool, intron canonicity indicator
            site: str, either donor or acceptor identifier

        Returns:
            None

        Raises:
            Error-exits if:
                * ambiguous site identifier was provided
                * input sequence is longer than consensus for the specified splice site
                * ambiguous symbol has been encountered in the sequence

        Note:
            Updates probability matrix per position for the indicated intron class and splice site
        """
        if site not in (DONOR, ACC):
            self._die(f"Unknown site type received as input: {site}")
        key: tuple[str, bool, str] = (intron_class, canon, site)
        ## CESAR2 operates fixed profile lengths for donor and acceptor
        exp_len: int = DONOR_PROFILE_LEN if site == DONOR else ACC_PROFILE_LEN
        ## sanity check: extracted sequence must not be longer than the profile length
        ## the reverse indicates a clear extraction error
        if len(seq) > exp_len:
            self._die(
                f"Sequence {seq} is longer than the expected {site} profile length of {exp_len}"
            )
        ## earlier TOGA2 code threw the same error if the profile was shorter than the consensus length
        ## however, for certain clades introns <22bp are quite common
        if len(seq) < exp_len:
            self._to_log(
                f"Sequence {seq} is shorter than the expected {site} profile length of {exp_len}",
                "warning",
            )
        ## position numeration differs
        start_pos: int = DONOR_START_POS if site == DONOR else ACC_START_POS
        for i in range(exp_len):
            ## WARNING: The trick was not tested in TOGA2 alpha since introns <30 bp
            ## were ignored in production runs and introns <70bp were excluded from profile preparation
            ## but the idea is to make position beyond the average intron length potentially equiprobable
            pos: int = i + start_pos
            if i >= len(seq):
                for nuc in NUCS:
                    self.profiles[key][pos][nuc] += 1
                continue
            nuc: str = seq[i]
            if nuc == N:
                continue
            if nuc not in NUCS:
                self._die(f"Ambiguous nucleotide encountered: {nuc}")
            self.profiles[key][pos][nuc] += 1

    def _infer_fusions_and_fissions(
        self, transcripts: list[tuple[str, Transcript]], nx: ModuleType
    ) -> tuple[FaultyModelReport, dict[tuple[str, int], list[str]]]:
        """
        An auxiliary method for gene model scanning: For a given list 
        of sorted transcripts, assess all contained gene models for 
        model fusions and fissions.

        Arguments:
            transcripts: list, (transcript_name, Transcript) tuples sorted by span coordinates
            nx: module, a (lazily imported) NetworkX module

        Returns:
            tuple, (FaultyModelReport, {unit: [transcripts]})
        """
        n: int = len(transcripts)
        intra_graph: nx.Graph = nx.Graph()
        intra_graph.add_nodes_from(tr for tr, _ in transcripts if tr not in self.rejection_reasons[TRANSCRIPT])
        inter_graph: nx.Graph = nx.Graph()
        inter_edges: list[tuple[str, str]] = []
        rank: dict[str, int] = {name: i for i, (name, _) in enumerate(transcripts)}

        ## start building the overlap graphs
        for i, (out_name, tr_out) in enumerate(transcripts):
            if out_name in self.rejection_reasons[TRANSCRIPT]:
                continue
            gene_out: str = self.tr2gene[out_name]
            for j in range(i + 1, n):
                in_name, tr_in = transcripts[j]
                if in_name in self.rejection_reasons[TRANSCRIPT]:
                    continue
                if tr_out.end <= tr_in.start:
                    break
                if cds_exon_overlap(tr_out, tr_in):
                    same_gene: bool = self.tr2gene[in_name] == gene_out
                    if same_gene: 
                        ## add the edge to the unit graph 
                        ## for further unit-level analysis
                        intra_graph.add_edge(out_name, in_name)
                    else:
                        ## record the edge for further unit overlap graph construction
                        inter_edges.append((out_name, in_name))
        ## collapse overlapping transcripts of the same gene into gene units:
        ## each collapsed component is represented as (gene_name, counter),
        ## where `counter` is the number of occurrences of this gene 
        ## up to this point minus one
        unit: dict[str, tuple[str, int]] = {} ## {transcript -> (gene, gene_counter)}
        members: dict[tuple[str, int], list[str]] = {} ## {unit -> transcripts} mapping
        counter: dict[str, int] = defaultdict(int) ## gene counter
        for i_component in nx.connected_components(intra_graph):
            gene: str = self.tr2gene[next(iter(i_component))]
            u: tuple[str, int] = (gene, counter[gene])
            counter[gene] += 1
            members[u] = sorted(i_component, key=rank.__getitem__)
            for tr in i_component:
                unit[tr] = u
        ## create an overlap graph between collapsed units
        cross_unit: dict[str, set[tuple[str, int]]] = defaultdict(set)
        for u, v in inter_edges:
            a, b = unit[u], unit[v]
            inter_graph.add_edge(a, b)
            cross_unit[u].add(b)
            cross_unit[v].add(a)
        ## perform gene model check
        report: FaultyModelReport = self._detect_faulty_models(
            members, 
            cross_unit, 
            set(nx.articulation_points(inter_graph)),
            partial(self._sides, adj=inter_graph)
        )
        return (report, members)


    def _process_fusions(
        self, 
        fusions: list[FusionReport], 
        members: dict[tuple[str, int], list[str]]
    ) -> tuple[set[str], set[str]]:
        """
        Processes a list of FusionReport elements:
        1. Collects all fusion isoforms;
        2. Marks gene from multi-isoform gene units;
        3. Stages genes consisting entirely of fusion isoforms, as well as 
            individual isoforms, for rejection

        Arguments:
            fusions: list, list of FusionReport objects
            members: dict, {unit: [transcripts]} mapping
        
        Returns:
            set, a collection of isoforms scheduled for removal
        """
        iso, genes = set(), set()
        for f in fusions:
            ## ignore self-fusions
            if f.bridge_type == BridgeType.SELF:
                continue
            iso.update(f.fusion_map.keys())
            self.fusions.update(f.fusion_map)
            ## if fusion concerns an entire unit rather than separate isoforms,
            ## schedule the progenitor gene for rejection
            if f.bridge_type != BridgeType.FUS_ISO:
                genes.add(self.tr2gene[members[f.node][0]])
        ## reject the genes consisting entirely of fusion isoforms
        genes: set[str] = {
            g for g in genes if g in self.gene2trs and 
            all(t in iso for t in self.gene2trs[g] if t not in self.rejection_reasons[TRANSCRIPT])
        }
        if self.remove_fusions:
            for g in genes:
                self.rejection_reasons[GENE][g].append(RejEnum.FUSION_GENE)
                self.gene2trs.pop(g)
            ## reject the fusion isoforms
            for tr in sorted(iso):
                g: str = self.tr2gene.pop(tr)
                if g in self.gene2trs:
                    self.gene2trs[g].remove(tr)
                self.rejection_reasons[TRANSCRIPT][tr].append(RejEnum.FUSION_TR)
                self.tr2annot.pop(tr)
        return iso, genes

    def _process_fissions(
        self, 
        fissions: list[FissionReport], 
        members: dict[tuple[str, int], list[str]]
    ) -> list[str]:
        """
        Processes a list of FissionReport items:
        * if fission removal was requested, schedule all units 
        and their isoforms for deletion
        * otherwise, assign new name to each split subunit, 
        rewire gene->transcripts and transcript->gene mapping, 
        and saved the mapping between new and old names  
        to ``gene2fissions`` 

        Arguments:
            fissions: list, a collection of FissionReport items
            members: dict, {unit: [transcripts]} mapping

        Returns:
            list, a list of potentially split gene models
        """
        split: list[str] = []
        for fission in fissions:
            split.append(fission.gene)
            if self.remove_fissions:
                self.rejection_reasons[GENE][fission.gene].append(RejEnum.FISSION_GENE)
                for tr in self.gene2trs.pop(fission.gene):
                    self.rejection_reasons[TRANSCRIPT][tr].append(RejEnum.FISSION_TR)
                    self.tr2gene.pop(tr)
                    self.tr2annot.pop(tr)
                continue
            self.gene2trs.pop(fission.gene)
            names: list[str] = []
            for c, u in enumerate(fission.units):
                sub: str = f"{fission.gene}{cluster_suffix(c)}"
                self.gene2trs[sub] = members[u]
                for t in members[u]:
                    self.tr2gene[t] = sub
                names.append(sub)
            self.gene2fissions[fission.gene] = tuple(names)
        return split

    def _build_rejection_line(self, name: str, level: str) -> str:
        """
        For the given item, builds a rejection log entry 
        for each rejection reason recorded

        Arguments:
            name: str, transcript/gene name

        Returns:
            str, one or more rejection log entries
        """
        reasons: list[RejEnum] = self.rejection_reasons[level][name]
        output: list[str] = []
        for reason in reasons:
            match reason:
                case RejEnum.NON_CODING:
                    output.append(RejectionReasons.NON_CODING_REJ_REASON.format(name))
                case RejEnum.FRAME:
                    output.append(RejectionReasons.FRAME_REJ_REASON.format(name))
                case RejEnum.CONTIG:
                    output.append(RejectionReasons.CONTIG_REJ_REASON.format(name))
                case RejEnum.NAME:
                    output.append(RejectionReasons.NAME_REJ_REASON.format(level, name))
                case RejEnum.FISSION_GENE:
                    output.append(RejectionReasons.FISSION.format(level, name))
                case RejEnum.FISSION_TR:
                    output.append(RejectionReasons.FISSION.format(level, name))
                case RejEnum.FUSION_GENE:
                    output.append(RejectionReasons.FUSION.format(level, name))
                case RejEnum.FUSION_TR:
                    output.append(RejectionReasons.FUSION.format(level, name))
                case RejEnum.ORPHAN:
                    output.append(RejectionReasons.ORPHAN_TR.format(name))
                case RejEnum.SHORT_EXON:
                    output.append(
                        RejectionReasons.SHORT_EXON_REASON.format(name, self.min_int_exon_length)
                    )
                case RejEnum.SHORT_INTRON:
                    output.append(
                        RejectionReasons.SHORT_INTRON_REASON.format(name, self.min_intron_length)
                    )
                case RejEnum.NMD:
                    output.append(RejectionReasons.NMD_REASON.format(name))
                case RejEnum.LENGTH_FILTER:
                    output.append(RejectionReasons.LENGTH_FILTER.format(name, self.rel_thresh))
                case RejEnum.GENE:
                    output.append(RejectionReasons.REJ_GENE.format(name))
        return "\n".join(output) + "\n"

    def _log_item_list(self, msg: str, items: list[str]) -> None:
        """
        Logs items subjected to removal or otherwise flagged by the code:

        Arguments:
            msg: str, message header
            items: list, a collection of transcripts (genes, etc.) in need of logging

        Returns:
            None

        Note:
            Logs ``items`` in the following format:
            ```
            msg:
            \titem1
            \titem2
            ...
            \titemN
            ```
        """
        self._to_log(msg + NEWTAB + NEWTAB.join(items), "warning")

    def set_logging(self) -> None:
        """Sets up logging system for a InputProducer instance"""
        super().set_logging(toga_module="prepare-input")
        self.logger.propagate = False
