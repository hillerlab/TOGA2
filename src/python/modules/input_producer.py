"""
Reference input annotation preparation functionality
"""

__author__ = "Yury V. Malovichko"
__credits__ = ("Michael Hiller", "Bernhard Bein")
__year__ = "2025"

import logging
import os
from collections import defaultdict
from contextlib import nullcontext
from enum import IntEnum
from shutil import which
from typing import Final

import click

from .constants import Headers, RejectionReasons
from .gxf import Gxf
from .shared import (
    CommandLineManager,
    consistent_name,
    dir_name_by_date,
    get_upper_dir,
    hex_dir_name,
    intersection,
    read_tab,
)

logging.root.handlers = []

DEFAULT_PREFIX: Final[str] = "TOGA2_ref_annotation"
LEVEL: Final[str] = "level"
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


class NmdConstants:
    STOP_JUNCTION_DIST: int = 55
    WEAK_STOP_JUNCTION_DIST: int = 80
    START_STOP_DIST: int = 100
    BIG_EXON_JUNCTION_DIST: int = 400


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


class IsNMD(IntEnum):
    NO = 0
    WEAK = 1
    STRONG = 2


class InputProducer(CommandLineManager):
    """
    Main class for input annotation preparation
    """

    __slots__ = (
        'all_intron_bed', 
        'annot', 
        'annot_format', 
        'bed2fraction_binary', 
        'bed6_exons', 
        'contigs', 
        'disable_cesar_profiles', 
        'disable_intron_classification', 
        'excluded_contigs', 
        'fasta_for_sleasy', 
        'fatotwobit_binary', 
        'filtered_annotation', 
        'filtered_isoforms', 
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
        'output', 
        'profile_dir', 
        'profiles', 
        'provisional_bed', 
        'provisional_isoforms', 
        'rejected_lines', 
        'rejected_transcripts', 
        'rejection_log', 
        'sleasy_2bit', 
        'strip_versions', 
        'tmp_dir', 
        'tr2annot', 
        'tr2coords', 
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
            the distance between the end of the coding sequence and the first UTR
            exon junction is more than 55 bp, the exon is considered a potential NMD target;

            4. The following conditions might hinder NMD, so fulfilling either of them
            downgrades the evidence for NMD to the "weak evidence" category:
                * CDS-UTR junction distance is below 80 nucleotides (relaxed "55 bp rule");
                * Stop-containing exon consists of both CDS and UTR, and the UTR
                sequence is at least 400 bases long;
                * Coding sequence is shorter than 100 bases.

        Arguments:
            * thin_start: int - sequence start in the genome, including UTRs;
            * thin_end: int - sequence end in the genome, including UTRs;
            * cds_start: int - coding sequence start in the genome;
            * cds_end: int - coding sequence end in the genome;
            * sizes: List[int] - list of exon sizes (BED12 field 11);
            * sizes: List[int] - list of exon relative start coordinates (BED12 field 12);
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

    def __init__(
        self,
        ref_2bit: click.Path,
        ref_annot: click.Path,
        annot_format: str = Gxf.BED,
        ref_isoforms: click.Path = None,
        output: click.Path | None = None,
        prefix: str = "",
        disable_transcript_filtering: bool = False,
        contigs: str | None = None,
        excluded_contigs: str | None = None,
        min_internal_exon_length: int = 0,
        min_intron_length: int = 0,
        nmd_filter_level: bool = 0,
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
        self.min_int_exon_length: int = min_internal_exon_length
        self.min_intron_length: int = min_intron_length
        self.nmd_level: bool = nmd_filter_level
        self.annot_format: str = annot_format
        self.disable_intron_classification: bool = disable_intron_classification
        self.disable_cesar_profiles: bool = disable_cesar_profiles
        self.min_intron_length_intronic: int = min_intron_length_intronic
        self.min_intron_length_cesar: int = min_intron_length_cesar
        self.output: str = (
            output if output is not None else hex_dir_name(DEFAULT_PREFIX)
        )
        self.strip_versions: bool = strip_versions
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

        self.tr2annot: dict[str, str] = {}
        self.rejected_transcripts: set[str] = set()
        self.immune: set[str] = set()
        self.rejected_lines: list[str] = []
        self.tr2coords: dict[str, tuple[str, int, int, bool]] = {}

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
        if self.isoforms is not None or self.annot_format != Gxf.BED:
            self._to_log("Refining the input isoform file")
            self.check_isoforms()
        ## write the results for steps 1 and 2
        self._to_log("Writing the annotation results")
        self.write_annotation()
        self.write_rejection_log()
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
            chrom: str = data[0]
            level: str = data[2]
            strand: bool = data[6] == "+"
            attrs: dict[str, str] = parse_gxf_attrs(
                data[8], gff3=(self.annot_format == Gxf.GFF3)
            )
            biotype: str | None = Gxf.get_biotype(attrs)
            tr_id: str = (
                attrs.get(Gxf.TR_ID, "")
                if self.annot_format == Gxf.GTF
                else attrs.get(Gxf.ID, "")
            )
            tr_id = Gxf.remove_tr_prefices(tr_id)
            if not tr_id:
                continue
            parent: str = Gxf.get_exon_parent(attrs)
            ## check the item's  biotype
            if (
                biotype in Gxf.CODING_TAGS
                or level in Gxf.CODING_TAGS
                or level == Gxf.CDS
            ):
                ## messenger RNA transcript
                is_coding: bool = True
            elif biotype in Gxf.VDJ_TAGS:
                ## IG fragment
                is_coding: bool = False
            else:
                ## neither a protein-coding entity nor an immune gene fragment;
                ## proceed further
                continue
            start: int = int(data[3]) - 1  ## GXF is 1-based, BED is 0-based
            end: int = int(data[4])
            # if level == Gxf.TRANSCRIPT:
            if level in Gxf.CODING_TAGS:
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
                continue
            if level not in Gxf.EXON_TYPES:
                continue
            parent: str = Gxf.get_exon_parent(attrs)
            ## an orphan exon encountered; continue
            if not parent or parent not in tr_objects:
                continue
            ## update transcript coding sequence boundaries
            if level == Gxf.CDS:
                if tr_objects[parent].cds_start is None:
                    tr_objects[parent].cds_start = start
                else:
                    tr_objects[parent].cds_start = min(
                        tr_objects[parent].cds_start, start
                    )
                if tr_objects[parent].cds_end is None:
                    tr_objects[parent].cds_end = end
                else:
                    tr_objects[parent].cds_end = max(tr_objects[parent].cds_end, end)
            ## add the exon to the coding sequence
            e: int = 0
            while e < len(tr_objects[parent].exons):
                prev_start, prev_end = tr_objects[parent].exons[e]
                if intersection(start, end, prev_start, prev_end) >= 0:
                    tr_objects[parent].exons.pop(e)
                    start, end = min(start, prev_start), max(end, prev_end)
                    continue
                e += 1
            tr_objects[parent].exons.append((start, end))
            tr_objects[parent].exons.sort(key=lambda x: (x[0], x[1]))
            ## the outcommented block below can significantly speed the exon addition up
            ## but it's better to assume that the user might provide an unsorted GXF file (why tho?)
            # prev_start, prev_end = tr_objects[parent].exons.pop(-1)
            # if intersection(start, end, prev_start, prev_end) >= 0:
            #     tr_objects[parent].exons.append(
            #         (min(start, prev_start), max(end, prev_end))
            #     )
            # else:
            #     tr_objects[parent].exons.extend(((prev_start, prev_end), (start, end)))
        with (
            open(self.provisional_bed, "w") as bed,
            open(self.provisional_isoforms, "w") as iso,
        ):
            for tr_obj in tr_objects.values():
                sizes: list[str] = []
                starts: list[str] = []
                count: int = 0
                for e_start, e_end in tr_obj.exons:
                    start: str = str(e_start - tr_obj.start)
                    size: str = str(e_end - e_start)
                    starts.append(start)
                    sizes.append(size)
                    count += 1
                sizes: str = ",".join(sizes) + ","
                starts: str = ",".join(starts) + ","
                strand: str = "+" if tr_obj.strand else "-"
                if self.strip_versions:
                    tr_obj.name = self.strip_v(tr_obj.name)
                    # tr_obj.parent = self.strip_v(tr_obj.parent)
                name: str = f"{tr_obj.name}#{tr_obj.parent}"
                if not tr_obj.is_coding:
                    cds_start: int = (
                        tr_obj.start if tr_obj.cds_start is None else tr_obj.cds_start
                    )
                    cds_end: int = (
                        tr_obj.start if tr_obj.cds_end is None else tr_obj.cds_end
                    )
                    self.immune.add(name)
                else:
                    cds_start: int = tr_obj.cds_start
                    cds_end: int = tr_obj.cds_end
                    if cds_start is None or cds_end is None:
                        continue
                bed_line: str = BED_LINE.format(
                    chrom=tr_obj.chrom,
                    start=tr_obj.start,
                    end=tr_obj.end,
                    name=name,
                    strand=strand,
                    cds_start=cds_start,
                    cds_end=cds_end,
                    exon_num=count,
                    sizes=sizes,
                    starts=starts,
                )
                bed.write(bed_line)
                iso_line: str = ISOFORM_LINE.format(gene=tr_obj.parent, tr=name)
                iso.write(iso_line)

    def check_annotation(self) -> None:
        """
        Filters reference annotation by the following criteria:
        * All transcripts in the final annotation must be

        Arguments:
            None

        Returns:
            None
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
                continue
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
                    "coding sequence start coordinate greated than the start coordinate"
                )
            if cds_start == cds_end:
                non_coding.add(name)
                continue
            ## iterate over exon entries to infer the CDS length
            frame_length: int = 0
            sizes: list[int] = [int(x) for x in data[10].split(",") if x]
            starts: list[int] = [int(x) for x in data[11].split(",") if x]
            flawed: bool = False
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
                    flawed = True
            has_short_exons: bool = False
            prev_end: int | None = None
            has_short_introns: bool = False
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
                if prev_end is not None:
                    intron_size: int = start - prev_end
                    if intron_size < self.min_intron_length:
                        has_short_introns = True
                prev_end = end
            if has_short_exons and is_coding:
                short_exons.add(name)
                flawed = True
            if has_short_introns and is_coding:
                short_introns.add(name)
                flawed = True
            if frame_length % 3:  # and not self.no_frame_filter:
                out_of_frame.add(name)
                flawed = True
            if flawed:
                continue
            self.tr2annot[name] = line
            self.tr2coords[name] = (chrom, cds_start, cds_end, strand)

        ## report the discarded items
        if illegal_name:
            log_msg: str = "\n\t".join(illegal_name)
            self._to_log(
                (
                    "The following transcripts were filtered out "
                    f"due to illegal symbols used in their names:\n\t{log_msg}"
                ),
                "warning",
            )
            self.rejected_transcripts.update(illegal_name)
            self.rejected_lines.extend(
                [RejectionReasons.NAME_REJ_REASON.format(x) for x in illegal_name]
            )
        if rejected_contigs:
            log_msg: str = "\n\t".join(rejected_contigs)
            self._to_log(
                (
                    "The following transcripts were filtered out "
                    f"due to their location in deprecated contigs:\n\t{log_msg}"
                ),
                "warning",
            )
            self.rejected_transcripts.update(rejected_contigs)
            self.rejected_lines.extend(
                [RejectionReasons.CONTIG_REJ_REASON.format(x) for x in rejected_contigs]
            )
        if non_coding:
            log_msg: str = "\n\t".join(non_coding)
            self._to_log(
                (
                    "The following transcripts were filtered out "
                    f"due to absence of coding sequence:\n\t{log_msg}"
                ),
                "warning",
            )
            self.rejected_transcripts.update(non_coding)
            self.rejected_lines.extend(
                [RejectionReasons.NON_CODING_REJ_REASON.format(x) for x in non_coding]
            )
        if short_exons:
            log_msg: str = "\n\t".join(short_exons)
            self._to_log(
                (
                    "The following transcripts were filtered out because they have "
                    f"short (<{self.min_int_exon_length} bp) internal exons: {log_msg}"
                ),
                "warning",
            )
            self.rejected_transcripts.update(short_exons)
            self.rejected_lines.extend(
                [RejectionReasons.SHORT_EXON_REASON.format(x) for x in short_exons]
            )
        if short_introns:
            log_msg: str = "\n\t".join(short_introns)
            self._to_log(
                (
                    "The following transcripts were filtered out because they have "
                    f"short (<{self.min_intron_length} bp) CDS introns:\n\t{log_msg}"
                ),
                "warning",
            )
            self.rejected_transcripts.update(short_introns)
            self.rejected_lines.extend(
                [
                    RejectionReasons.SHORT_INTRON_REASON.format(
                        x, self.min_intron_length
                    )
                    for x in short_introns
                ]
            )
        if out_of_frame:
            log_msg: str = "\n\t".join(out_of_frame)
            self._to_log(
                (
                    "The following transcripts were filtered out "
                    f"due to shifted reading frame:\n\t{log_msg}"
                ),
                "warning",
            )
            self.rejected_transcripts.update(out_of_frame)
            self.rejected_lines.extend(
                [RejectionReasons.FRAME_REJ_REASON.format(x) for x in out_of_frame]
            )
        if nmd_targets:
            log_msg: str = "\n\t".join(nmd_targets)
            self._to_log(
                (
                    "The following transcripts were filtered out "
                    f"as potential nonsense-mediated decay targets:\n\t{log_msg}"
                ),
                "warning",
            )
            ## TODO: Replace `self.rejected_lines` with dict to enable fast ejection mechanism
            ## (see `check_isoforms()`)
            self.rejected_transcripts.update(nmd_targets)
            self.rejected_lines.extend(
                [RejectionReasons.NMD_REASON.format(x) for x in nmd_targets]
            )
        ## proceed further

    def report_fusions(self) -> None:
        """
        Finds isoforms fusing two or more otherwise disjoint genes:
        """

    def report_fissions(self) -> None:
        """
        Finds isoforms attributed to the same gene but not overlapping by absolute coordinates:
        """

    def check_isoforms(self) -> None:
        """
        Filters the reference isoform (gene-to-transcript) mapping file
        by removing genes whose transcripts (isoforms) were discarded
        at the annotation filter step.
        In parallel, all the transcripts recorded at the annotation filter step
        and missing gene mapping in the isoforms file are further excluded
        from the final annotation.
        """
        gene2trs: dict[str, list[str]] = {}
        trs_found: set[str] = set()
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
            if gene not in gene2trs:
                gene2trs[gene] = []
            if tr in self.rejected_transcripts:
                continue
            gene2trs[gene].append(tr)
            trs_found.add(tr)
        rejected_genes: list[str] = []

        ## write the remaining isoforms to the output file
        with (
            open(self.filtered_isoforms, "w") as iso,
            open(self.gene_bed, "w") as bed,
            open(self.immune_isoforms, "w") as i_iso,
            open(self.immune_gene_bed, "w") as i_bed,
        ):
            for gene, trs in gene2trs.items():
                if not trs:
                    rejected_genes.append(gene)
                    continue
                chrom: str | None = None
                start: int | None = None
                end: int | None = None
                strand: bool | None = None
                immune: bool = False
                for tr in trs:
                    if tr in self.immune:
                        i_iso.write(f"{gene}\t{tr}\n")
                        immune = True
                    else:
                        iso.write(f"{gene}\t{tr}\n")
                    _chrom, _start, _end, _strand = self.tr2coords[tr]
                    if chrom is None:
                        chrom = _chrom
                    elif chrom != _chrom:
                        self._die(
                            f"Conflicting chromosomes encountered for gene {gene}: {chrom} and {_chrom}"
                        )
                    if start is None:
                        start = _start
                    else:
                        start = min(start, _start)
                    if end is None:
                        end = _end
                    else:
                        end = max(end, _end)
                    if strand is None:
                        strand = _strand
                    else:
                        if strand != _strand:
                            self._die(f"Conflicting strand encountered for gene {gene}")
                gene_bed_line: str = GENE_BED_LINE.format(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=gene,
                    strand=("+" if strand else "-"),
                )
                if immune:
                    i_bed.write(gene_bed_line)
                else:
                    bed.write(gene_bed_line)

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
            self.rejected_lines.extend(
                [RejectionReasons.REJ_GENE.format(x) for x in rejected_genes]
            )
        ## report transcripts for which genes were not found in the isoform file
        rejected_transcripts: list[str] = [
            x for x in self.tr2annot if x not in trs_found
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
            self.rejected_lines.extend(
                [RejectionReasons.ORPHAN_TR.format(x) for x in rejected_transcripts]
            )
            self.tr2annot = {
                k: v for k, v in self.tr2annot.items() if k not in rejected_transcripts
            }

    def write_annotation(self) -> None:
        """Writes the filtered reference annotation to the file"""
        if not self.tr2annot:
            self._die("All transcripts were filtered out for various reasons")
        self._to_log(f"Writing the filtered annotation to {self.filtered_annotation}")
        with (
            open(self.filtered_annotation, "w") as h,
            open(self.immune_annotation, "w") as ih,
        ):
            for name, line in self.tr2annot.items():
                if name in self.immune:
                    ih.write(line + "\n")
                else:
                    h.write(line + "\n")

    def write_rejection_log(self) -> None:
        """
        Writes the rejected items to the file

        Arguments:
            None

        Returns:
            None
        """
        if not self.rejected_lines:
            return
        self._to_log(f"Writing the rejected items to {self.rejection_log}")
        with open(self.rejection_log, "w") as h:
            h.write(Headers.REJ_LOG_HEADER)
            h.writelines(self.rejected_lines)

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
            h.writelines([equi_line] * ACC_PROFILE_LEN)
            # for _ in range(ACC_PROFILE_LEN):
            #     h.write(equi_line)
        equi_donor: str = os.path.join(self.profile_dir, EQUI_DONOR)
        with open(equi_donor, "w") as h:
            h.write(header + "\n")
            h.writeline([equi_line] * DONOR_PROFILE_LEN)
            # for _ in range(DONOR_PROFILE_LEN):
            #     h.write(equi_line)

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

    def set_logging(self) -> None:
        """Sets up logging system for a InputProducer instance"""
        super().set_logging(toga_module="prepare-input")
        self.logger.propagate = False
