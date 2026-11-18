"""
Given the extended BED files produced at the TOGA alignment step,
prepares a complete UCSC BigBed report file
"""

import os
from collections import defaultdict
from typing import Any, TextIO

import click

from .get_names_from_bed import BedNameRetriever
from .shared import (
    CONTEXT_SETTINGS,
    CommandLineManager,
    base_proj_name,
    get_proj2trans,
    parse_single_column,
)
from .ucsc_report import REF_LINK_PLACEHOLDER

__author__ = "Yury V. Malovichko"
__credits__ = ("Bogdan Kirilenko", "Björn Langer", "Michael Hiller")
__year__ = "2024"

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
BR: str = "<BR>"
DEFAULT_UCSC_PREFIX: str = "HLTOGAannot"
FRAGM_PROJ_MAIN_STUB: str = "Chain {}: {} <I>(this locus)</I>"
INTACT_COLORS: tuple[str, str] = ("0,0,100", "0,0,200")
# FRAGM_PROJ_ADD_STUB: str = (
#     'Chain {}: <a href="{}" target="_blank">{}</a><BR>'
# )
FRAGM_PROJ_ADD_STUB: str = "Chain {}: {}"


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument("stub", type=click.File("r", lazy=True), metavar="BED_STUB")
@click.argument("ref_bed", type=click.File("r", lazy=True), metavar="REF_BED")
@click.argument(
    "projection_features", type=click.File("r", lazy=True), metavar="PROJ_FEATURES_FILE"
)
@click.argument(
    "orthology_scores", type=click.File("r", lazy=True), metavar="ORTHOLOGY_SCORES_FILE"
)
@click.argument(
    "chrom_sizes", type=click.Path(exists=True), metavar="CHROM_SIZES_TABLE"
)
@click.argument("schema", type=click.Path(exists=True), metavar="SCHEMA_FILE")
@click.argument("output_dir", type=click.Path(exists=False), metavar="OUTPUT_DIR")
@click.option(
    "--alternative_bed_track",
    "-a",
    type=click.File("r", lazy=True),
    default=None,
    show_default=True,
    help=(
        "A path to an alternative annotation track, in Bed12 format. Contents for columns 1-12 "
        "will be supplanted with those from this file"
    ),
)
@click.option(
    "--link_file",
    "-i",
    type=click.File("r", lazy=True),
    default=None,
    show_default=True,
    help=(
        "A path to the two-column tab-separated file containing "
        "external links for reference transcripts"
    ),
)
@click.option(
    "--deprecated_projections",
    "-d",
    type=click.File("r", lazy=True),
    default=None,
    show_default=None,
    help=(
        "A path to a single-column file containing projections "
        "to be removed from the final BigBed file"
    ),
)
@click.option(
    "--paralogs",
    "-p",
    type=click.File("r", lazy=True),
    metavar="PARALOG_LIST_FILE",
    default=None,
    show_default=True,
    help=(
        "A path to a single-column file containing names "
        "of annotated paralogous projections; those will get "
        'the "#paralog" suffix in the BigBed file'
    ),
)
@click.option(
    "--processed_pseudogenes",
    "-pp",
    type=click.File("r", lazy=True),
    metavar="PPGENE_LIST_FILE",
    default=None,
    show_default=None,
    help=(
        "A path to a single-column file containing names  "
        "of processed pseudogene projections; those with loss statuses of FI/I "
        'will get the "#retro" suffix in the BigBed file, '
        "with the rest being excluded from the final annotation"
    ),
)
# @click.option(
#     '--make_processed_pseudogene_track',
#     '-pptrack',
#     if_flag=True,
#     default=False,
#     show_default=True,
#     help=(
#         'If set, creates a separate set of UCSC files for processed pseudogenes'
#     )
# )
@click.option(
    "--bedtobigbed_binary",
    type=click.Path(exists=True),
    metavar="BEDTOBIGBED_PATH",
    default=None,
    show_default=True,
    help=("A path to UCSC bedToBigBed executable. "),
)
@click.option(
    "--ixixx_binary",
    type=click.Path(exists=True),
    metavar="IXIXX_PATH",
    default=None,
    show_default=True,
    help=("A path to UCSC ixIxx executable. "),
)
@click.option(
    "--prefix",
    "-p",
    type=str,
    metavar="STR",
    default=DEFAULT_UCSC_PREFIX,
    show_default=True,
    help="A prefix to use in the output file names",
)
@click.option(
    "--log_file",
    "-l",
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help="A path to log file to write progress log in",
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
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls execution verbosity",
)
class BigBedProducer(CommandLineManager):
    """
    Given the extended BED files produced at the TOGA alignment step,
    prepares a complete UCSC BigBed report file.\n
    Arguments are:\n
    * BED_STUB \n
    * REF_BED \n
    * PROJ_FEATURES_FILE \n
    * ORTHOLOGY_SCORES_FILE \n
    * SCHEMA_FILE\n
    """

    __slots__ = (
        "alt_bed_track",
        "bed_index",
        "bedtobigbed_binary",
        "bigbed_file",
        "bigbed_index",
        "bigbed_ixx",
        "chrom_sizes",
        "deprecated_projs",
        "ixixx_binary",
        "log_file",
        "log_name",
        "longest_word",
        "out_bed_file",
        "output_dir",
        "paralogs",
        "ppgenes",
        "prefix",
        "proj2backbone",
        "proj2features",
        "proj2score",
        "ref2coords",
        "ref2link",
        "schema",
    )

    def __init__(
        self,
        stub: click.File,
        ref_bed: click.File,
        projection_features: click.File,
        orthology_scores: click.File,
        chrom_sizes: click.Path,
        schema: click.Path,
        output_dir: click.Path,
        alternative_bed_track: click.File | None,
        link_file: click.File | None,
        deprecated_projections: click.File | None,
        paralogs: click.File | None,
        processed_pseudogenes: click.File | None,
        bedtobigbed_binary: click.Path | None,
        ixixx_binary: click.Path | None,
        prefix: str | None,
        log_file: click.Path | None,
        log_name: str | None,
        verbose: bool = False,
    ) -> None:
        self.v: bool = verbose
        self.log_file: click.Path = log_file
        self.log_name: str = log_name
        self.set_logging(name=self.log_name, toga_module="ucsc_track")

        # self.prefix: str = prefix
        self.longest_word: int = 0
        self.proj2backbone: dict[str, dict[str, list[str]]] = defaultdict(dict)
        self._to_log("Parsing UCSC BED file stub")
        self.parse_backbone_bed(stub)
        self.ref2coords: dict[str, str] = {}
        self._to_log("Extracting reference coordinates from the reference BED file")
        self.parse_ref_bed(ref_bed)
        self.proj2features: dict[str, tuple[float]] = {}
        self._to_log("Extracting the necessary projection features")
        self.parse_projection_features(projection_features)
        self.proj2score: dict[str, str] = {}
        self._to_log("Extracting orthology probabilities")
        self.parse_scores(orthology_scores)
        self.ref2link: dict[str, str] = {}
        if link_file is not None:
            self._to_log("Parsing reference link file")
            self.parse_links(link_file)
        if deprecated_projections is not None:
            self._to_log("Parsing the list of deprecated projections")
        self.alt_bed_track: dict[str, list[str]] = {}
        if alternative_bed_track is not None:
            self._to_log("Parsing alternative Bed track")
            self.parse_alt_bed(alternative_bed_track)
        self.deprecated_projs: set[str] = parse_single_column(deprecated_projections)
        self.paralogs: set[str] = parse_single_column(paralogs)
        self.ppgenes: set[str] = parse_single_column(processed_pseudogenes)
        self.chrom_sizes: click.Path = chrom_sizes
        self.schema: click.Path = schema
        self.output_dir: click.Path = output_dir
        self.bedtobigbed_binary: click.Path = bedtobigbed_binary
        self.ixixx_binary: click.Path = ixixx_binary
        # self.make_pp_track: bool = make_processed_pseudogene_track
        self.out_bed_file: str = os.path.join(
            self.output_dir, "query_annotation.for_ucsc.bed34"
        )
        self.bigbed_file: str = os.path.join(self.output_dir, f"{prefix}.bb")
        self.bed_index: str = os.path.join(self.output_dir, "query_annotation.ix.txt")
        self.bigbed_index: str = os.path.join(self.output_dir, f"{prefix}.ix")
        self.bigbed_ixx: str = os.path.join(self.output_dir, f"{prefix}.ixx")
        # self.pp_bed_file: str = os.path.join(
        #     self.output_dir, 'query_annotation.processed_pseudogenes.bed34'
        # )
        # self.pp_bigbed_file: str = os.path.join(self.output_dir, f'{pseudogene_prefix}.bb')
        # self.pp_ix: str = os.path.join(self.output_dir, f'{pseudogene_prefix.}')

        self.run()

        ## fields 1-12 are original BED12 fields
        ## fields 13 and 14 are extracted from the reference annotation file
        ## field 15 is basically fields 1:2-3(1:7-8) of the original BED file merged
        ## field 16 is the orthology probability score
        ## fields 17-22 come from the projection feature file (columns 21 and 22 need additional shenanigans)
        ## fields 23-29 come from transcript meta (should be added to the extended BED)
        ## field 30 contains predicted protein sequence (should be added to the extended BED)
        ## field 31 contains SVG plot
        ## field 32 contains reference link
        ## field 33 contains mutation table
        ## field 34 contains exon alignment table
        ## field 35 contains nucleotide CDS

    def run(self) -> None:
        """Main method"""
        self._to_log("Creating UCSC output directory")
        self._mkdir(self.output_dir)
        sorted_projs: list[tuple[str, str, list[str]]] = sorted(
            [y for x in self.proj2backbone.values() for y in x.items()],
            key=lambda x: (x[1][0], int(x[1][1])),
        )
        self._to_log("Preparing UCSC input BED12+ file")
        with open(self.out_bed_file, "w") as h:
            for chain, backbone in sorted_projs:
                proj: str = backbone[3]
                basename: str = base_proj_name(proj)
                if self.deprecated_projs and basename in self.deprecated_projs:
                    continue
                if basename in self.ppgenes and backbone[8] not in INTACT_COLORS:
                    continue
                if "," not in proj and proj in self.alt_bed_track:
                    backbone[:12] = self.alt_bed_track[proj]
                # tr: str = '#'.join(proj.split('#')[:-1])
                tr: str = get_proj2trans(proj)[0]
                ref_coords: str = self.ref2coords.get(tr, "")
                if not ref_coords:
                    self._die(f"Transcript {tr} is missing from the reference BED file")
                if "," in proj:
                    query_coords: list[str] = [
                        FRAGM_PROJ_MAIN_STUB.format(
                            chain, f"{backbone[0]}:{backbone[1]}-{backbone[2]}"
                        )
                    ]
                    for other_chain in self.proj2backbone[proj]:
                        if other_chain == chain:
                            continue
                        _chrom, _start, _end = self.proj2backbone[proj][other_chain][:3]
                        other_locus_coords: str = f"{_chrom}:{_start}-{_end}"
                        query_coords.append(
                            FRAGM_PROJ_ADD_STUB.format(other_chain, other_locus_coords)
                        )
                    query_coords: str = BR.join(query_coords)
                    _proj: str = f"{tr}#{chain}"
                    orth_score: str = self.proj2score.get(_proj, "")
                    proj_features: list[str] = self.proj2features.get(_proj, [])
                else:
                    orth_score: str = self.proj2score.get(basename, "")
                    proj_features: list[str] = self.proj2features.get(basename, [])
                    query_coords: str = f"{backbone[0]}:{backbone[1]}-{backbone[2]}"
                if not orth_score:
                    self._die(f"Projection {proj} is missing from the orthology score file")
                if not proj_features:
                    self._die(f"Projection {proj} is missing from the projection feature file")
                link: str = self.ref2link.get(tr, REF_LINK_PLACEHOLDER)
                if basename in self.ppgenes and "#retro" not in proj:
                    backbone[3] = backbone[3] + "#retro"
                elif basename in self.paralogs and "#paralog" not in proj:
                    backbone[3] = backbone[3] + "#paralog"
                out_line: list[Any] = [
                    *backbone[:4],
                    "1000",
                    *backbone[5:12],
                    tr,
                    ref_coords,
                    query_coords,
                    orth_score,
                    *proj_features,
                    *backbone[12:21],
                    link,
                    *backbone[21:],  # *backbone[21:24], backbone[23], *backbone[24:]
                ]
                h.write("\t".join(map(str, out_line)) + "\n")
        self._to_log("UCSC Bed input prepared; creating BigBed file")
        if self.alt_bed_track:
            sort_cmd: str = (
                f"sort -k1,1 -k2,2n -o {self.out_bed_file} {self.out_bed_file}"
            )
            _ = self._exec(sort_cmd, "Browser report Bed file sorting failed")
        bigbed_cmd: str = (
            f"{self.bedtobigbed_binary} -type=bed12+26 {self.out_bed_file} {self.chrom_sizes} {self.bigbed_file} "
            f"-tab -extraIndex=name -as={self.schema}"
        )
        _ = self._exec(bigbed_cmd, "bedToBigBed failed")
        self._to_log("BigBed file successfully created")
        self._to_log("Retrieving names from the BED file")
        BedNameRetriever(
            (
                "--input", 
                self.out_bed_file, 
                "--output", 
                self.bed_index, 
                "--log_name", 
                self.log_name
            ),
            standalone_mode=False
        )
        self._to_log("Creating BigBed indices")
        bigbed_ix_cmd: str = (
            f"{self.ixixx_binary} {self.bed_index} {self.bigbed_index} {self.bigbed_ixx} "
            f"-maxWordLength={self.longest_word}"
        )
        self._exec(bigbed_ix_cmd, "BigBED file indexing failed")
        self._to_log("UCSC browser input preparation finished")

    def parse_backbone_bed(self, file: TextIO) -> None:
        """Parses extended BED output from the CESAR alignment step"""
        for i, line in enumerate(file, 1):
            data: list[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 26:
                self._die(
                    f"Line {i} in extended BED file does not contain expected number of columns (24)"
                )
            proj: str = data[3]
            chain: str = data[4]
            self.proj2backbone[proj][chain] = data
            self.longest_word = max(self.longest_word, len(proj))

    def parse_ref_bed(self, file: TextIO) -> None:
        """Extracts reference transcript coordinates from the BED file"""
        for i, line in enumerate(file, 1):
            data: list[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 4:
                self._die(
                    f"Line {i} in the reference BED file does not contain expected number of columns (4)"
                )
            chrom: str = data[0]
            if len(data) >= 7:
                start: int = data[6]
                end: int = data[7]
            else:
                start: int = data[1]
                end: int = data[2]
            name: str = data[3]
            self.ref2coords[name] = f"{chrom}:{start}-{end}"

    def parse_alt_bed(self, file: TextIO | None) -> None:
        """If provided, parses an alternative Bed12 track file"""
        if file is None:
            return
        for i, line in enumerate(file, 1):
            data: list[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 12:
                self._die(
                    f"Line {i} in the alternative Bed12 file is improperly formatted; "
                    f"expected 12 lines, got {len(data)}"
                )
            name: str = base_proj_name(data[3])
            self.alt_bed_track[name] = data

    def parse_projection_features(self, file: TextIO) -> None:
        """Extracts selected projection features from the feature table"""
        for line in file:
            data: list[str] = line.rstrip().split("\t")
            if len(data) < 16:
                self._die(
                    "Classification feature file differs from the expected format"
                )
            if data[0] == "transcript":
                continue
            trans: str = data[0]
            chain: str = data[2]
            proj: str = f"{trans}#{chain}"
            # if proj_list and proj not in proj_list:
            #     continue
            synt_: str = data[3]
            gl_exo_: str = data[5]
            loc_exon_: str = data[8]
            exon_cover_: str = data[9]
            intron_cover_: str = data[10]
            exon_fract_: str = data[13]
            intron_fract_: str = data[14]
            flank_cov_: str = data[15]

            exon_cov: str = (
                str(float(exon_cover_) / float(exon_fract_))
                if float(exon_fract_) != 0
                else "0"
            )
            intron_cov: str = (
                str(float(intron_cover_) / float(intron_fract_))
                if float(intron_fract_) != 0
                else "0"
            )
            self.proj2features[proj] = (
                synt_,
                flank_cov_,
                gl_exo_,
                loc_exon_,
                exon_cov,
                intron_cov,
            )

    def parse_scores(self, file: TextIO) -> None:
        """Extracts projection orthology probabilities"""
        for i, line in enumerate(file, 1):
            data: list[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 3:
                self._die(
                    f"Line {i} in the orthology score file does not contain expected number of columns (3)"
                )
            tr: str = data[0]
            if tr == "transcript":
                continue
            chain: int = data[1]
            proj: str = f"{tr}#{chain}"
            prob: str = data[2]
            try:
                float(prob)
            except ValueError:
                self._die(
                    f"Column 3 at line {i} does not contain a valid probability value"
                )
            self.proj2score[proj] = prob

    def parse_links(self, file: TextIO) -> None:
        """Extracts links to external sources for reference transcripts"""
        for i, line in enumerate(file, 1):
            data: list[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 2:
                self._die(
                    f"Line {i} in the reference link file does not contain expected number of columns (2)"
                )
            name: str = data[0]
            link: str = data[1]
            self.ref2link[name] = link


if __name__ == "__main__":
    BigBedProducer()
