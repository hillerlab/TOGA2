#!/usr/bin/env python3

"""
Given the extended BED files produced at the TOGA alignment step, 
prepares a complete UCSC BigBed report file
"""

from collections import defaultdict
from ucsc_report import REF_LINK_PLACEHOLDER
from shared import CommandLineManager, CONTEXT_SETTINGS, parse_single_column
from typing import Any, Dict, List, Optional, Set, TextIO, Tuple, Union

import click
import os

__author__ = 'Yury V. Malovichko'
__credits__ = ('Bogdan Kirilenko', 'Björn Langer', 'Michael Hiller')
__year__ = '2024'

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
MAKE_IX_SCRIPT: str = os.path.join(LOCATION, "get_names_from_bed.py")
BR: str = '<BR>'
DEFAULT_UCSC_PREFIX: str = 'HLTOGAannot'
FRAGM_PROJ_MAIN_STUB: str = (
    'Chain {}: {} <I>(this locus)</I>'
)
INTACT_COLORS: Tuple[str, str] = ('0,0,100', '0,0,200')
# FRAGM_PROJ_ADD_STUB: str = (
#     'Chain {}: <a href="{}" target="_blank">{}</a><BR>'
# )
FRAGM_PROJ_ADD_STUB: str = (
    'Chain {}: {}'
)

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'stub',
    type=click.File('r', lazy=True),
    metavar='BED_STUB'
)
@click.argument(
    'ref_bed',
    type=click.File('r', lazy=True),
    metavar='REF_BED'
)
@click.argument(
    'projection_features',
    type=click.File('r', lazy=True),
    metavar='PROJ_FEATURES_FILE'
)
@click.argument(
    'orthology_scores',
    type=click.File('r', lazy=True),
    metavar='ORTHOLOGY_SCORES_FILE'
)
@click.argument(
    'chrom_sizes',
    type=click.Path(exists=True),
    metavar='CHROM_SIZES_TABLE'
)
@click.argument(
    'schema',
    type=click.Path(exists=True),
    metavar='SCHEMA_FILE'
)
@click.argument(
    'output_dir',
    type=click.Path(exists=False),
    metavar='OUTPUT_DIR'
)
@click.option(
    '--alternative_bed_track',
    '-a',
    type=click.File('r', lazy=True),
    default=None,
    show_default=True,
    help=(
        'A path to an alternative annotation track, in Bed12 format. Contents for columns 1-12 '
        'will be supplanted with those from this file'
    )
)
@click.option(
    '--link_file',
    '-i',
    type=click.File('r', lazy=True),
    default=None,
    show_default=True,
    help=(
        'A path to the two-column tab-separated file containing '
        'external links for reference transcripts'
    )
)
@click.option(
    '--deprecated_projections',
    '-d',
    type=click.File('r', lazy=True),
    default=None,
    show_default=None,
    help=(
        'A path to a single-column file containing projections '
        'to be removed from the final BigBed file'
    )   
)
@click.option(
    '--processed_pseudogenes',
    '-pp',
    type=click.File('r', lazy=True),
    default=None,
    show_default=None,
    help=(
        'A path to a single-column file containing names  '
        'of processed pseudogene projections; those with loss statuses lower than FI/I '
        'will be excluded from the final annotation'
    )   
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
    '--prefix',
    '-p',
    type=str,
    metavar='STR',
    default=DEFAULT_UCSC_PREFIX,
    show_default=True,
    help='A prefix to use in the output file names'
)

@click.option(
    '--log_file',
    '-l',
    type=click.Path(exists=False),
    default=None,
    show_default=True,
    help='A path to log file to write progress log in'
)
@click.option(
    '--log_name',
    '-ln',
    type=str,
    metavar='STR',
    default=None,
    show_default=True,
    help='Logger name to use; relevant only upon main class import'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
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
        'prefix', 'proj2backbone', 'ref2coords', 'proj2features', 'proj2score', 
        'ref2link', 'alt_bed_track', 'chrom_sizes', 'schema', 'output_dir', 
        'out_bed_file', 'bigbed_file', 
        'bed_index', 'bigbed_index', 'bigbed_ixx', 'longest_word',
        # 'make_pp_track', '', '', '',
        'deprecated_projs', 'ppgenes',
        'log_file'
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
        alternative_bed_track: Optional[Union[click.File, None]],
        link_file: Optional[Union[click.File, None]],
        deprecated_projections: Optional[Union[click.File, None]],
        processed_pseudogenes: Optional[Union[click.File, None]],
        # make_processed_pseudogene_track: Optional[bool],
        prefix: Optional[str],
        # pseudogene_prefix: Optional[str],
        log_file: Optional[click.Path],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.log_file: click.Path = log_file
        self.set_logging(log_name)

        # self.prefix: str = prefix
        self.longest_word: int = 0
        self.proj2backbone: Dict[str, Dict[str, List[str]]] = defaultdict(dict)
        self._to_log('Parsing UCSC BED file stub')
        self.parse_backbone_bed(stub)
        self.ref2coords: Dict[str, str] = {}
        self._to_log('Extracting reference coordinates from the reference BED file')
        self.parse_ref_bed(ref_bed)
        self.proj2features: Dict[str, Tuple[float]] = {}
        self._to_log('Extracting the necessary projection features')
        self.parse_projection_features(projection_features)
        self.proj2score: Dict[str, str] = {}
        self._to_log('Extracting orthology probabilities')
        self.parse_scores(orthology_scores)
        self.ref2link: Dict[str, str] = {}
        if link_file is not None:
            self._to_log('Parsing reference link file')
            self.parse_links(link_file)
        if deprecated_projections is not None:
            self._to_log('Parsing the list of deprecated projections')
        self.alt_bed_track: Dict[str, List[str]] = {}
        if alternative_bed_track is not None:
            self._to_log('Parsing alternative Bed track')
            self.parse_alt_bed(alternative_bed_track)
        self.deprecated_projs: Set[str] = parse_single_column(deprecated_projections)
        self.ppgenes: Set[str] = parse_single_column(processed_pseudogenes)
        self.chrom_sizes: click.Path = chrom_sizes
        self.schema: click.Path = schema
        self.output_dir: click.Path = output_dir
        # self.make_pp_track: bool = make_processed_pseudogene_track
        self.out_bed_file: str = os.path.join(self.output_dir, 'query_annotation.for_ucsc.bed34')
        self.bigbed_file: str = os.path.join(self.output_dir, f'{prefix}.bb')
        self.bed_index: str = os.path.join(self.output_dir, 'query_annotation.ix.txt')
        self.bigbed_index: str = os.path.join(self.output_dir, f'{prefix}.ix')
        self.bigbed_ixx: str = os.path.join(self.output_dir, f'{prefix}.ixx')
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
        self._to_log('Creating UCSC output directory')
        self._mkdir(self.output_dir)
        sorted_projs: List[Tuple[str, str, List[str]]] = sorted(
            [y for x in self.proj2backbone.values() for y in x.items()], 
            key=lambda x: (x[1][0], int(x[1][1]))
        )
        self._to_log('Preparing UCSC input BED12+ file')
        with open(self.out_bed_file, 'w') as h:
            for chain, backbone in sorted_projs:
                proj: str = backbone[3]
                if self.deprecated_projs and proj in self.deprecated_projs:
                    continue
                if proj in self.ppgenes and backbone[8] not in INTACT_COLORS:
                    continue
                if ',' not in proj and proj in self.alt_bed_track:
                    backbone[:12] = self.alt_bed_track[proj]
                tr: str = '#'.join(proj.split('#')[:-1])
                ref_coords: str = self.ref2coords.get(tr, '')
                if not ref_coords:
                    self._die('Transcript %s is missing from the reference BED file' % tr)
                if ',' in proj:
                    query_coords: List[str] = [
                        FRAGM_PROJ_MAIN_STUB.format(
                            chain, f'{backbone[0]}:{backbone[1]}-{backbone[2]}'
                        )
                    ]
                    for other_chain in self.proj2backbone[proj]:
                        if other_chain == chain:
                            continue
                        _chrom, _start, _end =  self.proj2backbone[proj][other_chain][:3]
                        other_locus_coords: str = f'{_chrom}:{_start}-{_end}'
                        query_coords.append(
                            FRAGM_PROJ_ADD_STUB.format(other_chain, other_locus_coords)
                        )
                    query_coords: str = BR.join(query_coords)
                    _proj: str = f'{tr}#{chain}'
                    orth_score: str = self.proj2score.get(_proj, '')
                    proj_features: List[str] = self.proj2features.get(_proj, [])
                else:
                    orth_score: str = self.proj2score.get(proj, '')
                    proj_features: List[str] = self.proj2features.get(proj, [])
                    query_coords: str = f'{backbone[0]}:{backbone[1]}-{backbone[2]}'
                if not orth_score:
                    self._die(
                        'Projection %s is missing from the orthology score file' % proj
                    )
                if not proj_features:
                    self._die(
                        'Projection %s is missing from the projection feature file' % proj
                    )
                link: str = self.ref2link.get(tr, REF_LINK_PLACEHOLDER)
                if proj in self.ppgenes:
                    backbone[3] = backbone[3] + "#retro"
                out_line: List[Any] = [
                    *backbone[:4], '1000', *backbone[5:12], tr, 
                    ref_coords, query_coords, orth_score, *proj_features,
                    *backbone[12:21], link, *backbone[21:]#*backbone[21:24], backbone[23], *backbone[24:] 
                ]
                h.write('\t'.join(map(str, out_line)) + '\n')
        self._to_log('UCSC Bed input prepared; creating BigBed file')
        if self.alt_bed_track:
            sort_cmd: str = f'sort -k1,1 -k2,2n -o {self.out_bed_file} {self.out_bed_file}'
            _ = self._exec(sort_cmd, 'Browser report Bed file sorting failed')
        bigbed_cmd: str = (
            f'bedToBigBed -type=bed12+26 {self.out_bed_file} {self.chrom_sizes} {self.bigbed_file} '
            f'-tab -extraIndex=name -as={self.schema}'
        )
        _ = self._exec(bigbed_cmd, 'bedToBigBed failed')
        self._to_log('BigBed file successfully created')
        self._to_log('Creating BigBed indices')
        bed_ix_cmd: str = f'{MAKE_IX_SCRIPT} {self.out_bed_file} | sort -u > {self.bed_index}'
        _ = self._exec(bed_ix_cmd, 'BED file indexing failed')
        bigbed_ix_cmd: str = (
            f'ixIxx {self.bed_index} {self.bigbed_index} {self.bigbed_ixx} '
            f'-maxWordLength={self.longest_word}'
        )
        self._exec(bigbed_ix_cmd, 'BigBED file indexing failed')
        self._to_log('UCSC browser input preparation finished')
        

    def parse_backbone_bed(self, file: TextIO) -> None:
        """Parses extended BED output from the CESAR alignment step"""
        for i, line in enumerate(file, 1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 26:
                self._die(
                    'Line %i in extended BED file does not contain expected number of columns (24)' % i
                )
            proj: str = data[3]
            chain: str = data[4]
            self.proj2backbone[proj][chain] = data
            self.longest_word = max(self.longest_word, len(proj))

    def parse_ref_bed(self, file: TextIO) -> None:
        """Extracts reference transcript coordinates from the BED file"""
        for i, line in enumerate(file, 1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 4:
                self._die(
                    'Line %i in the reference BED file does not contain expected number of columns (4)' %i
                )
            chrom: str = data[0]
            if len(data) >= 7:
                start: int = data[6]
                end: int = data[7]
            else:
                start: int = data[1]
                end: int = data[2]
            name: str = data[3]
            self.ref2coords[name] = f'{chrom}:{start}-{end}'

    def parse_alt_bed(self, file: Union[TextIO, None]) -> None:
        """If provided, parses an alternative Bed12 track file"""
        if file is None:
            return
        for i, line in enumerate(file, 1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) != 12:
                self._die(
                    (
                        'Line %i in the alternative Bed12 file is improperly formatted; '
                        'expected 12 lines, got %i'
                    ) % (i, len(data))
                )
            name: str = data[3].replace('#retro', '')
            self.alt_bed_track[name] = data

    def parse_projection_features(self, file: TextIO) -> None:
        """Extracts selected projection features from the feature table"""
        for line in file:
            data: List[str] = line.rstrip().split('\t')
            if len(data) < 16:
                self._die(
                    'Classification feature file differs from the expected format'
                )
            if data[0] == 'transcript':
                continue
            trans: str = data[0]
            chain: str = data[2]
            proj: str= f'{trans}#{chain}'
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
                synt_, flank_cov_, gl_exo_, loc_exon_, exon_cov, intron_cov
            )

    def parse_scores(self, file: TextIO) -> None:
        """Extracts projection orthology probabilities"""
        for i, line in enumerate(file, 1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 3:
                self._die(
                    'Line %i in the orthology score file does not contain expected number of columns (3)' % i
                )
            tr: str = data[0]
            if tr == 'transcript':
                continue
            chain: int = data[1]
            proj: str = f'{tr}#{chain}'
            prob: str = data[2]
            try:
                float(prob)
            except ValueError:
                self._die('Column 3 at line %i does not contain a valid probability value' %i)
            self.proj2score[proj] = prob

    def parse_links(self, file: TextIO) -> None:
        """Extracts links to external sources for reference transcripts"""
        for i, line in enumerate(file, 1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 2:
                self._die(
                    'Line %i in the reference link file does not contain expected number of columns (2)' % i
                )
            name: str = data[0]
            link: str = data[1]
            self.ref2link[name] = link


if __name__ == '__main__':
    BigBedProducer()