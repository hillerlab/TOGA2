#!/usr/bin/env python3

"""
Infers proxies for query genes based on coding sequence intersection between
the projections
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from collections import defaultdict
from dataclasses import dataclass
from modules.cesar_wrapper_executables import AnnotationEntry, Exon, ExonDict
from modules.shared import CommandLineManager, CONTEXT_SETTINGS, intersection
from typing import Dict, List, Optional, Set, TextIO, Tuple, Union

import click
import networkx as nx

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = ('Bogdan V. Kirilenko')
__all__ = (None)

## check the NetworkX version
nx_v: str = nx.__version__
v_split: List[str] = [x for x in nx_v.split(".") if x.isnumeric()]
if len(v_split) > 1:
    NX_VERSION: float = float(f"{v_split[0]}.{v_split[1]}")
else:
    NX_VERSION: float = float(v_split[0])
PROJECTION: str = 'PROJECTION'
MISSING: Tuple[str, str] = ('M', 'N')
EXTENDED_HIGH_CONFIDENCE: Tuple[str, str] = ('FI', 'PI')
MIN_RELIABLE_EXON_COV: float = 0.6

REJ_ORTH_REASON: str = '\t'.join(
    (
        'PROJECTION', '{}', '0', 
        'Insufficiently covered exons in second-best projection', 
        'SECOND_BEST', '{}'
    )
)
REJ_PARA_REASON: str = '\t'.join(
    (
        'PROJECTION', '{}', '0', 
        'Redundant paralog overlapping orthologous projections', 
        'REDUNDANT_PARALOG', '{}'
    )
)
REJ_PPGENE_REASON: str = '\t'.join(
    (
        'PROJECTION', '{}', '0',
        'Processed pseudogene overlapping ortholog or paralog',
        'REDUNDANT_PPGENE', '{}'
    )
)

def parse_single_column(file: Union[TextIO, None]) -> Set[str]:
    """Parses a file as a newline-separated list of strings"""
    result: Set[str] = set()
    if file is None:
        return result
    for line in file:
        line = line.rstrip()
        if not line:
            continue
        result.add(line)
    return result


def get_orig_transcript(proj: str) -> str:
    """Extracts reference transcript name from projection identifier"""
    return '#'.join(proj.split('#')[:-1])


@dataclass
class Coords:
    __slots__ = ('name', 'chrom', 'start', 'end', 'strand')
    name: str
    chrom: str
    start: int
    end: int
    strand: bool


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
# @click.argument(
#     'query_bed',
#     type=click.File('r', lazy=True),
#     metavar='QUERY_BED'
# )
@click.argument(
    'query_exon_meta',
    type=click.File('r', lazy=True),
    metavar='QUERY_EXON_META'
)
@click.argument(
    'output',
    type=click.File('w', lazy=True),
    metavar='OUTPUT'
)
@click.option(
    '--bed6_output',
    '-b',
    type=click.File('w', lazy=True),
    metavar='OUTPUT_BED6',
    default=None,
    show_default=True,
    help='A path to the BED6 file containing query locus coordinates'
)
@click.option(
    '--ref_isoform_file',
    '-i',
    type=click.File('r', lazy=True),
    metavar='ISOFORM_FILE',
    default=None,
    show_default=True,
    help=(
        'A tab-separated two-column file containing reference gene-to-transcript mapping; '
        'if provided, projections from overlapping reference genes are not merged into one gene. '
        'WARNING: Valid only if --ref_bed is also provided'
    )
)
@click.option(
    '--ref_bed',
    '-r',
    type=click.File('r', lazy=True),
    metavar='BED_FILE',
    default=None,
    show_default=True,
    help=(
        'A reference annotation BED file containing at least 8 columns; '
        'if provided, projections from overlapping reference genes are not merged into one gene. '
        'WARNING: Valid only if --ref_isoform_file is also provided'
    )
)
@click.option(
    '--loss_summary_file',
    '-l',
    type=click.File('r', lazy=True),
    metavar='LOSS_SUMMARY_FILE',
    default=None,
    show_default=True,
    help=(
        'TOGA2 loss summary file containing at least three columns; '
        'if provided, missing transcripts will be excluded from gene inference'
    )
)
@click.option( ## TODO: Move to query gene inference code
    '--feature_file',
    '-pf',
    type=click.File('r', lazy=True),
    metavar='PROJECTION_FEATURES',
    default=None,
    show_default=True,
    help=(
        'Projection feature file produced by TOGA2 for projection classificaiton. '
        'Certain features can be used by orthology resolver module to filter out '
        'confounding orthology predictions.'
    )
)
@click.option(
    '--orthology_probabilities',
    '-op',
    type=click.File('r', lazy=True),
    metavar='ORTHOLOGY_PROB_FILE',
    default=None,
    show_default=True,
    help=(
        'Orthology probability table produced by TOGA2 at the projection classification step. '
        'If provided together with the feature file, probabilities '
        'will be used for overextended projection fitlering.'
    )
)
@click.option(
    '--orthology_threshold',
    '-ot',
    type=float,
    metavar='FLOAT',
    default=0.5,
    show_default=True,
    help='Probability threshold for considering projections as orthologous'
)
@click.option(
    '--paralog_list',
    '-p',
    type=click.File('r', lazy=True),
    metavar='PARALOG_LIST',
    default=None,
    show_default=True,
    help=(
        'A single-column file containing the list of paralogous projections included '
        'at the alignment step'
    )
)
@click.option(
    '--processed_pseudogene_list',
    '-pp',
    type=click.File('r', lazy=True),
    metavar='PROC_PSEUDOGENE_LIST',
    default=None,
    show_default=True,
    help=(
        'A single-column file containing the list of processed pseudogene projections included '
        'at the alignment step'
    )
)
@click.option(
    '--redundant_paralogs',
    '-d',
    type=click.File('w', lazy=True),
    metavar='DISCARDED_PARALOG_FILE',
    default=None,
    show_default=True,
    help=(
        'A path to write the names of "redundant" paralogous projections. '
        'A paralog is considered redundant if it corresponds to a query locus '
        'harbouring orthologous projections'
    )
)
@click.option(
    '--redundant_processed_pseudogenes',
    '-dpp',
    type=click.File('w', lazy=True),
    metavar='DISCARDED_PPGENE_FILE',
    default=None,
    show_default=True,
    help=(
        'A path to write the names of "redundant" retrogene/processed pseudogene projections. '
        'A processed pseudogene is considered redundant if it corresponds to a query locus '
        'harbouring orthologous or paralogous projections'
    )
)
@click.option(
    '--insufficiently_covered_orthologs',
    type=click.File('w', lazy=True),
    metavar='DISCARDED_ORTHOLOG_FILE',
    default=None,
    show_default=True,
    help=(
        'A path to write the names of orthologous projections '
        'with insufficient initial chain coverage which were discarded due to '
        'locus occupance by more trustworthy projections'
    )
)
@click.option(
    '--rejection_log',
    '-rl',
    type=click.File('w', lazy=True),
    default=None,
    show_default=True,
    help=(
        'A path to to write the discarded entries to'
    )
)
@click.option(
    '--log_file',
    '-log',
    type=click.Path(exists=False),
    metavar='LOG_FILE',
    default=None,
    show_default=True,
    help='A path to write execution log to'
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
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)

class QueryGeneCollapser(CommandLineManager):
    __slots__ = [
        'query_transcripts', 'output', 'bed6_output', 'tr2chrom', 'tr2exons',
        'proj2exon_cov', 'proj2status', 'proj2prob', 'tr2max_prob',
        'lost_projections', 'sorted_components', 'component_coords', 'component2trs',
        'ref_isoform2gene', 'intersecting_ref_genes', 'log_file', 'orthology_threshold',
        'paralog_list', 'proc_pseudogene_list', 
        'discarded_paralogs', 'discarded_paralogs_file',
        'discarded_ppgenes', 'discarded_ppgenes_file',
        'discarded_extensions_file', 'discarded_overextensions',
        'rejected_items_file'
    ]

    """
    Infers proxies for query genes based on coding sequence intersection between
    the projections.\nArguments are:\n
    * QUERY_BED - a BED12 output of the TOGA's CESAR step;\n
    * OUTPUT - a destination for the two-column file containing query region-to-
    projection mapping
    """

    def __init__(
        self,
        # query_bed: click.File,
        query_exon_meta: click.File,
        output: click.File,
        ref_isoform_file: Optional[Union[click.File, None]],
        bed6_output: Optional[Union[click.File, None]],
        ref_bed: Optional[Union[click.File, None]],
        loss_summary_file: Optional[Union[click.File, None]],
        feature_file: Optional[Union[click.File, None]],
        orthology_probabilities: Optional[Union[click.File, None]],
        orthology_threshold: Optional[float],
        paralog_list: Optional[Union[click.File, None]],
        processed_pseudogene_list: Optional[Union[click.File, None]],
        redundant_paralogs: Optional[Union[click.File, None]],
        redundant_processed_pseudogenes: Optional[Union[click.File, None]],
        insufficiently_covered_orthologs: Optional[Union[click.File, None]],
        rejection_log: Optional[Union[click.File, None]],
        log_file: Optional[Union[click.Path, None]],
        log_name: Optional[str],
        verbose: bool
    ) -> None:
        self.v: bool = verbose
        self.log_file: str = log_file
        self.set_logging(log_name)

        self._to_log('Initialising QueryGeneCollapser')
        if (feature_file is None) != (orthology_probabilities is None):
            self._die(
                'Options --feature_file and --orthology_probabilities cannot be used separately; '
                'please provide both files or remove the provided option'
            )
        # self.query_transcripts: Dict[str, List[AnnotationEntry]] = defaultdict(dict)
        self.query_transcripts: Dict[str, List[Coords]] = defaultdict(dict)
        self.output: click.File = output
        self.bed6_output: click.File = bed6_output
        self.tr2chrom: Dict[str, Tuple[str]] = defaultdict(tuple)
        self.tr2exons: Dict[str, Dict[int, Coords]] = defaultdict(dict)
        self.sorted_components: List[int] = []
        self.component_coords: Dict[int, List[Tuple[str, int, bool]]] = defaultdict(list)
        self.component2trs: Dict[int, List[str]] = defaultdict(list)
        # self._to_log('Parsing query projections BED file')
        # self.parse_bed(query_bed)
        self._to_log('Parsing query exon metadata file')
        self.discarded_overextensions: Set[str] = set()
        self.parse_exon_meta(query_exon_meta)
        self.intersecting_ref_genes: Dict[str, Set[str]] = defaultdict(set)
        self.ref_isoform2gene: Dict[str, str] = {}
        if ref_isoform_file is None and ref_bed is None:
            pass
        elif ref_isoform_file is not None and ref_bed is not None:
            self._to_log('Parsing reference annotation and reference isoforms files')
            self.get_intersections_in_ref(ref_isoform_file, ref_bed)
        else:
            self._die(
                '--ref_isoform_file and --ref_bed options do not work separately; '
                'please provide both files if you want the script to consider reference nested genes'
            )
        self.proj2exon_cov: Dict[str, float] = {}
        if feature_file:
            self._to_log('Extracting chain exon coverage from the projection feature file')
            self.parse_feature_file(feature_file)
        self.proj2prob: Dict[str, float] = {}
        self.tr2max_prob: Dict[str, float] = {}
        self.paralog_list: Set[str] = parse_single_column(paralog_list)
        self.proc_pseudogene_list: Set[str] = parse_single_column(processed_pseudogene_list)
        self.lost_projections: Set[str] = set()
        self.proj2status: Dict[str, str] = {}
        self.parse_loss_file(loss_summary_file)
        self.orthology_threshold: float = orthology_threshold
        self.parse_orthology_file(orthology_probabilities)
        self.discarded_paralogs: Set[str] = set()
        self.discarded_paralogs_file: str = redundant_paralogs
        self.discarded_ppgenes: Set[str] = set()
        self.discarded_ppgenes_file: str = redundant_processed_pseudogenes
        self.discarded_extensions_file: click.File = insufficiently_covered_orthologs
        self.rejected_items_file: click.File = rejection_log

        self.run()

    def run(self) -> None:
        """
        Major execution method
        """
        self.build_intersection_graph()
        self.write_output()
        self.write_bed6_output()
        self.write_discarded_items()
        # self.write_redundant_paralogs()
        # self.write_unreliable_projections()

    def parse_bed(self, bed_handle: TextIO) -> None:
        """
        Parses query BED file; each entry is turned into an AnnotationEntry object.
        Results are stored chromosome-wise and sorted by their start coordinate.
        """
        for line in bed_handle.readlines():
            data: List[str] = line.rstrip().split('\t')
            chrom: str = data[0]
            start: int = int(data[1])
            stop: int = int(data[2])
            name: str = data[3]
            strand: bool = data[5] == '+'
            sizes: List[int] = list(map(int, data[10].split(',')[:-1]))
            starts: List[int] = list(map(int, data[11].split(',')[:-1]))
            max_ex_num: int = len(sizes)
            exon_dict: ExonDict = ExonDict()
            for i in range(max_ex_num):
                ex_num: int = i + 1 if strand else max_ex_num - i
                exon_start: int = start + starts[i]
                exon_stop: int = exon_start + sizes[i]
                exon: Exon = Exon(chrom, exon_start, exon_stop)
                exon_dict[ex_num] = exon
            annot_entry: AnnotationEntry = AnnotationEntry(
                name, chrom,
                start, stop, strand,
                max_ex_num,
                exon_dict
            )
            self.tr2chrom[name] = (*self.tr2chrom.get(name, ()), chrom)
            self.query_transcripts[chrom][name] = annot_entry
        # for chrom in self.query_transcripts:
        #     self.query_transcripts[chrom].sort(key=lambda x: x.start)

    def parse_exon_meta(self, file: TextIO) -> None:
        """
        Parses query meta. The following data are further used for gene inference:
        * exon coordinates;
        * exon loss status;
        * exon presence in the initial chains 
        """
        annotated_projections: Set[str] = set()
        for line in file:
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if data[0] == 'projection':
                continue
            proj: str = data[0]
            annotated_projections.add(proj)
            exon: int = int(data[1])
            chrom: str = data[3]
            start: int = int(data[4])
            end: int = int(data[5])
            strand: bool = data[6] == '+'
            status: str = data[7]
            if status != 'I':
                continue
            gap_supported: bool = data[-3] == 'CHAIN_SUPPORTED'
            if not gap_supported:
                continue
            self.tr2exons[proj][exon] = Coords(str(exon), chrom, start, end, strand)
            if proj in self.query_transcripts[chrom]:
                # _start, _end = self.query_transcripts[chrom][proj][1:3]
                prev_coords: Coords = self.query_transcripts[chrom][proj]
                _start: int = prev_coords.start
                _end: int = prev_coords.end
                start = min(start, _start)
                end = max(end, _end)
            self.query_transcripts[chrom][proj] = Coords(proj, chrom, start, end, strand)
            self.tr2chrom[proj] = (*self.tr2chrom.get(proj, ()), chrom)
        for proj in annotated_projections:
            if proj not in self.tr2chrom.keys():
                self.discarded_overextensions.add(proj)

    def parse_loss_file(self, file: TextIO) -> None:
        """Parses loss summary file, recording data on missing projections"""
        if file is None:
            return
        for i, line in enumerate(file):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if data[0] == 'level':
                continue
            if len(data) < 3:
                self._die(
                    'Encountered less than three columns in the loss summary file at line %i' %i
                )
            if data[0] != PROJECTION:
                continue
            proj: str = data[1]
            status: str = data[2]
            self.proj2status[proj] = status
            if status not in MISSING:
                continue
            self.lost_projections.add(proj)

    def parse_feature_file(self, file: TextIO) -> None:
        """Extracts exon coverage from the projection feature file"""
        if file is None:
            return
        self._to_log('Parsing projection classification features data')
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 16:
                self._die(
                    'Memory report file contains less than 16 fields at line % i' % i
                )
            if data[0] == 'transcript':
                continue
            tr: str = data[0]
            chain: str = data[2]
            proj: str = f'{tr}#{chain}'
            # synteny: int = int(data[3])
            # exon_qlen: float = float(data[7])
            exon_cover: int = int(data[9])
            # intr_cover: int = int(data[10])
            ex_fract: int = int(data[13])
            exon_fraction: float = exon_cover / ex_fract if ex_fract else 0.0
            # intr_fract: int = int(data[14])
            # intron_fraction: int = intr_cover / intr_fract if intr_fract else 0.0
            # features: FilteringFeatures = FilteringFeatures(
            #     synteny, exon_qlen, exon_fraction, intron_fraction
            # )
            self.proj2exon_cov[proj] = exon_fraction

    def parse_orthology_file(self, file: TextIO) -> None:
        """
        Extracts orthology probabilities for projections and 
        estimates maximal probability per transcript
        """
        if file is None:
            return
        self._to_log('Parsing orthology probabilies table')
        for line in file:
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if data[0] == 'transcript':
                continue
            prob: float = float(data[2])
            if prob < self.orthology_threshold and prob >= 0:
                continue
            tr: str = data[0]
            chain: str = data[1]
            proj: str = f'{tr}#{chain}'
            if proj in self.lost_projections:
                continue
            self.proj2prob[proj] = prob
            self.tr2max_prob[tr] = max(prob, self.tr2max_prob.get(tr, 0.0))

    def _overextended_projection(self, proj: str) -> bool:
        """
        Estimates whether projection is an overextended unreliable ortholog if:
        1) it is not the most probable ortholog according to the TOGA2 orthology classsifier, AND
        2) it has initial exon coverage less than 60%
        """
        if ',' in proj:
            return False
        if proj in self.paralog_list or proj in self.proc_pseudogene_list:
            return False
        tr: str = '#'.join(proj.split('#')[:-1])
        max_prob: float = self.tr2max_prob.get(tr, -1)
        if max_prob < 0:
            return False
        second_probable: bool = self.proj2prob.get(proj, 0.0) < max_prob
        overextended: bool = self.proj2exon_cov[proj] < MIN_RELIABLE_EXON_COV
        return second_probable and overextended

    def get_intersections_in_ref(self, ref_isoform_file: TextIO, ref_bed_file: TextIO) -> None:
        """
        Infers nested and intersected genes in the reference genome:
        1) At the initial step, parses reference gene-to-transcript mapping;
        2) Then, parses reference annotation BED file, extracting coding sequence coordinates;
        3) Finally, estimates sequence overlaps between isoforms coming from different genes
        """
        chrom2tr_coords: Dict[str, List[Tuple[str, int, int]]] = defaultdict(list)
        for line in ref_isoform_file:
            data: List[str] = line.rstrip().split()
            if not data or not data[0]:
                continue
            gene, tr = data
            self.ref_isoform2gene[tr] = gene
        for i, line in enumerate(ref_bed_file, 1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if not len(data) >= 8:
                self._die(
                    'ERROR: %i in the reference annotation file does not contain at least 8 columns' % i
                )
            chrom: str = data[0]
            cds_start: int = int(data[6])
            cds_end: int = int(data[7])
            name: str = data[3]
            if name not in self.ref_isoform2gene:
                self._to_log(
                    'Isoform %s is absent from the reference isoform file' % name
                )
            chrom2tr_coords[chrom].append((name, cds_start, cds_end))
        for chrom, trs in chrom2tr_coords.items():
            trs = sorted(trs, key=lambda x: (x[1], x[2]))
            for i, tr1 in enumerate(trs):
                gene1: str = self.ref_isoform2gene.get(tr1[0], '')
                if not gene1:
                    continue
                for tr2 in trs[i:]:
                    if tr2[1] >= tr1[2]:
                        break
                    gene2: str = self.ref_isoform2gene.get(tr2[0], '')
                    if not gene2:
                        continue
                    if gene1 == gene2:
                        continue
                    inter: int = intersection(*tr1[1:], *tr2[1:]) ## TODO: Redundant unless we impose overlap threshold
                    if inter <= 0:
                        continue
                    self.intersecting_ref_genes[gene1].add(gene2)
                    self.intersecting_ref_genes[gene2].add(gene1)


    def build_intersection_graph(self) -> None:
        """
        Creates a graph of exon intersections between the projections.
        An edge is traversed between two projection nodes if any of their exons
        intersect by 90% of their length.
        After that, infers all connected components from the projection
        intersection graph. The resulting components are treated as proxies
        for query genes.
        """
        graph: nx.Graph = nx.Graph()
        self._to_log('Building intersection graph for query transcripts')
        for chrom in self.query_transcripts:
            chrom_trs: List[AnnotationEntry] = sorted(
                self.query_transcripts[chrom].values(), key=lambda x: x.start
            )
            for i, proj_out in enumerate(chrom_trs):
                # if proj_out.name in self.lost_projections:
                #     continue
                if proj_out.name in self.discarded_paralogs:
                    continue
                if proj_out.name in self.discarded_ppgenes:
                    continue
                out_is_paralog: bool = proj_out.name in self.paralog_list
                out_is_pseudo: bool = proj_out.name in self.proc_pseudogene_list
                out_is_lost: bool = proj_out.name in self.lost_projections
                out_non_orth: bool = out_is_paralog or out_is_pseudo
                if ',' in proj_out.name:
                    sufficient_exon_cov_out: bool = True
                else:
                    sufficient_exon_cov_out: bool = not self._overextended_projection(proj_out.name)
                ## store the putative edges in a temporary collection
                edges: Set[Tuple[str, str]] = set()
                edges_added: bool = False
                ## under certain circumstances, this projection wil be deemed redundant
                ## once this is established, proceed with the next projection
                was_discarded: bool = False
                for j, proj_in in enumerate(chrom_trs[i+1:]):
                    if was_discarded:
                        break
                    if proj_in.name in self.discarded_paralogs:
                        continue
                    if proj_in.name in self.discarded_ppgenes:
                        continue
                    ## do not overlap lost projections and (putatively) intact projections
                    if proj_in.name in self.lost_projections and not out_is_lost:
                        continue
                    if proj_out.end < proj_in.start:
                        # print('Not reached')
                        continue
                    ## if the `out` projection has not been discarded by that point,
                    ## edges starting in it can be safely added to the graph
                    if proj_out.start > proj_in.end:
                        ## retrogenes/processed pseudogenes do not participate in gene inference
                        if not out_is_pseudo:# and not was_discarded:
                            if proj_out.name not in graph.nodes:
                                graph.add_node(proj_out.name)
                            for _out, _in in edges:
                                graph.add_edge(_out, _in)
                        edges_added: bool = True
                        break
                    if proj_out.strand != proj_in.strand:
                        continue
                    in_is_paralog: bool = proj_in.name in self.paralog_list
                    in_is_pseudo: bool = proj_in.name in self.proc_pseudogene_list
                    in_non_orth: bool = in_is_paralog or in_is_pseudo
                    # if in_is_pseudo:
                    #     continue
                    # if out_non_orth != in_non_orth:
                    #     self.discarded_paralogs.add(proj_out.name if out_non_orth else proj_in.name)
                    #     continue
                    if ',' in proj_in.name:
                        sufficient_exon_cov_in: bool = True
                    else:
                        sufficient_exon_cov_in: bool = not self._overextended_projection(proj_in.name)
                    tr_in: str = get_orig_transcript(proj_in.name)
                    tr_out: str = get_orig_transcript(proj_out.name)
                    gene_in: str = self.ref_isoform2gene.get(tr_in, '')
                    gene_out: str = self.ref_isoform2gene.get(tr_out, '')
                    if gene_in and gene_out:
                        if gene_in in self.intersecting_ref_genes.get(gene_out, []) or gene_out in self.intersecting_ref_genes.get(gene_in, []):
                            self._to_log(
                                (
                                    'Projections %s and %s will not be merged '
                                    'since genes %s and %s overlap in the reference'
                                ) % (proj_out.name, proj_in.name, gene_out, gene_in),
                                'warning'
                            )
                            continue
                    ## do not overlap insufficiently covered second-best projections
                    # if sufficient_exon_cov_out != sufficient_exon_cov_in and gene_out != gene_in:
                    #     continue
                    has_intersection: bool = False
                    # for exon1 in proj_out.exons.values():
                    #     e_start1, e_stop1 = exon1.start, exon1.stop
                    for exon1 in sorted(self.tr2exons[proj_out.name].values(), key=lambda x: x.start):
                        # print(f'{exon1=}')
                        e_start1, e_stop1 = exon1.start, exon1.end
                        # e1_50p: int = int((e_stop1 - e_start1) * 0.5)
                        # for exon2 in proj_in.exons.values():
                        #     e_start2, e_stop2 = exon2.start, exon2.stop
                        for exon2 in sorted(self.tr2exons[proj_in.name].values(), key=lambda x: x.start):
                            # print(f'{exon2=}')
                            if exon1.chrom != exon2.chrom:
                                continue
                            e_start2, e_stop2 = exon2.start, exon2.end
                            if e_start1 > e_stop2:
                                continue
                            if e_stop1 < e_start2:
                                break
                            # e2_50p: int = int((e_stop2 - e_start2) * 0.5)
                            inter_size: int = intersection(
                                e_start1, e_stop1, e_start2, e_stop2
                            )
                            ## at least two exons intersect plausibly,
                            ## no need to check for others
                            # if inter_size >= e1_50p and inter_size >= e2_50p:
                            if inter_size > 0:
                                has_intersection = True
                                break
                        if has_intersection:
                            break
                    # print(f'{has_intersection=}')
                    if has_intersection:
                        ## discarded any of the two current projections if the following applies
                        ## processed pseudogenes cannot overlap orthologs/paralogs
                        if out_is_pseudo and not in_is_pseudo:
                            self.discarded_ppgenes.add(proj_out.name)
                            was_discarded = True
                            continue
                        if in_is_pseudo and not out_is_pseudo:
                            self.discarded_ppgenes.add(proj_in.name)
                            continue
                        ## paralogs cannot overlap orthologs
                        if out_is_paralog and not in_is_paralog:
                            self.discarded_paralogs.add(proj_out.name)
                            was_discarded = True
                            continue
                        if in_is_paralog and not out_is_paralog:
                            self.discarded_paralogs.add(proj_in.name)
                            continue
                        edges.add((proj_out.name, proj_in.name))
                        # graph.add_edge(proj_out.name, proj_in.name)
                if not was_discarded and not out_is_pseudo:
                    if proj_out.name not in graph.nodes:
                        graph.add_node(proj_out.name)
                    for _out, _in in edges:
                        graph.add_edge(_out, _in)

        ## extract the connected components
        if NX_VERSION < 2.4:
            raw_components = list(nx.connected_component_subgraphs(graph))
        else:
            raw_components = [
                graph.subgraph(c) for c in nx.connected_components(graph)
            ]

        ## for each component, estimate their coordinates in the query
        coords_for_sorting: Dict[int, Tuple[str, int]] = {}
        self._to_log('Inferring query genes from transcript intersection graph')
        curr: int = 0
        for i, c in enumerate(raw_components, start=1):
            c = [x for x in c if x not in self.discarded_paralogs]
            if not c:
                self._die('Projection clique %i consists entirely of redundant entities' % i)
            visited: List[str] = []
            starts: Dict[str, int] = {}
            stops: Dict[str, int] = {}
            strands: Dict[str, bool] = {}
            ## extended (=insufficiently covered) projection filter
            if self.proj2exon_cov:
                # sufficiently_covered: List[str] = [
                #     x for x in c if self.proj2exon_cov.get(x, MIN_RELIABLE_EXON_COV) >= MIN_RELIABLE_EXON_COV
                # ]
                # insufficiently_covered: List[str] = [
                #     x for x in c if x not in sufficiently_covered
                # ]
                insufficiently_covered: List[str] = [
                    x for x in c if self._overextended_projection(x)
                ]
                sufficiently_covered: List[str] = [
                    x for x in c if x not in insufficiently_covered
                ]
                if insufficiently_covered:
                    ## if the locus contains both regular and extended projections, 
                    ## leave only properly covered ones
                    if sufficiently_covered:
                        c = sufficiently_covered
                        self.discarded_overextensions.update(insufficiently_covered)
                    ## otherwise, check how many genes were projected to this locus
                    else:
                        reliable_projs: List[str] = [
                            x for x in insufficiently_covered if 
                            self.proj2status[x] in EXTENDED_HIGH_CONFIDENCE
                        ]
                        ## no high-quality projections -> trash with potential to confound orthology inference
                        if not reliable_projs:
                            self.discarded_overextensions.update(insufficiently_covered)
                            continue
                        reliable_genes_in_locus: Set[str] = {
                            self.ref_isoform2gene['#'.join(x.split('#')[:-1])] for x in reliable_projs
                        }
                        ## overextended projections from at least two reference genes -> bad call;
                        ## drop this locus, it is most likely a false call
                        if len(reliable_genes_in_locus) > 1:
                            self.discarded_overextensions.update(insufficiently_covered)
                            continue
                        ## only one genes ends up with reliable projections -> all good,
                        ## but leave only those high quality projections
                        unreliable_projs: List[str] = [
                            x for x in c if x not in reliable_projs
                        ]
                        self.discarded_overextensions.update(unreliable_projs)
                        c = reliable_projs
            ## TODO:
            ## 1) think of how to store transcript-to-chromosome mapping
            ## 2) think of how to resolve multiple chromosomes in
            for tr in c:
                if tr in visited:
                    continue
                visited.append(tr)
                self.component2trs[curr].append(tr)
                chroms: Tuple[str] = self.tr2chrom[tr]
                for chrom in chroms:
                    tr_obj: AnnotationEntry = self.query_transcripts[chrom][tr]
                    strand: bool = self.query_transcripts[chrom][tr].strand
                    strands[chrom] = strand
                    starts[chrom] = (
                        tr_obj.start if chrom not in starts else
                        min(starts[chrom], tr_obj.start)
                    )
                    # stops[chrom] = (
                    #     tr_obj.stop if chrom not in stops else
                    #     max(stops[chrom], tr_obj.stop)
                    # )
                    stops[chrom] = (
                        tr_obj.end if chrom not in stops else
                        max(stops[chrom], tr_obj.end)
                    )
            min_chrom: int = min(starts.keys())
            for chrom in starts:
                locus_start: int = starts[chrom]
                locus_stop: int = stops[chrom]
                locus_strand: bool = strands[chrom]
                self.component_coords[curr].append(
                    (chrom, locus_start, locus_stop, locus_strand)
                )
                if chrom == min_chrom:
                    coords_for_sorting[curr] = (chrom, locus_start, locus_stop)
            curr += 1
        ## sort components by their coordinates
        self.sorted_components = sorted(
            coords_for_sorting.keys(), key=lambda x: coords_for_sorting[x]
        )

    def write_output(self) -> None:
        """
        """
        self._to_log('Writing query genes data')
        for i, c in enumerate(self.sorted_components, start=1):
            name: str = f'reg_{i}'
            for tr in self.component2trs[c]:
                self.output.write(f'{name}\t{tr}\n')

    def write_bed6_output(self) -> None:
        """
        """
        if not self.bed6_output:
            return
        self._to_log('Writing query genes coordinates')
        for i, c in enumerate(self.sorted_components, start=1):
            name: str = f'reg_{i}'
            for locus in self.component_coords[c]:
                chrom, start, stop, strand = locus
                strand = '+' if strand else '-'
                out_str: str = f'{chrom}\t{start}\t{stop}\t{name}\t0\t{strand}'
                self.bed6_output.write(out_str + '\n')

    def write_discarded_items(self) -> None:
        """Save discarded items to respective files"""
        for rej_orth in self.discarded_overextensions:
            if self.discarded_extensions_file is not None:
                self.discarded_extensions_file.write(rej_orth + '\n')
            if self.rejected_items_file is not None:
                status: str = self.proj2status[rej_orth]
                orth_rej_line: str = REJ_ORTH_REASON.format(rej_orth, status)
                self.rejected_items_file.write(orth_rej_line + '\n')
        for rej_par in self.discarded_paralogs:
            if self.discarded_paralogs_file is not None:
                self.discarded_paralogs_file.write(rej_par + '\n')
            if self.rejected_items_file is not None:
                status: str = self.proj2status[rej_par]
                par_rej_line: str = REJ_PARA_REASON.format(rej_par, status)
                self.rejected_items_file.write(par_rej_line + '\n')
        for rej_ppgene in self.discarded_ppgenes:
            if self.discarded_ppgenes_file is not None:
                self.discarded_ppgenes_file.write(rej_ppgene + '\n')
            if self.rejected_items_file is not None:
                status: str = self.proj2status[rej_ppgene]
                ppgene_rej_line: str = REJ_PPGENE_REASON.format(rej_ppgene, status)
                self.rejected_items_file.write(ppgene_rej_line + '\n')

    def write_redundant_paralogs(self) -> None:
        """Write the names of redundant paralogous projections to a file"""
        if self.discarded_paralogs_file is None:
            return
        for par in self.discarded_paralogs:
            self.discarded_paralogs_file.write(par + '\n')

    def write_unreliable_projections(self) -> None:
        """Write the names of projections which were ignored when inferring orthology"""
        if self.discarded_extensions_file is None:
            return
        for proj in self.discarded_overextensions:
            self.discarded_extensions_file.write(proj + '\n')


if __name__ == '__main__':
    QueryGeneCollapser()
