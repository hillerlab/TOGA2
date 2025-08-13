#!/usr/bin/env python3

"""
Infers orthology relationships between reference genes and collapsed projections
of respective transcripts in the query
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from cesar_wrapper_constants import (
    FI, I, PI, UL, M, L, PG, PP, N, 
    MAX_QLEN_FOR_ORTH, MIN_COV_FOR_ORTH, 
    MIN_INTRON_COV_FOR_ORTH, CLASS_TO_NUM
)
from collections import defaultdict
from constants import PHYLO_NOT_FOUND, PRE_CLEANUP_LINE
from heapq import heappop, heappush
from shared import (
    CommandLineManager, CONTEXT_SETTINGS, flatten, get_proj2trans,
    get_connected_components, SPLIT_JOB_HEADER
)
from shutil import which
from typing import (
    Any, Dict, Iterable, List, Optional, Set, TextIO, Tuple, Union
)

import click
import h5py
import networkx as nx

__author__ = 'Yury V. Malovichko'
__year__  = '2024'
__credits__ = 'Bogdan M. Kirilenko'

FINE_RESOLVER: str = os.path.join(LOCATION, 'fine_orthology_resolver.py')

Q_PREFIX: str = "#Q#"
R_PREFIX: str = "#R#"
ABS_EDGE_THRESHOLD: float = 0.75
REL_EDGE_THRESHOLD: float = 0.9
ALL_LOSS_SYMBOLS: Tuple[str] = (FI, I, PI, UL, M, L, PG, N)
DEFAULT_LOSS_SYMBOLS: Tuple[str] = (FI, I, PI, UL)
ONE2ZERO: str = 'one2zero'
ONE2ONE: str = 'one2one'
ONE2MANY: str = 'one2many'
MANY2ONE: str = 'many2one'
MANY2MANY: str = 'many2many'
HEADER: str = 't_gene\tt_transcript\tq_gene\tq_transcript\torthology_class\n'
REMOVED_ORTH_REASON: str = 'TRANSCRIPT\t{}\t0\tNo projections reached the orthology step\tALL_ORTHS_DISCARDED\t{}'
WEAK_EDGE_REASON: str = 'PROJECTION\t{}\t0\tWeak orthology graph edge\tWEAK_EDGE\t{}'
TOUCH: str = 'touch {}'


class FilteringFeatures:
    """
    An auaxiliary data classs for storing projection features used for orthology resolution
    """
    __slots__ = (
        'synteny', 'exon_cov_ratio', 'exon_coverage', 'intron_coverage'
    )

    def __init__(
        self, 
        synteny: int, 
        exon_cov_ratio: float, 
        exon_coverage: float, 
        intron_coverage: float
    ) -> None:
        self.synteny: int = synteny
        self.exon_cov_ratio: float = exon_cov_ratio
        self.exon_coverage: float = exon_coverage
        self.intron_coverage: float = intron_coverage


def extract_names_from_bed(file: TextIO) -> List[str]:
    """
    """
    output: List[str] = []
    for line in file:
        data: List[str] = line.rstrip().split('\t')
        if not data:
            continue
        try:
            output.append(data[3])
        except IndexError:
            raise Exception('BED file provided was improperly formatted')
    return output

def parse_isoforms(
    file: TextIO,
    allowed_names: List[str] = [],
    is_ref: bool = True
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    """
    gene2tr: Dict[str, List[str]] = defaultdict(list)
    tr2gene: Dict[str, str] = {}
    prefix: str = R_PREFIX if is_ref else Q_PREFIX
    if file is None:
        for tr in allowed_names:
            gene_subst: str = f'{prefix}{tr}'
            gene2tr[gene_subst] = [tr]
            tr2gene[tr] = gene_subst
        return gene2tr, tr2gene
    for line in file:
        data: List[str] = line.rstrip().split('\t')
        gene: str = f'{prefix}{data[0]}'
        tr: str = data[1]
        if allowed_names and tr not in allowed_names:
            continue
        gene2tr[gene].append(tr)
        tr2gene[tr] = gene
    return gene2tr, tr2gene

def ordered_edges(edges: Iterable[Tuple[Any]]) -> Iterable[Tuple[Any]]:
    """
    """
    ordered: List[Tuple[Any]] = []
    for i, edge in enumerate(edges):
        if R_PREFIX in edge[0]:
            ordered.append(edge)
            continue
        temp: str = edge[1]
        edge = (edge[1], edge[0], *edge[2:])
        ordered.append(edge)
    return ordered

def get_orth_class(ref_genes: Iterable[str], query_genes: Iterable[str]) -> str:
    """
    For two sets of genes in the subgraph, returns the orthology relationship class
    """
    ref_num: int = len(ref_genes)
    que_num: int = len(query_genes)
    if not ref_num:
        raise ValueError('Empty reference gene list was provided')
    if not que_num:
        return ONE2ZERO
    if ref_num == que_num == 1:
        return ONE2ONE
    if ref_num == 1:
        return ONE2MANY
    if que_num == 1:
        return MANY2ONE
    return MANY2MANY

def is_complete_bipartite(clique: nx.Graph) -> bool:
    """
    Returns whether the (sub)graph is complete bipartite, i.e. whether each
    node from one part is connected to each node from the counterpart
    """
    edges: Iterable[Tuple[str]] = ordered_edges(clique.edges())
    if not edges:
        return False
    ref_genes_cons: Dict[str, Set[str]] = defaultdict(set)
    query_genes_cons: Dict[str, Set[str]] = defaultdict(set)
    for edge in edges:
        r, q = edge
        ref_genes_cons[r].add(q)
        query_genes_cons[q].add(r)
    ref_gene_num: int = len(ref_genes_cons)
    query_gene_num: int = len(query_genes_cons)
    if any(len(y) != query_gene_num for y in ref_genes_cons.values()):
        return False
    if any(len(y) != ref_gene_num for y in query_genes_cons.values()):
        return False
    return True

def strongly_connected_subgraphs(graph: nx.Graph) -> Dict[str, List[str]]:
    """
    Selects reference gene nodes sharing more than one neighbor in the query
    """
    r2q: Dict[str, Set[str]] = defaultdict(set)
    q2r: Dict[str, Set[str]] = defaultdict(set)
    strong_connections: Dict[str, Set[str]] = {}
    for r, q in ordered_edges(graph.edges()):
        r2q[r].add(r)
        q2r[q].add(r)
    for r_gene, q_genes in r2q.items():
        counter: Dict[str, int] = defaultdict(int)
        for q_gene in q_genes:
            ref_neighbors: Set[str] = q2r[q_gene]
            for ref_n in ref_neighbors:
                counter[ref_n] += 1
        strong_connections[r_gene] = {k for k, v in counter.items() if v > 1}
    return strong_connections


def strong_connections_disrupted(
    graph: nx.Graph, strong_connections: Dict[str, Set[str]]
) -> bool:
    """
    Checks if graph connects all strongly connected nodes
    """
    if not strong_connections:
        return False
    nodes_in_graph: Dict[str, Set[str]] = {
        k:v for k,v in strong_connections.items() if k in graph.nodes()
    }
    if not nodes_in_graph:
        return False
    for v, conn in nodes_in_graph.items():
        if any(x not in graph.nodes() for x in conn):
            return True
    return False


def split_graph(clique: nx.Graph) -> Tuple[List[nx.Graph], List[Tuple[str, str]]]:
    """
    Splits a connected graph into subgraphs according to a set of rules
    """
    ## find 'leaves' (nodes with degree of 1)
    leaves: List[str] = [x for x in clique.nodes() if clique.degree(x) == 1]
    ## graphs with no leaves are deemed to complicated to solve
    if not leaves:
        return [clique], []
    ## get nodes adjacent (or equal) to leaves
    leaf_edges: List[Tuple[str]] = [
        x for x in clique.edges() if any(y in leaves for y in x)
    ]
    leaf_adjacent: Set[str] = set(flatten(leaf_edges))
    non_leaf_adjacent: Set[str] = set(clique.nodes()).difference(leaf_adjacent)
    ## get strongly connected reference nodes (nodes with >1 shared neighbor)
    strongly_connected: Dict[str, Set[str]] = strongly_connected_subgraphs(clique)
    ## now, create a local graph copy
    clique_copy: nx.Graph = clique.copy()
    ## subtract the leaf-adjacent nodes from a total graph
    leafless_graph: nx.Graph = clique_copy.subgraph(non_leaf_adjacent)
    ## get its counterpart and split it into connected
    components: List[nx.Graph] = get_connected_components(
        clique_copy.edge_subgraph(leaf_edges)
    )
    ## check whether leaf-trimmed subgraph has no isolated vertices
    if any(nx.isolates(leafless_graph)):
        return [clique], []
    ## add leafless subgraph as an extra component if it is non-empty
    if len(leafless_graph.nodes()):
        components.append(leafless_graph)
    ## check that no strong connections between reference genes are disrupted
    if any(strong_connections_disrupted(g, strongly_connected) for g in components):
        return [clique], []
    ## finally, check that all removed edges weigh consistently less
    ## than the preserved ones
    preserved_edges: List[Tuple[str, str, Dict[str, Any]]] = ordered_edges(
        [x for y in components for x in y.edges(data=True)]
    )
    removed_edges: List[Tuple[str, str, Dict[str, Any]]]= [
        x for x in ordered_edges(clique.edges(data=True))
        if x not in preserved_edges
    ]
    preserved_edge_scores: List[float] = [x[2]['weight'] for x in preserved_edges]
    # removed_edge_scores: List[float] = [
    #     x[2]['weight'] for x in ordered_edges(clique.edges(data=True))
    #     if x not in preserved_edges
    # ]
    removed_edge_scores: List[float] = [
        x[2]['weight'] for x in removed_edges
    ]
    removed_edges = [x[:2] for x in removed_edges]
    ## if all removed thresholds weigh less than 0.75, graph reduction is likely adequate
    if all(x < ABS_EDGE_THRESHOLD for x in removed_edge_scores):
        return components, removed_edges
    ## calculated a relative threshold and compare all the removed edges to it
    max_removed_score: float = max(removed_edge_scores)
    min_preserved_score: float = min(preserved_edge_scores)
    preservation_threshold: float = REL_EDGE_THRESHOLD * min_preserved_score
    ## removed edges do not differ significantly (sic) from the preserved ones;
    ## return the original graph
    if max_removed_score > preservation_threshold:
        return [clique], []
    ## reduction omitted adequately low-scoring edges; remove the extracted subgraphs
    return components, removed_edges


def resolve_many2many(clique: nx.Graph) -> Tuple[List[Tuple[Any]], str]:
    """
    """
    edges: List[Tuple[str, str]] = ordered_edges(clique.edges())
    ref_genes: List[str] = list({x[0] for x in edges})
    query_genes: List[str] = list({x[1] for x in edges})
    ## TODO: Should be a method???
    is_complete: bool = is_complete_bipartite(clique)
    if is_complete:
        return [(ref_genes, query_genes, MANY2MANY)], []
    output: List[Tuple[Any]] = []
    # removed_edges: List[str] = []
    subgraphs, removed_edges = split_graph(clique)
    # for sub in split_graph(clique):
    for sub in subgraphs:
        ref_in_sub: List[str] = [x for x in sub.nodes() if x in ref_genes]
        query_in_sub: List[str] = [x for x in sub.nodes() if x in query_genes]
        orthology_status: str = get_orth_class(ref_in_sub, query_in_sub)
        output.append((ref_in_sub, query_in_sub, orthology_status))
        # removed_edges.extend(_edges)
    return output, removed_edges


def undefined_only(record: str) -> bool:
    """
    Checks if sequence field in an amino acid FASTA record 
    contains only undefined symbols
    """
    seq: str = record.split('\n')[-1]
    return all(x in ('X', 'x', '*') for x in seq)


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'ref_bed',
    type=click.File('r', lazy=True),
    metavar='REF_BED'
)
@click.argument(
    'query_bed',
    type=click.File('r', lazy=True),
    metavar='QUERY_BED'
)
@click.argument(
    'loss_data',
    type=click.File('r', lazy=True),
    metavar='LOSS_SUMMARY'
)
@click.argument(
    'orthology_probs',
    type=click.File('r', lazy=True),
    metavar='ORTHOLOGY_SCORES'
)
@click.argument(
    'output',
    type=click.Path(exists=False),
    metavar='OUTPUT_DIR'
)
@click.option(
    '--ref_isoforms',
    '-ri',
    type=click.File('r', lazy=True),
    metavar='REF_ISOFORMS',
    default=None,
    show_default=True,
    help='A file containing reference gene-to-isoforms mapping'
)
@click.option(
    '--query_isoforms',
    '-qi',
    type=click.File('r', lazy=True),
    metavar='QUERY_ISOFORMS',
    default=None,
    show_default=True,
    help=(
        'A file containing query gene-to-isoforms mapping '
        '(query_genes.tsv file from TOGA output)'
    )
)
@click.option(
    '--paralogs',
    '-p',
    type=click.File('r', lazy=True),
    metavar='PARALOG_FILE',
    default=None,
    show_default=True,
    help='A single-column file containing paralogous projections\' identifiers'
)
@click.option(
    '--processed_pseudogenes',
    '-pp',
    type=click.File('r', lazy=True),
    metavar='PROCESSED_PSEUDOGENE_FILE',
    default=None,
    show_default=True,
    help='A single-column file containing processed pseudogene projections\' identifiers'
)
@click.option(
    '--accepted_loss_symbols',
    '-l',
    type=str,
    metavar='LOSS_SYMBOLS',
    default=','.join(DEFAULT_LOSS_SYMBOLS),
    show_default=True,
    help=(
        'A comma-separated list of loss status symbols; only projections of '
        'respective statuses will be considered when creating a connection graph. '
        'Keyword ALL lets all possible statuses in.'
    )
)
@click.option( ## TODO: Move to query gene inference code
    '--memory_report',
    '-mr',
    type=click.File('r', lazy=True),
    metavar='MEMORY_REPORT_FILE',
    default=None,
    show_default=False,
    help=(
        'Memory report file produced by TOGA2 after the data preprocessing step. '
        'If provided, chain-covered fraction for each projection will be used '
        'for orthology refinement'
    )
)
@click.option(
    '--projection_features', ## IDEA: Move it to gene inference codess
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
    '--schedule_tree_resolver_jobs',
    '-st',
    is_flag=True,
    default=False,
    show_default=True,
    help='If set, schedules fine resolution jobs for convoluted many2many clades'
)
@click.option(
    '--max_clique_size',
    '-mcs',
    type=int,
    metavar='INT',
    default=100,
    show_default=True,
    help=(
        'A maximum number of sequences in many:many cliques to be resolved with '
        'PRANK+RAxML pipeline'
    )
)
@click.option(
    '--fasta_file',
    '-f',
    type=click.Path(exists=True),
    metavar='FASTA_FILE',
    default=None,
    show_default=True,
    help=(
        'A protein (?) FASTA file produced by TOGA used for fine orthology resolution'
    )
)
@click.option(
    '--fasta_is_hdf5',
    '-fh',
    is_flag=True,
    default=False,
    show_default=False,
    help=(
        'If set, the code expectes the FASTA file provided to be an HDF5 storage'
    )
)
@click.option(
    '--job_number',
    '-j',
    type=int,
    metavar='INT',
    default=50,
    show_default=True,
    help=(
        'A number of jobs to split orthology fine resolution commands into'
    )
)
@click.option(
    '--use_raxml',
    '-raxml',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'Use RAxML (raxmlHPC-PTHREADS-AVX) instead of IqTree2 for tree inference'
    )
)
@click.option(
    '--prank_binary',
    '-pb',
    type=click.Path(exists=True),
    metavar='PRANK_BINARY',
    default=None,
    show_default=True,
    help=(
        'A path to the PRANK executable to be used at fine resolution step. '
        'If not provided, the program will try to infer its location from the PATH'
    )
)
@click.option(
    '--tree_binary',
    '-rb',
    type=click.Path(exists=True),
    metavar='TREE_BINARY',
    default=None,
    show_default=True,
    help=(
        'A path to the phylogeny tree executable to be used at fine resolution step. '
        'If not provided, the program will try to infer its location from the PATH.'
    )
)
@click.option(
    '--tree_cpus',
    '-rc',
    type=int,
    metavar='INT',
    default=1,
    show_default=True,
    help='A number of CPUs to parallel RAxML run onto'
)
@click.option(
    '--tree_bootstrap',
    '-rs',
    type=int,
    metavar='INT',
    default=100,
    show_default=True,
    help='A number of bootstrap replications for RAxML run'
)
@click.option(
    '--job_directory',
    '-jd',
    type=click.Path(exists=False),
    metavar='JOB_DIR',
    default=None,
    show_default=False,
    help=(
        'A path to write the fine resolution job files to '
        '[default: OUTPUT_DIR/jobs]'
    )
)
@click.option(
    '--fasta_directory',
    '-fd',
    type=click.Path(exists=False),
    metavar='FASTA_DIR',
    default=None,
    show_default=False,
    help=(
        'A path to write the fine resolution input FASTA files to '
        '[default: OUTPUT_DIR/fasta]'
    )
)
@click.option(
    '--results_directory',
    '-rd',
    type=click.Path(exists=False),
    default=None,
    show_default=False,
    help=(
        'A path to write the fine resolution results files to '
        '[default: OUTPUT_DIR/res]'
    )
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
    help='Control logging verbosity'
)

class InitialOrthologyResolver(CommandLineManager):
    """
    Resolves orthology relationships from the raw TOGA output; if specified,
    also schedules jobs for more elaborated PRANK+RAxML resolution.\n
    Arguments are:\n
    * REF_BED is a reference BED12 (BED4+) file used as TOGA input;\n
    * QUERY_BED is a query BED12 (BED4+) file produced by TOGA;\n
    * LOSS_SUMMARY is a TOGA output file containing transcript loss statuses;\n
    * ORTHOLOGY_SCORES is a three-column TOGA output file containing projection
    orthology probability scores;\n
    * OUTPUT is a name of directory to contain output files
    """

    __slots__ = (
        'gene2tr_ref', 'tr2gene_ref', 'gene2tr_que', 'tr2gene_que', 'tr2proj', 'proj2tr',
        'loss_status', 'proj2prob', 'output', 'paralogs', 'processed_pseudogenes', 
        'accepted_losses', 'proj2cov', 'proj2filt_features',
        'tree_resolver', 'max_clique_size', 'fasta_file', 'hdf5_fasta', 'job_num',
        'use_raxml', 'prank_bin', 'tree_bin', 'tree_cpus', 'tree_bootnum',
        'job_dir', 'fasta_dir', 'res_dir',
        'graph', 'orthology_report', 'removed_edges', 'removed_projections', 
        'cliques_to_resolve',
        'orthology_file', 'missing_transcripts', 'weak_orthology_projections',
        'rejection_file', 'jobs2cliques', 'jobfile'
    )

    def __init__(
        self,
        ref_bed: click.File,
        query_bed: click.File,
        loss_data: click.File,
        orthology_probs: click.File,
        output: click.Path,
        ref_isoforms: Optional[Union[click.File, None]],
        query_isoforms: Optional[Union[click.File, None]],
        paralogs: Optional[Union[click.File, None]],
        processed_pseudogenes: Optional[Union[click.File, None]],
        accepted_loss_symbols: Optional[str],
        memory_report: Optional[Union[click.File, None]],
        projection_features: Optional[Union[click.File, None]],
        schedule_tree_resolver_jobs: Optional[bool],
        max_clique_size: Optional[int],
        fasta_file: Optional[click.Path],
        fasta_is_hdf5: Optional[bool],
        job_number: Optional[int],
        use_raxml: Optional[bool],
        prank_binary: Optional[Union[click.Path, None]],
        tree_binary: Optional[Union[click.Path, None]],
        tree_cpus: Optional[int],
        tree_bootstrap: Optional[int],
        job_directory: Optional[Union[click.Path, None]],
        fasta_directory: Optional[Union[click.Path, None]],
        results_directory: Optional[Union[click.Path, None]],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)

        self._to_log('Extracting reference transcripts names for orthology resolution')
        ref_bed_names: List[str] = extract_names_from_bed(ref_bed)
        self._to_log('Extracting query projection names for orthology resolution')
        query_bed_names: List[str] = extract_names_from_bed(query_bed)
        self._to_log('Inferring reference gene-to-transcripts mapping for orthology resolution')
        self.gene2tr_ref, self.tr2gene_ref = parse_isoforms(
            ref_isoforms, ref_bed_names
        )
        self._to_log('Inferring query gene-to-transcripts mapping for orthology resolution')
        self.gene2tr_que, self.tr2gene_que = parse_isoforms(
            query_isoforms, query_bed_names, is_ref=False
        )
        self.tr2proj: Dict[str, List[str]] = defaultdict(list)
        self.proj2tr: Dict[str, str] = {}
        self._to_log('Inferring projection-to-transcript mapping for orthology resolution')
        self.get_tr2proj_mapping()


        if accepted_loss_symbols == 'ALL':
            self.accepted_losses: List[str] = ALL_LOSS_SYMBOLS
        else:
            self.accepted_losses: List[str] = [
                x for x in accepted_loss_symbols.split(',') if x
            ]
            invalid_symbols: List[str] = [
                x for x in self.accepted_losses if x not in ALL_LOSS_SYMBOLS
            ]
            if invalid_symbols:
                self._die(
                    f'ERROR: Invalid loss status symbols provided: {",".join(invalid_symbols)}'
                )

        self.loss_status: Dict[str, str] = {}
        self._to_log('Parsing gene loss summary')
        self.parse_loss_file(loss_data)
        self.proj2prob: Dict[str, float] = {}
        self._to_log('Parsing orthology score file')
        self.parse_score_file(orthology_probs)
        self.output: click.Path = output
        self.paralogs: List[str] = []
        self.parse_paralog_file(paralogs)
        self.processed_pseudogenes: List[str] = []
        self.parse_processed_pseudogene_file(processed_pseudogenes)
        self.proj2cov: Dict[str, float] = {}
        self.parse_memory_report(memory_report)
        # self.remove_overextended_projections()
        self.proj2filt_features: Dict[str, FilteringFeatures] = {}
        self.parse_feature_file(projection_features)
        self.remove_dubious_orthologs()

        self.tree_resolver: bool = schedule_tree_resolver_jobs
        self.max_clique_size: int = max_clique_size
        if self.tree_resolver and fasta_file is None:
            self._die(
                'ERROR: Fine orthology resolution was prompted with no FASTA '
                'file provided'
            )
        self.fasta_file: Union[click.Path, None] = fasta_file
        self.hdf5_fasta: bool = fasta_is_hdf5
        self.job_num: int = job_number
        self.use_raxml:  bool = use_raxml
        self.prank_bin: Union[click.Path, None] = prank_binary
        self.tree_bin: Union[click.Path, None] = tree_binary
        self.tree_cpus: int = tree_cpus
        self.tree_bootnum: int = tree_bootstrap

        self.job_dir: str = (
            job_directory if job_directory is not None
            else os.path.join(self.output, 'jobs')
        )
        self.fasta_dir: str = (
            fasta_directory if fasta_directory is not None
            else os.path.join(self.output, 'fasta')
        )
        self.res_dir: str = (
            results_directory if results_directory is not None
            else os.path.join(self.output, 'res')
        )

        self.graph: nx.Graph = nx.Graph()
        self.orthology_report: List[Tuple[Any]] = []
        self.removed_edges: List[Tuple[str, str]] = []
        self.removed_projections: List[str] = []
        self.cliques_to_resolve: List[List[str]] = []

        self.orthology_file: str = os.path.join(
            output, 'orthology_classification.tsv'
        )
        self.missing_transcripts: str = os.path.join(
            output, 'missing_transcripts.txt'
        )
        self.weak_orthology_projections: str = os.path.join(
            output, 'rejected_projection_names.txt'
        )
        self.rejection_file: str = os.path.join(output, 'rejected_by_graph_reduction.tsv')
        self.jobs2cliques: Dict[int, List[int]] = defaultdict(list)

        self.jobfile: str = os.path.join(self.job_dir, 'joblist')

        self.run()

    def run(self) -> None:
        """
        """
        self._mkdir(self.output)
        self._to_log('Out dir for orthology resolution successfully created')
        if self.tree_resolver:
            self.check_executables()
            self._to_log('Executables for gene tree-based orthology resolution checked')
            self._mkdir(self.job_dir)
            self._mkdir(self.fasta_dir)
            self._mkdir(self.res_dir)
        self.create_gene_graph()
        self.extract_connected_components()
        self.write_resolved_orthologies()
        if self.tree_resolver:
            self.schedule_tree_jobs()
            # self.write_job_files()
            self.extract_seqs_and_write_jobs()

    def check_executables(self) -> None:
        """
        Checks whether valid PRANK and RAxML/IQTree2 executables are available
        """
        if self.prank_bin is None:
            prank_path: Union[str, None] = which('prank')
            if prank_path is None:
                self._die(
                'ERROR: PRANK executable is missing from both PATH and user input'
            )
            self.prank_bin = os.path.abspath(prank_path)
        if self.tree_bin is None:
            tree_path: Union[str, None] = which(
                'iqtree2' if not self.use_raxml else 'raxmlHPC-PTHREADS-AVX'
            )
            if tree_path is None:
                err_msg: str = PHYLO_NOT_FOUND.format(
                    'IQTree2' if not self.use_raxml else 'raxmlHPC-PTHREADS-AVX'
                )
                self._die(err_msg)
            self.tree_bin = os.path.abspath(tree_path)

    def parse_loss_file(self, file: TextIO) -> None:
        """
        """
        for line in file:
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if data[0] != 'PROJECTION':
                continue
            proj: str = data[1]
            status: str = data[2]
            if status not in self.accepted_losses and proj in self.tr2gene_que:
                gene: str = self.tr2gene_que[proj]
                if proj in self.gene2tr_que[gene]:
                    self.gene2tr_que[gene].remove(proj)
                del self.tr2gene_que[proj]
                continue
            self.loss_status[proj] = status

    def parse_score_file(self, file: TextIO) -> None:
        """
        Parses the orthology score file, populating the {projection: prob} dictionary
        """
        for line in file:
            data: List[str] = line.rstrip().split('\t')
            if data[0] == 'transcript':
                continue
            proj: str = f'{data[0]}#{data[1]}'
            prob: float = float(data[2])
            self.proj2prob[proj] = prob

    def parse_paralog_file(self, file: Union[TextIO, None]) -> None:
        """
        Simply parses a single-column paralog file
        """
        if file is None:
            return
        self._to_log('Parsing paralogous projection file')
        for line in file:
            line = line.rstrip()
            if not line:
                continue
            self.paralogs.append(line)

    def parse_processed_pseudogene_file(self, file: Union[TextIO, None]) -> None:
        if file is None:
            return
        self._to_log('Parsing processed pseudogene projection file')
        for line in file:
            line = line.rstrip()
            if not line:
                continue
            self.processed_pseudogenes.append(line)

    def parse_memory_report(self, file: Union[TextIO, None]) -> None:
        if file is None:
            return
        self._to_log('Parsing projection metadata report')
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if data[0] == 'transcript':
                continue
            if len(data) < 11:
                self._die(
                    'Memory report file contains less than 11 fields at line % i' % i
                )
            tr: str = data[0]
            chain: str = data[1]
            proj: str = f'{tr}#{chain}'
            coverage: float = float(data[6])
            self.proj2cov[proj] = max(coverage, self.proj2cov.get(proj, 0.0))

    def parse_feature_file(self, file: Union[TextIO, None]) -> None:
        """
        Extracts data on exon length ratio in query and reference intron coverage 
        from the projection feature file
        """
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
            synteny: int = int(data[3])
            exon_qlen: float = float(data[7])
            exon_cover: int = int(data[9])
            intr_cover: int = int(data[10])
            ex_fract: int = int(data[13])
            exon_fraction: float = exon_cover / ex_fract if ex_fract else 0.0
            intr_fract: int = int(data[14])
            intron_fraction: int = intr_cover / intr_fract if intr_fract else 0.0
            features: FilteringFeatures = FilteringFeatures(
                synteny, exon_qlen, exon_fraction, intron_fraction
            )
            self.proj2filt_features[proj] = features

    ## DEPRECATED
    def remove_overextended_projections(self) -> None:
        """
        DEPRECATED
        Adjusts the transcript-to-projection relationships by original reference transcript coverage 
        before sequence extrapolation in the following fashion:
        * If any projection covered more than 50% of the reference transcript coding sequence , 
        remove all projections that covered less than 50% from the orthology graph;
        * Otherwise, leave the projection with the highest original coverage
        """
        if not self.proj2cov:
            return
        self._to_log('Removing projections with insufficient initial coverage')
        ## iterate over all reference genes
        for ref_g, ref_trs in self.gene2tr_ref.items():
            ## iterate over all isoforms for a given gene
            for ref_tr in ref_trs:
                ## and iterate over all the recorded projections
                projs: List[str] = [
                    x for x in self.tr2proj.get(ref_tr, []) if self.loss_status.get(x, N) in self.accepted_losses
                ]
                # print(f'{ref_tr=}, {projs=}')
                if not projs:
                    continue
                exceed_min_cov: List[str] = [x for x in projs if self.proj2cov.get(x, 0.0) >= MIN_COV_FOR_ORTH]
                if not exceed_min_cov:
                    max_cov: float = max(self.proj2cov.get(x, 0.0) for x in projs)
                    exceed_min_cov: List[str] = [x for x in projs if self.proj2cov.get(x, 0.0) == max_cov]
                for j, proj in enumerate(projs):
                    if proj in exceed_min_cov:
                        continue
                    if proj in self.paralogs:
                        continue
                    if proj in self.processed_pseudogenes:
                        continue
                    if ',' in proj.split('#')[-1]:
                        continue
                    self._to_log(
                        (
                            'Projection %s is excluded from gene orthology graph '
                            'due to insufficient original coverage'
                        ) % proj,
                        'warning'
                    )
                    self.tr2proj[ref_tr].remove(proj)
                    # rej_line: str = COVERAGE_REJ_REASON.format(proj)
                    # self.rejected_projections.append()

    def remove_dubious_orthologs(self) -> None:
        """
        Removes untrustworthy orthologs suspected to be false calls 
        from the classification step:
        * Projections with 
        """
        if not self.proj2filt_features:
            return
        ## iterate over all reference genes
        for ref_g, ref_trs in self.gene2tr_ref.items():
            ## iterate over all isoforms for a given gene
            for ref_tr in ref_trs:
                ## and iterate over all the recorded projections
                projs: List[str] = [
                    x for x in self.tr2proj[ref_tr]
                ]
                if not projs:
                    continue
                projs.sort(key=lambda x: -self.proj2prob.get(x, 0.0))
                top_proj: str = projs[0]
                top_proj_status: str = self.loss_status.get(top_proj[0], N)
                top_proj_prob: float = self.proj2prob.get(top_proj, 0.0)
                if not top_proj_prob:
                    continue
                top_proj_features: FilteringFeatures = self.proj2filt_features[top_proj]
                for other_proj in projs[1:]:
                    if other_proj in self.paralogs:
                        continue
                    if other_proj in self.processed_pseudogenes:
                        continue
                    if ',' in other_proj:
                        continue
                    other_proj_status: str = self.loss_status.get(other_proj, N)
                    if CLASS_TO_NUM[other_proj_status] > CLASS_TO_NUM[top_proj_status]:
                        continue
                    other_proj_prob: float = self.proj2prob[other_proj]
                    if other_proj_prob == top_proj_prob:
                        continue
                    other_proj_features: FilteringFeatures = self.proj2filt_features[other_proj]
                    to_remove: bool = False
                    if other_proj_features.exon_cov_ratio > MAX_QLEN_FOR_ORTH:
                        to_remove = True
                    if other_proj_features.intron_coverage < MIN_INTRON_COV_FOR_ORTH:
                        to_remove = True
                    if to_remove:
                        self._to_log(
                            'Removing projection %s as an uncertain ortholog' % other_proj,
                            'warning'
                        )
                        self.tr2proj[ref_tr].remove(other_proj)

    def get_tr2proj_mapping(self) -> None:
        """
        Populates a {reference_transcript:[projections]}
        """
        for proj in self.tr2gene_que:
            tr: str = '#'.join(proj.split('#')[:-1])
            self.tr2proj[tr].append(proj)
            self.proj2tr[proj] = tr

    def create_gene_graph(self) -> None:
        """
        """
        ## initialize all reference genes as graph nodes
        # print(f'{self.gene2tr_ref.keys()=}')
        self.graph.add_nodes_from(self.gene2tr_ref.keys())
        ## here come the loops!
        ## first, iterate over all reference genes
        for ref_g, ref_trs in self.gene2tr_ref.items():
            ## iterate over all of its isoforms
            query_g_scores: Dict[str, float] = defaultdict(float)
            for ref_tr in ref_trs:
                ## and iterate over all the recorded projections
                # trusted_projections: List[str] = []
                for proj in self.tr2proj.get(ref_tr, []):
                    # v = proj == 'ENST00000322861.12#MYL11#325045'
                    if self.loss_status.get(proj, N) not in self.accepted_losses:
                        # if v:
                        #     print('ENST00000322861.12#MYL11#325045 is lost')
                        continue
                    if proj in self.paralogs:
                        # if v:
                        #     print('ENST00000322861.12#MYL11#325045 is a paralog')
                        continue
                    if proj in self.processed_pseudogenes:
                        # if v:
                        #     print('ENST00000322861.12#MYL11#325045 is a ppgene')
                        continue
                    # trusted_projections.append(proj)
                    # if v:
                    #     print('ENST00000322861.12#MYL11#325045 is in the graph')
                    query_g: str = self.tr2gene_que[proj]
                    prob: float = self.proj2prob.get(proj, 0.0)
                    query_g_scores[query_g] = max(query_g_scores[query_g], prob)
            ## once all possible query genes have been captured,
            ## add the respective edges to the graph
            for query_g, edge_score in query_g_scores.items():
                # edge_key: Tuple[str, str] = (ref_g, query_g)
                self.graph.add_edge(ref_g, query_g, weight=edge_score)

    def extract_connected_components(self) -> None:
        """
        """
        ## get connected components from the resulting graph
        raw_components: List[nx.Graph] = get_connected_components(self.graph)
        ## and try to resolve every connected component found
        for component in raw_components:
            ## get the nodes for this component
            nodes: Iterable[str] = component.nodes()
            ref_nodes: List[str] = [x for x in nodes if x in self.gene2tr_ref]
            query_nodes: List[str] = [x for x in nodes if x in self.gene2tr_que]
            ## get the orthology class
            orthology_class: str = get_orth_class(ref_nodes, query_nodes)
            # if '#R#ENSG00000180209' in ref_nodes:
            #     print(f'Clique with ENSG00000180209: {ref_nodes=}, {query_nodes=}, {orthology_class=}')
            if orthology_class == MANY2MANY:
                resolved_components, removed_edges = resolve_many2many(component)
                # for resolved_component, removed_edges in resolve_many2many(component):
                for resolved_component in resolved_components:
                    ref_in_comp, query_in_comp, comp_stat = resolved_component
                    if comp_stat == MANY2MANY:
                        clique_size: int = len(ref_in_comp) + len(query_in_comp)
                        if self.tree_resolver and clique_size <= self.max_clique_size:
                            # self.schedule_resolver_jobs
                            self.cliques_to_resolve.append(ref_in_comp + query_in_comp)
                            # continue ## Add all the cliques to output for now
                    self.orthology_report.append(resolved_component)
                self.removed_edges.extend(removed_edges)
                continue
            self.orthology_report.append((ref_nodes, query_nodes, orthology_class))

    def write_resolved_orthologies(self) -> None:
        """
        """
        projected_transcripts: Set[str] = set()
        with open(self.orthology_file, 'w') as h:
            h.write(HEADER)
            for (ref_genes, query_genes, status) in self.orthology_report:
                # for r_gene in ref_genes:
                #     _r_gene: str = r_gene[3:]
                #     for q_gene in query_genes:
                #         _q_gene: str = q_gene[3:]
                #         projections: List[str] = self.gene2tr_que[q_gene]
                #         for proj in projections:
                #             tr, _ = get_proj2trans(proj)
                #             projected_transcripts.add(tr)
                #             l: str = f'{_r_gene}\t{tr}\t{_q_gene}\t{proj}\t{status}\n'
                #             h.write(l)
                #     if status == ONE2ZERO:
                #         for tr in self.gene2tr_ref[r_gene]:
                #             l: str = f'{_r_gene}\t{tr}\tNone\tNone\t{status}\n'
                #             h.write(l)
                for r_gene in ref_genes:
                    _r_gene: str = r_gene[3:]
                    ref_trs: List[str] = self.gene2tr_ref[r_gene]
                    for tr in ref_trs:
                        if status == ONE2ZERO:
                            projected_transcripts.add(tr)
                            l: str = f'{_r_gene}\t{tr}\tNone\tNone\t{status}\n'
                            h.write(l)
                        else:
                            projections: List[str] = self.tr2proj[tr]
                            for proj in projections:
                                # v = proj == 'ENST00000322861.12#MYL11#325045'
                                if proj not in self.tr2gene_que:
                                    # if v:
                                    #     print('ENST00000322861.12#MYL11#325045 not in tr2gene_que')
                                    continue
                                q_gene: str = self.tr2gene_que[proj]
                                if q_gene not in query_genes:
                                    # if v:
                                    #     print(f'Gene {q_gene} for ENST00000322861.12#MYL11#325045 is not in query_genes ({query_genes})')
                                        # if (r_gene, q_gene) in self.removed_edges:
                                        #     print('Yep, the edge was removed')
                                    continue
                                projected_transcripts.add(tr)
                                _q_gene: str = q_gene[3:]
                                l: str = f'{_r_gene}\t{tr}\t{_q_gene}\t{proj}\t{status}\n'
                                h.write(l)
        non_projected_transcripts: Set[str] = set(self.tr2gene_ref.keys()).difference(
            projected_transcripts
        )
        if non_projected_transcripts:
            with open(self.missing_transcripts, 'w') as h:
                for tr in non_projected_transcripts:
                    h.write(tr + '\n')
        with open(self.rejection_file, 'w') as h, open(self.weak_orthology_projections, 'w') as h1:
            for tr in non_projected_transcripts:
                if tr in self.tr2proj and self.tr2proj[tr]:
                    status: str = max(
                        [self.loss_status.get(x, 'N') for x in self.tr2proj[tr]],
                        key=lambda x: CLASS_TO_NUM[x]
                    )
                else:
                    status: str = 'N'
                line: str = REMOVED_ORTH_REASON.format(tr, status)
                h.write(line + '\n')
            for (r_gene, q_gene) in self.removed_edges:
                if r_gene not in self.gene2tr_ref:
                    self._die('Reference gene %s is missing from the mapping' % r_gene)
                r_trs: str = ','.join(self.gene2tr_ref[r_gene])
                r_trs: List[str] = self.gene2tr_ref[r_gene]
                all_projs: List[str] = [
                    x for y in r_trs for x in self.tr2proj[y]
                ]
                if q_gene not in self.gene2tr_que:
                    self._die('Query gene %s is missing from the mapping' % q_gene)
                removed_projections: List[str] = [
                    x for x in all_projs if x in self.gene2tr_que[q_gene]
                ]
                self.removed_projections.extend(removed_projections)
                for proj in removed_projections:
                    status: str = self.loss_status[proj] 
                    line: str = WEAK_EDGE_REASON.format(proj, status)
                    h.write(line + '\n')
                    h1.write(proj + '\n')

    def schedule_tree_jobs(self) -> None:
        """
        """
        jobs: List[Tuple[int, int]] = [(0, i) for i in range(self.job_num)]
        for i, clique in enumerate(self.cliques_to_resolve):
            clique_size: int = len(clique)
            job_size, job = heappop(jobs)
            self.jobs2cliques[job].append(i)
            heappush(jobs, (job_size + clique_size, job))

    def extract_from_fasta(self, names: List[str]) -> Dict[str, str]:
        """
        Extracts entries from the TOGA output FASTA file for the provided sequences
        """
        ## TODO: Add longest isoform selection!!!
        # names = [x[3:] if x[:3] in (R_PREFIX, Q_PREFIX) else x for x in names]
        output_dict: Dict[str, str] = {}
        gene2longest_isoform: Dict[str, Tuple[str, int]] = {}
        gene2status: Dict[str, str] = {}
        header: str = ''
        seq: str = ''
        prev_gene: str = ''
        gene: str = ''
        proj: str = ''
        tr: str = ''
        prefix: str = ''
        with open(self.fasta_file, 'r') as h:
            for line in h:
                line = line.rstrip()
                if not line:
                    continue
                if line[0] == '>': ## fasta entry header's encountered
                    ## finish the prervious entry
                    if header:
                        seq_len: int = len(seq.replace('-', ''))
                        loss_status: str = self.loss_status[tr]
                        ## check if any isoform has been already considered for this gene
                        if gene in gene2longest_isoform:
                            prev_isoform, prev_status = gene2status[gene]
                            prev_isoform, prev_seq_len = gene2longest_isoform[gene]
                            better_conservation: bool = CLASS_TO_NUM[loss_status] > CLASS_TO_NUM[prev_status]
                            same_conservation: bool = CLASS_TO_NUM[loss_status] ==  CLASS_TO_NUM[prev_status]
                            new_is_longer: bool = seq_len > prev_seq_len 
                            # if CLASS_TO_NUM[prev_status] CLASS_TO_NUM[loss_status]
                            ## if current isoform is the longest one, delete the previous entries
                            ## and record data for this one
                            if better_conservation or same_conservation and new_is_longer:
                            # if seq_len > prev_seq_len:
                                del output_dict[f'>{prev_isoform}']
                                output_dict[header] = seq
                                gene2longest_isoform[gene] = (f'{prefix}{tr}', seq_len)
                            ## otherwise, we are not interested in this entry, so don't add it
                        else:
                            ## instance of this gene has not been yet encountered;
                            ## add this isoform for now
                            output_dict[header] = seq
                            gene2longest_isoform[gene] = (f'{prefix}{tr}', seq_len)
                            gene2status[gene] = (f'{prefix}{tr}', loss_status)
                        seq = ''
                        header = ''
                    ## now, process the new header
                    split_header: List[str] = line.split(' | ')
                    seq_id: str = split_header[0]
                    tr = seq_id[1:]
                    if tr not in self.loss_status or self.loss_status[tr] not in self.accepted_losses:
                        continue
                    if tr in self.paralogs:
                        continue
                    ref_tr: str = '#'.join(tr.split('#')[:-1])
                    if tr not in self.tr2proj[ref_tr]:
                        continue
                    source: str = split_header[2]
                    prefix: str = R_PREFIX if source == 'REFERENCE' else Q_PREFIX
                    ## check if gene is in the gene name list
                    if source == 'REFERENCE':
                        ref_tr = '#'.join(tr.split('#')[:-1])
                        gene = self.tr2gene_ref.get(ref_tr, None)
                    elif source == 'QUERY':
                        gene = self.tr2gene_que.get(tr, None)
                    if gene not in names:
                        gene = ''
                        tr = ''
                        source = ''
                        prefix = ''
                        continue
                    seq_id = seq_id.replace('>', f'>{prefix}')
                    header = seq_id
                    continue
                ## if this is not a header but a valid one has been encountered
                ## earlier, extend the sequence line
                if header:
                    seq += line.replace('-', '').replace('*', '')
            if seq:
                output_dict[header] = seq
        return output_dict

    def extract_from_hdf(self, names: List[str]) -> Dict[str, str]:
        """
        Extracts entries from the HDF5-converted pairwise FASTA
        """
        output_dict: Dict[str, str] = {}
        # gene2longest_isoform: Dict[str, Tuple[str, int]] = {}
        # gene2status: Dict[str, str] = {}
        for gene in names:
            is_ref: bool = False
            if gene in self.gene2tr_ref:
                is_ref = True
                trs: List[str] = self.gene2tr_ref[gene]
                postfix: str = '_ref'
            elif gene in self.gene2tr_que:
                trs: List[str] = self.gene2tr_que[gene]
                postfix: str = '_query'
            else:
                self._die(
                    'Gene %s does not have any transcripts in either reference or query' % gene,
                )
            with h5py.File(self.fasta_file, 'r') as f:
                if is_ref:
                    projs: List[str] = [
                        x.rstrip(postfix) for x in f.keys() 
                        if postfix in x and '#'.join(x.split('#')[:-1]) in trs
                    ]
                    projs = [
                        x for x in projs if x in self.tr2proj['#'.join(x.split('#')[:-1])]
                    ]
                    best_status: str = max(
                        [
                            self.loss_status.get(x, N) for x in projs if 
                            x not in self.processed_pseudogenes and x not in self.paralogs
                        ],
                        key=lambda y: CLASS_TO_NUM[y]
                    )
                else:
                    projs: List[str] = trs
                    best_status: str = max(
                        [
                            self.loss_status.get(x, N) for x in projs if 
                            x not in self.processed_pseudogenes and x not in self.paralogs
                        ],
                        key=lambda y: CLASS_TO_NUM[y]
                    )
                best_projs: List[str] = [x for x in projs if self.loss_status.get(x, N) == best_status]
                header, seq = '', ''
                for proj in best_projs:
                    entry: str = f[f'{proj}{postfix}'][()].decode('utf8')
                    _header, _seq = entry.split('\n')
                    if len(_seq) > len(seq):
                        header, seq = _header, _seq
                if not header or not seq:
                    self._die(
                        'Failed to find the longest relevant isoform for gene %s' % gene
                    )
                output_dict[header] = seq
        return output_dict
                    

    def write_job_files(self) -> None:
        """
        Writes input FASTA files, job files, and a job list
        for fine orthology resolution step
        """
        with open(self.jobfile, 'w') as jl:
            for j, cliques in self.jobs2cliques.items():
                job_path: str = os.path.join(self.job_dir, f'batch{j}.ex')
                table_path: str = os.path.join(self.fasta_dir, f'batch{j}.txt')
                res_path: str = os.path.join(self.res_dir, f'batch{j}')
                all_fasta_files: List[str] = []
                for c in cliques:
                    clique: List[str] = self.cliques_to_resolve[c]
                    self._to_log(
                        f'Writing FASTA input for clique {c} '
                        f'({len(clique)} sequences)'
                    )
                    if self.hdf5_fasta:
                        fasta_seqs: Dict[str, str] = self.extract_from_hdf(clique)
                    else:
                        fasta_seqs: Dict[str, str] = self.extract_from_fasta(clique)
                    # print(f'{fasta_seqs=}')
                    fasta_path: str = os.path.join(
                        self.fasta_dir, f'batch{j}_clique{c}.fa'
                    )
                    with open(fasta_path, 'w') as fp:
                        for header, seq in fasta_seqs.items():
                            fp.write(header + '\n' + seq + '\n')
                    fasta_path = os.path.abspath(fasta_path)
                    all_fasta_files.append(fasta_path)
                with open(table_path, 'w') as t:
                    for fasta_file in all_fasta_files:
                        t.write(fasta_file + '\n')
                table_path = os.path.abspath(table_path)
                res_path = os.path.abspath(res_path)
                cmd: str = (
                    f'{FINE_RESOLVER} {table_path} {res_path} -t '
                    f'-pb {self.prank_bin} -rb {self.tree_bin} '
                    f'-rc {self.tree_cpus} -rs {self.tree_bootnum}'
                )
                job_path = os.path.abspath(job_path)
                with open(job_path, 'w') as jf:
                    jf.write('\n'.join(SPLIT_JOB_HEADER) + '\n')
                    jf.write(cmd + '\n')
                file_mode: bytes = os.stat(job_path).st_mode
                file_mode |= (file_mode & 0o444) >> 2
                os.chmod(job_path, file_mode)
                jl.write(job_path + '\n')

    def pick_representatives(self, clique: List[str], storage: Any) -> List[str]:
        """
        """
        output_list: List[str] = []
        query_genes: List[str] = [x for x in clique if x in self.gene2tr_que]
        gene2longest: Dict[str, Tuple[str, int]] = {}
        gene2best_status: Dict[str, Tuple[str, str]] = {}
        all_projections: List[str] = []
        for query_gene in query_genes:
            projections: List[str] = self.gene2tr_que[query_gene]
            for proj in projections:
                if proj in self.paralogs: 
                    continue
                if proj in self.processed_pseudogenes:
                    continue
                if proj in self.removed_projections:
                    continue
                tr: str = '#'.join(proj.split('#')[:-1])
                if proj not in self.tr2proj[tr]:
                    continue
                all_projections.append(proj)
                loss_status: str = self.loss_status.get(proj, N)
                seq_id: str = f'{proj}_query'
                seq: str = storage[seq_id][()].decode('utf8').split('\n')[1].replace('-', '')
                prev_header, prev_seq = gene2longest.get(query_gene, ('', ''))
                new_is_longer: bool = len(seq) > len(prev_seq)
                prev_status: str = gene2best_status.get(query_gene, N)
                better_conservation: bool = CLASS_TO_NUM[loss_status] > CLASS_TO_NUM[prev_status]
                same_conservation: bool = CLASS_TO_NUM[loss_status] ==  CLASS_TO_NUM[prev_status]
                if better_conservation or same_conservation and new_is_longer:
                    gene2longest[query_gene] = (proj, seq)
                    gene2best_status[query_gene] = loss_status
        # selected_isoforms: List[str] = [
        #     v[0] for v in gene2longest.values()
        # ]
        ref_genes: List[str] = [x for x in clique if x not in query_genes]
        ref_gene2longest: Dict[str, Tuple[str, str]] = {}
        ref_gene2best_status: Dict[str, str] = {}
        for ref_gene in ref_genes:
            # ref_tr: List[str] = next(x for x in self.gene2tr_ref[ref_gene] if x in
            ref_trs: List[str] = self.gene2tr_ref[ref_gene]
            for proj in all_projections:
                progenitor: str = '#'.join(proj.split('#')[:-1])
                if progenitor not in ref_trs:
                    continue
                loss_status: str = self.loss_status.get(proj, N)
                seq_id: str = f'{proj}_ref'
                seq: str = storage[seq_id][()].decode('utf8').split('\n')[1].replace('-','')
                prev_header, prev_seq = ref_gene2longest.get(ref_gene, ('', ''))
                new_is_longer: bool = len(seq) > len(prev_seq)
                prev_status: str = ref_gene2best_status.get(ref_gene, N)
                better_conservation: bool = CLASS_TO_NUM[loss_status] > CLASS_TO_NUM[prev_status]
                same_conservation: bool = CLASS_TO_NUM[loss_status] ==  CLASS_TO_NUM[prev_status]
                if better_conservation or same_conservation and new_is_longer:
                    ref_gene2longest[ref_gene] = (proj, seq)
                    ref_gene2best_status[ref_gene] = loss_status
        for name, seq in gene2longest.values():
            fasta_record: str = f'>{Q_PREFIX}{name}\n{seq}'
            output_list.append(fasta_record)
        for name, seq in ref_gene2longest.values():
            fasta_record: str = f'>{R_PREFIX}{name}\n{seq}'
            output_list.append(fasta_record)
        return output_list

    def extract_seqs_and_write_jobs(self) -> None:
        """
        Writes the tree-based orthology jobs and the respective input,
        extracting the sequences on the fly
        """
        if not self.jobs2cliques:
            return
        with open(self.jobfile, 'w') as jl, h5py.File(self.fasta_file, 'r') as f:
            for j, cliques in self.jobs2cliques.items():
                job_path: str = os.path.join(self.job_dir, f'batch{j}.ex')
                table_path: str = os.path.join(self.fasta_dir, f'batch{j}.txt')
                res_path: str = os.path.join(self.res_dir, f'batch{j}')
                all_fasta_files: List[str] = []
                for n, c in enumerate(cliques):
                    clique: List[str] = self.cliques_to_resolve[c]
                    self._to_log(
                        f'Writing FASTA input for clique {c} '
                        f'({len(clique)} sequences)'
                    )
                    # if self.hdf5_fasta:
                    #     fasta_seqs: Dict[str, str] = self.extract_from_hdf(clique)
                    # else:
                    #     fasta_seqs: Dict[str, str] = self.extract_from_fasta(clique)
                    # print(f'{fasta_seqs=}')
                    fasta_seqs: List[str] = self.pick_representatives(clique, f)
                    if all(len(x.split('\n')[1]) < 50 for x in fasta_seqs):
                        self._to_log(
                            (
                                'Clique %i in job %i has all of its representatives shorter than 50 amino acids; '
                                'full representation of all 20 amino acid states required by RAxML cannot '
                                'be ascertained; skipping'
                            ) % (n, j),
                            'warning'
                        )
                        continue
                    if any(undefined_only(x) for x in fasta_seqs):
                        self._to_log(
                            (
                                'Clique %i in job %i has at least one sequence containing '
                                'only undefined symbols; skipping'
                            ) % (n, j), 
                            'warning'
                        )
                        continue
                    fasta_path: str = os.path.join(
                        self.fasta_dir, f'batch{j}_clique{c}.fa'
                    )
                    with open(fasta_path, 'w') as fp:
                        for entry in fasta_seqs:
                            fp.write(entry + '\n')
                    fasta_path = os.path.abspath(fasta_path)
                    all_fasta_files.append(fasta_path)
                with open(table_path, 'w') as t:
                    for fasta_file in all_fasta_files:
                        t.write(fasta_file + '\n')
                table_path = os.path.abspath(table_path)
                res_path = os.path.abspath(res_path)
                cmd: str = (
                    f'{FINE_RESOLVER} {table_path} {res_path} -t '
                    f'-pb {self.prank_bin} -rb {self.tree_bin} '
                    f'-rc {self.tree_cpus} -rs {self.tree_bootnum}'
                )
                if self.use_raxml:
                    cmd += ' -raxml'
                job_path = os.path.abspath(job_path)
                ok_file: str = os.path.join(res_path, '.ok')
                with open(job_path, 'w') as jf:
                    jf.write('\n'.join(SPLIT_JOB_HEADER) + '\n')
                    jf.write(PRE_CLEANUP_LINE.format(res_path) + '\n')
                    jf.write(cmd + '\n')
                    jf.write(TOUCH.format(ok_file) + '\n')
                file_mode: bytes = os.stat(job_path).st_mode
                file_mode |= (file_mode & 0o444) >> 2
                os.chmod(job_path, file_mode)
                jl.write(job_path + '\n')


if __name__ == '__main__':
    InitialOrthologyResolver()
