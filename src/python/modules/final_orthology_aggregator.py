#!/usr/bin/env python3

"""
Aggregates the results from the initial (graph-based) and fine (tree-based)
orthology resolution steps
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from cesar_wrapper_constants import I, PI, UL, M, L, PG, N
from constants import Headers
from collections import defaultdict
from shared import (
    CommandLineManager, CONTEXT_SETTINGS, get_proj2trans
)
from typing import Dict, List, Optional, Set, TextIO, Tuple, Union

import click

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = ('Amy Stephen', 'Michael Hiller',  'Bogdan M. Kirilenko')
__all__ = (None)

Q_PREFIX: str = "#Q#"
R_PREFIX: str = "#R#"
# ABS_EDGE_THRESHOLD: float = 0.75
# REL_EDGE_THRESHOLD: float = 0.9
# ALL_LOSS_SYMBOLS: Tuple[str] = (I, PI, UL, M, L, PG, N)
# DEFAULT_LOSS_SYMBOLS: Tuple[str] = (I, PI, UL)
ONE2ZERO: str = 'one2zero'
ONE2ONE: str = 'one2one'
ONE2MANY: str = 'one2many'
MANY2ONE: str = 'many2one'
MANY2MANY: str = 'many2many'
ORTH_REJ_TEMPLATE: str = 'TRANSCRIPT\t{}\t0\tRejected after the gene resolution step\tGENE_TREE_REJECTION\t{}'

def restore_fragmented_proj_id(proj: str) -> str:
    """
    Restores commas in the fragmented projection names replaced with underscores by RAxML
    """
    components: List[str] = proj.split('#')
    chain_id: str = components[-1]
    if chain_id.isdigit():
        return proj
    if '_' in chain_id and ',' not in chain_id:
        chain_id = chain_id.replace('_', ',')
        components[-1] = chain_id
        return '#'.join(components)
    else:
        return proj


def get_tr(proj: str) -> str:
    """Strips the projection name of the chain identifier(s) and returns the progenitor transcript's name"""
    return '#'.join(proj.split('#')[:-1])


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'init_results',
    type=click.File('r'),
    metavar='INIT_ORTH_RESULTS'
)
# @click.argument(
#     'ref_genes',
#     type=click.File('r'),
#     metavar='REF_GENES'
# )
# @click.argument(
#     'query_genes',
#     type=click.File('r'),
#     metavar='QUERY_GENES'
# )
@click.argument(
    'resolved_leaves',
    type=click.File('r'),
    metavar='RESOLVED_LEAVES'
)
@click.option(
    '--output',
    '-o',
    type=click.File('w', lazy=True),
    metavar='OUTPUT_FILE',
    default=sys.stdout,
    show_default=False,
    help='A path to write the updated orthology data to [default: stdout]'
)
@click.option(
    '--one2zero_file',
    '-o2z',
    type=click.File('w', lazy=True),
    metavar='ONE2ZERO_FILE',
    default=None,
    show_default=True,
    help='A path to write the genes rendered 1:0 after the tree reconciliation step'
)
@click.option(
    '--rejection_log',
    type=click.File('w', lazy=True),
    metavar='REJECTED_TRANSCRIPTS',
    default=None,
    show_default=True,
    help='A path to write the rejected transcript items to'
)
@click.option(
    '--loss_summary',
    type=click.File('r', lazy=True),
    metavar='LOSS_SUMMARY_FILE',
    default=None,
    show_default=True,
    help=(
        'TOGA-produced conservation/loss summary file. If provided, loss statuses '
        'from the file will be used in the rejected items report'
    )
)
@click.option(
    '--log_file',
    '-l',
    type=click.Path(exists=False),
    metavar='FILE',
    default=None,
    help='A file to write the execution progress log to'
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
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls the execution verbosity'
)

class FinalOrthologyResolver(CommandLineManager):
    __slots__ = (
        'ref_gene2tr', 'ref_tr2gene', 'query_gene2tr', 'query_tr2gene',
        'tr2proj', 'proj2tr', 'r2q', 'q2r', 'removed_genes', 'one2zero_genes', 
        'rejected_items', 'rejection_log',
        'loss_file', 'proj2loss',
        'out_lines', 'output', 'one2zero_file', 'log_file'
    )

    def __init__(
        self,
        init_results: click.File,
        # ref_genes: click.File,
        # query_genes: click.File,
        resolved_leaves: click.File,
        output: Optional[click.File],
        one2zero_file: Optional[Union[click.File, None]],
        rejection_log: Optional[Union[click.File, None]],
        loss_summary: Optional[Union[click.File, None]],
        log_file: Optional[click.Path],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.log_file: click.Path = log_file
        self.set_logging(log_name)

        self.ref_gene2tr: Dict[str, List[str]] = defaultdict(list) ## TODO: Redundant???
        self.ref_tr2gene: Dict[str, str] = {}
        self.query_gene2tr: Dict[str, List[str]] = defaultdict(list)
        self.query_tr2gene: Dict[str, str] = {}

        self.tr2proj: Dict[str, List[str]] = defaultdict(list)
        self.proj2tr: Dict[str, List[str]] = defaultdict(list)
        self.r2q: Dict[str, Set[str]] = defaultdict(set)
        self.q2r: Dict[str, Set[str]] = defaultdict(set)
        self.removed_genes: Set[str] = set() ## genes fully resolved by gene trees
        self.one2zero_genes: List[str] = [] ## genes rendered one2zero after gene tree reconciliation
        self.out_lines: List[str] = []
        self.one2zero_file: Union[click.File, None] = one2zero_file
        self.rejected_items: List[str] = []
        self.rejection_log: Union[click.File, None] = rejection_log
        self.proj2loss: Dict[str, str] = {}
        self.loss_file: Union[click.File, None] = loss_summary
        self.output: click.File = output

        # self.parse_isoforms(ref_genes, is_ref=True)
        # self.parse_isoforms(query_genes, is_ref=False)
        self.parse_init_results(init_results)
        self.parse_loss_summary()
        self.parse_resolved_leaves(resolved_leaves)
        self.process_remaining_clades()
        self.write_output()

    def parse_isoforms(self, file: TextIO, is_ref: bool) -> None:
        """
        Parses gene-to-isoform two-column file
        """
        for line in file:
            line = line.strip()
            if not line:
                continue
            data: List[str] = line.split('\t')
            if len(data) != 2:
                self._die('ERROR: Isoform file provided is improperly formatted')
            gene: str = data[0]
            tr: str = data[1]
            if is_ref:
                self.ref_gene2tr[gene].append(tr)
                self.ref_tr2gene[tr] = gene
            else:
                self.query_gene2tr[gene].append(tr)
                self.query_tr2gene[tr] = gene

    def parse_init_results(self, file: TextIO) -> None:
        """
        Parses initial orthology resolution results. The logic goes as follows:
        * Lines referring to orthology relationships other than 'many2many' are
          added to the output as-is;
        * Lines referring to 'many2many' cliques are processed
          to save the following information:
          * reference gene-to-query gene mapping;
          * gene-to-transcript mapping for both reference and query;
          * transcript-to-gene mapping for both reference and query
        """
        for line in file:
            if line == Headers.ORTHOLOGY_TABLE_HEADER:
                continue
            line = line.rstrip()
            if not line:
                continue
            data: List[str] = line.split('\t')
            ## TODO: Sanity check for column number
            if data[-1] != MANY2MANY:
                self.out_lines.append(line)
                continue
            ref_gene: str = data[0]
            query_gene: str = data[2]
            self.r2q[ref_gene].add(query_gene)
            self.q2r[query_gene].add(ref_gene)
            ref_tr: str = data[1]
            query_tr: str = data[3]
            self.ref_gene2tr[ref_gene].append(ref_tr)
            self.ref_tr2gene[ref_tr] = ref_gene
            self.query_gene2tr[query_gene].append(query_tr)
            self.query_tr2gene[query_tr] = query_gene

    def parse_loss_summary(self) -> None:
        """
        If loss summary files was provided, 
        extracts query projection loss statuses
        """
        if self.loss_file is None:
            return
        for i, line in enumerate(self.loss_file):
            data: List[str] = line.rstrip().split('\t')
            if not data:
                continue
            if len(data) != 3:
                self._die(
                    (
                        'Wrong data format encountered in the loss summary file at line %i; '
                        'expected 3 fields, got %i'
                    ) % (i + 1, len(data))
                )
            if data[0] != 'PROJECTION':
                continue
            self.proj2loss[data[1]] = data[2]

    def parse_resolved_leaves(self, file: TextIO) -> None:
        """
        Parses a two-column file containing resolved single ortholog pairs
        """
        transcript_pairs: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
        deprecated_cliques: Set[str] = set()
        for line in file:
            line = line.rstrip()
            if not line:
                continue
            data: List[str] = line.split('\t')
            if not data or not data[0]:
                continue
            if data[0] == 'reference':
                continue
            ## TODO: Sanity check for column number
            try:
                ref_tr: str = next(
                    x for x in data if R_PREFIX in x
                ).replace(R_PREFIX, '')
            except StopIteration:
                ref_tr: str = data[0]
            try:
                query_tr: str = next(
                    x for x in data if Q_PREFIX in x
                ).replace(Q_PREFIX, '')
            except StopIteration:
                query_tr: str = data[1]
            query_tr = restore_fragmented_proj_id(query_tr)
            ref_tr = '#'.join(ref_tr.split('#')[:-1])
            ref_gene: str = self.ref_tr2gene[ref_tr]
            query_gene: str = self.query_tr2gene[query_tr]
            ### CURRENT IDEA: First, check whether all the query genes 
            ### actually have projections from the reference orthologs
            ### does not seem trivial since the input does not contain any indications of the original clique
            if query_gene not in self.r2q[ref_gene] and ref_gene not in self.q2r[query_gene]:
                self._to_log(
                    'Ortholog pair resolved in the gene tree has no connection in the original graph: %s and %s' % (ref_tr, query_tr),
                    'warning'
                )
                clique: List[str] = self._restore_original_clique(ref_gene)
                deprecated_cliques = deprecated_cliques.union(clique)
            transcript_pairs[(ref_gene, query_gene)].append((ref_tr, query_tr))
        
        for (ref_gene, query_gene), transcripts in transcript_pairs.items():
            if ref_gene in deprecated_cliques or query_gene in deprecated_cliques:
                continue
            ## remove the edges connecting reference and query genes
            if query_gene in self.r2q[ref_gene]:
                self.r2q[ref_gene].remove(query_gene)
            if ref_gene in self.q2r[query_gene]:
                self.q2r[query_gene].remove(ref_gene)
            ## remove all connections to the reference gene
            for q in self.r2q[ref_gene]:
                self.q2r[q].remove(ref_gene)
            ## remove all connections to the query gene
            for r in self.q2r[query_gene]:
                self.r2q[r].remove(query_gene)
            del self.r2q[ref_gene]
            del self.q2r[query_gene]
            self.removed_genes.add(ref_gene)
            self.removed_genes.add(query_gene)
            for ref_tr, query_tr in transcripts:
                ## CAVEAT: Query transcript used  for gene tree reconstruction 
                ## might have come from a gene other than the newly established ortholog;
                ## pick the one used for the tree reconstruction only if the established 
                ## ortholog has no projection in the query gene
                progenitor_tr: str = get_tr(query_tr)
                if progenitor_tr not in self.ref_tr2gene:
                    self._die('Transcript %s is missing from the reference gene-to-transcript mapping' % progenitor_tr)
                progenitor_gene: str = self.ref_tr2gene[progenitor_tr]
                if progenitor_gene != ref_gene:
                    self._to_log(
                        (
                            'Projection %s was used for gene tree reconstruction according to which the actual '
                            'ortholog of %s is %s; looking for representative transcripts in the original graph'
                        ) % (query_tr, query_gene, ref_gene), 'warning'
                    )
                    recorded_lines: bool = False
                else:
                    # out_line: str = '\t'.join(
                    #     (ref_gene, ref_tr, query_gene, query_tr, ONE2ONE)
                    # )
                    # self.out_lines.append(out_line)
                    recorded_lines: bool = True
                for other_query_tr in self.query_gene2tr[query_gene]:
                    other_ref_tr: str = '#'.join(other_query_tr.split('#')[:-1])
                    # if other_ref_tr == ref_tr and recorded_lines:
                    #     continue
                    ## projections from other genes are counted as rejected
                    if self.ref_tr2gene[other_ref_tr] != ref_gene:
                        self._to_log(f'Skipping {other_query_tr} since it does not belong to the original reference gene')
                        loss_status: str = self.proj2loss.get(other_query_tr, 'N')
                        rej_line: str = ORTH_REJ_TEMPLATE.format(other_query_tr, loss_status)
                        self.rejected_items.append(rej_line)
                        continue
                    out_line: str = '\t'.join(
                        (ref_gene, other_ref_tr, query_gene, other_query_tr, ONE2ONE)
                    )
                    self.out_lines.append(out_line)
                    recorded_lines = True
                if not recorded_lines:
                    self._die(
                        (
                            'No transcripts from gene %s were projected to query locus %s; '
                            'transcript %s will be used for query gene annotation instead'
                        ) % (ref_gene, query_gene, query_tr)
                    )
                    # self._to_log(
                    #     (
                    #         'No transcripts from gene %s were projected to query locus %s; '
                    #         'transcript %s will be used for query gene annotation instead'
                    #     ) % (ref_gene, query_gene, query_tr), 'warning'
                    # )
                    # mod_query_tr: str = f'[{ref_gene}]{query_tr}'
                    # out_line: str = '\t'.join(
                    #     (ref_gene, ref_tr, query_gene, mod_query_tr, ONE2ONE)
                    # )
                    # self.out_lines.append(out_line)
            # self.ref_gene2tr[ref_gene].remove(ref_tr)
            # self.query_gene2tr[query_gene].remove(query_tr)

    def process_remaining_clades(self) -> None:
        """
        Adds the remaining components of (partially) resolved orthology cliques
        to the output
        """
        for ref_gene, query_genes in self.r2q.items():
            ## sanity check
            if ref_gene in self.removed_genes:
                continue
            query_genes = {
                x for x in query_genes if x not in self.removed_genes
            }
            ref_genes: Set[str] = {
                x for q in query_genes for x in self.q2r[q]
                if x not in self.removed_genes
            }
            ref_gene_num: int = len(ref_genes)
            other_query_genes: Set[str] = {
                x for r in ref_genes for x in self.r2q[r]
                if x not in self.removed_genes
            }
            query_gene_num: int = len(other_query_genes)
            # if any(not x for x in (ref_gene_num, query_gene_num)):
            if not query_gene_num:
                self._to_log(
                    'Reference gene %s was reduced to ONE2ZERO after the tree-based resolution step' % ref_gene,
                    'warning'
                )
                self.one2zero_genes.append(ref_gene) if ref_gene not in self.one2zero_genes else None
                for ref_tr in self.ref_gene2tr[ref_gene]:
                    out_line: str = '\t'.join(
                            (ref_gene, ref_tr, 'None', 'None', ONE2ZERO)
                        )
                    self.out_lines.append(out_line)
            # if any(not x for x in (ref_gene_num, query_gene_num)):
            # if not ref_gene_num:
            #     self._die(
            #         'ERROR: Certain clades were reduced to ONE2ONE but missing '
            #         'from the resolved pairs file: %s|%s' % (ref_gene, ','.join(query_genes)) 
            #     )
            status: str = (
                ONE2MANY if ref_gene_num == 1 else 
                MANY2ONE if query_gene_num == 1 else 
                ONE2ZERO if not query_gene_num else 
                MANY2MANY
            )
            for query_gene in query_genes:
                query_trs: List[str] = self.query_gene2tr[query_gene]
                for query_tr in query_trs:
                    ref_tr: str = '#'.join(query_tr.split('#')[:-1])
                    if self.ref_tr2gene[ref_tr] != ref_gene:
                        continue
                    if status == ONE2ZERO:
                        query_gene = 'None'
                        query_tr = 'None'
                    out_line: str = '\t'.join(
                        (ref_gene, ref_tr, query_gene, query_tr, status)
                    )
                    self.out_lines.append(out_line)

    def write_output(self) -> None:
        """
        Writes the updated orthology inference results to the output file
        """
        self.output.write(Headers.ORTHOLOGY_TABLE_HEADER + '\n')
        for line in self.out_lines:
            self.output.write(line + '\n')
        if self.one2zero_file is not None and self.one2zero_genes:
            for gene in self.one2zero_genes:
                self.one2zero_file.write(gene + '\n')
        if self.rejection_log is not None and self.rejected_items:
            for line in self.rejected_items:
                self.rejection_log.write(line + '\n')

    def _restore_original_clique(self, start: str) -> List[str]:
        """
        **
        :start: is the gene name to start the search with. Must refer to a reference gene
        """
        queue: List[str] = [start]
        visited: List[str] = []
        while queue:
            node: str = queue.pop()
            visited.append(node)
            if node in self.r2q:
                neighbours: List[str] = self.r2q[node]
            elif node in self.q2r:
                neighbours: List[str] = self.q2r[node]
            else:
                self._die('Gene %s not found in either reference or query gene list', node)
            for neighbour in neighbours:
                if neighbour in visited:
                    continue
                queue.append(neighbour)
        return visited

if __name__ == '__main__':
    FinalOrthologyResolver()
