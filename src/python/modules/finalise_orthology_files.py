#!/usr/bin/env python3

"""
Assign the names based on referene gene IDs to query genes, finalising the orthology step 
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from collections import defaultdict
from modules.constants import Headers
from modules.shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Dict, List, Optional, Set, TextIO, Tuple, Union

import click

MANY2MANY: str = 'many2many'
NONE: str = 'None'
ONE2MANY: str = 'one2many'
T_GENE: str = 't_gene'

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'orthology_results',
    type=click.File('r', lazy=True),
    metavar='ORTHOLOGY_CLASSIFICATION'
)
@click.argument(
    'query_gene_file',
    type=click.File('r', lazy=True),
    metavar='QUERY_GENE_TABLE'
)
@click.argument(
    'output_dir',
    type=click.Path(exists=False),
    metavar='OUTPUT_DIR',
    default='orthology_output'
)
@click.option(
    '--query_gene_bed_file',
    '-qb',
    type=click.File('r', lazy=True),
    metavar='QUERY_GENE_BED_FILE',
    help='A path to query gene BED file'
)
@click.option(
    '--log_file',
    '-l',
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


class QueryGeneNamer(CommandLineManager):
    """
    Renames the query genes from the reg_{x} notation after the reference genes they correspond to:\n
        * Genes with a sole orthologous gene in the query gain its name after it;\n
        * Genes with up to three orthologs in the reference get a composite name;\n
        * Genes with more than three orthologs are named after 
        the highest scoring orthology with plus symbol at the end\n
        Arguments are:\n
        * ORTHOLOGY_CLASSIFICATION is the final orthology file created by TOGA2 
        (orthology_classification.tsv by default);\n
        * QUERY_GENE_TABLE is a query gene-to-transcript (projection) mapping file created by TOGA2 
        (query_genes.tsv by default);\n
        * OUTPUT_DIR is an output directory to store the results in
    """
    __slots__: Tuple[str] = (
        'log_file', 'orthology_file', 'query_gene_table', 
        'query_gene_bed', 'gene2new_name'
    )

    def __init__(
        self,
        orthology_results: click.File,
        query_gene_file: click.File,
        output_dir: click.Path,
        query_gene_bed_file: Optional[click.File],
        log_file: Optional[click.Path],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.log_file: Union[str, None] = log_file
        self.set_logging(log_name)

        self._to_log('Creating output directory for orthology wrap-up')
        self._mkdir(output_dir)
        self.orthology_file: str = os.path.join(output_dir, 'orthology_classification.tsv')
        self.query_gene_table: str = os.path.join(output_dir, 'query_genes.tsv')
        self.query_gene_bed: str = os.path.join(output_dir, 'query_genes.bed')

        self._to_log(
            'Inferring orthology-based query gene names '
            'and adding them to the orthology classification file'
        )
        self.gene2new_name: Dict[str, str] = {}
        self.modify_orthology_classification(orthology_results)
        self._to_log('Renaming query genes in the gene-to-transcript mapping file')
        self.modify_gene_mapping(query_gene_file)
        if query_gene_bed_file:
            self._to_log('Renaming query genes in the gene BED file')
            self.modify_gene_bed(query_gene_bed_file)

    def modify_orthology_classification(self, file: TextIO) -> None:
        """
        Parses the original orthology file, inferres query gene names, 
        and writes the updated orthology file 
        """
        one2many_naming: Dict[str, int] = defaultdict(int)
        many2many_naming: Dict[str, int] = defaultdict(int)
        query_gene2orth_data: Dict[str, List[str]] = defaultdict(list)
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 5:
                self._die(
                    'Line %i in the orthology classification file has less than 5 fields' % i
                )
            if data[0] == T_GENE:
                continue
            query_gene: str = data[2]
            query_gene2orth_data[query_gene].append(data)
        with open(self.orthology_file, 'w') as h:
            h.write(Headers.ORTHOLOGY_TABLE_HEADER)
            for gene, lines in query_gene2orth_data.items():
                if gene == NONE:
                    new_query_name: str = NONE
                else:
                    ref_gene_names: List[str] = []
                    lines = sorted(lines, key=lambda x: int(x[3].split('#')[-1].split(',')[0]))
                    for line in lines:
                        if line[0] not in ref_gene_names:
                            ref_gene_names.append(line[0])
                    status: str = lines[0][4]
                    # ref_gene_names: Set[str] = {x[0] for x in lines}
                    upd_ref_gene_names: List[str] = []
                    for ref_gene in ref_gene_names:
                        if status == ONE2MANY:
                            one2many_naming[ref_gene] += 1
                            ref_gene += f'_{chr(96 + one2many_naming[ref_gene])}'
                        elif status == MANY2MANY:
                            many2many_naming[ref_gene] += 1
                            ref_gene += f'_{many2many_naming[ref_gene]}'
                        upd_ref_gene_names.append(ref_gene)
                    if len(upd_ref_gene_names) == 1:
                        new_query_name: str = upd_ref_gene_names.pop()
                    elif len(upd_ref_gene_names) <= 3:
                        # ref_gene_names = []
                        # for line in sorted(lines, key=lambda x: int(x[3].split('#')[-1].split(',')[0])):
                        #     if line[0] in ref_gene_names:
                        #         continue
                        #     ref_gene_names.append(line[0])
                        new_query_name: str = ','.join(upd_ref_gene_names)
                    else:
                        # ref_gene_names = [
                        #     x[0] for x in sorted(lines, key=lambda x: int(x[3].split('#')[-1].split(',')[0]))
                        # ]
                        new_query_name: str = upd_ref_gene_names[0] + '+'
                    # print(f'{status=}')
                    # if status == ONE2MANY:
                    #     one2many_naming[new_query_name] += 1
                    #     new_query_name += f'_{chr(96 + one2many_naming[new_query_name])}'
                    # elif status == MANY2MANY:
                    #     many2many_naming[new_query_name] += 1
                    #     new_query_name += f'_{many2many_naming[new_query_name]}'
                    self.gene2new_name[gene] = new_query_name
                for line in lines:
                    line[2] = new_query_name
                    h.write('\t'.join(line) + '\n')

    def modify_gene_mapping(self, file: TextIO) -> None:
        """
        Modifies the query gene-to-transcript mapping file, substituting reg_{x} 
        symbols with the corresponding orthology-based names
        """
        with open(self.query_gene_table, 'w') as h:
            h.write(Headers.QUERY_GENE_HEADER)
            for i, line in enumerate(file, start=1):
                data: List[str] = line.rstrip().split('\t')
                if not data or not data[0]:
                    continue
                if len(data) < 2:
                    self._die(
                        'Line %i in the query gene mapping file has less than two columns' %i
                    )
                old_name: str = data[0]
                if old_name not in self.gene2new_name:
                    self._to_log(
                        'Gene %s at line %i in the gene table has no proven orthologs in the reference' % (
                            old_name, i
                        ),
                        'warning'
                    )
                    gene_name: str = old_name
                else:
                    gene_name: str = self.gene2new_name[old_name]
                data[0] = gene_name
                h.write('\t'.join(data) + '\n')

    def modify_gene_bed(self, file: TextIO) -> None:
        """
        Modifies the query gene BED file, substituting reg_{x} 
        symbols with the corresponding orthology-based names
        """
        with open(self.query_gene_bed, 'w') as h:
            for i, line in enumerate(file, start=1):
                data: List[str] = line.rstrip().split('\t')
                if not data or not data[0]:
                    continue
                if len(data) < 4:
                    self._die(
                        'Line %i in the query gene BED file has less than four columns' %i
                    )
                old_name: str = data[3]
                if old_name not in self.gene2new_name:
                    self._to_log(
                        'Gene %s at line %i in the gene table has no proven orthologs in the reference' % (
                            old_name, i
                        ),
                        'warning'
                    )
                    gene_name: str = old_name
                else:
                    gene_name: str = self.gene2new_name[old_name]
                data[3] = gene_name
                h.write('\t'.join(data) + '\n')


if __name__ == '__main__':
    QueryGeneNamer()