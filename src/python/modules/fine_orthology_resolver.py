#!/usr/bin/env python3

"""
Resolves many2many orthologies with PRANK+RAxML pipeline
"""

import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from Bio import Phylo
from collections import defaultdict
from constants import IQTREE_ACCEPTED_MODELS, PHYLO_NOT_FOUND
from shared import CommandLineManager, CONTEXT_SETTINGS
from shutil import which
from tree_analysis import can_resolve, make_cat_tree
from typing import Dict, List, Optional, Tuple, Union

import click

__author__ = ('Amy Stephen', 'Yury V. Malovichko')
__year__ = '2024'

R_PREFIX: str = "#R#"
Q_PREFIX: str = "#Q#"
IQTREE_R_PREFIX: str = '_R_'
IQTREE_Q_PREFIX: str = '_Q_'
PARENT_REF: str = '..'
BOOTSTRAPPED_TREE: str = 'RAxML_bipartitionsBranchLabels'
IQTREE_CMD_STUB: str = (
    '{} -s {} --seqtype AA -T AUTO --threads-max {} --seed 12345 --subsample-seed 12345   --prefix {} --alrt {} -B {} --keep-ident --mset {} -redo'
)
RAXML_CMD_STUB: str = '{} -T {} -f a -s {} -n {} -w {} -m PROTGAMMAAUTO -x 12345 -p 12345 -# {}'

def format_abspath(path: str) -> str:
    """
    Formats absolute paths into RAxML-friendly format
    """
    path_components: List[str] = path.split(os.sep)
    first_ref_to_parent: bool = -1
    parent_ref_num: int = 0
    for i, comp in enumerate(path_components):
        if comp == PARENT_REF:
            first_ref_to_parent = i if first_ref_to_parent < 0 else first_ref_to_parent
            parent_ref_num += 1
    if first_ref_to_parent < 0:
        return path
    last_pre_ref: int = first_ref_to_parent - parent_ref_num
    first_post_ref: int = first_ref_to_parent + parent_ref_num + 1
    filt_components: List[str] = path_components[:last_pre_ref] + path_components[first_post_ref:]
    return os.sep.join(filt_components)


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'input',
    type=click.Path(exists=True),
    metavar='INPUT'
)
@click.argument(
    'output',
    type=click.Path(exists=False),
    metavar='OUTPUT'
)
@click.option(
    '--table_input',
    '-t',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, expects INPUT to be a single-column text file containing paths '
        'to cliquewise FASTA files'
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
    '--tree_binary',
    '-rb',
    type=click.Path(exists=True),
    metavar='TREE_BINARY',
    default=None,
    show_default=True,
    help=(
        'A path to the RAxML/IQTree2 executable to be used at fine resolution step. '
        'If not provided, the program will try to infer its location from the PATH'
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
    help='A number of bootstrap replications for RAxML/IQTree2 run'
)
@click.option(
    '--bootstrap_threshold',
    '-bt',
    type=click.FloatRange(min=0.0, max=100.0),
    metavar='FLOAT',
    default=90.0,
    help='Minimum bootstrap support value for a tree clade to be considered resolvable'
)
@click.option(
    '--tmp_dir',
    '-d',
    type=click.Path(exists=False),
    metavar='PATH',
    default=None,
    show_default=False,
    help='A temporary directory for PRANK and RAxML results [default: OUTPUT/tmp]'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)

class FineOrthologyResolver(CommandLineManager):
    """
    A module for tree-based complex orthology resolution. Given a FASTA file
    or a list of those, runs the PRANK+RAxML pipeline, roots the resulting
    phylogenetic tree, and resolves terminal clades corresponding to orthologous
    pairs between reference and query species.\n
    Arguments are:\n
    * INPUT is either a FASTA file (default input format) or a single-column
    text file contaiing paths to one or more FASTA files
    (expected if --table_input flag is enabled);\n
    * OUTPUT is path to directory containing orthology resolution results.
    Resolved orthologous pairs are saved in the "resolved_pairs.tsv" file, and
    unresolved homologs are written to "unresolved_clades.txt". Note that, depending
    on the orthology resolution results, either of the files could be missing from
    the output directory.
    """
    __slots__ = (
        'input_files', 'output', 'tree_cpus', 'tree_bootnum',
        'prank_bin', 'use_raxml', 'tree_bin', 'boot_thresh', 'tmp',
        'original_names', 'prank_dir', 'tree_dir',
        'resolved_leaves', 'unresolved_clades',
        'resolved_leaves_file', 'unresolved_clades_file'
    )

    ## TODO: Add the initial orthology relationships as an input file
    @staticmethod
    def get_clique_name(file: str) -> str:
        """Returns the clique name for the given cliquewise input file"""
        filename: str = file.split(os.sep)[-1]
        clique: str = filename.split('_')[1].split('.')[0]
        return clique

    def __init__(
        self,
        input: click.Path,
        output: click.Path,
        table_input: Optional[bool],
        prank_binary: Optional[Union[click.Path, None]],
        use_raxml: Optional[bool],
        tree_binary: Optional[Union[click.Path, None]],
        tree_cpus: Optional[int],
        tree_bootstrap: Optional[int],
        bootstrap_threshold: Optional[float],
        tmp_dir: Optional[Union[click.Path, None]],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging()

        self.input_files: List[str] = []
        if table_input:
            self.parse_input_table(input)
        else:
            self.input_files.append(input)
        self.output: click.Path = output
        self.tree_cpus: int = tree_cpus
        self.tree_bootnum: int = tree_bootstrap
        self.prank_bin: str = prank_binary
        self.use_raxml: bool = use_raxml
        self.tree_bin: str = tree_binary
        self.boot_thresh: float = bootstrap_threshold
        self.tmp: str = (
            tmp_dir if tmp_dir is not None else os.path.join(self.output, 'tmp')
        )

        self.original_names: Dict[str, str] = {}

        self.prank_dir: str = os.path.join(self.tmp, 'prank')
        self.tree_dir: str = os.path.join(self.tmp, 'tree')

        self.resolved_leaves: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
        self.unresolved_clades: List[List[str]] = []

        self.resolved_leaves_file: str = os.path.join(
            self.output, 'resolved_pairs.tsv'
        )
        self.unresolved_clades_file: str = os.path.join(
            self.output, 'unresolved_clades.txt'
        )

        self.run()

    def run(self) -> None:
        """
        Major executor method
        """
        self.check_executables()
        self._mkdir(self.output)
        self._mkdir(self.tmp)
        self._mkdir(self.prank_dir)
        self._mkdir(self.tree_dir)
        """
        For each FASTA file in input:
        * run PRANK
        * run RAXML
        * midpoint-root the tree
        * the structure
        """
        for c, file in enumerate(self.input_files):
            if not self.use_raxml:
                self.save_fasta_names(file)
            tree_pref: str = self.get_clique_name(file)
            prank_pref: str = os.path.join(self.prank_dir, f'{tree_pref}.msa')
            prank_out: str = f'{prank_pref}.best.fas'
            tree_tmp: str = os.path.join(self.tree_dir, f'{tree_pref}_tmp')
            self._mkdir(tree_tmp)
            self.run_prank(file, prank_pref)
            pref_for_tree: str = tree_pref if self.use_raxml else os.path.join(tree_tmp, tree_pref)
            self.run_tree(prank_out, pref_for_tree, tree_tmp)
            tree_out: str = (
                # os.path.join(
                #     tree_tmp, f'{tree_pref}.contree'
                # ) 
                f'{BOOTSTRAPPED_TREE}.{tree_pref}' if self.use_raxml else f'{tree_pref}.contree' 
            )
            tree_out_path: str = os.path.join(
                tree_tmp, tree_out
            )
            rooted_tree: str = os.path.join(self.tree_dir, f'{tree_pref}_rooted.nwk')
            self.midpoint_root(tree_out_path, rooted_tree)
            self.resolve_tree_orthology(rooted_tree, self.output, file) ## TODO: Format output
        self.write_output()

    def check_executables(self) -> None:
        """
        Checks whether valid PRANK and RAxML executables are available
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

    def parse_input_table(self, table: str) -> None:
        """
        """
        with open(table, 'r') as h:
            for line in h:
                path: str = line.strip()
                if not path:
                    continue
                if not os.path.exists(path):
                    self._die(
                        f'ERROR: path {path} listed in the input file {table} '
                        'does not exist. Please check that all paths in the '
                        'input table lead to valid FASTA files'
                    )
                self.input_files.append(path)

    def mktmp(self) -> str:
        """
        Creates a temporary directory for PRANK runs
        """
        # prank -d=MS_unaligned.fa -o=output_aligned_AA.fasta
        cmd: str = 'mktemp -d TEMPtreeResolveXXXXX'

        # run prank time
        stdout: str = self._exec(
            cmd, 'ERROR: temporary directory creation failed'
        )
        return stdout.strip()

    def save_fasta_names(self, file: str) -> None:
        """
        Extracts sequence names from the input Fasta file and saves them to 
        restore the original naming when parsing IQTree2 results
        """
        with open(file, 'r') as h:
            for line in h:
                line = line.rstrip()
                if line[0] != '>':
                    continue
                orig_name = line[1:]
                iqtree_compatible: str = orig_name.replace('#', '_').replace(',', '_')
                self.original_names[iqtree_compatible] = orig_name

    def run_prank(self, in_file: str, out_pref: str) -> None:
        """
        Runs PRANK alignment command
        """
        cmd: str = f'{self.prank_bin} -d={in_file} -o={out_pref} -protein '#-nomafft'
        _ = self._exec(cmd, 'ERROR: PRANK alignment failed with the following error:')

    def run_tree(self, aln_file: str, tree_file: str, tmp_dir: str) -> None:
        """
        """
        abs_tmp_dir: str = os.path.abspath(tmp_dir)
        formatted_dir: str = format_abspath(abs_tmp_dir)
        if not self.use_raxml:
            cmd: str = IQTREE_CMD_STUB.format(
                self.tree_bin, aln_file, self.tree_cpus, tree_file, 
                self.tree_bootnum, self.tree_bootnum, IQTREE_ACCEPTED_MODELS
            )
        else:
            cmd: str = RAXML_CMD_STUB.format(
                self.tree_bin, self.tree_cpus, aln_file, 
                tree_file, formatted_dir, self.tree_bootnum
            )
        # cmd: str = (
        #     f'{self.tree_bin} -T {self.tree_cpus} -f a '
        #     # f'{self.tree_bin} --all --threads {self.tree_cpus}'
        #     f'-s {aln_file} -n {tree_file} -w {formatted_dir} '
        #     f'-m PROTGAMMAAUTO -x 12345 -p 12345 -# {self.tree_bootnum}'
        # )
        print(f'{cmd=}')
        _ = self._exec(cmd, 'ERROR: RAxML run failed with the following error:')

    def midpoint_root(self, unroot_path: str, root_path: str) -> None:
        """
        Roots the RAxML-resulting tree by the midpoint
        """
        with open(unroot_path, 'r') as ur, open(root_path, 'w') as r:
            tree = Phylo.read(ur, 'newick')
            tree.root_at_midpoint()
            Phylo.write(tree, r, 'newick')

    def resolve_tree_orthology(self, tree_path: str, out_path: str, source: str) -> None:
        """
        Resolves orthology relationships by the rooted phylogenetic tree
        """
        with open(tree_path, 'r') as h:
            tree: Phylo.BaseTree.Tree = Phylo.read(h, 'newick')
        make_cat_tree(tree)
        resolved_leaves: List[Tuple[str]] = can_resolve(tree, self.boot_thresh, not self.use_raxml)
        if resolved_leaves:
            self.resolved_leaves[source].extend(resolved_leaves)
            # self.resolved_leaves.extend(resolved_leaves)
        resolved_flattened: List[str] = [x for y in resolved_leaves for x in y]
        unresolved_leaves: List[str] = [
            x.name for x in tree.get_terminals()
            if x.name not in resolved_flattened
        ]
        if unresolved_leaves:
            self.unresolved_clades.append(unresolved_leaves)

    def write_output(self) -> None:
        """
        Writes resolution results to files in the output directory:
        * Resolved one2one pairs are written to a two-column file
          in the {ref}\t{query} format;
        * Unresolved (sub)cliques are written in a single-column file
          as comma-separated lists
        """
        with open(self.resolved_leaves_file, 'w') as h:
            for source, results in self.resolved_leaves.items():
                for orth_pair in results:
                    if not self.use_raxml:
                        r_: str = next(
                            x for x in orth_pair if IQTREE_R_PREFIX in x
                        )#.replace(IQTREE_R_PREFIX, R_PREFIX)
                        q_: str = next(
                            x for x in orth_pair if IQTREE_Q_PREFIX in x
                        )#.replace(IQTREE_Q_PREFIX, Q_PREFIX)
                        # r_components: List[str] = r_.split('_')
                        # r: str = '_'.join(r_components[:-1]) + '#' + r_components[-1]
                        # q_components: List[str] = q_.split('_')
                        # q: str = '_'.join(q_components[:-1]) + '#' + q_components[-1]
                        r: str = self.original_names[r_]
                        q: str = self.original_names[q_]
                    else:
                        r: str = next(x for x in orth_pair if R_PREFIX in x)
                        q: str = next(x for x in orth_pair if Q_PREFIX in x)
                    h.write(f'{r}\t{q}\t{source}\n')
        with open(self.unresolved_clades_file, 'w') as h:
            for clade in self.unresolved_clades:
                ## IQTree2 replaces sharps with underscores; restore the original naming first
                if not self.use_raxml:
                    out_line: str = ''
                    last_el: int = len(clade) - 1
                    for i, el in enumerate(clade):
                        # if IQTREE_R_PREFIX in el:
                        #     el = el.replace(IQTREE_R_PREFIX, R_PREFIX)
                        # elif IQTREE_Q_PREFIX in el:
                        #     el = el.replace(IQTREE_Q_PREFIX, Q_PREFIX)
                        # components: List[str] = el.split('_')
                        # tr: str = '_'.join(components[:-1])
                        # el = f'{tr}#{components[-1]}'
                        el = self.original_names[el]
                        out_line += el
                        if i < last_el:
                            out_line += ','
                else:
                    out_line = ','.join(clade)
                h.write(out_line + '\n')


if __name__ == '__main__':
    FineOrthologyResolver()
