#!/usr/bin/env python3

"""
Creates a summary table for the gene tree-based orthology refinement step
"""

from constants import Headers
from shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Dict, List, Optional, Union, TextIO

import click
import os

RESOLVED_LEAVES: str = 'resolved_pairs.tsv'
IQTREE_PATH_COMPONENTS: str = os.path.join('tmp', 'tree', '{}_tmp', '{}.log')
RAXML_PATH_COMPONENTS: str = os.path.join('tmp', 'tree', '{}_tmp', 'RAxML_info.{}')

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'input_dir',
    type=click.Path(exists=True),
    metavar='INPUT_DIR'
)
@click.argument(
    'results_dir',
    type=click.Path(exists=True),
    metavar='RESULT_DIR'
)
@click.argument(
    'output',
    type=click.File('w', lazy=True),
    metavar='OUTPUT'
)
@click.option(
    '--raxml',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, assumes output directory to comply with RAxML (raxmlHPC-PTHREADS-AVX) output structure; expects '
        'IqTree2 output otherwise'
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
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)


class GeneTreeSummary(CommandLineManager):

    @staticmethod
    def line_num(file: TextIO, start_pattern: str = '') -> int:
        start_byte: bytes = start_pattern.encode('utf8')
        counter: int = 0
        with open(file, 'rb') as h:
            for line in h:
                if not line.startswith(start_byte):
                    continue
                counter += 1
        return counter

    @staticmethod
    def get_clique_dict_stub(file: str) -> Dict[str, int]:
        """
        Parses a single-column gene tree batch input file, 
        returns a {line:0} dictionary stub
        """
        output: Dict[str, int] = {}
        with open(file, 'r') as h:
            for line in h:
                line = line.rstrip()
                if not line:
                    continue
                output[line] = 0
        return output

    @staticmethod
    def get_pairs_per_clique(file: str) -> Dict[str, int]:
        output: Dict[str, int] = {}
        with open(file, 'r') as h:
            for line in h:
                data: List[str] = line.rstrip().split('\t')
                source: str = data[2]
                if source in output:
                    output[source] += 1
                else:
                    output[source] = 1
        return output

    @staticmethod
    def get_model_raxml(file: str) -> str:
        """Retrieves the best-fit model from a RAxML log"""
        model: str = 'NA'
        likelihood: float = 0.0
        with open(file, 'r') as h:
            for line in h:
                line = line.rstrip()
                if not line:
                    continue
                if 'best-scoring AA model:' not in line:
                    continue
                data: List[str] = line.split(': ')[-1].split(' ')
                ml: str=  float(data[2])
                if ml <= likelihood:
                    model = data[0]
                    likelihood = ml
        return model

    @staticmethod
    def get_model_iqtree(file: str) -> str:
        """Retrieves the best-fit model from an IQTree2 log"""
        model: str = 'NA'
        with open(file, 'r') as h:
            for line in h:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('Best-fit model:'):
                    model: str = line.split(' ')[2]
                    break
        return model

    @staticmethod
    def get_clique_name(file: str) -> str:
        """Strips clique name from the filename"""
        filename: str = file.split(os.sep)[-1]
        clique: str = filename.split('_')[1].split('.')[0]
        return clique

    def __init__(
        self,
        input_dir: click.Path,
        results_dir: click.Path,
        output: click.File,
        raxml: Optional[bool],
        log_name: Optional[Union[str, None]],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)

        output.write(Headers.TREE_SUMMARY_HEADER)
        batches: List[str] = [x for x in os.listdir(results_dir) if x.startswith('batch')]
        batches.sort(
            key=lambda x: int(x.replace('batch', ''))
        )
        # clique_inputs: List[str] = [
        #     x for x in os.listdir(input_dir) if x.startswith('batch') and  x.endswith('.fa')
        # ]
        for batch in batches:
            resolved_leaves_path: str = os.path.join(results_dir, batch, RESOLVED_LEAVES)
            if not os.path.exists(resolved_leaves_path):
                self._to_log(
                    'No resolved pair file found for batch %s' % batch,
                    'warning'
                )
                continue
            config_path: str = os.path.join(input_dir, f'{batch}.txt')
            clique2pairs: Dict[str, int] = self.get_clique_dict_stub(config_path)
            # clique2pairs: Dict[str, int] = self.get_pairs_per_clique(resolved_leaves_path)
            clique2pairs.update(self.get_pairs_per_clique(resolved_leaves_path))
            # clique_files: List[str] = [
            #     x for x in clique_inputs if x.startswith(f'{batch}_')
            # ]
            for clique_file, pair_num in clique2pairs.items():
                clique_pref: str = self.get_clique_name(clique_file)
                print(f'{clique_file=}, {clique_pref=}, {pair_num=}')
                # batch_num: str = batch.replace('batch', '')
                # input_file: str = f'{batch}_clique{batch_num}.fa'
                # input_file: str = f'{batch}_clique{batch_num}.fa'
                # input_path: str = os.path.join(input_dir, input_file)
                input_path: str = os.path.join(input_dir, clique_file)
                if not os.path.exists(input_path):
                    self._to_log(
                        'No input FASTA file found for batch %s' % batch,
                        'warning'
                    )
                    continue
                seq_num: int = self.line_num(input_path, start_pattern='>')
                if raxml:
                    log_path: str = os.path.join(
                        results_dir, batch, RAXML_PATH_COMPONENTS.format(clique_pref, clique_pref)
                    )
                    model: str = self.get_model_raxml(log_path)
                else:
                    log_path: str = os.path.join(
                            results_dir, batch, IQTREE_PATH_COMPONENTS.format(clique_pref, clique_pref)
                        )
                    model: str = self.get_model_iqtree(log_path)
                # pair_num: int = self.line_num(resolved_leaves_path)
                # pair_num: int = resolved_leaves_path[clique_pref]
                out_line: str = '\t'.join(map(str, (batch, clique_pref, seq_num, pair_num, model)))
                output.write(out_line + '\n')


if __name__ == '__main__':
    GeneTreeSummary()