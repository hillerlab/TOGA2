#!/usr/bin/env python3

"""
Given a projection feature table, classifies projections in terms of their orthology
"""

## TODO: Move to 'modules'?

from collections import Counter, defaultdict
from modules.shared import CONTEXT_SETTINGS, CommandLineManager
from typing import Any, Dict, Iterable, List, Optional, Set

import click
import joblib
import os
import numpy as np
import pandas as pd
import xgboost as xgb

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__credits__ = ['Bogdan M. Kirilenko']
__all__ = [None]

class Constants:
    # __slots__ = (
    #     '', 'FINAL_COLUMNS',
    #     'ORTH', 'PARA', 'SPAN', 'P_PGENE',
    #     'TR2CHAIN_HEADER', 'UNCLASS_TEMPLATE'
    # )
    SE_MODEL_FEATURES: List[str] = [
        'gl_exo', 'flank_cov', 'exon_perc', 'synt_log'
    ]
    ME_MODEL_FEATURES: List[str] = [
        'gl_exo', 'loc_exo', 'flank_cov', 'synt_log', 'intr_perc'
    ]
    LD_MODEL_FEATURES: List[str] = [
        'gl_exo', 'flank_cov', 'exon_perc', 'synt_log', 'loc_exo',
        'intr_perc', 'score', 'single_exon'
    ]
    PP_FEATURES: List[str] = [
        'clipped_exon_qlen', 'clipped_intr_cover'
    ]
    PP_CLIPPED_EXON_QLEN: float = 0.3
    PP_CLIPPED_INTRON_QLEN: float = 0.1
    ORTH: str = 'ORTH'
    PARA: str = 'PARA'
    SPAN: str = 'SPAN'
    P_PGENE: str = 'P_PGENE'
    TR2CHAIN_HEADER: str = 'TRANSCRIPT\tORTH\tPARA\tSPAN\tP_PGENE'
    FINAL_COLUMNS: List[str] = ['transcript', 'chain', 'pred']
    UNCLASS_TEMPLATE: str = (
        'TRANSCRIPT\t{}\t0\tNo classifiable projections found\tNO_PROJ\tM'
    )
    UNDERSCORED_TEMPLATE: str = (
        'PROJECTION\t{}\t0\tChain score below set threshold ({})\tINSUFFICIENT_CHAIN_SCORE\tL'
    )


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'feature_table',
    type=click.Path(exists=True),
    metavar='FEATURE_TABLE'
)
@click.argument(
    'output_dir',
    type=click.Path(exists=False),
    metavar='OUTPUT_DIR'
)
@click.argument(
    'single_exon_model',
    type=click.Path(exists=True),
    metavar='SINGLE_EXON_MODEL'
)
@click.argument(
    'multi_exon_model',
    type=click.Path(exists=True),
    metavar='MULTI_EXON_MODEL'
)
@click.option(
    '--orthology_threshold',
    '-t',
    type=click.FloatRange(min=0.0, max=1.0),
    metavar='FLOAT',
    default=0.5,
    show_default=True,
    help='Probability threshold for classifying projections as orthologous'
)
@click.option(
    '--long_distance_model',
    '-ld',
    type=click.Path(exists=True),
    metavar='LONG_DISTANCE_MODEL',
    default=None,
    show_default=True,
    help=(
        'If set, applies extra classifier for distantly related species. '
        'By default, expects TOGA long-distance model relevant '
        'at molecular distances ~1sps'
    )
)
@click.option(
    '--min_orthologous_chain_score',
    '-minscore',
    type=click.IntRange(min=0, max=None),
    metavar='INT',
    default=0,
    show_default=True,
    help=(
        'Minimal score for chains to be potentially classified as orthologous. Chains with '
        'score less than that are discarded unless they are classified as retrogenes/processed pseudogenes'
    )
)
@click.option(
    '--legacy',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, expects the results of the legacy (TOGA1) feature extraction step, '
        'and applies legacy spanning chain definition'
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
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)

class ChainClassifier(CommandLineManager):
    __slots__ = [
        'feature_table', 'output', 'se_model', 'me_model', 'orthology_threshold',
        'ld_model', 'min_orth_chain_score', 'legacy', 'v',
        'orthology_class_table', 'tr2chain_classes', 'rejection_log', 
        'df', 'underscored_chain_projections'
    ]

    def __init__(
        self,
        feature_table: click.Path,
        output_dir: click.Path,
        single_exon_model: click.Path,
        multi_exon_model: click.Path,
        orthology_threshold: float,
        long_distance_model: Optional[click.Path],
        min_orthologous_chain_score: Optional[int],
        legacy: Optional[bool],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        """
        Runs XGBoost classifier for projection orthology status prediction.\n
        This is a ripoff of TOGA 1.0 classification master script.\n
        Arguments are:\n
        * FEATURE_TABLE is a path to projection (transcript-chain pair) feature
        table used for projection classification;\n
        * OUTPUT_DIR is a path to directory to store results in; this directory
        does not have to necessarily exist prior to code execution;\n
        * SINGLE_EXON_MODEL is a path to model used for classification of single
        exon transcripts' projections;\n
        * MULTI_EXON_MODEL is path to multi-exon transcript model counterpart
        """
        self.v: bool = verbose
        self.set_logging(log_name)

        self._to_log('Initializing ChainClassifier')
        self._to_log('Uploading feature dataset')
        self.df: pd.core.frame.DataFrame = pd.read_csv(
            feature_table, header=0, sep='\t'
        )
        self.underscored_chain_projections: List[str] = []
        self.output: click.Path = output_dir
        self._to_log('Uploading classification models')

        self.se_model = self._load_model(single_exon_model)
        self.me_model = self._load_model(multi_exon_model)
        self.ld_model = (
            self._load_model(long_distance_model)
            if long_distance_model is not None else None
        )

        self.orthology_threshold: float = orthology_threshold
        self.min_orth_chain_score: int = min_orthologous_chain_score

        self.orthology_class_table: str = os.path.join(
            self.output, 'orthology_scores.tsv'
        )
        self.tr2chain_classes: str = os.path.join(
            self.output, 'trans_to_chain_classes.tsv'
        )
        self.rejection_log: str = os.path.join(
            self.output, 'rejection_report.tsv'
        )
        self.legacy: bool = legacy

        self.run()

    def _load_model(self, model_path: click.Path) -> Any:
        """
        Uploads XGBoost classification model
        """
        try:
            return joblib.load(model_path)
        except (xgb.core.XGBoostError, AttributeError):
            xgboost_version: str = xgb.__version__
            err_msg: str = (
                f'Cannot load models located at {model_path}. '
                'Probably, models were trained with a different version of '
                f'XGBoost. You used XBGoost version: {xgboost_version}; '
                'Please make sure you called train_model.py with the same version.'
            )
            self._die(err_msg)

    def run(self) -> None:
        """
        Main executing method
        """
        ## create an output directory
        self._mkdir(self.output)

        ## extract unique names
        init_tr_set: Set[str] = set(self.df['transcript'])

        ## get indices of processed pseudogene projections
        ## if a projection has synteny = 1, introns are deleted, and exon number > 1,
        ## this is likely a processed pseudogene

        ## extract spanning projections: chains do not cover any coding exons
        ## but have high synteny due to flanking sequence alignment
        if self.legacy:
            spanning_lines: pd.core.frame.DataFrame = self.df[
                (self.df['exon_cover'] == 0) & (self.df['synt'] > 1)
            ]
        else:
            spanning_lines: pd.core.frame.DataFrame = self.df[self.df['exon_cover'] == 0]

        ## filter the dataframe of spanning and non-syntenic projections;
        ## this dataframe will be used for classification
        self.df = self.df[(self.df['exon_cover'] > 0) & (self.df['synt'] > 0)]

        ## compute the derived features
        self.df["exon_perc"] = self.df["exon_cover"] / self.df["ex_fract"]
        self.df["chain_len_log"] = np.log10(self.df["chain_len"])
        self.df["synt_log"] = np.log10(self.df["synt"])
        self.df["intr_perc"] = self.df["intr_cover"] / self.df["intr_fract"]
        self.df = self.df.fillna(0.0)  ## fill NA values with zeros

        ## split df into two: for single and multi exon models
        df_se: pd.core.frame.DataFrame = self.df[self.df['ex_num'] == 1]
        df_me: pd.core.frame.DataFrame = self.df[self.df['ex_num'] > 1]

        ## extract predictor data for both single- and multi-exon projections
        X_se: pd.core.frame.DataFrame = df_se[Constants.SE_MODEL_FEATURES]
        X_me: pd.core.frame.DataFrame = df_me[Constants.ME_MODEL_FEATURES]
        df_me_pp: pd.core.Frame.DataFrame = df_me[Constants.PP_FEATURES]

        ## run prediction model
        se_pred: Iterable[float] = (
            self.se_model.predict_proba(X_se)[:, 1] if len(X_se) > 0 else np.array([])
        )
        me_pred: Iterable[float] = (
            self.me_model.predict_proba(X_me)[:, 1] if len(X_me) > 0 else np.array([])
        )

        ## add predicted values to initial data frames
        df_se['pred'] = se_pred
        df_me['pred'] = me_pred

        ## if long distance model is specified, convert predictions
        if self.ld_model is not None:
            df_se = df_se.rename({'pred': 'score'}, axis='columns').assign(single_exon = 1)
            df_me = df_me.rename({'pred': 'score'}, axis='columns').assign(single_exon = 0)
            X_se = df_se[Constants.LD_MODEL_FEATURES]
            X_me = df_me[Constants.LD_MODEL_FEATURES]
            se_ld_pred = (
                self.ld_model.predict_proba(X_se)[:, 1] if len(X_se) > 0 else np.array([])
            )
            me_ld_pred = (
                self.ld_model.predict_proba(X_me)[:, 1] if len(X_me) > 0 else np.array([])
            )
            df_se['pred'] = se_ld_pred
            df_me['pred'] = me_ld_pred

        ## assign a probability placeholder of -1 to spanning projections
        spanning_lines['pred'] = -1

        ## identify processed pseudogene projections:
        ## those are multi-exon projections with orthology probability below
        ## the threshold, minimal syntenty, and high exonic fraction;
        ## those get a probability placeholder of -2
        df_me.loc[
            (df_me['synt'] == 1)
            & (df_me['exon_qlen'] > 0.95)
            & (df_me['pred'] < self.orthology_threshold)
            & (df_me['exon_perc'] > 0.65),
            'pred',
        ] = -2
        ## TOGA2 speciial: identify processed pseudogenes 
        ## by CDS-clipped alignment-to-query span 
        ## adn CDS-clipped intron coverage
        all_orthologs: List[str] = df_me[
            df_me['pred'] > self.orthology_threshold
        ]['transcript'].to_list()
        ortholog_counter: Dict[str, int] = Counter(all_orthologs)
        # max_prob_per_multiorth: Dict[str, float] = {
        #     x: max(df_me[df_me['gene'] == x]['pred'].to_list()) 
        #     for x,y in ortholog_counter.items() if y > 1
        # }
        max_prob_per_multiorth: Dict[str, float] = df_me[
            df_me.apply(lambda x: ortholog_counter[x['transcript']] > 1, axis=1)
        ].groupby(
            ['transcript']
        )['pred'].max().to_dict()
        # print(max_prob_per_multiorth)
        df_me.loc[
            (
                (
                    df_me.apply(
                        lambda x: x['transcript'] in max_prob_per_multiorth and 
                        x['pred'] < max_prob_per_multiorth[x['transcript']], 
                        axis=1
                    )
                ) & (
                    df_me['pred'] >= self.orthology_threshold
                ) | (
                    df_me['pred'] < self.orthology_threshold
                )
            ) & (
                df_me_pp['clipped_exon_qlen'] > Constants.PP_CLIPPED_EXON_QLEN
            ) & (
            df_me_pp['clipped_intr_cover'] < Constants.PP_CLIPPED_INTRON_QLEN
            ) & (
                df_me_pp['clipped_intr_cover'] >= 0
            ),
            'pred'
        ] = -2

        ## if minimal chain score was set, 
        ## override predictions for chains with scores less than that
        ## unless they correspond to retrogenes/processed pseudogenes
        if self.min_orth_chain_score > 0:
            # underscored_chain_projections: List[str] = []
            if df_se.shape[0]:
                deprecated_se_names: pd.core.frame.DataFrame = df_se.loc[
                    df_se['gl_score'] < self.min_orth_chain_score
                ].apply(
                    lambda x: f'{x["transcript"]}#{x["chain"]}',
                    axis=1  
                )
                deprecated_se: List[str] = (
                    deprecated_se_names.to_list() if deprecated_se_names.shape[0] else []
                )
            else:
                deprecated_se: List[str] = []
            self.underscored_chain_projections.extend(deprecated_se)
            df_se = df_se[df_se['gl_score'] >= self.min_orth_chain_score]
            if df_me.shape[0]:
                deprecated_me_names: pd.core.frame.DataFrame = df_me[
                    (
                        df_me['gl_score'] < self.min_orth_chain_score
                    ) & (
                        df_me['pred'] != -2.0
                    )
                ].apply(
                    lambda x: f'{x["transcript"]}#{x["chain"]}', 
                    axis=1
                )
                deprecated_me: List[str] = (
                    deprecated_me_names.to_list() if deprecated_me_names.shape[0] else []
                )
            else:
                deprecated_me: List[str] = []
            self.underscored_chain_projections.extend(deprecated_me)
            df_me = df_me[
                (df_me['gl_score'] >= self.min_orth_chain_score) | (df_me['pred'] == -2.0)
            ]
        ## concatenate the results
        results: pd.core.Frame.DataFrame = pd.concat(
            [
                df_se.loc[:, ~df_se.columns.isin(['score', 'single_exon'])],
                df_me.loc[:, ~df_me.columns.isin(['score', 'single_exon'])],
                spanning_lines
            ]
        )
        results = results.loc[:, Constants.FINAL_COLUMNS]

        ## and write those to a file
        with open(self.orthology_class_table, 'w') as h:
            results.to_csv(h, header=True, index=False, sep='\t')

        ## there are two more files to write
        ## first, create a transcript-to-classified-chains table
        summary_dict: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
        for line in results.itertuples():
            gene: str = line.transcript
            chain: str = str(line.chain)
            prob: float = line.pred
            if prob == -1:
                summary_dict[gene][Constants.SPAN].append(chain)
            elif prob == -2:
                summary_dict[gene][Constants.P_PGENE].append(chain)
            elif prob < self.orthology_threshold:
                summary_dict[gene][Constants.PARA].append(chain)
            else:
                summary_dict[gene][Constants.ORTH].append(chain)

        with open(self.tr2chain_classes, 'w') as h:
            h.write(Constants.TR2CHAIN_HEADER + '\n')
            for g, data in summary_dict.items():
                orthologs: str = ','.join(data.get(Constants.ORTH, ['0']))
                paralogs: str = ','.join(data.get(Constants.PARA, ['0']))
                spanning: str = ','.join(data.get(Constants.SPAN, ['0']))
                pseudo: str = ','.join(data.get(Constants.P_PGENE, ['0']))
                h.write(
                    '\t'.join([g, orthologs, paralogs, spanning, pseudo]) + '\n'
                )

        ## second, detect and record unclassified genes
        rejected_transcripts: Set[str] = init_tr_set.difference(set(results.transcript))
        if not rejected_transcripts and not self.underscored_chain_projections:
            return
        with open(self.rejection_log, 'w') as h:
            for tr in rejected_transcripts:
                rej_line: str = Constants.UNCLASS_TEMPLATE.format(tr)
                h.write(rej_line + '\n')
            for proj in self.underscored_chain_projections:
                rej_line: str = Constants.UNDERSCORED_TEMPLATE.format(proj, self.min_orth_chain_score)
                h.write(rej_line + '\n')


if __name__ == '__main__':
    ChainClassifier()
