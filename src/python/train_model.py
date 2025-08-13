#!/usr/bin/env python3
"""Script to train XGBoost models."""

from numpy import log10, ndarray
from modules.constants import Constants
from modules.shared import CommandLineManager, CONTEXT_SETTINGS, get_upper_dir
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from typing import Optional

import click
import joblib
import os
import pandas as pd
import xgboost as xgb
# from version import __version__

__author__ = ('Bogdan Kirilenko', 'Yury V. Malovichko')
__year__ = '2025'
__credits__ = (
    'Michael Hiller', 'Virag Sharma', 'David Jebb'
)

"""
The script is expected to reside at TOGA2/src/python/ . 
Unless you're using custom model input and intend to save the models 
at a user-defined directory, please do not move the script around!
"""

DEFAULT_MODEL_DIR: str = os.path.join(get_upper_dir(__file__, 3), 'models')

DEFAULT_DF: str = os.path.join(DEFAULT_MODEL_DIR, 'train.tsv')
LD_DF: str = os.path.join(DEFAULT_MODEL_DIR, 'train_long_distance.tsv')

LD_MODEL: str = 'ld_model.dat'
ME_MODEL: str = 'me_model.dat'
SE_MODEL: str = 'se_model.dat'

RANDOM_STATE: str = 777
TREE_NUM: int = 50
MAX_DEPTH: int = 3
LD_MAX_DEPTH: int = 4
LEARNING_RATE: float = 0.1
N_SPLITS: int = 5

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '--default_dataset',
    type=click.File('r', lazy=True),
    metavar='TSV',
    default=DEFAULT_DF,
    show_default=True,
    help=(
        'An alternative path to the default train set. '
        'These data will be used for training default multi-exon ("me_model") '
        'and single-exon ("se_model") transcript orthology classification models.'
    )
)
@click.option(
    '--long_distance_dataset',
    type=click.File('r', lazy=True),
    metavar='TSV',
    default=LD_DF,
    show_default=True,
    help=(
        'An alternative path to the long-distance model ("ld_model") train set. '
        'This model is used optionally for orthology prediction refinement '
        'for evolutionary distant reference-query pairs.'
    )
)
@click.option(
    '--output',
    '-o',
    type=click.Path(exists=False),
    metavar='OUTPUT_DIR',
    default=DEFAULT_MODEL_DIR,
    show_default=True,
    help='A path to directory to store the trained models at'
)

class ModelTrainer(CommandLineManager):

    """
    Trains the multi-exon, single-exon, and long-distance models for TOGA2. By default, this script 
    is run once when building TOGA2.\n

    Note that the model joblib files are sensitive to Python3, joblib, and XGBoost versions and 
    might not be compatible with different versions of these utilities. If you want to run TOGA2 with environment 
    different from the one it was installed with, you will most likely have to train models anew.\n

    Also note that training metaparameters are hardcoded for the sake of reproducibility. 
    """

    __slots__ = (
        'df_def', 'df_ld', 'se_path', 'me_path', 'ld_path'
    )

    def __init__(
        self, 
        default_dataset: Optional[click.Path], 
        long_distance_dataset: Optional[click.Path],
        output: Optional[click.Path]
    ) -> None:
        """Entry point"""
        self.v: bool = True
        self.set_logging()

        self._to_log('Initializing model training module')
        self.df_def: click.Path = default_dataset
        self.df_ld: click.Path = long_distance_dataset

        self._mkdir(output)
        self.se_path: str = os.path.join(output, SE_MODEL)
        self.me_path: str = os.path.join(output, ME_MODEL)
        self.ld_path: str = os.path.join(output, LD_MODEL)

        self.train_default_models()
        self.train_ld_model()

    def train_default_models(self) -> None:
        """
        Trains default orthology clasification models for 
        single-exon ("se_model") and multi-exon ("me_model") reference transcripts
        """
        self._to_log('Training default models for single- and multi-exon transcripts')
        df: pd.core.frame.DataFrame = pd.read_csv(self.df_def, header=0, sep='\t')

        ## default train set supplied with TOGA2 contains data for both single- 
        ## and multi-exon models; split those into two separate dataset
        df_se: pd.core.frame.DataFrame = df[df['single_exon'] == 1]
        df_me: pd.core.frame.DataFrame = df[df['single_exon'] == 0]

        ## get predictors and target for single-exon model
        X_se: pd.core.frame.DataFrame = df_se.copy()
        X_se = X_se[Constants.SE_MODEL_FEATURES]
        y_se: pd.core.series.Series = df_se['y']
        ## and train the respective model
        self._to_log('Training default model for single exon transcripts')
        self._train(X_se, y_se, self.se_path, name='single-exon transcript', long_distance=False)

        ## rinse and repeat for multi-exon transcripts
        X_me: pd.core.frame.DataFrame = df_me.copy()
        X_me = X_me[Constants.ME_MODEL_FEATURES]
        y_me: pd.core.series.Series = df_me['y']
        self._to_log('Training default model for multi-exon transcripts')
        self._train(X_me, y_me, self.me_path, name='multi-exon transcript', long_distance=False)

    def train_ld_model(self) -> None:
        """
        Trains long-distance model ("ld_model") for orthology classification 
        in evolutionarily distant reference-query species pairs.

        NOTE: The code was taken from a standalone Jupyter notebook from TOGA1 
        dev branch. The model has not been revised for TOGA2 release, hence both the model 
        predictors and model training code can change in future
        """
        self._to_log('Training long-distance orthology classification model')
        ## upload the long-distance model train set
        df: pd.core.frame.DataFrame = pd.read_csv(self.df_ld, header=0, sep='\t')
        ## calculate the derived features
        df['exon_perc'] = df['exon_cover'] / df['ex_fract']
        df['synt_log'] = log10(df['synt'])
        df['intr_perc'] = df["intr_cover"] / df["intr_fract"]
        df = df.fillna(0.0)
        df['single_exon'] = df['ex_num'] == 1
        X: pd.core.frame.DataFrame = df[Constants.LD_MODEL_FEATURES].copy()
        y: pd.core.series.Series = df['y'].copy()
        ## liftoff!
        self._train(X, y, self.ld_path, name='long-distance', long_distance=True)

    def _train(
            self,
            x: ndarray, 
            y: ndarray, 
            save_to: str, 
            name: str,
            long_distance: bool
        ):
        """Model training method"""
        # create, fit, and cross-validate the model
        model = xgb.XGBClassifier(
            n_estimators=TREE_NUM, 
            max_depth=LD_MAX_DEPTH if long_distance else MAX_DEPTH, 
            learning_rate=LEARNING_RATE,
            random_state=RANDOM_STATE
        )
        kfold = StratifiedKFold(n_splits=N_SPLITS, random_state=RANDOM_STATE, shuffle=True)
        results = cross_val_score(model, x, y, cv=kfold)
        model.fit(x, y)
        y_lst = list(y)
        self._to_log('Summary for %s model: ' % name)
        self._to_log('Number of training samples: %i' % len(x))
        self._to_log('Features used: %s' % ','.join(x.columns))
        self._to_log(
            '#positives: %i; #negatives: %i' % (y_lst.count(1), y_lst.count(0))
        )
        self._to_log(
            'Accuracy: {0:.3f} {1:.3f}'.format(results.mean() * 100, results.std() * 100)
        )
        joblib.dump(model, save_to)  # save the model
        self._to_log('Model saved to: %s' % save_to)

if __name__ == '__main__':
    ModelTrainer()
