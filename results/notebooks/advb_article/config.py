# Configurations for adverbial article notebooks
# includes file paths and data loading

# standard & data science packages
import collections
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.serif'] = ['SBL Biblit']
pd.set_option('display.max_rows', 100)
import seaborn as sns
from bidi.algorithm import get_display # bi-directional text support for plotting

# custom packages (see /tools)
import paths
import stats.significance as mystats
import stats.pca as my_pca
from tf_tools.load import load_tf
import tf_tools.formatting as form

# launch Text-Fabric with custom data
TF, API, A = load_tf('funct_assoc nhead top_assoc', silent='deep')
A.displaySetup(condenseType='phrase')
F, E, T, L = A.api.F, A.api.E, A.api.T, A.api.L # corpus analysis methods

# Configure article data paths and loads
article_data = paths.data.joinpath('advb_article')

# import narrow dataset
functions_data = article_data.joinpath('function_data.csv')
functs_df = pd.read_csv(functions_data)
functs_df.set_index(['node'], inplace=True)

# make broad dataset path
broad_dataset = article_data.joinpath('broad_dataset.csv')

def savefig(name):
    """Formats filename and save a figure"""
    file_name = paths.figs.joinpath(f'advb_article/{name}.svg')
    plt.savefig(file_name, format='svg', bbox_inches='tight')
