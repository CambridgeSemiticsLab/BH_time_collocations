"""
Standard imports for the analysis notebooks.
"""
import re
import json
import matplotlib.pyplot as plt

import seaborn as sns
import pandas as pd
pd.set_option('display.max_rows', 200)
idx = pd.IndexSlice
from adjustText import adjust_text

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
scaler = StandardScaler()
from scripts.stats.pca import apply_pca
from scripts.stats import significance as sig
from scripts.df_styles import TextShower, df_highlighter
from scripts.counting import pivot_ct, join_ct_pr
import numpy as np
from bidi.algorithm import get_display as bidi_get_display

def remove_shindots(string):
    """Remove dots from ש"""
    return(
        string
            .replace('\u05c1', '') 
            .replace('\u05c2', '') 
    )  

# mod get display to remove shin dots
def get_display(string):
    rem_accents = ''.join(re.findall('[\u05D0-\u05EA]', string))
    return bidi_get_display(rem_accents)

# latex tags
def textheb(string):
    return '\\texthebrew{%s}'%string

# custom modules
from .paths import paths
from .export import Exporter
from .plotting import heatmap
from .string_refs import get_verserefs

# load the data
df = pd.read_csv(paths['time_dataset'], index_col='node')
df_sg = df.query("(n_times == 1) and (is_advb == False)").copy()

# pretty-show Hebrew text from a df
ts = TextShower(
    default=['verse', 'clause'],
    stylize=['clause']
)
