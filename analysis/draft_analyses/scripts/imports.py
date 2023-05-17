"""
Standard imports for the analysis notebooks.
"""
import re

import pandas as pd
pd.set_option('display.max_rows', 200)
idx = pd.IndexSlice

from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
from scripts.df_styles import TextShower
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

# thesis tags
def textheb(string):
    return '\\texthebrew{%s}'%string

# custom modules
from analysis_tools.src.analysis_tools.paths import paths

# load the source_data
df = pd.read_csv(paths['time_dataset'], index_col='node')
df_sg = df.query("(n_times == 1) and (is_advb == False)").copy()

# pretty-show Hebrew text from a df
ts = TextShower(
    default=['verse', 'clause'],
    stylize=['clause']
)
