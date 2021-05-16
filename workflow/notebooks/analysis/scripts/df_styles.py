"""
This module contains functions for stylizing dataframes.

On showing text:
Dataframes in this project contain columns of text. 
When analyzing various categories, we often want to
access these  text columns and browse particular
samples without looking at all of the columns. 
We also want the Hebrew text to show up big and 
readable with a serif font. The TextShower (as in "show"er) 
does this and can also provide an even subsample spread 
across the dataframe for getting a sense of the dataset.

On highlighting max values in a table:
The highlighter functions can highlight in red those values which
are the max across their rows (for example).
""" 

import numpy as np

def pretty_hebrew(val):
    """Render Hebrew in a dataframe."""
    return 'font-size:20px; font-family: Times New Roman; text-align: right; max-width: 500px'

def get_spread(array, n): 
    """Retrieve an even spread of indices for an array/Series.
    
    This allows us to access representative samples
    from across the corpus.
    see https://stackoverflow.com/a/50685454/8351428
    
    Args:
        array: either numpy array or Pandas Series 
            (with multiple indices allowed)
        n: number of indices to return
    Returns:
        indexed array or series
    """
    end = len(array) - 1 
    spread = np.ceil(np.linspace(0, end, n)).astype(int)
    indices = np.unique(spread)
    try:
        return array[indices]
    except KeyError:
        return array.iloc[indices]

class TextShower:
    """Show stylized text examples from a DataFrame"""
    
    def __init__(self, default=['ref', 'sentence', 'text_full'], 
                 stylize=['sentence', 'text_full']):
        """Initialize a TS object.

        Objects can be set up in advance to handle showing
        text examples in various ways.

        Args:
            default: default columns to show always
            stylize: columns with Hebrew text to stylize
        """
        self.default = default
        self.stylize = stylize

    def show(self, df, extra=[], spread=0):
        """Display text from pandas dataframe in a readable way."""
        original_shape = df.shape
        df = df[self.default + extra]
        if spread > 0:
            spread_i = get_spread(df.index, spread)
            df = df.loc[spread_i]
        print(f'showing {df.shape[0]} of {original_shape[0]}')
        return df.style.applymap(pretty_hebrew, subset=self.stylize)

def highlight_max(s):
    """Highlight max value in a df column."""
    is_max = s == s.max()
    return ['color: red' if v else '' for v in is_max]

def highlight_sig(s, sig_up=1.3, sig_down=-1.3):
    """Highlights values of significance (Fishers > or <  1.3)"""
    if s > sig_up:
        color = 'red'
    elif sig_down and s < sig_down:
        color = 'blue'
    else:
        color = ''
    return f'color: {color}'

def df_highlighter(pr_df, rule='max'):
    """Show proportion dataframe with highlighting."""
    if rule == 'max':
        return pr_df.style.apply(highlight_max, 1)
    elif rule == 'fishers':
        return pr_df.style.applymap(highlight_sig)
    else:
        return pr_df.style
