'''
This module contains scripts for visualizing
data in the analysis notebooks. The visualizations
are built with Seaborn and Matplotlib.
'''

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

class Figures:
    """Auto config figure titles"""

    def __init__(self, chapter=1):
        self.chapter = chapter
        self.fig2num = {}

    def title(self, title_str, chapter=None):
        """Configure a title string"""
        
        chapter = chapter or self.chapter

        # configure fig number
        nums = self.fig2num
        num = nums.setdefault(
            title_str,
            max(nums.values(), default=0) + 1
        )

        return f'fig.{chapter}.{num}_{title_str}'

def convert2pandas(counterdict):
    ''' 
    Converts a counter dict to a sorted Pandas DF
    '''
    return pd.DataFrame.from_dict(counterdict, orient='index', columns=['Total']).sort_values(by='Total', ascending=False)

def get_index(thislist, index, default=None):
    ''' 
    A safe way to get index from 
    a list/tuple. If indexError returns default.
    '''
    try:
        return thislist[index]
    except IndexError:
        return default

def reverse_hb(hb_text):
    '''
    Reverses a string of Hebrew text 
    in order to display right-to-left
    correctly in Matplotlib.
    '''
    return ''.join(reversed(hb_text))

## TODO! Decide which barplot counter to keep

def barplot_counts(count_dict, title='', reverse_labels=False, size=(8, 6), text_size=14, rotation=None, limit=None, save=''):
    '''
    Makes simple barplot of counts contained
    in a collections.Counter dictionary.
    '''
    count_df = pd.DataFrame.from_dict(count_dict, orient='index', columns=['count']).sort_values(ascending=False, by='count')
    plotme = count_df.head(limit) if limit else count_df
    
    n_bars = list(range(0, plotme.shape[0]))
    x_labels = [''.join(reversed(prep)) for prep in plotme.index] if reverse_labels else plotme.index
    plt.figure(figsize=size)
    sns.barplot(n_bars, plotme['count'], color='darkblue')
    plt.xticks(n_bars, x_labels, size=text_size, rotation=rotation)
    plt.yticks(size=text_size)
    plt.title(title, size=text_size)
    plt.ylabel('count', size=text_size)
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    plt.show()
    return count_df

def countBarplot(count_df, 
                 title='', 
                 column='Total', 
                 reverse_labels=False, 
                 size=(8, 6),
                 xlab_rotation=None,
                 ylim=None,
                 save=None,
                 xlabel=None,
                ):
    '''
    Makes simple barplot from collections.Counter type objects.
    '''
    n_bars = list(range(0, count_df.shape[0]))
    x_labels = [''.join(reversed(prep)) for prep in count_df.index] if reverse_labels else count_df.index
    plt.figure(figsize=size)
    sns.barplot(n_bars, count_df[column], color='darkblue')
    plt.xticks(n_bars, x_labels, size=18, rotation=xlab_rotation)
    plt.yticks(size=18)
    if ylim:
        plt.ylim(top=ylim[0], bottom=ylim[1])
    if xlabel:
        plt.xlabel(xlabel,size=18)
    plt.ylabel(column, size=18)    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    plt.title(title, size=18,  y=1.05)
    plt.show()
