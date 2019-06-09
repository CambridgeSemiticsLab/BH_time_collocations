'''
This module contains scripts for visualizing
data in the analysis notebooks. The visualizations
are built with Seaborn and Matplotlib.
'''

import seaborn as sns
import matplotlib.pyplot as plt

def reverse_hb(hb_text):
    '''
    Reverses a string of Hebrew text 
    in order to display right-to-left
    correctly in Matplotlib.
    '''
    return ''.join(reversed(hb_text))

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