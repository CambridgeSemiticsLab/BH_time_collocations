'''
Functions that perform
various helper tasks for the analyses.
'''

import pandas as pd

def convert2pandas(counterdict):
    '''
    Converts a counter dict to a sorted Pandas DF
    '''
    return pd.DataFrame.from_dict(counterdict, orient='index', columns=['Total']).sort_values(by='Total', ascending=False)