'''
This module contains scripts for testing statistical associations.
'''

import collections
import numpy as np
import pandas as pd
import scipy.stats as stats

def contingency_table(df, sample_axis, feature_axis):
    """Build 2x2 contingency table for calculating association measures.

    A 2x2 contingency table is defined as:
    ----------------------------------
   |              feature     ¬feature |
   |  sample         A            B    |
   |  ¬sample        C            D    |
    -----------------------------------
    
    Where sample is, e.g., a given word  and feature 
    is a given co-occurrence construction (Levshina 2015, 224); 
    ¬sample  (i.e. "not" sample) is every sample besides 
    a given word and ¬feature is every feature besides a 
    given feature. "A", "B", "C", "D" are the frequency 
    integers between the given categories. I follow Levshina's 
    2015 explanation (224) on setting up 2x2 contingency tables 
    for testing collocations of linguistic constructions.  

    This method calculates A, B, C, D and "expected frequency"
    for a supplied dataset. The values are built into matrices
    of the same dimensions as the input matrix, so that for co-
    occurrence of sample*feature, one can access, e.g., the A
    or B value for that individual example, e.g.:

        >> df_b[wordX][featureY] = B value for wordX*featureY 
    
    Given a sample and a feature count in a dataset, the math 
    for finding A, B, C, D (see Levshina 2015, 223) is:

        >> A = frequency of sample w/ feature (in dataset)
        >> B = sum(sample) - A
        >> C = sum(feature) - A
        >> D = sum(dataset) - (A+B+C)

    And the expected frequency (ibid., 211) is:

        >> E = sum(sample) * sum(feature) / sum(dataset)

    Args:
        df: a dataframe with co-occurrence frequencies in shape
            of samples*features or feature*samples
        sample_axis: 0 (row) or 1 (column); axis that contains 
            the sample population
        feature_axis: 0 (row) or 1 (column); axis that contains
            the collocating features on samples

    Returns: 
        5-tuple of dataframes as (a, b, c, d, e)
    """    

    # put data in sample * feature format for calculations
    # will flip it back at end if needed
    if sample_axis == 1 and feature_axis == 0:
        df = df.T
    elif sample_axis != 0 and feature_axis != 1:
        raise Exception('Invalid axis! Should be 0 or 1')

    # get observation sums across samples / features
    # fill each row in a column with the sum across the whole column
    # and do same for columns
    samp_margins = df.apply(
        lambda row: row.sum(), 
        axis=1, 
        result_type='broadcast' # keeps same shape
    ) 
    feat_margins = df.apply(
        lambda col: col.sum(), 
        axis=0, 
        result_type='broadcast'
    ) 
    total_margin = df.sum().sum()
    b = samp_margins.sub(df) # NB "df" == a
    c = feat_margins.sub(df)
    # for d, make table where every cell is total margin
    # use that table to make subtractions:
    d = pd.DataFrame.copy(df, deep=True)
    d[:] = total_margin 
    d = d.sub(df+b+c)
    e = samp_margins * feat_margins / total_margin
    # flip axes back if needed:
    if sample_axis == 1:
        df,b,c,d,e = df.T, b.T, c.T, d.T, e.T
    return (df, b, c, d, e)
    
def apply_fishers(df, sample_axis, feature_axis, logtransform=True):
    """Calculate Fisher's Exact Test with optional log10 transform.

    This function applies Fisher's Exact test to every 
    value in a co-occurrence matrix. It includes default 
    option to log-transform the results based on log10 
    and expected frequency condition. This is based on 
    the method of Stefanowitsch and Gries 2003, "Collostructions". 
    The resulting values "range from - infinitity (mutual repulsion) 
    to + infinity (mutual attraction)" (Levshina 2015, 232). 

    Args:
        df: a dataframe with co-occurrence frequencies in shape
            of samples*features or feature*samples
        sample_axis: 0 (row) or 1 (column); axis that contains 
            the sample population
        feature_axis: 0 (row) or 1 (column); axis that contains
            the collocating features on samples

    Returns:
        2-tuple of (p-values, odds_ratios) in DataFrames
    """

    # put data in sample * feature format for calculations
    # will flip it back at end if needed
    if sample_axis == 1 and feature_axis == 0:
        df = df.T
    elif sample_axis != 0 and feature_axis != 1:
        raise Exception('Invalid axis! Should be 0 or 1')
    a_df, b_df, c_df, d_df, e_df = contingency_table(df, 0, 1)
    ps = collections.defaultdict(lambda: collections.defaultdict())
    odds = collections.defaultdict(lambda: collections.defaultdict())

    # Calculate Fisher's value-by-value
    # I'm not yet sure if there's a better way to do this
    for sample in df.index:
        for feature in df.columns: 

            # exctract contingencies 
            a = df[feature][sample]
            b = b_df[feature][sample]
            c = c_df[feature][sample]
            d = d_df[feature][sample]
            
            # run Fisher's
            contingency = np.matrix([[a, b], [c, d]])
            oddsratio, p_value = stats.fisher_exact(contingency)
            
            # save and transform? scores
            odds[feature][sample] = oddsratio
            if not logtransform:
                ps[feature][sample] = p_value
            else:
                expected_freq = e_df[feature][sample]
                if a < expected_freq:
                    with np.errstate(divide='ignore'):
                        strength = np.log10(p_value)
                else:
                    with np.errstate(divide='ignore'):
                        strength = -np.log10(p_value)
                ps[feature][sample] = strength

    # package into dfs, flip axis back if needed
    orient = 'columns' if sample_axis == 0 else 'index' 
    ps = pd.DataFrame.from_dict(ps, orient=orient)
    odds = pd.DataFrame.from_dict(odds, orient=orient)
    return (ps, odds)
