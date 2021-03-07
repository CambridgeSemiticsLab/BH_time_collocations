"""
This module and function calculates the
statistical association between a given
semantic head-word in a phrase and its
statistical association with the phrase's 
function in the clause (e.g. Subj, Time, Objc, etc.)
"""

import collections
import pandas as pd
import numpy as np
from stats.significance import contingency_table, apply_fishers

def calculate(ph_functions, api):
    """Run significance tests on all phrase heads.

    The output is a node feature on a phrase's head word. 
    The feature is an integer value which is a rounded 
    significance score which is log10 transformed with a 
    negative/positive sign applied to indicate kind of 
    association. 
    
    See Gries and Stefanowitsch, "Collostructions," 2003.

    Args:
        dictionary containing node:function string
    Returns:
        dictionary of node:float where float is a measure
            of attraction.
    """

    print('Running association features...')        

    # set up methods
    F, E, L = api.F, api.E, api.L
    
    # combine BHSA functions to prevent over-splitting
    funct_maps = {'PreO': 'Pred', 'PreS': 'Pred', 'PtcO': 'Pred',
                  'IntS': 'Intj', 'NCoS': 'NCop','ModS': 'Modi',
                  'ExsS': 'Exst'}

    # make head lex co-occurrence counts
    print('Beginning head lexeme // phrase function co-occurrence counts...')

    # store counts below:
    # functions is for top assoc
    # lex2funct2word is for funct_assoc, and is thus contextual 
    functions = collections.defaultdict(
        lambda: collections.Counter()
    )
    lex2funct2word = collections.defaultdict(
        lambda: collections.defaultdict(list)
    ) 

    # count head lexemes per function
    for phrase, funct in ph_functions.items(): 

        # skip complex phrases for accuracy
        # nhead is less reliable on complex phrases
        if len(L.d(phrase,'phrase_atom')) > 1:
            continue

        # make mappings and counts
        function = funct_maps.get(funct, funct)
        for head in E.nhead.t(phrase):
            head_lex = L.u(head,'lex')[0]
            functions[function][head_lex] += 1
            lex2funct2word[head_lex][function].append(head)

    functions = pd.DataFrame(functions).fillna(0)
    print(f'\t{functions.shape[0]} unique head lexemes counted...')
    print(f'\t{functions.shape[1]} unique phrase functions counted...')
            
    # apply association 
    # the math is done by apply_fishers
    size = functions.shape[0] * functions.shape[1]
    print(f'Applying Fisher\'s tests to {size} pairwise relations...')
    functions, odds = apply_fishers(functions, 0, 1)
    print('\tcalculations DONE')
    
    print(f'Substituting infinite scores with max and min associations...')
    ds_max = functions[functions != np.inf].max().max()
    ds_min = functions[functions != -np.inf].min().min()
    print(f'\tmin association: {round(ds_min)}')
    print(f'\tmax association: {round(ds_max)}')
    # replace with min/max scores
    for funct in functions.columns:
        for lex in functions.index:
            if functions[funct][lex] == np.inf:
                functions[funct][lex] = ds_max
            elif functions[funct][lex] == -np.inf:
                functions[funct][lex] = ds_min
    

    # package scores for TF
    nodeFeatures = {
        'top_assoc': {},
        'funct_assoc': {},
    }    
    for lexnode in functions.index:
        
        # get top associations 
        top_assoc = functions.loc[lexnode].sort_values(ascending=False).index[0]
        for word in L.d(lexnode, 'word'):
            nodeFeatures['top_assoc'][word] = top_assoc
            
        # get head association scores
        for function in functions.columns:
            assoc_score = round(functions[function][lexnode])
            for word in lex2funct2word[lexnode][function]:
                nodeFeatures['funct_assoc'][word] = int(assoc_score)
    
    len_assocs = len(nodeFeatures['funct_assoc']) 
    len_topassocs = len(nodeFeatures['top_assoc'])
    print('Exporting TF data...')
    print(f'\t{len_assocs} association scores...')
    print(f'\t{len_topassocs} top associations...')
  
    return nodeFeatures
