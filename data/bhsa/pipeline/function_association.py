'''
This module and class calculates the
statistical association between a given
semantic head-word in a phrase and its
statistical association with the phrase's 
function in the clause (e.g. Subj, Time, Objc, etc.)
'''
import collections, sys
from tf.fabric import Fabric
import pandas as pd
import numpy as np
import scipy.stats as stats
sys.path.append('../../analysis/tools')
from significance import contingency_table, apply_fishers

class FunctAssoc:
    
    def __init__(self, locs, metadata):
        
        # initialize Text-Fabric methods
        locas = [locs['bhsa'], locs['heads'], locs['output']]
        TF = Fabric(locations=locas, silent=True)
        load_features = ['nhead','function']
        self.bhsa = TF.load(' '.join(load_features), silent=True)    
        self.output_dir = locs['output']
        self.meta = metadata
        self.nodeFeatures = collections.defaultdict(lambda:collections.defaultdict()) # put new features here
        
    def execute(self):
        '''
        Calculates the statistical significance
        scores and exports TF data. The data 
        is a node feature on a phrase's head word. 
        The feature is an integer value which
        is a rounded significance score
        which is log10 transformed with a negative/positive
        sign applied to indicate kind of association.
        See Gries and Stefanowitsch, "Collostructions," 2003.
        '''
        
        F, E, L = self.bhsa.F, self.bhsa.E, self.bhsa.L
        
        # mappings to strings to prevent unnecessary splitting
        funct_maps = {'PreO': 'Pred', 'PreS': 'Pred', 'PtcO': 'Pred',
                      'IntS': 'Intj', 'NCoS': 'NCop','ModS': 'Modi',
                      'ExsS': 'Exst'}

        # make head lex co-occurrence counts
        print('Beginning head lexeme // phrase function co-occurrence counts...')
        functions = collections.defaultdict(lambda: collections.Counter())
        lex2funct2word = collections.defaultdict(lambda: collections.defaultdict(list)) # for node feature assignment
        for phrase in F.otype.s('phrase'):
            funct = F.function.v(phrase)
            function = funct_maps.get(funct, funct)
            # map to head lex
            for head in E.nhead.t(phrase):
                head_lex = L.u(head, 'lex')[0]
                functions[function][head_lex] += 1
                lex2funct2word[head_lex][function].append(head)
        functions = pd.DataFrame(functions).fillna(0)
        print(f'\t{functions.shape[0]} unique head lexemes counted...')
        print(f'\t{functions.shape[1]} unique phrase functions counted...')
                
        # apply association 
        # the math is done by apply_fishers
        size = functions.shape[0] * functions.shape[1]
        print(f'Applying Fisher\'s tests to {size} pairwise relations...')
        functions = apply_fishers(functions)
        print('\tcalculations DONE')
        
        print(f'Substituting infinite scores with max and min associations...')
        ds_max = functions[functions != np.inf].max().max()
        ds_min = functions[functions != -np.inf].min().min()
        print(f'\tmin association: {round(ds_min)}')
        print(f'\tmax association: {round(ds_max)}')
        # replace with min/max scores
        for funct in functions:
            for lex in functions.index:
                if functions[funct][lex] == np.inf:
                    functions[funct][lex] = ds_max
                elif functions[funct][lex] == -np.inf:
                    functions[funct][lex] = ds_min
        
        # get association scores
        for lexnode in functions.index:
            
            # get top associations 
            top_assoc = functions.loc[lexnode].sort_values(ascending=False).index[0]
            for word in L.d(lexnode, 'word'):
                self.nodeFeatures['top_assoc'][word] = top_assoc
                
            # get head association scores
            for function in functions.columns:
                assoc_score = round(functions[function][lexnode])
                for word in lex2funct2word[lexnode][function]:
                    self.nodeFeatures['funct_assoc'][word] = int(assoc_score)
        
        len_assocs = len(self.nodeFeatures['funct_assoc']) 
        len_topassocs = len(self.nodeFeatures['top_assoc'])
        print('Exporting TF data...')
        print(f'\t{len_assocs} association scores...')
        print(f'\t{len_topassocs} top associations...')
        
        meta = {'':self.meta,
                'funct_assoc': {'valueType':'int',
                             'description':'a feature on words that function as a head in their enclosing phrase; integer tells how attracted the head word is to its phrase\'s functions',
                             'interpreting scores':'score > 1.3 is significantly attracted; score < -1.3 is significantly repelled'},
                'top_assoc': {'valueType':'str',
                             'description':'top associated function to this word',
                             'interpreting scores':'score > 1.3 is significantly attracted; score < -1.3 is significantly repelled'}
                }
        
        TF = Fabric(locations=self.output_dir, silent=True)
        TF.save(nodeFeatures=self.nodeFeatures, metaData=meta)
        print('SUCCESS')