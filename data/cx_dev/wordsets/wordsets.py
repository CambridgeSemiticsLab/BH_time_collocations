'''
This module delivers a series of sets
which are needed to process word relations
and semantic heads. This includes:

    • accent classification sets
    • custom nominal sets
    • custom quantifier sets
    • custom preposition sets
    • attested coordinate pair sets
    • attested construct pair sets
    
The sets are built by querying the corpus
for matching patterns. 
'''

import pickle
from tf.app import use
from accents import Accents
from nominals import Nominals
from quantifiers import Quants
from prepositions import Preps
from pairs import Conjunction, Construct

output = 'wsets.pickle'

class WordSets:
    '''
    Delivers word sets by executing
    their respective scripts.
    '''
    def __init__(self, tf, silent=True):
        
        self.silent = silent
        
        self.report('processing accents...')
        self.accents = Accents(tf)
        self.report('\tdone')
        
        self.report('processing quants...')
        self.quants = Quants(tf).quants
        self.report('\tdone')
        
        self.report('processing preps...')
        self.preps = Preps(tf).preps
        self.report('\tdone')
        
        self.report('processing nominals...')
        self.noms = Nominals(tf, preps=self.preps).nominals
        self.report('\tdone')
        
        base_sets = {
            'quants':self.quants,
            'preps':self.preps,
            'noms': self.noms
        }
        
        self.report('processing conjunction pairs...')
        self.conj = Conjunction(tf, **base_sets)
        self.report('\tdone')
        
        self.report('processing construct pairs...')
        self.cons = Construct(tf, **base_sets)
        self.report('\tdone')
         
    def report(self, mssg):
        if not self.silent:
            print(mssg)


# set up TF 
print('Setting up Text-Fabric...')
A = use('bhsa', hoist=globals(), silent=True)
print('\tdone...')

print('\n-- RUNNING WORDSETS --\n')
wsets = WordSets(A, silent=False)
print('\n-- WSETS COMPLETE --')

print('\npickleing word sets...')
export = {
    'noms': wsets.noms,
    'preps': wsets.preps,
    'quants': wsets.quants,
    'accent_type': wsets.accents.accenttype,
    'mwords': wsets.accents.mwords,
    'conj_pairs': wsets.conj.pairs,
    'cons_pairs': wsets.cons.pairs,
}

pickle.dump(export, open(output, 'wb'))

print('\n!*!*!*! DONE !*!*!*!')
