'''
Defines a series of observed word-pairs in the corpus 
which are used for distinguishing ambiguous cases.
These functions require pre-processed word sets 
of quantifiers and prepositions.
'''

import collections
from positions import Positions
from context import Mom

class Conjunction:
    '''
    Assembles attested conjunction word-pairs in BHSA.
    Returns a dict mapping from a word to all its pairs.
    '''
    
    def __init__(self, tf, **wsets):
        
        self.tf = tf
        F = tf.api.F
        covered = set() # skip items already matched
        self.wsets = wsets
        self.pairs = collections.defaultdict(set)
        self.pairresults = collections.defaultdict(lambda:collections.defaultdict(list))
    
        # gather pairs
        for w in F.otype.s('word'):

            # skip words already visited in a chain
            # or those that are not nominal parts of speech
            if w in covered or w not in wsets['noms']:
                continue

            # check for chain
            chain = list(self.conj_climber(w))
            if not chain:
                continue

            # add pairs
            for i in chain:
                for j in chain:
                    if i == j:
                        continue
                    self.pairs[F.lex.v(i)].add(F.lex.v(j))
                    self.pairresults[F.lex.v(i)][F.lex.v(j)].append((i, j))

    def conj_climber(self, a):
        '''
        Climbs down conjunction chains recursively
        and yields the connected words. Start with first word.
        '''
        F = self.tf.api.F

        yield a

        b = Mom(a, self.tf, **self.wsets).coordinate()

        if b:
            yield from self.conj_climber(b)
            
class Construct:
    '''
    Collects all potential goverened
    construct pairs, where key is in
    construct state and value is a set
    of valid goverened genitives.
    '''
    
    def __init__(self, tf, **wsets):
        
        F = tf.api.F
        
        self.pairs = collections.defaultdict(set)
        self.pairresults = collections.defaultdict(lambda:collections.defaultdict(list))
        
        for word in set(F.otype.s('word')) & set(F.st.s('c')):
        
            # skip non-nominal words
            if word not in wsets['noms']:
                continue
        
            construct = Mom(word, tf, **wsets).construct()
            
            if not construct:
                continue
                
            self.pairs[F.lex.v(word)].add(F.lex.v(construct))
            self.pairresults[F.lex.v(word)][F.lex.v(construct)].append((word, construct))