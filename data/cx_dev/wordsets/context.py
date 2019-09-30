'''
Classes that collect data on the grammatical 
environment surrounding a given node.
'''

import collections

from positions import Positions, Getter, Evaluator, getnext
   
class Mom:
    '''
    Identifies syntactic "relations"*
    from a word to its goverened elements.
    This class thus returns nodes that are 
    DEPENDENT on the mother (this node).
    
    Conditions are stored as strings and converted
    to booleans in order to easily diagnose why certain 
    items are validated.
    
    * the term "relation" is used loosely to
    refer to recurring noun patterns of modification.
    These patterns reflect productive constructions.
    '''
    
    def __init__(self, n, tf, **kwargs):
        quants = kwargs['quants'] # word sets
        preps = kwargs['preps']
        nominals = kwargs['noms']
        self.context = kwargs.get('context', 'phrase_atom')
        
        # set up variables needed for processing / storing
        P = Positions(n, self.context, tf).get
        quants = quants
        preps = preps
        self.P = P
        self.kids = {}
        self.explain = {}
        
        # set up evaluator with namespace as is
        # NB: Evaluator will rely on variables in
        # __init__, including P, quants, and preps
        self.conddict = Evaluator(locals()).conddict
        
    def analyze(self):
        '''
        Store ALL relationships.
        '''
        self.construct()
        self.adjacency()
        self.prep_mod()
        self.coordinate()
        self.quantifiers()

    def construct(self):
        '''
        Finds construct relations.
        '''
        P = self.P
        conddict = self.conddict

        const = (
                (P(1), conddict(
                    "P(0,'st') == 'c'",
                    "P(1,'sp') != 'art'",
                )),

                (P(2), conddict(
                    "P(0,'st') == 'c'",
                    "P(1,'sp') == 'art'",
                )),
        )

        self.kids['const'] = getnext(const)
        self.explain['const'] = const
        return getnext(const)

    def adjacency(self):
        '''
        Find adjacent nominals to this word.
        NB: a subset of these are so-called adjectives, but 
        this pattern does not yet parse down to that level.
        '''
        P = self.P
        conddict = self.conddict
        
        adja = (
            (P(1), conddict(
                "P(1) in nominals",
                "P(0,'st') in {'a', 'NA'}",
            )),

            (P(2), conddict(
                "P(1,'sp') == 'art'",
                "P(2) in nominals",
                "P(0,'st') in {'a', 'NA'}",
            )),
        )

        self.kids['adja'] = getnext(adja)
        self.explain['adja'] = adja
        return getnext(adja)

    def prep_mod(self):
        '''
        NB: noun + prep NOT prep + noun
        '''
        P = self.P
        conddict = self.conddict

        prep = (
            (P(1), conddict(
                "P(1) in preps"
            )),
        )

        self.kids['prep'] = getnext(prep)
        self.explain['prep'] = prep
        return getnext(prep)

    def coordinate(self):
        '''
        gets coordinate patterns
        NB: kid == thismom AND thatkid
        '''
        P = self.P
        conddict = self.conddict

        coord = (

            (P(2), conddict(
                "P(0) not in preps",
                "P(1,'sp') == 'conj'",
                "P(2) in nominals",
                "P(2) not in preps"
            )),

            (P(3), conddict(
                "P(0) not in preps",
                "P(1,'sp') == 'conj'",
                "P(2,'sp') == 'art'",
            )),

            (P(3), conddict(
                "P(1,'sp') == 'conj'",
                "P(2) in preps",
                "P(3) in nominals",
                "P(-1) in preps",
            )),

            (P(4), conddict(
                "P(1,'sp') == 'conj'",
                "P(2) in preps",
                "P(3,'sp') == 'art'",
                "P(-1,'sp') == 'art'",
                "P(-2) in preps",
            )),
        )

        self.kids['coord'] = getnext(coord)
        self.explain['coord'] = coord
        return getnext(coord)

    def quantifiers(self):
        '''
        -- quantifier patterns --
        '''
        P = self.P
        conddict = self.conddict

        quant = (

            (P(1), conddict(
                "P(0) not in quants",
                "P(1) in quants",
            )),

            (P(2), conddict(
                "P(0) not in quants",
                "P(2) in quants",
                "P(1,'sp') == 'art'",
            )),

            (P(-1), conddict(
                "P(0) not in quants",
                "P(-1) in quants",
            )),

            (P(-2), conddict(
                "P(0) not in quants",
                "P(-2) in quants",
                "P(-1,'sp') == 'art'",
            )),
        )

        self.kids['quant'] = getnext(quant)
        self.explain['quant'] = quant
        return getnext(quant)
    
class Relas:
    '''
    Provides mom-kid relations using the Mom class.
    '''
    def __init__(self, tf, **wsets):
        self.mom = collections.defaultdict(dict)
        self.kid = collections.defaultdict(dict)
        
        for w in tf.api.F.otype.s('word'):
            
            if w not in wsets['noms']:
                continue
            
            momma = Mom(w, tf, **wsets)
            momma.analyze()
            self.kid[w].update(momma.kids)
            for rela, kid in momma.kids.items():
                self.mom[kid].update({rela:w})
