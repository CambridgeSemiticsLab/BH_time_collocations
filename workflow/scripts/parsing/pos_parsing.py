import json
from tf.fabric import Fabric
from .positions import PositionsTF

# a class to parse BHSA word nodes into parts of speech
class posParser:

    def __init__(self, tf_api):
        self._tfapi = tf_api
        self._F = tf_api.F
        self._E = tf_api.E
        self._T = tf_api.T
        self._L = tf_api.L

        # precedence assignments for when multiple matches are found
        # higher value == the function will be matched first
        self._prec = {
            'CONJPP': 1, # match first to eliminate some eroneous CONJGP matches
        }

    def _getprec(self, name):
        """Get precedence for parse values."""
        return self._prec.get(name, 0)

    def _getparsers(self):
        """Retrieve defined parsers in this class.

        Parsers are identified by those methods which
        are named without a prefixed underscore.
        """

        # retrieve defined parsers
        parsenames = [
            m for m in dir(self) 
                if not m.startswith('_')
        ]

        # sort by precedence
        parsenames = sorted(
            parsenames,
            key=self._getprec,
            reverse=True,
        )

        # remap method names to actual functions for calls
        parsers = {
            p:getattr(self, p) for p in parsenames
        }
        
        return parsers

    def _parse(self, wordn):
        """Parse a supplied BHSA word node."""

        # get parsers sorted by precedence
        parsers = self._getparsers()

        # return first match in order of precedence
        for value, parser in parsers.items():
            if parser(wordn):
                return value

    def _getP(self, node, context='phrase_atom'):
        """Get Positions object for a TF node."""
        return PositionsTF(node, context, self._tfapi).get

    # -- PARSERS defined below --
    # The tag for each parser method is the name
    # of the method itself. Parser methods return
    # a boolean on whether the given tag is a match 
    # Text-Fabric methods are heavily utilized
    # to access the text-graph and its features

    def PREP(self, w):
        """A preposition word."""
        
        P = self._getP(w)
        F, L = self._F, self._L
        phn = L.u(w, 'phrase')[0] # containing phrase node

        return any([

            # prepositional דרך
            (   
                F.ls.v(w) == 'ppre'
                and F.lex.v(w) != 'DRK/'
            ),

            # prepositional ראשׁ
            (
                F.lex.v(w) == 'R>C/'
                and F.st.v(w) == 'c'
                and P(-1,'pdp') == 'prep'
                and F.function.v(phn) in {
                        'Time', 'Adju', 
                        'Cmpl', 'Loca',
                }
            ),
            (
                F.lex.v(w) == 'R>C/'
                and F.st.v(w) == 'c'
                and F.function.v(phn) == 'Time'
            ),

            # key lexemes in construct
            (
                F.lex.v(w) in {
                    'PNH/','TWK/', 
                    'QY/', 'QYH=/', 
                    'QYT/',
                }
                and F.prs.v(w) == 'absent'
                and F.st.v(w) == 'c'
            ),
            
            # key lexemes with certain prepositions
            (
                F.lex.v(w) == 'PNH/'
                and P(-1, 'lex') in {'L', 'MN'}
            ),
            (
                F.lex.v(w) == 'BD/'
                and P(-1,'lex') == 'L'
            ),
            (
                F.lex.v(w) == '>XRJT/'
                and F.st.v(w) == 'c'
                and not ({P(1,'lex'), P(2,'lex')} & {'>JWB/', 'RC</'})
            ),

            # fix: this should be an adverb
            #(
            #    F.lex.v(w) == '<YM/',
            #    and F.st.v(w) == 'c',
            #    and F.function.v(phn) == 'Time'
            #)
        ])
    
    def PROPN(self, w):
        """A proper noun."""
        return any([
            self._F.pdp.v(w) == 'nmpr'
        ])
    
    def NOUN(self, w):
        """
        Substantives
        """
        F = self._F
        return any([
            # potentially substantive participles
            (
                F.sp.v(w) == 'verb'
                and F.vt.v(w) in {'ptcp', 'ptca'}
                and F.pdp.v(w) not in {'advb', 'adjv'}
            )
        ])    
    
    def CARD(self, w):
        """A cardinal number followed by other cardinals."""
        P = self._getP(w)
        return any([
            (
                self._F.ls.v(w) == 'card'
                and (
                    P(-1,'ls') == 'card'
                    or (P(1,'ls') == 'card')
                    or (P(-1,'lex') == 'W' and P(-2,'ls') == 'card')
                    or (P(1,'lex') == 'W' and P(2,'ls') == 'card')
                )
            ),
        ])

    def CARD1(self, w):
        """A standalone cardinal number."""
        P = self._getP(w)
        return any([
            (
                self._F.ls.v(w) == 'card'
                and not self.CARD(w)
            ),
        ])
    
    def ORDN(self, w):
        """An ordinal word."""
        return any([
            self._F.ls.v(w) == 'ordn',
        ])
    
    def QUANT(self, w):
        """A qualitative quantifier word."""
        F = self._F
        return any([
            F.lex.v(w) in {
                'KL/', 'M<V/', 'JTR/',
                'XYJ/', 'C>R=/', 'MSPR/', 
                'RB/', 'RB=/',
            },
            F.lex.v(w) in {
                'M<FR/', '<FRWN/',
                'XMJCJT/',
            },
        ])

    def CONJPP(self, w):
        """A conjunction followed by a preposition.

        This aids the parser which cannot look ahead 
        more than 1 step. By this method we provide a
        token that has already looked ahead and prevents
        the parser from unecessarily shifting tokens over
        instead of reducing."""
        P = self._getP(w)
        nextw = P(1)
        if nextw:
            nw_is_prep = (
                self.PREP(nextw)
                or self._F.pdp.v(nextw) == 'prep'
            )
        else:
            nw_is_prep = False
        return any([
            (
                self._F.pdp.v(w) == 'conj' 
                and nw_is_prep
            )
        ])

    # -- Conjunction Types --
    # These aid the phrase parser later on, which 
    # cannot look ahead more than 1 step. By these POS tags, 
    # we provide a token that has already looked around,
    # preventing the parser from unecessarily shifting tokens 
    # over instead of reducing.

    def CONJGP(self, w):
        """A conjunction followed by a construct."""
        P = self._getP(w)
        return any([
            (
                P(-2,'st') == 'c'
                and self._F.pdp.v(w) == 'conj'
                and P(1,'st') == 'c'
                and P(1,'ls') != 'card'
            ),
            (
                P(-3,'st') == 'c'
                and P(-2,'lex') == 'H'
                and self._F.pdp.v(w) == 'conj'
                and P(1,'st') == 'c'
                and P(1,'ls') != 'card'
            )
        ])

    def CONJADJV(self, w):
        """Conj connecting 2 adjv phrases"""
        P = self._getP(w)
        return any([
            (
                P(-1,'pdp') == 'adjv'
                and self._F.pdp.v(w) == 'conj'
                and P(2,'pdp') == 'adjv'
            )
        ])

    def CONJCARD(self, w):
        """
        Conj surrounded by cardinals.
        """
        P = self._getP(w)
        return any([
            (
                P(-1,'ls') == 'card'
                and self._F.pdp.v(w) == 'conj'
                and P(1,'ls') == 'card'
            )
        ])

    def ADVB(self, w):
        """Identify adverbs.
    
        Currently only identifies עוד
        """
        return self._F.lex.v(w) == '<WD/'
            

def parse_pos(data_locs, slot2pos_path, uniquepos_path):
    """Apply parsing using Text-Fabric."""

    # initialize Text-Fabric
    tf = Fabric(data_locs)
    tf_api = tf.load('''
        pdp sp ls lex st 
        prs function vt
    ''')

    F = tf_api.F

    # assign parts of speech
    # if posParser does not return a value, we 
    # assign it the default BHSA tag in uppercase
    parser = posParser(tf_api)
    pos_map = {'subs':'NOUN'}
    slot2pos = {}
    for slot in F.otype.s('word'):
        pdp = F.pdp.v(slot)
        pos = parser._parse(slot)
        pos = pos or pos_map.get(pdp, pdp.upper())
        slot2pos[slot] = pos 

    uniquepos = sorted(set(slot2pos.values()))

    # export to JSON
    with open(slot2pos_path, 'w') as outfile:
        json.dump(slot2pos, outfile, indent=1)

    # export unique values, which are needed to configure the phrase parser
    with open(uniquepos_path, 'w') as outfile:
        outfile.write(',\n'.join(uniquepos))
