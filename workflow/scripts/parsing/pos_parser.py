from positions import PositionsTF

# a class to parse BHSA word nodes into parts of speech
class posParser:

    def __init__(self, tf_api):
        self._tfapi = tf_api
        self._F = tf_api.F
        self._E = tf_api.E
        self._T = tf_api.T
        self._L = tf_api.L

        # precedence assignments for when multiple matches are found
        # higher value == higher precedence
        self._prec = {}

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

    def _getP(self, node, context='phrase'):
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

            # vanilla subs
            F.pdp.v(w) == 'subs',

            # potentially substantive participles
            (
                F.sp.v(w) == 'verb'
                and F.vt.v(w) in {'ptcp', 'ptca'}
            )
        ])    
    
    def CARD(self, w):
        """A cardinal number."""
        return any([
            self._F.ls.v(w) == 'card',
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

def tag_word(w, tf_api):
    """Tag grammatical word features on a given BHSA word node."""

    F = tf_api.F
    gender = F.gn.v(w)
    number = F.nu.v(w)
    state = F.st.v(w)
    sffx = F.prs.v(w)
