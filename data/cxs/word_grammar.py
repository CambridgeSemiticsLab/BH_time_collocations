from cx_analysis.cx import Construction
from cx_analysis.build import CXbuilder, CXbuilderTF

class Words(CXbuilderTF):
    """Build word constructions."""
    
    def __init__(self, tf, **kwargs):
        
        """Initialize with Constructions attribs/methods."""
        CXbuilderTF.__init__(self, tf, **kwargs)
        
        # Order matters! More specific meanings last
        self.cxs = (
            self.pos,
            self.prep,
            self.qual_quant,
            self.card,
            self.ordn,
            self.name,
            self.cont_ptcp,
        )
        
        self.kind = 'word_cx'
    
    def cxdict(self, slotlist):
        """Map all TF word slots to a construction.
        
        Method returns a dictionary of slot:cx
        mappings.
        """
        
        slot2cx = {}
        for w in slotlist:
            for cx in self.findall(w):
                slot2cx[w] = cx
    
        return slot2cx
    
    def pos(self, w):
        """A drip-bucket part of speech CX.
        
        The standard ETCBC feature is pdp,
        which is "phrase-dependent part of
        speech." I.e. it is a contextually
        sensive pos label.
        """
        
        F = self.F
        
        # map
        pdplabel = {
            'subs': 'cont',
            'adjv': 'cont',
            'advb': 'cont',
        }
        pdp = F.pdp.v(w)
        
        return self.test(
            {
                'element': w,
                'name': f'{pdplabel.get(pdp, pdp)}',
                'kind': self.kind,
                'roles': {'head': w},
                'conds': {
                    f'bool(F.pdp.v({w}))':
                        bool(F.pdp.v(w)),
                }
            }
        )
    
    def prep(self, w):
        """A preposition word."""
        
        P = self.getP(w)
        F, L = self.F, self.L
        name = 'prep'
        roles = {'head': w}
        return self.test(
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'ETCBC pdp',
                'roles': roles,
                'conds': {
                    'F.pdp.v(w) == prep':
                        F.pdp.v(w) == 'prep',
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'ETCBC ppre words',
                'roles': roles,
                'conds': {
                    'F.ls.v(w) == ppre':
                        F.ls.v(w) == 'ppre',
                    'F.lex.v(w) != DRK/':
                        F.lex.v(w) != 'DRK/',
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'prep+R>C/',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) == R>C/':
                        F.lex.v(w) == 'R>C/',
                    'F.st.v(w) == c':
                        F.st.v(w) == 'c',
                    'P(-1,pdp) == prep':
                        P(-1,'pdp') == 'prep',
                    'phrase is adverbial':
                        F.function.v(
                            L.u(w,'phrase')[0]
                        ) in {
                            'Time', 'Adju', 
                            'Cmpl', 'Loca',
                        },
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'R>C/ no prep',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) == R>C/':
                        F.lex.v(w) == 'R>C/',
                    'F.st.v(w) == c':
                        F.st.v(w) == 'c',
                    'phrase is adverbial':
                        F.function.v(L.u(w,'phrase')[0]) == 'Time',
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'construct lexs',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) in lexset':
                        F.lex.v(w) in {
                            'PNH/','TWK/', 
                            'QY/', 'QYH=/', 
                            'QYT/',
                        },
                    'F.prs.v(w) == absent':
                        F.prs.v(w) == 'absent',
                    'F.st.v(w) == c':
                        F.st.v(w) == 'c'
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'L/M PNH/',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) == PNH/':
                        F.lex.v(w) == 'PNH/',
                    'P(-1, lex) in {L, MN}':
                        P(-1, 'lex') in {'L', 'MN'},
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'L+BD',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) == BD/':
                        F.lex.v(w) == 'BD/',
                    'P(-1,lex) == L':
                        P(-1,'lex') == 'L',
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': '>XRJT/',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) == >XRJT/':
                        F.lex.v(w) == '>XRJT/',
                    'F.st.v(w) == c':
                        F.st.v(w) == 'c',
                    'P(1,lex) or P(2,lex) not >JWB|RC</':
                        not {
                            P(1,'lex'), P(2,'lex')
                        } & {
                            '>JWB/', 'RC</'
                        }
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': '<YM/ time',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) == <YM/':
                        F.lex.v(w) == '<YM/',
                    'F.st.v(w) == c':
                        F.st.v(w) == 'c',
                    'F.function.v(phrase) == Time':
                        F.function.v(
                            L.u(w,'phrase')[0]
                        ) == 'Time',
                }
            }
        )
    
    def name(self, w):
        """A name word (i.e. proper noun)."""
        return self.test(
            {
                'element': w,
                'name': 'name',
                'kind': self.kind,
                'roles': {'head': w},
                'conds': {
                    'F.pdp.v(w) == nmpr':
                        self.F.pdp.v(w) == 'nmpr'
                }
            }
        )
    
    def cont_ptcp(self, w):
        """A content word participle.
        
        A participle which can potentially
        function like a "noun" i.e. a content word.
        """
        
        F = self.F
        
        return self.test(
            {
                'element': w,
                'name': 'cont',
                'kind': self.kind,
                'pattern': 'participle',
                'roles': {'head': w},
                'conds': {
                    'F.sp.v(w) == verb':
                        F.sp.v(w) == 'verb',
                    'F.vt.v(w) in {ptcp, ptca}':
                        F.vt.v(w) in {'ptcp', 'ptca'},
                }
            },
        )    
    
    def card(self, w):
        """A cardinal number."""
        
        F = self.F
        P = self.getP(w)
        name = 'card'
        roles = {'head': w}
        
        return self.test(
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'roles': roles,
                'conds': {
                    'F.ls.v(w) == card':
                        F.ls.v(w) == 'card',
                }
            },
        )
    
    def ordn(self, w):
        """An ordinal word."""
        
        F = self.F
        P = self.getP(w)
        roles = {'head': w}
        
        return self.test(
            {
                'element': w,
                'name': 'ordn',
                'kind': self.kind,
                'pattern': 'ETCBC ls',
                'roles': roles,
                'conds': {
                    'F.ls.v(w) == ordn':
                        F.ls.v(w) == 'ordn',
                }
            },
        )
    
    def qual_quant(self, w):
        """A qualitative quantifier word."""
        
        F = self.F
        P = self.getP(w)
        name = 'qquant'
        roles = {'head': w}
        
        return self.test(
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'qualitative',
                'roles': roles,
                'conds': {
                    f'{F.lex.v(w)} in lexset':
                        F.lex.v(w) in {
                            'KL/', 'M<V/', 'JTR/',
                            'XYJ/', 'C>R=/', 'MSPR/', 
                            'RB/', 'RB=/',
                        },
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'portion',
                'roles': roles,
                'conds': {
                    f'{F.lex.v(w)} in lexset':
                        F.lex.v(w) in {
                            'M<FR/', '<FRWN/',
                            'XMJCJT/',
                        },
                }
            },
        )
