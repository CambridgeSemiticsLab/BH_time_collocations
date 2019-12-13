import re
import pickle
import networkx as nx
from Levenshtein import distance as lev_dist
from .cx import Construction
from .build import CXbuilder, CXbuilderTF

class Subphrases(CXbuilderTF):
    """Class for building time phrase constructions."""
    
    def __init__(self, wordcxs, semdist, tf, **kwargs):
        
        """Initialize with Constructions attribs/methods."""
        CXbuilderTF.__init__(self, tf, **kwargs)
        
        self.words = wordcxs
        self.sdist = semdist
        # get maximum semantic distance value from vector space
        self.max_dist = max((
            semdist[lex1][lex2] for lex1, lexs in semdist.items()
                for lex2 in lexs
        ))
        
        # map cx searches for full analyses
        self.cxs = (
            self.defi,
            self.card_chain,
            self.demon,
            self.adjv,
            self.advb,
            self.attrib,
            self.geni,
            self.numb,
            self.prep,
            self.appo,
            self.appo_name,
        )
        
        self.dripbucket = (
            self.wordphrase,
        )
        
        self.kind = 'subphrase'
        
        appo_yield = {
            'numb_ph', 'attrib_ph', 
            'appo', 'appo_name',
            'geni_ph',
        }
        
        # submit these cxs to cx in set 
        self.yieldsto = {
            'card_chain': {'numb_ph'},
            'word_cx': {self.kind},
            'appo': appo_yield,
            'appo_name': appo_yield,
        }
        
    def word(self, w):
        """Safely get word CX"""
        return self.words.get(w, Construction())
        
    def wordphrase(self, w):
        """A phrase construction for one word.
        
        Returns first matching word cx for a word.
        """
        return self.word(w)
        
    def getindex(self, indexable, index, default=None):
        """Safely get an index on an item"""
        try:
            return indexable[index]
        except IndexError:
            return default
        
    def defi(self, w):
        """Matches a definite construction."""
        
        P = self.getP(w)
        
        return self.test( 
            {
                'element': w,
                'name': 'defi_ph',
                'kind': self.kind,
                'roles': {'art': self.word(w), 'head': self.word(P(1))},
                'conds': {

                    f'F.sp.v({w}) == art':
                        self.F.sp.v(w) == 'art',

                    'bool(P(1))':
                        bool(P(1))
                }
            }
        )
    
    def prep(self, w):
        """Matches a preposition with a modified element."""
                
        P = self.getP(w)
        Wk =  self.getWk(w)
        F = self.F
        
        return self.test(
            {
                'element': w,
                'name': 'prep_ph',
                'kind': self.kind,
                'roles': {'prep':self.word(w), 'head':self.word(P(1))},
                'conds': {

                    f'({w}).name == prep':
                        self.word(w).name == 'prep',

                    f'F.prs.v({w}) == absent':
                        self.F.prs.v(w) == 'absent',
                    
                    'bool(P(1))':
                        bool(P(1)),
                }
            },
            {
                'element': w,
                'name': 'prep_ph',
                'pattern': 'suffix',
                'kind': self.kind,
                'roles': {'prep': self.word(w), 'head': self.word(w)},
                'conds': {
                    
                    f'({w}).name == prep':
                        self.word(w).name == 'prep',
                    
                    'F.prs.v(w) not in {absent, NA}':
                        F.prs.v(w) not in {'absent', 'NA'},
                }
                
            },
            {
                'element': w,
                'name': 'prep_ph',
                'pattern': 'prep...on',
                'kind': self.kind,
                'roles': {'prep': self.word(w), 'head': self.word(w)},
                'conds': {
                    f'{F.lex.v(w)} in lexset':
                        F.lex.v(w) in {'M<L/', 'HL>H'},
                    f'Wk.back(({w}).name == prep)':
                        bool(Wk.back(lambda n: self.word(n).name=='prep'))
                }
                
            }
        )
        
    def geni_cl(self, wordnode):
        """Match genitive clauses.
        
        To be used in conjunction with self.geni
        """
        E, F, L, = self.E, self.F, self.L
        geni_cl = next(
            (n for n in E.mother.t(wordnode)
                if F.rela.v(n) in {'RgRc'}),
            0
        )
        return self.test(
            { 
                'element': geni_cl,
                'name': 'clause',
                'kind': 'tf_node',
                'roles': {f'{w}':w for w in L.d(geni_cl,'word')},
                'conds': {
                    'genitive edge relation found':
                        bool(geni_cl)
                 }
            }
        )


    def geni(self, w):
        """Queries for "genitive" relations on a word."""
        
        P = self.getP(w)
        word = self.word
        geni_cl = self.geni_cl(w)    
    
        return self.test(
            {
                'element': w,
                'name': 'geni_ph',
                'kind': self.kind,
                'roles': {'geni': self.word(w), 'head': self.word(P(-1))},
                'conds': {

                    'P(-1, st) == c': 
                        P(-1,'st') == 'c',

                    'P(-1).name not in {qquant,card}':
                        word(P(-1)).name not in {'qquant','card'},
                    
                    'P(-1).name != prep':
                        word(P(-1)).name != 'prep',
                }
            },
            {
                'element': w,
                'name': 'geni_ph',
                'pattern': 'clause',
                'kind': self.kind,
                'roles': {'geni':geni_cl, 'head': self.word(w)},
                'conds': {
                    'genitive edge found':
                        bool(geni_cl)
                }
            },
        )

#     def advb_cl(self, w):
#         """Match adverb clauses.
        
#         To be used in conjunction with self.advb.
#         """
#         P = self.getP(w)
#         F = self.F
#         P_cl = self.getP(w, 'clause')
#         Wk_cl = self.getWk(w, 'clause')
#         clause = self.L.u(w,'clause')[0]
        
#         # collect subsequent words,
#         # separate verb to mark it as the head, if present
#         # else, the "head" will be the first word
#         is_verb = lambda n: self.F.pdp.v(n) == 'verb'
#         ahd_vb = Wk_cl.ahead(
#             is_verb,
#         )
#         ahd_wd  = Wk_cl.ahead(
#             lambda n: not is_verb(n),
#             every=True,
#             default=[]
#         )
#         # assign roles
#         roles = {}
#         if ahd_vb:
#             roles['head'] = ahd_vb
#         else:
#             roles['head'] = ahd_wd.pop(0) if ahd_wd else None
#         for aword in ahd_wd:
#             roles[aword] = aword
            
#         # set of lexemes which produce a 
#         # dependent clause
#         dep_advbs = {
#             '>Z', '>XR/',
#             'VRM/'
#         }
        
#         return self.test(
#             {
#                 'element': clause,
#                 'name': 'clause',
#                 'kind': 'tf_node',
#                 'roles': roles,
#                 'conds': {
#                     'lex(word) in lex set':
#                         F.lex.v(w) in dep_advbs,
#                     'no P1 in phrase':
#                         not P(1),
#                     'no P-1 in phrase':
#                         not P(-1),
#                     'there are words ahead in clause':
#                         bool(ahd_vb) or bool(ahd_wd),
#                     'no suffix on adverb':
#                         self.F.prs.v(w) in {'n/a', 'absent'},
#                 }
#             }
#         )
        
    def advb(self, w):
        """Match and adverb and its mod."""
        
        P = self.getP(w)
        F = self.F
        word = self.word
        name = 'advb'
        # advb_cl = self.advb_cl(w)
        
        return self.test(
           {
                'element': w,
                'name': name,
                'kind': self.kind,
                'roles': {'advb': word(w), 'head': word(P(1))},
                'conds': {
                    f'F.sp.v({w}) == advb':
                        self.F.sp.v(w) == 'advb',
                    'P(-1,sp) != art':
                        P(-1,'sp') != 'art',
                    'bool(P(1))':
                        bool(P(1)),
                    'P(1,sp) != conj': # ensure not a nominal use
                        P(1,'sp') != 'conj',
                    'P(-1).name != prep': # ensure not nominal
                        word(P(-1)).name != 'prep',
                    f'F.lex.v({F.lex.v(w)}) not in noadvb_set':
                        F.lex.v(w) not in {'JWMM'},
                }
            },               
           {
                'element': w,
                'name': name,
                'pattern': 'advb_lex',
                'kind': self.kind,
                'roles': {'advb': word(w), 'head': word(P(1))},
                'conds': {
                    'F.lex.v(w) in lex set':
                        F.lex.v(w) in {'L>', '>Z', '<WD/'},
                    'name(P1) in {cont, art}':
                        word(P(1)).name in {'cont','art'},
                }
            },
#            {
#                 'element': w,
#                 'name': name,
#                 'pattern': 'advb_clause',
#                 'kind': self.kind,
#                 'roles': {'advb': word(w), 'head': advb_cl},
#                 'conds': { 
#                     'an adverb clause has been found (self.advb_cl)':
#                         bool(advb_cl)
#                 }
#             },
        )
    
    def adjv(self, w):
        """Matches a word serving as an adjective."""
        
        P = self.getP(w)
        F = self.F
        word = self.word
        name = 'adjv_ph'
        
        # check for recursive adjective matches 
        a2match = self.adjv(P(-1)) if P(-1) else Construction()
        a2match_head = int(a2match.getrole('head', 0))
        
        common = {
            
            'w.name not in {qquant,card}':
                word(w).name not in {'qquant','card'},
            
            'P(-1).name == cont':
                word(P(-1)).name == 'cont',
                        
            'P(-1, st) & {NA, a}': 
                P(-1,'st') in {'NA', 'a'},   
            
            'P(-1).name != quant':
                word(P(-1)).name != 'quant',
            
            'P(-1).name != prep':
                word(P(-1)).name != 'prep',
        }
                
        return self.test(
            
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'adjv (1x)',
                'roles': {'adjv':word(w), 'head': word(P(-1))},
                'conds': dict(common, **{
                    'F.sp.v(w) in {adjv, verb}':
                        F.sp.v(w) in {'adjv', 'verb'},
                })
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'adjv (2x)',
                'roles': {'adjv': word(w), 'head': word(a2match_head)},
                'conds': dict(common, **{
                    
                    'F.sp.v(w) in {adjv, verb}':
                        F.sp.v(w) in {'adjv', 'verb'},
                    
                     'self.adjv(P(-1)) and target != P(0)':
                        bool(a2match) and a2match_head != P(0)
                })
            },
        )
     
    def attrib(self, w):
        """Identify elements in a attrib construction.
        
        In Hebrew this construction typically consists of four slots:
            > ה + A + ה + B
        Attrib identifies each of these elements and labels them.
        A is assumed to be the head, or modified, element and B
        is assumed to be an adjectival element.
        """
                
        # CX consists of two constituent cxs
        # start walk from head of first match
        P = self.getP(w)
        defi1 = self.defi(w)
        d1head = int(defi1.getrole('head', 0))
        Wk = self.getWk(d1head)

        # walk to next valid defi match
        # and allow adjectives to intervene:
        defi2 = Wk.ahead(
            lambda n: self.defi(n),
            go=lambda n: self.F.sp.v(n)=='adjv',
            output=True
        ) if Wk else Construction()
        defi2 = defi2 or Construction()

        # check for single_defi (only two cases)
        defi_p1 = self.defi(P(1))
        
        return self.test(
            {
                'element': w,
                'name': 'attrib_ph',
                'pattern': 'double_defi',
                'kind': self.kind,
                'roles': {'head': defi1, 'attrib': defi2},
                'conds': {
                    'bool(defi1)':
                        bool(defi1),
                    'bool(defi2)':
                        bool(defi2), 
                }
            },
            {
                'element': w,
                'name': 'attrib_ph',
                'pattern': 'single_defi',
                'kind': self.kind,
                'roles': {'head': self.word(w), 'attrib': defi_p1},
                'conds': {
                    'name(w) == cont':
                        self.word(w).name == 'cont',
                    'F.st.v(w) == a':
                        self.F.st.v(w) == 'a',
                    'P(-1,lex) != H':
                        P(-1,'lex') != 'H',
                    'bool(defi_p1)':
                        bool(defi_p1),
                }
            }
        )
        
    def numb(self, w):
        """Defines numerical relations with an non-quant word.
        
        Often but not always indicates quantification as other
        semantic relations are possible.
        """

        P = self.getP(w)
        Wk = self.getWk(w)
        word = self.word
        is_nom = (
            lambda n: word(n).name == 'cont'
        )
        
        # for the quant ahead check
        # should stop at a preposition or another quantifier
        stop_ahead = (
            lambda n: (word(n).name == 'prep'
                or word(n).name in {'card', 'qquant'} and word(n).name != word(w).name)
        )
        
        behind_nom = Wk.back(is_nom, stop=lambda n: not is_nom(n)) 
        
        return self.test(
        
            {
                'element': w,
                'name': 'numb_ph',
                'kind': self.kind,
                'pattern': 'numbered forward',
                'roles': {'numb': word(w), 'head': word(P(1))},
                'conds': {
                    
                    'w.name in {qquant,card}':
                     word(w).name in {'qquant', 'card'},
                    
                    'bool(P(1))':
                        bool(P(1)),
                    
                    'P(1,sp) != conj':
                        P(1,'sp') != 'conj',
                    
                    'P(1).name not in {qquant,card,prep}':
                        word(P(1)).name not in {'qquant','card','prep'},
        
                    'P(-1,sp) != art':
                        P(-1,'sp') != 'art',
                },
            },  
            {
                'element': w,
                'name': 'numb_ph',
                'kind': self.kind,
                'pattern': 'numbered backward',
                'roles': {'numb': word(w), 'head': word(behind_nom)},
                'conds': {
                    
                    'w.name in {qquant,card}':
                        word(w).name in {'qquant','card'},
                    
                    'not Wk.ahead(is_nominal)':
                        not Wk.ahead(is_nom, stop=stop_ahead),
                    
                    'bool(Wk.back(is_nominal))':
                        bool(behind_nom),
                    
                    'F.st.v(behind_nom) in {a, NA}':
                        self.F.st.v(behind_nom) in {'a', 'NA'},
                }
            }
        )
        
    def card_chain(self, w):
        """Defines cardinal number chain constructions"""
        
        P = self.getP(w)
        F = self.F
        word = self.word
        
        return self.test(
            {
                'element': w,
                'name': 'card_chain',
                'kind': self.kind,
                'pattern': 'adjacent',
                'roles': {'card':word(w), 'head':word(P(-1))},
                'conds': {
                    
                    'F.ls.v(w) == card':
                        F.ls.v(w) == 'card',
                    'P(-1,ls) == card':
                        P(-1,'ls') == 'card',                    
                }
            },
            {
                'element': w,
                'name': 'card_chain',
                'kind': self.kind,
                'pattern': 'conjunctive',
                'roles': {'card': word(w), 'head': word(P(-2)), 'conj': word(P(-1))},
                'conds': {
                    'F.ls.v(w) == card':
                        F.ls.v(w) == 'card',
                    'P(-1,lex) == W':
                        P(-1,'lex') == 'W',
                    'P(-2,ls) == card':
                        P(-2,'ls') == 'card',   
                }
            }
        )
    
    def demon(self, w):
        """Defines an adjacent demonstrative construction."""
        
        P = self.getP(w)
        word = self.word
        F = self.F
        name = 'demon_ph'
        
        return self.test(
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'adjacent forward',
                'roles': {'demon': word(w), 'head': word(P(1))},
                'conds': {
                    'prde in {F.pdp.v(w), F.sp.v(w)}':
                        'prde' in {F.pdp.v(w), F.sp.v(w)},
                    
                    'P(-1,sp) != art': # ensure not part of attrib pattern
                        P(-1,'sp') != 'art',
                    
                    'P(-1).name != prep':
                        word(P(-1)).name != 'prep',
                    
                    'bool(P(1))':
                        bool(P(1)),
                    
                    'P(1).name == cont':
                        word(P(1)).name == 'cont',
                }
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'adjacent back',
                'roles': {'demon':word(w), 'head':word(P(-1))},
                'conds': {
                    'prde in {F.pdp.v(w), F.sp.v(w)}':
                        'prde' in {F.pdp.v(w), F.sp.v(w)},
                    
                    'P(-1).name not in {prep,qquant,card}':
                        word(P(-1)).name not in {'prep','qquant','card'},
                    
                    'P(-1,sp) == subs':
                        P(-1,'sp') == 'subs',
                }
            }
        )
    
    def get_distance(self, w1, w2, default=None):
        """Retrieve semantic distance between two word nodes."""
        default = default or self.max_dist
        lex1, lex2 = self.F.lex.v(w1), self.F.lex.v(w2)
        return self.sdist.get(lex1,{}).get(lex2, default)
        
                'conds': {
                    
                    'name(w) == cont':
                        wd(w).name == 'cont',
                    
                    'not adjv(w)':
                        not self.adjv(w),
                    'not advb(w)':
                        not self.advb(w),
                    
                    'name(P-1) == cont':
                        wd(P(-1)).name == 'cont',
                    
                    'st(P-1) == a':
                        F.st.v(P(-1)) == 'a',
                    
                    f'{dist} < 1.2 or {ldist} < 2':
                        (dist < 1.2) or (ldist < 2)
                }
            }
        )
    
    def appo_name(self, w):
        """Match an apposition of name"""
    
        name = 'appo_name'
        F = self.F
        P = self.getP(w)
        Wk = self.getWk(w)
        wd = self.word
        
        # get word back with only intervention of article
        bk = Wk.back(
            lambda n: wd(n).name == 'name',
            go=lambda n: wd(n).name == 'art'
        )
        
        return self.test(
        
            {
                'element': w,
                'name': name,
                'pattern': 'name_entity',
                'kind': self.kind,
                'roles': {'head': wd(w), 'name': wd(bk)},
                'conds': {
                    
                    'name(w) == cont':
                        wd(w).name == 'cont',
                    
                    'name(back) == name':
                        wd(bk).name == 'name',
                    
                    'F.st.v(back) == a':
                        F.st.v(bk) == 'a',
                    
                    f'F.nu.v({w}) == F.nu.v({bk}) or lex exception':
                        (F.nu.v(w) == F.nu.v(bk)) or F.lex.v(w) in {'>LHJM/'},
                    
                    # NB:
                    # rule below reveals the need to be able to say
                    # what head_slot should be; i.e., the lexeme should
                    # be semantically consistent with the ID of the proper name
                    # of a person, head_slot should ~ person, etc.
                    # but for now I'll use a work-around solution
                    'F.lex.v(w) not in timeword set':
                        F.lex.v(w) not in {'CNH/'}
                },
            },
            {
                'element': w,
                'name': name,
                'kind': self.kind,
                'pattern': 'entity_name',
                'roles': {'head': wd(w), 'name': wd(P(1))},
                'conds': {
                    
                    'name(w) == cont':
                        wd(w).name == 'cont',
                    
                    'name(P1) == name':
                        wd(P(1)).name == 'name',
                    
                    'F.st.v(w) == a':
                        F.st.v(w) == 'a',
                    
                    f'F.nu.v({w}) == F.nu.v({P(1)}) or lex exception':
                        (F.nu.v(w) == F.nu.v(P(1))) or F.lex.v(w) in {'>LHJM/'},
                },
            },
        )

class Phrases(CXbuilder):
    """Build complete phrase constructions."""
    
    def __init__(self, phrase2cxs, semdist, tf):
        CXbuilder.__init__(self)
        
        # set up tf methods
        self.tf = tf
        self.F, self.T, self.L = tf.api.F, tf.api.T, tf.api.L
        
        # map cx to phrase node for context retrieval
        self.cx2phrase = {
            cx:ph 
                for ph in phrase2cxs
                    for cx in phrase2cxs[ph]
        }
        
        self.phrase2cxs = phrase2cxs
        self.semdists = semdist
        
        self.cxs = (        
            self.appo,
            self.coord,
        )
        self.dripbucket = (
            self.cxph,
        )
        
#         self.yieldsto = {
#             'appo': True,
#             'coord': True,
#         }
        
        self.kind = 'phrase'
        
    def cxph(self, cx):
        """Dripbucket function that returns cx as is."""
        return cx
        
    def get_context(self, cx):
        """Get context for a given cx."""
        phrase = self.cx2phrase.get(cx, None)
        if phrase:
            return self.phrase2cxs[phrase]
        else:
            return tuple()
        
    def getP(self, cx):
        """Index positions on phrase context"""
        positions = self.get_context(cx)
        if positions:
            return Positions(
                cx, positions, default=Construction()
            ).get
        else:
            return Dummy

    def getWk(self, cx):
        """Index walks on phrase context"""
        positions = self.get_context(cx)
        if positions:
            return Walker(cx, positions)
        else:
            return Dummy()
    
    def getindex(
        self, indexable, index, 
        default=Construction()
    ):
        """Safe index on iterables w/out IndexErrors."""
        try:
            return indexable[index]
        except:
            return default
    
    def getname(self, cx):
        """Get a cx name"""
        return cx.name
    
    def getkind(self, cx):
        """Get a cx kind."""
        return cx.kind
    
    def getsuccrole(self, cx, role, index=-1):
        """Get a cx role from a list of successive roles.
        
        e.g.
        [big_head, medium_head, small_head][-1] == small_head
        """
        cands = list(cx.getsuccroles(role))
        try:
            return cands[index]
        except IndexError:
            return Construction()
    
    def string_plus(self, cx, plus=1):
        """Stringifies a CX + N-slots for Levenshtein tests."""
        
        # get all slots in the context for plussing
        allslots = sorted(set(
            s for scx in self.get_context(cx)
                for s in scx.slots
        ))
        
        # get plus slots
        P = (Positions(self.getindex(cx.slots, -1), allslots).get
                 if cx.slots and allslots else Dummy)
        plusses = []
        for i in range(plus, plus+1):
            plusses.append(P(i,-1)) # -1 for null slots (== empty string in T.text)
        plusses = [p for p in plusses if type(p) == int]
        
        # format the text string for Levenshtein testing
        ptxt = T.text(
            cx.slots + tuple(plusses),
            fmt='text-orig-plain'
        ) if cx.slots else ''
        
        return ptxt

    def rank_candidates(self, cx, cx_patterns=[]):
        """Ranks preceding phrases on likelihood of a relationship
        
        TODO: Give a thorough explanation
        """
        
        F, T = self.F, self.T
        P = self.getP(cx)
        semdist = self.semdists
        Wk = self.getWk(cx)
                         
        # get all top-level cxs behind this one that match in name
        cx_behinds = Wk.back(
            lambda c: c.name == cx.name,
            every=True,
            stop=lambda c: (
                c.name == 'conj' and (c != P(-1))
            )
        )
        
        # if top level phrases produce no results,
        # use subphrases instead
        if not cx_behinds:
            topcontext = self.get_context(cx)
            
            # gather all valid subphrase candidates
            subcontext = []
            for topcx in topcontext:
                for subcx in topcx.subgraph():
                    if type(subcx) == int: # skip TF slots
                        continue
                    if (
                        subcx in topcontext or subcx.name != 'conj'
                        and subcx not in cx
                    ):
                        subcontext.append(subcx)        
            
            # walk the new candidates
            Wk2 = Walker(cx, subcontext)
            cx_behinds = Wk2.back(
                lambda c: c.name != 'conj', 
                default=[P(-2)],
                every=True,
                stop=lambda c: (
                    c.name == 'conj' and (c != P(-1))
                )
            )
        
        # map each back-cx to its last slot to make sure
        # every candidate is the last item in its phrase
        # check is made in next series of lines
        cx2last = {
            cxb:self.getindex(sorted(cxb.slots), -1, 0)
                for cxb in cx_behinds
        }
        
        # find coordinate candidate subphrases that stand
        # at the end of the phrase
        cx_subphrases = []
        
        for cx_back in cx_behinds:
            for cxsp in cx_back.subgraph():
                if type(cxsp) == int:
                    continue
                elif (
                    cx2last[cx_back] in cxsp.slots
                    and cxsp.getrole('head')
                ):
                    cx_subphrases.append(cxsp)
        
        # get subphrase heads for semantic tests
        cx2heads = [
            (cxsp, self.getsuccrole(cxsp,'head'))
                for cxsp in cx_behinds
        ]

        # get head of this cx
        head1 = self.getsuccrole(cx,'head')     
        head1lex = F.lex.v(head1)
        
        # sort on a set of priorities
        # the default sort behavior is used (least to greatest)
        # thus when a bigger value should be more important, 
        # a negative is added to the number
        stringp = self.string_plus
        
        # arrange candidates by priority
        cxpriority = []
        for cxsp, headsp in cx2heads:
            name_eq = 0 if cxsp.name == cx.name else 1
            semantic_dist = semdist.get(
                head1lex,{}
            ).get(F.lex.v(headsp), np.inf)
            size = -len(cxsp.slots)
            levenshtein = lev_dist(stringp(cx), stringp(cxsp))
            slot_dist = -next(iter(cxsp.slots), 0)
            heads = (head1, headsp) # for reporting purposes only
            
            cxpriority.append((
                name_eq,
                semantic_dist,
                size,
                levenshtein,
                slot_dist,
                heads,
                cxsp
            ))
            
        # make the sorting
        candidates = sorted(cxpriority, key=lambda k: k[:-1])
        
        # select the first priority candidate
        cand = next(iter(candidates), (0,0,Construction()))
        
        # add data for conds report / debugging
        stats = collections.defaultdict(str)
        for namescore,sdist,leng,ldist,lslot,heads,cxp in candidates:
            # name equality
            stats['namescore'] += f'\n\t{cxp} namescore: {namescore}'
            # semantic distance
            stats['semdists'] += (
                f'\n\t{round(sdist, 2)}, {F.lex.v(heads[0])} ~ {F.lex.v(heads[1])}, {cxp}'
            )
            # size of cx
            stats['size'] += f'\n\t{cxp} length: {abs(leng)}'
            
            # Levenstein distance
            stats['ldist'] += f'\n\t{cxp} dist: {ldist}'
            
            # dist of last slot
            stats['lslot'] += f'\n\t{cxp} last slot: {abs(lslot)}'
    
        return (candidates, cand, stats)
    
    def coord(self, cx):
        """A coordinate construction.
        
        In order to match a coordinate cx, we need to determine
        which item in the previous phrase this cx belongs with. 
        This is done using a semantic vector space, which can
        quantify the approximate semantic distance between the
        heads of this cx and a candidate cx.
        
        Criteria utilized in validating a coordinate cx between
        an origin cx and a candidate cx are the following:
            TODO: fill in
        """
        
        P = self.getP(cx)        
        cands, cand, stats = self.rank_candidates(cx)
        
        return self.test(
            {
                'element': cx,
                'name': 'coord',
                'kind': self.kind,
                'roles': {'part2':cx, 'conj': P(-1), 'part1': cand[-1]},
                'conds': {
                    'P(-1).name == conj':
                        P(-1).name == 'conj',
                    'bool(cand)':
                        bool(cand[-1]),
                    f'name matches {stats["namescore"]}\n':
                        bool(cands),
                    f'is shortest sem. distance of {stats["semdists"]}\n':
                        bool(cands),
                    f'is longest length of: {stats["size"]}\n':
                        bool(cands),
                    f'is shortest Levenshtein distance: {stats["ldist"]}\n':
                        bool(cands),
                    f'is closest last slot of: {stats["lslot"]}\n':
                        bool(cands)
                }
            },
        )
    
    
    def appo(self, cx):
        """Find appositional cxs"""
        
        P = self.getP(cx)
        cands, cand, stats = self.rank_candidates(cx)
                
        return self.test(
            {
                'element': cx,
                'name': 'appo',
                'pattern': 'NP',
                'kind': self.kind,
                'roles': {'appo':cx, 'head': cand[-1]},
                'conds': {
                    'name(cx) not in not_NPset':
                        cx.name not in {'prep_ph','conj'},
                    'P(-1).name != conj':
                        P(-1).name != 'conj',
                    'bool(cand)':
                        bool(cand[-1]),
                    f'name matches {stats["namescore"]}\n':
                        bool(cands),
                    f'is shortest sem. distance of {stats["semdists"]}\n':
                        bool(cands),
                    f'is longest length of: {stats["size"]}\n':
                        bool(cands),
                    f'is shortest Levenshtein distance: {stats["ldist"]}\n':
                        bool(cands),
                    f'is closest last slot of: {stats["lslot"]}\n':
                        bool(cands),
                }
            },
            {
                'element': cx,
                'name': 'appo',
                'pattern': 'PP',
                'kind': self.kind,
                'roles': {'appo':cx, 'head': cand[-1]},
                'conds': {
                    'name(cx) == prep':
                        cx.name == 'prep_ph',
                    'P(-1).name != conj':
                        P(-1).name != 'conj',
                    'bool(cand)':
                        bool(cand[-1]),
                    f'name matches {stats["namescore"]}\n':
                        bool(cands),
                    f'is shortest sem. distance of {stats["semdists"]}\n':
                        bool(cands),
                    f'is longest length of: {stats["size"]}\n':
                        bool(cands),
                    f'is shortest Levenshtein distance: {stats["ldist"]}\n':
                        bool(cands),
                    f'is closest last slot of: {stats["lslot"]}\n':
                        bool(cands)
                }
            }        
        )
    
    
    def adjacent(self, cx):
        """Find adjacent CXs"""
        
        P = self.getP(cx)
        
        return self.test(
            {
                'element': cx,
                'name': 'appo',
                'kind': self.kind,
                'roles': {'head':cx, 'appo':P(1)},
                'conds': {
                    'cx.name != conj':
                        cx.name != 'conj',
                    'P(1).name != prep':
                        P(1).name != 'prep',
                    'bool(P(1))':
                        bool(P(1)),
                    f'name({P(1).name}) not in (conj, prep_ph)':
                        P(1).name not in {'conj','prep_ph'},
                }
            }
        )
