"""
This module contains CXBuilders for
phrase classification of time adverbials.
"""

import collections
from cx_analysis.build import CXbuilder
from cx_analysis.cx import Construction

class SinglePhrase(CXbuilder):
    """Modify cx classifications for single phrase CXs
    
    Arguments:
        cxset: a set of Construction objects
        goodheads: A set of acceptable head lexemes for classifying.
            The subphrase tuples fed into the Builder at
            this point contain some subphrases that themselves do
            not function as a timephrase, but only in conjunction 
            with one. This includes phrases that relate events to 
            people (e.g. למלך), for example. At present, the best
            way to avoid classifying non-timephrases is to only
            accept those whose heads function in stand-alone phrases
            (i.e. tuples with only 1 subphrase), or whose heads have
            been manually added to the goodheads set.
    """
    
    def __init__(self, cxset, goodheads, tf):
        CXbuilder.__init__(self) # initialize with standard CXbuilder methods
        
        self.cxset = cxset
        self.api = tf
        self.F, self.E, self.L = tf.api.F, tf.api.E, tf.api.L
        
        # cx queries
        # NB: order matters!
        self.cxs = (
            self.prep,
            self.bare,
            self.definite,
            self.def_appo,
            self.genitive,
            self.quantified,
            self.adjective,
        )
        self.goodheads = goodheads
        self.prereq = self.single
        self.kind = 'time_class'
        self.class2cx = collections.defaultdict(set)
        
    def test_result(self, test, *cases):
        """Add class attributes to CX results.
        
        Unlike other test_result methods, this one
        does not return anything. It only operates
        on the existing CXs in-place.
        """
        if test:
            result = test[-1]
            cx = result['element']
            classi= result['class']
            cx.__dict__.setdefault('classification', []).extend(classi)
            cx.match = result
            cx.conds = result['conds']
            cx.cases = (result,) + cx.cases
            return cx # return to show conditions
        else:
            return Construction(cases=cases, **cases[0])
    
    def findall(self, element):
        """Find all results with prerequisite
        
        Iterates through a cx_data tuple and 
        yields results.
        """        
        # tag all CXs stored in the tuple element
        for cx in element:
            self.prereq(cx, element) # run prerequisite tagging
            if 'single' in cx.__dict__.get('classification', {}):
                for funct in self.cxs:
                    funct(cx) # tag the CX
            else:
                self.not_single(cx) # apply non-single categories
    
    def label_cxs(self):
        """Run all queries against dataset
        
        Nothing is returned since CXs are 
        tagged in-place.
        """
        
        # tag eligible CXs
        for cxtuple in self.cxset:
            self.findall(cxtuple)
            
        # organize classified CXs by tags
        for cxtuple in self.cxset:
            for cx in cxtuple:
                if cx.__dict__.get('classification'):
                    for tag in cx.classification:
                        self.class2cx[tag].add(cx)
    
    def geta(self, item, attrib, default=None):
        """Safely retrieve attribute from object
        
        Some objects in a CX graph are TF integer
        nodes, while most are CX objects. In order
        to safely call attributes on a given position,
        we need to handle attribute errors when called
        on an integer.
        """
        try:
            return item.__dict__[attrib]
        except AttributeError:
            return default
    
    def get_headword(self, cx):
        """Get a word that serves as head"""
        head = list(cx.getsuccroles('head'))[-1]
        return head
    
    def get_head_modi(self, head, cx, name, default=Construction()):
        """Retrieve a modifier on a particular head"""
        for c in cx.graph:
            if (self.geta(c,'name') == name) and (head in c):
                return c
        # unsuccessful search
        return default
    
    def single(self, cx, cxtuple):
        """Tag CXs as singles.
        
        NB: single phrase CXs can be 
        parts of more complex multi-phrase
        CXs. This method tags cases where that
        both is and is not true.
        """
        
        F, E, L = self.F, self.E, self.L
        relas = set(
            self.geta(c,'name') for c in cx
        )
        headpath = list(cx.getsuccroles('head'))
        head = headpath[-1]
        head_cx = next(iter(cx.graph.pred[head]))
        bhsa_phrase = L.u(cx.slots[0], 'phrase')[0]
        attr_cl = E.mother.t(bhsa_phrase)
        
        custom_goodhead = {
            'XG/',
        }
        
        return self.test(
            {
                'element': cx,
                'class': ['single'],
                'kind': self.kind,
                'conds': {
                    'len(cxtuple) == 1':
                        len(cxtuple) == 1,
                    'no apposition in cx':
                        not relas & {'appo'},
                    'no attributive clause on phrase':
                        not attr_cl,
                }
            },
            {
                'element': cx,
                'class': ['single', 'component'],
                'kind': self.kind,
                'conds': {
                    'len(cxtuple) > 1':
                        len(cxtuple) > 1,
                    'head(cx) is good':
                        F.lex.v(head) in self.goodheads|custom_goodhead,
                    'no apposition in cx':
                        not relas & {'appo'},
                    'no attributive clause on phrase':
                        not attr_cl,
                }
            },
            {
                'element': cx,
                'class': ['single', 'component'],
                'kind': self.kind,
                'conds': {
                    'len(cxtuple) > 1':
                        len(cxtuple) > 1,
                    'name(head_cx) in goodset':
                        head_cx.name in {'card',},
                    'no apposition in cx':
                        not relas & {'appo'},
                    'no attributive clause on phrase':
                        not attr_cl,
                }
            },
        )
    
    def not_single(self, cx):
        """Tag and track CXs that are not single timephrases.
        
        This is the drip-bucket category, only executed
        if self.single fails. See self.findall for code.
        """
        return self.test(
            {
                'element': cx,
                'class': ['not_single'],
                'kind': cx.kind, # preserve the kind
                'conds': {
                    'self.single failed':
                        True,
                },
            },
        )
    
    def prep(self, cx):
        """Tag prepositional cxs"""
                
        return self.test(
            {
                'element': cx,
                'class': ['prep'],
                'kind': self.kind,
                'conds': {
                    'cx.name == prep_ph':
                        cx.name == 'prep_ph',
                }
            },
            {
                'element': cx,
                'class': ['øprep'],
                'conds': {
                    'cx.name != prep_ph':
                        cx.name != 'prep_ph',
                }
            }
        )

    def bare(self, cx):
        """Tag bare, non-modified cxs"""
        F, E = self.F, self.E
        head_path = list(cx.getsuccroles('head'))
        head = head_path[-1]
        etcbc_phrase = self.L.u(int(head),'phrase')[0]
        
        # two types of units allowed in the path:
        # word cxs or prep_ph
        # trace path to head and collect relations along the way
        cx_name = cx.name if cx.kind != 'word_cx' else cx.kind
        head_phs = {cx_name}
        for c in head_path:
            if self.geta(c,'kind') == 'subphrase':
                head_phs.add(c.name)
            else:
                head_phs.add('word_cx')
        
        prereqs = {
            'head_phs is subset of {word_cx, prep_ph}':
                head_phs.issubset({'word_cx', 'prep_ph', 'advb'}),
            'F.st.v(head) != c':
                F.st.v(int(head)) != 'c',
            'not daughters(etcbc_phrase)':
                not E.mother.t(etcbc_phrase),
        }
        
        return self.test(
            {
                'element': cx,
                'class': ['bare'],
                'kind': self.kind,
                'conds': dict({
                    'F.prs.v(head) in {n/a, absent}':
                        F.prs.v(int(head)) in {'n/a', 'absent'},
                }, **prereqs)
            },
            {
                'element': cx,
                'class': ['suffix'],
                'kind': self.kind,
                'conds': dict({
                    'F.prs.v(head) not in {n/a, absent}':
                        F.prs.v(int(head)) not in {'n/a', 'absent'},
                }, **prereqs)
            },
        )
    
    def definite(self, cx):
        """A definite phrase"""
        head = self.get_headword(cx)
        def_ph = self.get_head_modi(head, cx, 'defi_ph')
        
        return self.test(
            {
                'element': cx,
                'class': ['definite'],
                'kind': self.kind,
                'conds': {
                    'cx contains defi phrase with head':
                        bool(def_ph)
                }
            }
        
        )
    
    def def_appo(self, cx):
        """Definite apposition"""
        
        F = self.F
        geta = self.geta
        head = self.get_headword(cx)
        
        # get attribute cx if it contains head word
        att_ph = self.get_head_modi(head, cx, 'attrib_ph')
        
        return self.test(
            {
                'element': cx,
                'class': ['def_apposition'],
                'kind': self.kind,
                'conds': {
                    f'cx contains attrib ph with head':
                        bool(att_ph)
                }
            },
            {
                'element': cx,
                'class': ['def_apposition', 'demonstrative'],
                'kind': self.kind,
                'conds': {
                    f'cx contains attrib ph with head':
                        bool(att_ph),
                    'apposition contains demonstrative':
                        {'prde', 'prps'} & set(
                            F.pdp.v(w) for w in att_ph.getrole('attrib', Construction()).slots
                        )
                }
            },
            {
                'element': cx,
                'class': ['def_apposition', 'ordinal'],
                'kind': self.kind,
                'conds': {
                    f'cx contains attrib ph with head':
                        bool(att_ph),
                    'apposition contains ordinal':
                        'ordn' in set(
                            geta(c,'name') for c in att_ph.graph
                        ),
                }
            },
        )
    
    def genitive(self, cx):
        """Genitive relation on head"""
        head = self.get_headword(cx)
        geni_ph = self.get_head_modi(head, cx, 'geni_ph')
        geni_items = set(
            self.geta(c, 'name') for c in geni_ph
        )
        return self.test(
            {
                'element': cx,
                'class': ['genitive'],
                'kind': self.kind,
                'conds': {
                    'cx contains geni phrase on head':
                        bool(geni_ph)
                }
            },
            {
                'element': cx,
                'class': ['geni_cardinal'],
                'kind': self.kind,
                'conds': {
                    'cx contains geni phrase on head':
                        bool(geni_ph),
                    
                    'a cardinal is genitive to this word':
                        'card' in geni_items,
                }
            }
        )
    
    def quantified(self, cx):
        """Find quantified time phrases"""
        geta = self.geta
        head = self.get_headword(cx)
        headpath = list(cx.getsuccroles('head'))
        headpath_cxs = set(
            geta(cx,'name') for cx in headpath
                if geta(cx,'kind') == 'subphrase'
                or geta(cx,'name') == 'card'
        )
        quant_ph = self.get_head_modi(head, cx, 'numb_ph')
        
        return self.test(
            {
                'element': cx,
                'class': ['quantified', 'cardinal'],
                'kind': self.kind,
                'conds': {
                    'cx contains numbered phrase on head':
                        bool(quant_ph),
                    
                    'does not contain qualitative quant':
                        'qquant' not in set(
                            geta(c,'name') for c in quant_ph
                        )
                }
            },
            {
                'element': cx,
                'class': ['quantified', 'qualitative'],
                'kind': self.kind,
                'conds': {
                    'cx contains numbered phrase on head':
                        bool(quant_ph),
                    
                    'contains qualitative quant':
                        'qquant' in set(
                            geta(c,'name') for c in quant_ph
                        )
                }
            },
            {
                'element': cx,
                'class': ['cardinal'],
                'kind': self.kind,
                'conds': {
                    '{card_chain, card} & headpath':
                        {'card_chain','card'} & headpath_cxs,
                    f'only (card_chain, prep_ph) in headpath: {headpath_cxs}':
                        headpath_cxs.issubset({'card_chain', 'card', 'prep_ph'}),
                }
            },
            {
                'element': cx,
                'class': ['cardinal'],
                'kind': self.kind,
                'conds': {
                    'cx.name in {card_chain,card}':
                        cx.name in {'card_chain','card'}
                }
            }
        )
    
    def adjective(self, cx):
        """Adjectival modifications via non-definite apposition"""
        head = self.get_headword(cx)
        adjv_ph = self.get_head_modi(head, cx, 'adjv_ph')
        
        return self.test(
            {
                'element': cx,
                'class': ['adjective'],
                'kind': self.kind,
                'conds': {
                    'cx contains adjectival phrase on head':
                        bool(adjv_ph),
                }
            },
            {
                'element': cx,
                'class': ['demonstrative'],
                'kind': self.kind,
                'conds': {
                    'cx is a demonstrative phrase':
                        cx.name == 'demon_ph',
                }
            }
        )
