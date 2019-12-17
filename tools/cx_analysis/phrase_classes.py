"""
This module contains CXBuilders for
phrase classification of time adverbials.
"""

import collections
from .build import CXbuilder
from .cx import Construction

class SinglePhrase(CXbuilder):
    """Modify cx classifications for single phrase CXs"""
    
    def __init__(self, cxset, tf):
        CXbuilder.__init__(self) # initialize with standard CXbuilder methods
        
        self.cxset = cxset
        self.api = tf
        self.F, self.L = tf.api.F, tf.api.L
        
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
        self.prereq = self.single
        self.kind = 'time_class'
        
        self.class2cx = collections.defaultdict(set)
        
    def test_result(self, test, *cases):
        """Add class attributes to CX results"""
        if test:
            result = test[-1]
            cx = result['element']
            classi= result['class']
            cx.__dict__.setdefault('classification', []).extend(classi)
            cx.match = result
            cx.conds = result['conds']
            cx.cases = (result,) + cx.cases
            return cx
        else:
            return Construction(cases=cases, **cases[0])
    
    def findall(self, element):
        """Find all results with prerequisite
        
        NB this version of findall only returns
        a single result: the construction object
        itself, since it is modified in-place.
        This version expects cx tuples with
        only one cx.
        """
        results = []
        if self.prereq(element):
            for funct in self.cxs:
                cx = funct(element)
                if cx:
                    results.append(cx)
        if results:
            return results[0] # NB, only 1st matters as all are same obj
        else:
            return None
    
    def label_cxs(self):
        """Run all queries against dataset"""
        for cxtuple in self.cxset:
            cx = self.findall(cxtuple)
            if cx:
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
    
    def single(self, cxtuple):
        """Tag CXs as singles"""
        cx1 = cxtuple[0]
        relas = set(
            self.geta(c,'name') for c in cx1
        )
        bhsa_phrase = L.u(cx1.slots[0], 'phrase')[0]
        attr_cl = E.mother.t(bhsa_phrase)
        
        return self.test(
            {
                'element': cxtuple[0],
                'class': ['single'],
                'kind': self.kind,
                'conds': {
                    'len(cxtuple) == 1':
                        len(cxtuple) == 1,
                    'no apposition in cx':
                        not relas & {'appo'},
                    'no attributive clause on phrase':
                        not attr_cl
                }
            }
        )
    
    def prep(self, cxtuple):
        """Tag prepositional cxs"""
        cx = cxtuple[0]
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

    def bare(self, cxtuple):
        """Tag bare, non-modified cxs"""
        F = self.F
        cx = cxtuple[0]
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
    
    def definite(self, cxtuple):
        """A definite phrase"""
        cx = cxtuple[0]
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
    
    def def_appo(self, cxtuple):
        """Definite apposition"""
        
        F = self.F
        geta = self.geta
        cx = cxtuple[0]
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
    
    def genitive(self, cxtuple):
        """Genitive relation on head"""
        cx = cxtuple[0]
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
    
    def quantified(self, cxtuple):
        """Find quantified time phrases"""
        cx = cxtuple[0]
        head = self.get_headword(cx)
        quant_ph = self.get_head_modi(head, cx, 'numb_ph')
        geta = self.geta
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
                    'cx.name == card_chain':
                        cx.name == 'card_chain'
                }
            }
        )
    
    def adjective(self, cxtuple):
        """Adjectival modifications via non-definite apposition"""
        cx = cxtuple[0]
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