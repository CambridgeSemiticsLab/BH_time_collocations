import collections
from tools.langtools import PositionsTF

class Bunch(object):
    """Stores variables for shorthand and safe access.
    
    Like a dot-dictionary.
    """
    def __init__(self, vardict):
        """Initialize variables object with dict."""
        self.dict = vardict
        for k,v in vardict.items():
            setattr(self, k, v)
    def __getattr__(self, name):
        return None
    def __deepcopy__(self, memo=None):
        """Handle deep copy errors by instancing a diff Bunch object."""
        return Bunch({k:v for k,v in self.dict.items()})
    
class Construction(object):
    """A linguistic construction and its attributes."""
    
    def __init__(self, **specs):
        """Initialize construction item.
        
        **specs:
            name: A name for the construction.
            roles: A dict which maps roles
                to either another Construction item
                or to a Text-Fabric word node.
            cases: A tuple containing condition dicts
                that were evaluated when processing this
                Construction. Key is string containing condition,
                value is Boolean.
            conds: A condition dict containing all of the
                conditions that evaluated to True to validate
                this Construction.
        """
        for k,v in specs.items():
            setattr(self, k, v)
        self.match = specs.get('match', {})
        self.name = specs.get('name', '')
        self.kind = specs.get('kind', '')
        self.pattern = specs.get('pattern', specs.get('name', ''))
        self.roles = Bunch(specs.get('roles', {}))
        self.conds = specs.get('conds', {})
        self.cases = specs.get('cases', tuple())
        self.indexslots()
        
    def __bool__(self):
        """Determine truth value of CX."""
        if self.match:
            return True
        else:
            return False
        
    def __repr__(self):
        """Display CX name with slots."""
        if self:
            return f'CX {self.name} {self.slots}'
        else:
            return '{CX EMPTY}'
            
    def __eq__(self, other):
        """Determine slot/role-based equality between CXs."""
        if (
            self.slots2role == other.slots2role
            and self.name == other.name
        ):
            return True
        else:
            return False
        
    def __hash__(self):
        return hash(
            (
                self.name, 
                 tuple(self.slots2role.items())
            )
        )
    
    def __int__(self):
        """Provide integers for first slot in cx.
        
        Most relevant for word-level CXs and for
        using TF methods on those objects.
        """
        return next(iter(sorted(self.slots)), 0)
        
    def __contains__(self, cx):
        """Determine whether certain CX is contained in this one."""
        return cx in list(self.unfoldcxs())

    def mapslots(self, rolesdict, rolename=None):
        """Recursively map all slots to top embedding role name.

        Match items contain a roles key which can contain
        any number of other match items. This function maps
        all constituent words (Text-Fabric "slots") to their
        top-level linguistic unit (linguistic role).
        """
        for role, item in rolesdict.items():
            if type(item) == Construction:
                self.mapslots(
                    item.roles.dict,
                    rolename=rolename or role
                )
            elif type(item) == int:
                self.role2slots[rolename or role].add(item)
                self.slots.add(item)   
            
    def indexslots(self):
        """Indexes slots contained in this CX."""
        self.role2slots = collections.defaultdict(set)
        self.slots = set()
        self.mapslots(self.roles.dict) # populates role2slots and slots
        self.slots = set(sorted(self.slots)) # sort slots
        self.slots2role = {
            tuple(sorted(slots)):role 
                for role, slots in self.role2slots.items()
        }  
      
    def unfoldroles(self, cx=None):
        """Return all contained construction roles as a dict.

        Recursively calls down into construction objects to convert
        to role.dict with TF slots.
        """
        cx = cx or self
        roledict={}
        roledict['__cx__'] = cx.name
        for role, item in cx.roles.dict.items():
            if type(item) == Construction:
                roledict[role] = self.unfoldroles(item)
            elif type(item) == int:
                roledict[role] = item
        return roledict
    
    def unfoldrole(self, role, cx=None):
        """Return a role that is recursively embedded.
        
        e.g.
        head
            head
                head
        """
        cx = cx or self
        for findrole, item in cx.roles.dict.items():
            if findrole == role:
                yield item
                if (
                    type(item) == Construction
                    and role in item.roles.dict
                ):
                    yield from self.unfoldrole(role, cx=item)

        
    def unfoldcxs(self, cx=None):
        """Return all contained constructions with flattened structure.
        
        Recursively calls down into construction objects and yields them.
        """
        cx = cx if cx is not None else self
        yield cx
        for role, item in cx.roles.dict.items():
            if type(item) == Construction:
                yield from self.unfoldcxs(item)  
    
    def unfoldcxpath(self, slots):
        """Return all contained constructions along a path.

        Recursively calls down into construction objects and yields them.
        """
        cx_path = []
        for cx in self.unfoldcxs():
            if set(slots).issubset(cx.slots):
                cx_path.append(cx)
        return cx_path
                
    def slots2cx(self, slottuple):
        """Return the embedded Construction to which a span of slots belong"""
        for cx in self.unfoldcxs():
            for slots, role in cx.slots2role.items():
                if slots == slottuple:
                    return cx
    
    def getslotrole(self, slot):
        """Returns the role to which a slot belongs to."""
        for role, slots in self.role2slots.items():
            if slot in slots:
                return role
                
    def updaterole(self, role, newitem):
        """Updates a role in the CX."""
        setattr(self.roles, role, newitem)
        self.roles.dict[role] = newitem
        self.indexslots() # remap slots

class CXbuilder(object):
    """Identifies and builds constructions using Text-Fabric nodes."""
    
    def __init__(self):
        """Initialize CXbuilder, giving methods for CX detection."""
        
        # cache matched constructions for backreferences
        self.cache = collections.defaultdict(
            lambda: collections.defaultdict()
        )
        
        # NB: objects below should be overwritten 
        # and configured for the particular cxs needed
        self.cxs = tuple()
        self.yieldsto = {} 
        
        # for drip-bucket categories
        self.dripbucket = tuple()
    
    def debugmess(self, msg, toggle):
        """Prints debugging messages if toggled."""
        if toggle:
            sys.stderr.write(msg+'\n')
    
    def cxcache(self, element, name, method):
        """Get cx from cache or run."""
        try:
            return self.cache[element][name]
        except KeyError:
            return method(element)
    
    def test(self, *cases):
        """Populate Construction obj based on a cases's all Truth value.
        
        The last-matching case will be used to populate
        a Construction object. This allows more complex
        cases to take precedence over simpler ones.
        
        Args:
            cases: an arbitrary number of dictionaries,
                each of which contains a string key that
                describes the test and a test that evals 
                to a Boolean.
        
        Returns:
            a populated or blank Construction object
        """
        
        # find cases where all cnds == True
        test = [
            case for case in cases
                if all(case['conds'].values())
                    and all(case['roles'].values())
        ]
        
        # return last test
        if test:
            cx = Construction(
                match=test[-1],
                cases=cases,
                **test[-1]
            )
            self.cache[cx.element][cx.name] = cx
            return cx
        else:
            return Construction(cases=cases, **cases[0])
        
    def findall(self, element):
        """Runs analysis for all constructions with an element.
        
        Returns as dict with test:result as key:value.
        """
        results = []
        
        # add cxs from this builder
        for funct in self.cxs:
            cx = funct(element)
            if cx:
                results.append(cx)
        
        # apply drip-bucket categories
        if not results:
            for funct in self.dripbucket:
                cx = funct(element)
                if cx:
                    results.append(cx)
        
        return results
                        
    def sortbyslot(self, cxlist):
        """Sort constructions by order of contained slots."""
        sort = sorted(
            ((sorted(cx.slots), cx) for cx in cxlist),
            key=lambda k: k[0]
        )
        return [cx[-1] for cx in sort]
    
    def clusterCXs(self, cxlist):
        """Cluster constructions which overlap in their slots/roles.

        Overlapping constructions form a graph wherein the constructions 
        are nodes and the overlaps are edges. This algorithm retrieves all 
        interconnected constructions. It does so with a recursive check 
        for overlapping slot sets. Merging the slot sets produces new 
        overlaps. The algorithm passes over all constructions until no 
        further overlaps are detected.

        Args:
            cxlist: list of Construction objects

        Returns:
            list of lists, where each embedded list 
            is a cluster of overlapping constructions.
        """

        clusters = []
        cxlist = [i for i in cxlist] # operate on copy

        # iterate until no more intersections found
        thiscluster = [cxlist.pop(0)]
        theseslots = set(s for s in thiscluster[0].slots)

        # loop continues as it snowballs and picks up slots
        # loop stops when a complete loop produces no other matches
        while cxlist:

            matched = False # whether loop was successful

            for cx in cxlist:
                if theseslots & cx.slots:
                    thiscluster.append(cx)
                    theseslots |= cx.slots
                    matched = True

            # cxlist shrinks; when empty, it stops loop
            cxlist = [
                cx for cx in cxlist 
                    if cx not in thiscluster
            ]

            # assemble loop
            if not matched:
                clusters.append(thiscluster)
                thiscluster = [cxlist.pop(0)]
                theseslots = set(s for s in thiscluster[0].slots)
        
        # add last cluster
        clusters.append(thiscluster)

        return clusters

    def test_yield(self, cx1, cx2):
        """Determine whether to submit a cx1 to cx2."""
        
        # get name or class yields
        cx1yields = self.yieldsto.get(
            cx1.name,
            self.yieldsto.get(cx1.kind, set())
        )
        # test yields
        if type(cx1yields) == set:
            return bool({cx2.name, cx2.kind} & cx1yields)
        elif type(cx1yields) == bool:
            return cx1yields
           
    def weaveCX(self, cxlist, cx=None, debug=False):
        """Weave together constructions on their intersections.

        Overlapping constructions form a graph wherein constructions 
        are nodes and the overlaps are edges. The graph indicates
        that the constructions function together as one single unit.
        weaveCX combines all constructions into a single one. Moving
        from right-to-left (Hebrew), the function consumes and subsumes
        subsequent constructions to previous ones. The result is a 
        single unit with embedding based on the order of consumption.
        Roles in previous constructions are thus expanded into the 
        constructions of their subsequent constituents.
        
        For instance, take the following phrase in English:
        
            >    "to the dog"
            
        Say a CXbuilder object contains basic noun patterns and can
        recognize the following contained constructions:
        
            >    cx Preposition: ('prep', to), ('obj', the),
            >    cx Definite: ('art', the), ('noun', dog)
        
        When the words of the constructions are compared, an overlap
        can be seen:
        
            >    cx Preposition:    to  the
            >    cx Definite:           the  dog
        
        The overlap in this case is "the". The overlap suggests that
        the slot filled by "the" in the Preposition construction 
        should be expanded. This can be done by remapping the role
        filled by "the" alone to the subsequent Definite construction.
        This results in embedding:
        
            >    cx Preposition: ('prep', to), 
                                 ('obj', cx Definite: ('art', the), 
                                                      ('noun', dog))
        
        weaveCX accomplishes this by calling the updaterole method native
        to Construction objects. The end result is a list of merged 
        constructions that contain embedding.
        
        Args: 
            cxlist: a list of constructions pre-sorted for word order;
                the list shrinks throughout recursive iteration until
                the job is finished
            cx: a construction object to begin/continue analysis on
            debug: an option to display debugging messages for when 
                things go wrong 🤪
                
        Prerequisites:
            self.yieldsto: A dictionary in CXbuilder that tells weaveCX
                to subsume one construction into another regardless of
                word order. Key is name of submissive construction, value
                is a set of dominating constructions. Important for, e.g., 
                cases of quantification where a head-noun might be preceded 
                by a chain of quantifiers but should still be at the top of 
                the structure since it is more semantically prominent.
                
        Returns:
            a list of composed constructions
        """
        debugmess = self.debugmess
        
        debugmess(f'\nReceived {cx} with cxlist {cxlist}', debug)

        # the search is complete, stop here
        if not cxlist:
            debugmess(f'\tSearch complete with {cx} with roles: {cx.roles.dict}', debug)
            return cx
        
        # or no search necessary, stop here
        elif cx is None and len(cxlist) == 1:
            debugmess(f'\tSearch complete with {cxlist[0]} with roles: {cxlist[0].roles.dict}', debug)
            return cxlist[0]

        # Copy constructions and operate on copies
        cx1 = cx or copy.deepcopy(cxlist.pop(0))
        cx2 = copy.deepcopy(cxlist.pop(0))
        debugmess(f'\t comparing {cx1} & {cx2}', debug)

        # replace cx1 if already contained in cx2
        if (cx1 in cx2):
            debugmess(f'\t Discarding cx1 because cx2 already contains it...', debug)
            return self.weaveCX(cxlist, cx2, debug=debug)

        # get first slot of intersection between cx1 and 2
        # get that slot's role in both cxs
        link = tuple(sorted(cx1.slots & cx2.slots))
        debugmess(f'\t link is {link}', debug)

        # retrieve lowest-contained construction with link
        cx1link = cx1.slots2cx(link)
        debugmess(f'\t\t cx1link is {cx1link}', debug)

        # submit cx1 to cx2 if cx2 is semantically dominant
        if self.test_yield(cx1link, cx2):

            debugmess(f'\t cx2 is semantically dominant over cx1...', debug)

            link1path = list(cx1.unfoldcxpath(cx1link.slots))[:-1]
            debugmess(f'\t\t searching {link1path}', debug)

            # submit until cx2's dominance ends
            while (
                link1path
                and self.test_yield(link1path[-1], cx2)
            ):
                cx1link = link1path.pop()

            debugmess(f'\t\t cx1link iterated upward to {cx1link}', debug)

            # subsume cx1 to cx2
            debugmess(f'\t\t cx2 role [{cx2.slots2role[link]}] remapping to {cx1link}...', debug)
            cx2.updaterole(cx2.slots2role[link], cx1link)
            debugmess(f'\t\t remapping done with {cx2} containing roles {cx2.roles.dict}', debug)     

            # subsume cx2 to enclosing cx
            bigcx = link1path[-1] if link1path else None 
            if bigcx:
                debugmess(f'\t assigning cx2 {cx2} to bigcx {bigcx}', debug)
                
                #return self.weaveCX([cx2]+cxlist, bigcx, debug=debug)
                
                biglink = tuple(sorted(cx1link.slots & bigcx.slots))
                bigrole = bigcx.slots2role[biglink]
                debugmess(f'\t\t big role [{bigrole}] remapping to {cx2}...', debug)
                bigcx.updaterole(bigrole, cx2)
                cx1.indexslots()
                return self.weaveCX(cxlist, cx1, debug=debug)

            # or continue with cx2
            else:
                debugmess(f'\t\t moving on with cx2 as new primary: {cx2}', debug)
                return self.weaveCX(cxlist, cx2, debug=debug)
        
        # submit cx2 to cx1
        else:
            linkcx1 = cx1.slots2cx(link)
            debugmess(f'\t submitting {cx2} to {linkcx1}...', debug)
            linkrole1 = linkcx1.slots2role[link]
            debugmess(f'\t\t role [{linkrole1}] remapping to {cx2}...', debug)
            linkcx1.updaterole(linkrole1, cx2)
            cx1.indexslots()
            debugmess(f'\t\t remapping done with {linkcx1} containing roles {linkcx1.roles.dict}', debug)
            return self.weaveCX(cxlist, cx1, debug=debug)
    
    def analyzestretch(self, stretch, debug=False):
        """Analyze an entire stretch of a linguistic unit.
        
        Applies construction tests for every constituent 
        and merges all overlapping constructions into a 
        single construction.
        
        Args:
            stretch: an iterable containing elements that
                are tested by construction tests to build
                Construction objects. e.g. stretch might be 
                a list of TF word nodes.
            debug: option to display debuggin messages
        
        Returns:
            list of merged constructions
        """
        
        # match elements to constructions based on tests
        rawcxs = [
            match for element in stretch
                for match in self.findall(element)
                    if match
        ]
        
        self.debugmess(f'rawcxs found: {rawcxs}...', debug)
        
        # return empty results
        if not rawcxs:
            self.debugmess(f'!no cx pattern matches! returning []', debug)
            return []
            
        # cluster and sort matched constructions
        clsort = [
            self.sortbyslot(cxlist)
                for cxlist in self.clusterCXs(rawcxs)    
        ]
    
        self.debugmess(f'cxs clustered into: {clsort}...', debug)
    
        self.debugmess(f'Beginning weaveCX method...', debug)
        # merge overlapping constructions
        cxs = [
            self.weaveCX(cluster, debug=debug)
                for cluster in clsort
        ]
        
        return cxs
    
class CXbuilderTF(CXbuilder):
    """Build Constructions with TF integration."""
    
    def __init__(self, tf, **kwargs):
        
        # set up TF data for tests
        self.tf = tf
        self.F, self.T, self.L = tf.api.F, tf.api.T, tf.api.L
        self.context = kwargs.get('context', 'timephrase')
        
        # set up CXbuilder
        CXbuilder.__init__(self)

    def getP(self, node):
        """Get Positions object for a TF node.
        
        Return Dummy object if not node.
        """
        if not node:
            return Dummy()
        return PositionsTF(node, self.context, self.tf).get
    
    def getWk(self, node):
        """Get Walker object for a TF word node.
        
        Return Dummy object if not node.
        """
        if not node:
            return Dummy()
        
        # format tf things to send
        thisotype = self.F.otype.v(node)
        context = self.L.u(node, self.context)[0]
        positions = self.L.d(context, thisotype)        
        return Walker(node, positions)
    
class wordConstructions(CXbuilderTF):
    """Build word constructions."""
    
    def __init__(self, tf, **kwargs):
        
        """Initialize with Constructions attribs/methods."""
        CXbuilderTF.__init__(self, tf, **kwargs)
        
        # Order matters! More specific meanings last
        self.cxs = (
            self.prep,
            self.qual_quant,
            self.card,
            self.ordn,
            self.name,
            self.cont_ptcp,
        )
        
        self.dripbucket = (
            self.pos,
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
        
        F, L = self.F, self.L
        
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
                'pattern': 'R>C/',
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
                'pattern': 'construct lexs',
                'roles': roles,
                'conds': {
                    'F.lex.v(w) in lexset':
                        F.lex.v(w) in {
                            'PNH/','TWK/', 
                            'QY/', 'QYH=/', 
                            'QYT/', '<WD/'
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