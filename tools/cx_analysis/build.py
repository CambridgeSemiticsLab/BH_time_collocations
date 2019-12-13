"""
Classes used to identify and build Construction objects
"""

import collections
import copy
import networkx as nx
from positions import Dummy, Positions, PositionsTF, Walker
from debugging import Debugger
from .cx import Construction

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
    
    def cxcache(self, element, name, method):
        """Get cx from cache or run."""
        try:
            return self.cache[element][name]
        except KeyError:
            return method(element)
    
    def test_result(self, test, *cases):
        """Return the result of a test as a new Construction object"""
        
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
        def get_roles(case):
            # get roles dict
            roles = case.get('roles', {'':True})
            return roles

        test = [
            case for case in cases
                if all(case['conds'].values())
                    and all(get_roles(case).values())
        ]
        return self.test_result(test, *cases) 
        
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
                if theseslots & set(cx.slots):
                    thiscluster.append(cx)
                    theseslots |= set(cx.slots)
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

    def yields(self, cx1, cx2):
        """Determine whether to submit cx1 to cx2."""
        
        # determine which yields dict to use
        # yielding can be configured generally
        # or specific to a pattern and its rules
        yieldsto = cx1.__dict__.get('yieldsto', self.yieldsto)
        
        # get name or class yields
        cx1yields = yieldsto.get(
            cx1.name,
            yieldsto.get(cx1.kind, set())
        )
        # test yields
        if type(cx1yields) == set:
            return bool({cx2.name, cx2.kind} & cx1yields)
        elif type(cx1yields) == bool:
            return cx1yields
        
    def interslots(self, cx1, cx2):
        """Get the intersecting slots of two CXs
        
        Return as sorted tuple.
        """
        return tuple(sorted(
            set(cx1.slots) & set(cx2.slots)
        ))
    
    def slots2node(self, cx, slots):
        """Get a CX node from a tuple of slots."""
        for node in nx.bfs_tree(cx.graph, cx):
            if cx.getslots(node) == slots:
                return node
    
    def intersect_node(self, cx1, cx2):
        """Get node from cx1 with slots common with cx2."""
        intersect = self.interslots(cx1, cx2)
        return self.slots2node(cx1, intersect)

    def weaveCX(self, cxlist, debug=False):
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
        
        db = Debugger(debug)
        
        db.say(f'\nReceived cxlist {cxlist}', 0)

        # compile all cxs to here
        root = copy.deepcopy(cxlist.pop(0))
        
        db.say(f'Beginning analysis with {root}')
        
        # begin matching and remapping
        while cxlist:
            
            # get next cx
            ncx = copy.deepcopy(cxlist.pop(0))
            
            # find root node with slots intersecting next cx
            db.say(f'comparing {root} with {ncx}', 1)
            node = self.intersect_node(root, ncx)
            db.say(f'intersect is at {node}')
            
            # remove cxs covered by larger version
            if root in ncx:
                db.say(f'root {root} in ncx {ncx}...replacing root with ncx')
                root = ncx
            
            # update yielded nodes
            elif self.yields(node, ncx):
                
                db.say(f'{node} being yielded to {ncx}')
                   
                # get top-most yielding node
                path = nx.shortest_path(root.graph, root, node)
                while path and self.yields(path[-1], ncx):
                    node = path.pop(-1)
                
                db.say(f'top-yielding node is {node}', 2)
                   
                # update ncx graph
                db.say(f'comparing {ncx} with {node}')
                ncxnode = self.intersect_node(ncx, node)
                db.say(f'intersect is at {ncxnode}')
                ncx.updategraph(ncxnode, node)
                db.say(f'ncx updated to {ncx}')
                
                # update root graph or remap root to ncx
                if root != node:
                    rnode = self.intersect_node(root, ncx)
                    db.say(f'replacing node {rnode} in root {root} with {ncx}')
                    root.updategraph(rnode, ncx)
                    
                else:
                    # switch root and ncx
                    db.say(f'switching {root} with {ncx}')
                    root = ncx
                 
            # update all non-yielding nodes
            else:
                db.say(f'\tupdating {node} in root with {ncx}')
                root.updategraph(node, ncx)
            
        return root
            
    def analyzestretch(self, stretch, duplicate=False, debug=False):
        """Analyze an entire stretch of a linguistic unit.
        
        Applies construction tests for every constituent 
        and merges all overlapping constructions into a 
        single construction.
        
        Args:
            stretch: an iterable containing elements that
                are tested by construction tests to build
                Construction objects. e.g. stretch might be 
                a list of TF word nodes.
            duplicate: whether to keep a copy of an analyzed
                cx
            debug: option to display debuggin messages
        
        Returns:
            list of merged constructions
        """
                   
        db = Debugger(debug)
        
        # match elements to constructions based on tests
        rawcxs = []
        covered = set()
        for element in stretch:
            matches = self.findall(element)
            if matches:
                rawcxs.extend(matches)
                covered |= set(
                    el for cx in matches 
                        for el in cx.graph
                )
                # keep copy of the cx
                if duplicate:
                    rawcxs.append(element)
                    covered.add(element)
        
        # apply drip-bucket categories
        for element in set(stretch) - covered:
            for funct in self.dripbucket:
                dripcx = funct(element)
                if dripcx:
                    rawcxs.append(dripcx)
        
        db.say(f'rawcxs found: {rawcxs}...')
        
        # return empty results
        if not rawcxs:
            db.say(f'!no cx pattern matches! returning []')
            return []
            
        # cluster and sort matched constructions
        clsort = [
            self.sortbyslot(cxlist)
                for cxlist in self.clusterCXs(rawcxs)    
        ]
    
        db.say(f'cxs clustered into: {clsort}...')
    
        db.say(f'Beginning weaveCX method...')
        # merge overlapping constructions
        cxs = [
            self.weaveCX(cluster, debug=debug)
                for cluster in clsort
        ]
        
        return self.sortbyslot(cxs)
    
class CXbuilderTF(CXbuilder):
    """Build Constructions with TF integration."""
    
    def __init__(self, tf, **kwargs):
        
        # set up TF data for tests
        self.tf = tf
        self.F, self.E, self.T, self.L = tf.api.F, tf.api.E, tf.api.T, tf.api.L
        self.context = kwargs.get('context', 'timephrase')
        
        # set up CXbuilder
        CXbuilder.__init__(self)

    def getP(self, node, context=None):
        """Get Positions object for a TF node.
        
        Return Dummy object if not node.
        """
        context = context or self.context
        if not node:
            return Dummy
        return PositionsTF(node, context, self.tf).get
    
    def getWk(self, node, context=None):
        """Get Walker object for a TF word node.
        
        Return Dummy object if not node.
        """
        if not node:
            return Dummy()
        
        # format tf things to send
        thisotype = self.F.otype.v(node)
        get_context = context or self.context
        context = self.L.u(node, get_context)[0]
        positions = self.L.d(context, thisotype)        
        return Walker(node, positions)
