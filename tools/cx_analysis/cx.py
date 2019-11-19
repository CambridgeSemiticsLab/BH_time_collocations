import copy
import uuid
import networkx as nx

class Construction(object):
    """A linguistic construction and its attributes.
    
    This is version 2, which utilizes NetworkX graphs
    instead of standard dictionaries.
    """
    
    def __init__(self, **specs):
        """Make a new construction.
        
        **specs:
            name: A name for the construction (CX).
            kind: A kind for the CX.
            pattern: Name of pattern that matched to license
                this CX.
            conds: A dictionary of conditions that all eval to
                True to license this CX. Keys are strings that
                describe what was tested; values are booleans.
            cases: A tuple containing all of the possible conds
                dicts that were tested and their results, including 
                non-matches. Useful for debugging.
                
        Key Attributes:
            slots: An ordered tuple of TF slot integers which
                describe what span of words in the corpus this
                CX represents.
            graph: A NetworkX graph object that contains the
                internal structure of this cx. Top node of
                the graph is this object; edges have values of
                "role" that give semantic role of each node.
            parent: a parent CX if this one is contained in 
                another's graph.
        """
        
        # map optional attributes
        for k,v in specs.items():
            setattr(self, k, v)
            
        # map obligatory attributes
        self.element = specs.get('element', str(uuid.uuid4()))
        self.match = specs.get('match', {})
        self.name = specs.get('name', '')
        self.kind = specs.get('kind', '')
        self.pattern = specs.get('pattern', specs.get('name', ''))
        self.conds = specs.get('conds', {})
        self.cases = specs.get('cases', tuple())
        
        # map roles and slots
        self.graph = nx.DiGraph()
        self.populate_graph(specs.get('roles', {}))
        self.slots = tuple()
        self.updateslots() # populates self.slots
    
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
        
    def _cx_att(self, attr, item):
        """Get an attribute on a cx or return int"""
        if type(item) == Construction:
            return item.__dict__[attr]
        elif type(item) == int:
            return item
            
    def _rolestuple(self):
        return tuple(
            (n1, n2, self.graph[n1][n2]['role'])
                 for n1, n2 in nx.bfs_edges(self.graph, self)
        )
            
    def __eq__(self, other):
        """Determine slot/role-based equality between CXs."""
        if (
            self.__class__ == other.__class__
            and self.name == other.name
            and str(self._rolestuple) == str(other._rolestuple)
        ):
            return True
        else:
            return False
        
    def __hash__(self):
        return hash(
            (self.name, self.element)
        )
    
    def __int__(self):
        """Provide integers for first slot in cx.
        
        Most relevant for word-level CXs and for
        using TF methods on those objects.
        """
        return next(iter(sorted(self.slots)), 0)
        
    def __contains__(self, cx):
        """Determine whether certain CX is contained in this one."""
        return cx in self.subgraph()
        
    def __deepcopy__(self, memo):
        """Return a copied version of this CX"""
        roles = {
            self.graph[self][node]['role']:node 
                for node in self.graph.succ[self]
        }
        attribs = {
            k:v for k,v in self.__dict__.items()
                if k != 'graph'
        }
        attribs['roles'] = roles
        return Construction(**attribs)
        
    def getslots(self, item):
        """Get TF integer slots as tuple."""
        slots = self._cx_att('slots', item)
        if type(slots) == tuple:
            return slots
        else:
            return (slots,)
            
    def populate_graph(self, rolesdict):
        """Populate the graph with the CX's structure"""
        
        # populate graph with roles
        self.graph.add_node(self)
        for role, child in rolesdict.items():
            
            # create unique copy of child 
            # esp. relevant for CX objects
            # that are shared between other CXs
            child = copy.deepcopy(child)
            
            # add child to graph
            self.graph.add_edge(self, child, role=role)
            
            # import child's graph structure
            if type(child) == Construction:
                self.graph.update(child.graph)
                child.graph = self.graph # assign graph to child
    
    def subgraph(self):
        """Return graph governed by this CX"""        
        # return subgraph
        return self.graph.subgraph(nx.bfs_tree(self.graph, self))
    
    def updategraph(self, oldnode, newnode):
        """Update the internal structure of CX graph.
        
        Change oldnode to newnode.
        """
        
        # get predecessor for reassignment
        pred = next(iter(self.graph.pred[oldnode]))
        
        # get replacement role 
        role = self.graph[pred][oldnode]['role']

        # remove old node
        self.graph.remove_node(oldnode)

        # make unique copy of newnode
        newnode = copy.deepcopy(newnode)

        # add new node
        self.graph.add_edge(pred, newnode, role=role)

        # add new nodes's constituents & roles to graph
        if type(newnode) == Construction:
            self.graph.update(newnode.graph)
            newnode.graph = self.graph # assign graph to child
            
        # remap slots to reflect new nodes
        self.updateslots()
        
        # remap slots for constituent cxs
        for node in self.graph:
            if type(node) == Construction:
                node.updateslots()
        
    def updateslots(self):
        """Update slots covered by this CX."""
        self.slots = tuple(sorted(set(
            slot for node in nx.bfs_tree(self.graph, self)
                for slot in self.getslots(node)
        )))
        
    def getrole(self, role, default=None):
        """Retrieves a succesor node of a specific role.
        
        If node is not present, return default.
        """
        for node in self.graph.succ[self]:
            if self.graph[self][node]['role'] == role:
                return node
        return default
    
    def getsuccroles(self, role, start=None):
        """Retrieve successive roles.
        
        Recursively calls down the graph looking
        for successive roles.
        E.g. 
        >    head -> head -> head
        but not
        >    head -> adjv -> head
        """
        start = start or self
        for adj_node in self.graph.adj[start]:
            if self.graph[start][adj_node]['role'] == role:
                yield adj_node
                yield from self.getsuccroles(role, start=adj_node)
                
    def unfoldroles(self, cx=None):
        """Return all contained construction roles as a dict.

        Recursively calls down into graph nodes to populate
        a recursive dict along with labels.
        """
        cx = cx if cx is not None else self
        roledict = {}
        roledict['__cx__'] = cx.name
        for child in self.graph.succ[cx]:
            role = self.graph[cx][child]['role']
            if type(child) == Construction:
                roledict[role] = self.unfoldroles(child)
            elif type(child) == int:
                roledict[role] = child
        return roledict