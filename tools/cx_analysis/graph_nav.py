"""
Retrieve data from Construction object graphs
"""

import networkx as nx
from itertools import cycle
from copy import deepcopy

class Navigator:

    def __init__(self):
        pass
    
    def get_attribute(self, node, attr, default=None):
        """Return a requested attribute from a graph node.

        Construction (CX) graphs contain other CX objects but also
        integers (Text-Fabric node numbers). This method
        enables attribute calls on CX objects without
        erroring out on integer objects.
            
        Args:
            node: a node in a CX graph
            attr: an attribute string to call on node
            default: a default to return if attr not found
        """      
        if type(node) == int:
            return default
        try:
            return node.__dict__[attr]
        except KeyError:
            return default

    def unfold_roles(self, cx):
        """Return all contained construction roles as a dict.

        Recursively calls down into graph nodes to populate
        a recursive dict along with labels.
        """
        roledict = {}
        roledict['__cx__'] = cx.name
        for child in self.graph.succ[cx]:
            role = self.graph[cx][child]['role']
            if type(child) == Construction:
                roledict[role] = self.unfoldroles(child)
            elif type(child) == int:
                roledict[role] = child
        return roledict

    def get_role(self, cx, role, default=None):
        """Retrieves a succesor node of a specific role.
        
        If node is not present, return default.
        """
        for node in cx.graph.succ[cx]:
            if cx.graph[cx][node]['role'] == role:
                return node
        return default

    def get_succroles(self, cx, role, graph=None):
        """Retrieve successive roles.
        
        Recursively calls down the graph looking
        for successive roles.
        E.g. 
        >    head -> head -> head
        but not
        >    head -> adjv -> head
        """ 
        # assign top-level graph
        graph = graph or cx.graph
        for node in graph.adj[cx]:
            if graph[cx][node]['role'] == role:
                yield node
                yield from self.get_succroles(node, role, graph)

    def find_paths(self, cx, paths, graph=None, edge={'role':'head'}):
        """Find nodes which match a validator function

        Nodes sit along paths of edges. Both nodes and 
        edges possess attributes. Sometimes one needs to 
        retrieve a specific path wherein each node and 
        edge matches a criteria. The nodes/edges might
        also need to match criteria in a specific order.
        This method selects such paths. 

        Args:
            cx: Construction object with graph
            paths: an iterable of functions which return
                True/False on whether to yield node. Function
                requires node and edge dict as args.
            graph: A NetworkX graph object. Start at None by
                default, reset by algorithm to top-level graph.
        """

        # begin new root level
        # create rotating carousel of validator functions
        if graph is None:
            graph = cx.graph
            paths = list(
                (i, cycle(path)) for i, path in enumerate(paths)
            )

        # search for matching paths
        for path in paths:
            n_path, vfunct = path[0], next(path[1])
            if vfunct(cx, edge, graph):
                yield cx 
                for node, edge in graph.adj[cx].items():
                    yield from self.find_paths(node, [deepcopy(path)], graph, edge)

#        for node, edge in graph.adj[cx].items():
#            print(f'val id: {id(paths)}')
#            if debug: print(f'testing node: {node} edge: {edge}')
#            if vfunct(node, edge, graph):
#                if debug: print(f'\tn_path {n_path}: True')
#                yield from self.find_paths(node, copy(paths), graph, debug=debug)
#            elif debug:
#                print(f'\tn_path {n_path}: False')
