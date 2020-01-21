"""
Retrieve data from Construction object graphs
"""

import networkx as nx

class Getter:

    def __init__(self):
        pass

    def getsuccroles(self, cx, role):
        """Retrieve successive roles.
        
        Recursively calls down the graph looking
        for successive roles.
        E.g. 
        >    head -> head -> head
        but not
        >    head -> adjv -> head
        """
        for node in cx.graph.adj[cx]:
            if cx.graph[cx][node]['role'] == role:
                yield node
                yield from self.getsuccroles(node, role)

    def unfoldroles(self, cx):
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
