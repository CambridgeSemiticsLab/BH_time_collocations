"""
Retrieve data from Construction object graphs
"""

import networkx as nx
from itertools import cycle
from copy import deepcopy

def get_headword(cx, default=None):
    """Get a word that serves as head"""
    try:
        return list(get_succroles(cx,'head'))[-1]
    except (IndexError, AttributeError):
        return default

def get_modifiers(head, cx, name, default=None):
    """Retrieve a modifier on a particular head"""
    success = False
    for c in cx.graph:
        if (get_attribute(c,'name') == name) and (head in c): 
            yield c
            success = True
    # unsuccessful search
    if not success:
        yield default

def get_attribute(node, attr, default=None):
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

def has_node(cx, attribute, value):
    """Determine whether a CX object contains node with value."""
    return value in set(get_attribute(c,attribute) for c in cx)

def unfold_roles(cx):
    """Return all contained construction roles as a dict.

    Recursively calls down into graph nodes to populate
    a recursive dict along with labels.
    """
    roledict = {}
    roledict['__cx__'] = cx.name
    for child in graph.succ[cx]:
        role = graph[cx][child]['role']
        if type(child) == Construction:
            roledict[role] = unfoldroles(child)
        elif type(child) == int:
            roledict[role] = child
    return roledict

def get_role(cx, role, default=None):
    """Retrieves a succesor node of a specific role.
    
    If node is not present, return default.
    """
    for node in cx.graph.succ[cx]:
        if cx.graph[cx][node]['role'] == role:
            return node
    return default

def get_succroles(cx, role, graph=None):
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
            yield from get_succroles(node, role, graph)

def get_predecessor(cx, graph, default=None):
    """Select first predecessor in graph.

    In Construction graphs, any node only 
    has one predecessor. This function returns
    that item.
    """
    return next(iter(graph.pred[cx]), default)

def find_paths(cx, paths, graph=None):
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

    ## NB THIS IS AN INCOMPLETE ALGORITHM!
    # after some fiddling with this, I've realized
    # for now it is not worth the time. If I need it
    # later I can come back to it

    # create rotating carousel of validator functions
    paths = list(cycle(path) for path in paths)
    graph = cx.graph
    walks = [
        walk for node in graph 
            for walk in list(nx.all_simple_paths(graph, cx, node))
    ]

    for path in paths:
        for walk in walks:
            tests = deepcopy(path)
            for node in walk:
                pred = get_predecessor(node, graph)
                edge = graph[pred][node] if pred else {'role':'head'}
                # continue walk or break it off
                if next(tests)(node, edge, graph):
                    yield node
                else:
                    break


