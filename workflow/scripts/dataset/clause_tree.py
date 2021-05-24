"""
Functions for retrieving clause relations up and down 
the BHSA clause relation tree.
"""

def get_predecessors(clause_atom, api, start=True):
    """Retrieve all of a clause_atom's predecessors in BHSA tree"""
    E = api.E
    if not start:
        yield clause_atom
    mother = E.mother.f(clause_atom)
    if mother:
        yield from get_predecessors(mother[0], api, start=False)
        
def get_successors(clause_atom, api, start=True):
    """Retrieve all of a clause_atom's successors in BHSA tree"""
    E = api.E
    if not start:
        yield clause_atom
    children = E.mother.t(clause_atom)
    for child in children:
        yield from get_successors(child, api, start=False)
