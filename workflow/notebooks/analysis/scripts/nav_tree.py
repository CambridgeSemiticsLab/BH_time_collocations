"""
Functions for navigating parse-trees as created in this project.

A parse tree consists of recursively embedded lists, where each
list represents a separate phrase that contains 3 items:

    source: a node or list representing the entity which
        is modifying the target entity with some relationship
    target: a node or list representing the entity which
        is being modified
    relation: the name of the relationship between source and
        target

All together these lists form a graph network of relationships.
"""

def get_slots(phrase):
    """Recursively retrieve slots from a phrase tree."""
    if type(phrase) == int:
        yield phrase
        return
    src, tgt, rela = phrase
    if type(src) == int:
        yield src
    else:
        yield from get_slots(src)
    if type(tgt) == int:
        yield tgt
    else:
        yield from get_slots(tgt)

def get_head(phrase):
    """Retrieve a phrase head."""
    src, tgt, rela = phrase
    if type(tgt) == int:
        return tgt
    else:
        return get_head(tgt)
    
def get_head_path(phrase):
    """Yield phrases all the way down to the right-most item (the head).
    
    This method ignores modifications of head modifiers.
    For instance, if an adverb modifies an adjective that 
    modifies a head, that adverb will not fall inside the 
    path of the head item and will therefore not be yielded."""
    yield phrase
    src, tgt, rela = phrase
    if type(tgt) != int:
        yield from get_head_path(tgt)
        
def traverse_tree(phrase):
    """Traversing down a phrase tree."""
    src, tgt, rela = phrase
    branches = []
    # NB: storing branches before yielding them allows
    # the tree to be changed during iteration
    # without causing an infinite loop;
    # this is because we only access the whole tree once
    # and return the results all at once
    if type(src) == list:
        branches.append(traverse_tree(src))
    if type(tgt) == list:
        branches.append(traverse_tree(tgt))
    yield phrase
    for branch in branches:
        yield from branch

def unfold_paras(phrase):
    """Take a parallel phrase and yield its standalone phrases.
    
    This is done recursively by picking which phrases to descend
    down to and yield. For instance, note that a rela==CONJ will 
    result in the src item (the conjunction) being bypassed whereas
    the rela==PARA results in descending to both src and tgt.
    
    NB: if used on a non-parallel phrase, 
    this will simply yield that phrase.
    """
    src, tgt, rela = phrase
    if rela == 'PARA':
        # NB: The order of execution matters below;
        # since the parser has placed conjunctive items first
        # in the parse list (src), we must first process
        # tgt to yield the items in their natural order
        if type(tgt) != int:
            yield from unfold_paras(tgt)
        else:
            yield [tgt]
        if type(src) != int:
            yield from unfold_paras(src)
        else:
            yield [src]
            
    elif rela == 'CONJ':
        if type(tgt) != int:
            yield from unfold_paras(tgt)
        else:
            yield [tgt]
    else:
        yield phrase

def show_relas(parse_tree, stringer, joiner='\n'):
    """Visualize relationship in a tree."""
    subphrases = list(traverse_tree(parse_tree))
    show_strings = []
    for src, tgt, rela in subphrases:
        src_slots = sorted(get_slots(src))
        tgt_slots = sorted(get_slots(tgt))
        src_str = stringer(src_slots)
        tgt_str = stringer(tgt_slots)
        show_strings.append(f'{src_str}  --{rela}-->  {tgt_str}')
    return joiner.join(show_strings)
