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
    yield phrase
    src, tgt, rela = phrase
    src_slots = sorted(get_slots(src))
    tgt_slots = sorted(get_slots(tgt))
    head = get_head(phrase)
    if type(src) == list:
        yield from traverse_tree(src)
    if type(tgt) == list:
        yield from traverse_tree(tgt)
        
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
        yield from unfold_paras(tgt)
        yield from unfold_paras(src)
    elif rela == 'CONJ':
        yield from unfold_paras(tgt)
    else:
        yield phrase