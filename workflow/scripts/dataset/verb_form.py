from .positions import PositionsTF

def get_verbform(node, api):
    """Remap BHSA verb tense values to custom values
    
    Args:
        node: int, representing BHSA node id
        api: an instance of Text-Fabric with BHSA loaded
    """
    
    tense_map = {
            'impf': 'yqtl',
            'perf': 'qtl',
            'ptca': 'ptcp',
    }
    
    bhsa_tense = api.F.vt.v(node)
    verb_form = tense_map.get(bhsa_tense, bhsa_tense)
    
    # adjust weqatal
    P = PositionsTF(node, 'clause', api)
    if verb_form == 'qtl' and P.get(-1, 'lex') == 'W':
        verb_form = 'wqtl'
    
    return verb_form

def get_cl_verbform(clause, api):
    """Tag a verbform for a supplied clause."""

    F, L = api.F, api.L

    # retrieve a list of verb candidates
    verb_cands = [word for word in L.d(clause, 'word')
                     if F.pdp.v(word) == 'verb']
    
    # filter out non-predicative participles
    if len(verb_cands) > 1:
        verb_cands = [word for word in verb_cands
                         if F.vt.v(word) not in {'ptcp', 'ptca'}]
        
    # return appropriate tag
    if verb_cands:
        return get_verbform(verb_cands[0], api)
    else:
        return 'Ø'
