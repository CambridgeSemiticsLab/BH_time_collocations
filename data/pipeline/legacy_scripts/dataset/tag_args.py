"""
Tag clause arguments
"""

import re

def clause_objects(verb, clause_atom, clause, api):
    """Search a given clause for any marked objects."""
    
    F, E, L = api.F, api.E, api.L

    clause_phrases = L.d(clause, 'phrase')
    daughters = E.mother.t(clause_atom)
    daught_relas = set(F.rela.v(d) for d in daughters)
    daught_codes = set(F.code.v(d) for d in daughters)
    phrase_functs = set(F.function.v(p) for p in clause_phrases)
    obj_phrases = list(p for p in clause_phrases
                            if F.function.v(p) == 'Objc')
    
    # we can count direct speech clauses following amar as direct objects
    amar_obj = F.lex.v(verb) == '>MR[' and 999 in daught_codes 
    
    data = {}

    # evaluate conditions for object presence
    data['has_objc'] = any([
        bool(obj_phrases),
        'PreO' in phrase_functs,
        'PtcO' in phrase_functs,
        'Objc' in daught_relas,
        amar_obj,
    ])

    # retrieve object positionality for phrasal objects
    if data['has_objc']:
        if obj_phrases:
            pred_phrase = L.u(verb, 'phrase')[0]
            phrases = sorted([pred_phrase] + obj_phrases)

            # build order label with allowance for double objects
            order_label = ''
            for ph in phrases:
                if ph == pred_phrase:
                    order_label += 'V'
                else:
                    order_label += 'O'
            data['objc_pos'] = order_label
        else: 
            data['objc_pos'] = 'VO?'
    else:
        data['objc_pos'] = ''

    return data


def get_loca_assocs(api):
    """Identifies lexemes statistically associated with location."""
    
    F, E, L = api.F, api.E, api.L
    lexs = set()
    for cmpl in F.function.s('Loca'):
        heads = E.nhead.t(cmpl)
        for head in heads:
            loc_assoc = F.funct_assoc.v(head) or 0
            if loc_assoc >= 1.5:
                lexs.add(F.lex.v(head))
    return lexs

def clause_locas(verb, loca_lexs, api):
    """Tag location source_data within a clause"""

    F, E, L = api.F, api.E, api.L

    data = {
        'has_loca': 0,
        'loca_type': [],
        'loca_heads': [],
    }

    clause = L.u(verb, 'clause')[0]
    clause_phrases = L.d(clause, 'phrase')
    
    # check for location
    loca_ph = [p for p in clause_phrases if F.function.v(p) == 'Loca']
    
    # attempt to find a locative complement
    # Remove for now
#     if not loca_ph:
#         cmpl_phrs = [p for p in clause_phrases if F.function.v(p) == 'Cmpl']
#         for cp in cmpl_phrs:
#             for head in E.nhead.t(cp):
#                 if F.lex.v(head) in loca_lexs:
#                     loca_ph.append(cp)

    # detect presence of location
    if loca_ph:
        data['has_loca'] = bool(loca_ph)

    # tag various features on the location phrases
    for lp in loca_ph:
        data['loca_type'].append(F.typ.v(lp))
        for head in E.nhead.t(lp):
            data['loca_heads'].append(F.lex.v(head))
        
    
    data['loca_type'] = '|'.join(data['loca_type'])
    data['loca_heads'] = '|'.join(data['loca_heads'])

    return data

def clause_time(verb, api):
    """Look for an adjunctive Time element in the clause."""
    L, F = api.L, api.F
    clause_atom  = L.u(verb, 'clause_atom')[0]
    time_phs = [
        p for p in L.d(clause_atom, 'phrase')
            if F.function.v(p) == 'Time'
    ]
    data = {'has_time': 0}
    if time_phs:
        data['has_time'] = 1
    return data
    

# -- Tag Arguments in General -- 

def tag_ph_arg(ph, api, functions):
    """Tag various phrase arguments of interest."""
    F, E, L, T = api.F, api.E, api.L, api.T
    function2tag = {
        'Objc': 'O',
        'Ques': 'Q',
        'Rela': 'R',
        'Subj': 'S',
        'Cmpl': 'C',
        'Adju': 'A',
        'Time': 'T',
        'Loca': 'L',
        'Modi': 'M',
        'PreC': 'PrC',
        'PrAd': 'PrA',
        'Intj': 'I',
        'PreO': 'V-O',
        'PtcO': 'V-O',
    }
    typ2tag = {
        'IPrP': 'Q',
        'InrP': 'Q',
    }
    function = functions[ph]
    typ = F.typ.v(ph)
    
    # NB: participles require a special approach    
    # normally Pred is ignored here, but in participial
    # clauses the ptcp takes PreC function (most often)
    # So we must remap the function so that it does not
    # get matched here
    clause = L.u(ph, 'clause')[0]
    is_ptcp = (
        F.typ.v(clause) == 'Ptcp'
        and function  == 'PreC'
        and typ == 'VP'
    )
    if is_ptcp:
        function = ''    

    ph_lexs = [F.lex.v(w) for w in L.d(ph, 'word')]
    if typ in typ2tag:
        return typ2tag[typ]
    elif function in function2tag:
        return function2tag[function]
    elif function == 'Conj' and 'W' not in ph_lexs:
        return 'C'

def tag_word_arg(word, api, functions):
    """Tag various word arguments of interest."""
    F = api.F
    pdp_set = {'conj'}
    tag = ''
    if F.pdp.v(word) in pdp_set:
        tag += f'{F.lex.v(word)}'
    elif F.sp.v(word) == 'verb':
        tag += 'V'
    return tag 

def tag_cl_arg(cl, api, functions):
    """Tag various daughter clause arguments of interest."""
    F = api.F
    rela2tag = {
        'Objc': 'O',
        'Subj': 'S',
        'Adju': 'A',
        'Cmpl': 'C',
        'Spec': 'Sp',
        'PrAd': 'PrA',
        'PreC': 'PrC',
        'Attr': 'Atr',
    }
    rela = F.rela.v(cl)
    return rela2tag.get(rela, '')

def get_args(node, slot_getter, arg_getter, covered_slots, api, functions):
    """Get args by running arg_getter function and add to a list."""

    L = api.L
    arg_tag = arg_getter(node, api, functions)
    slots = slot_getter(node)

    # skip non-tagged elements or covered elements
    if not arg_tag or set(slots) & covered_slots:
        return []
    # return listed element with tag and slots
    else:
        for slot in slots:
            covered_slots.add(slot)
        return [(slots, arg_tag)] 

def clause_args(verb, api, functions):
    """Tag key clause arguments."""
    
    F, E, L = api.F, api.E, api.L

    clause = L.u(verb, 'clause')[0]
    clause_atom = L.u(verb, 'clause_atom')[0]
    cl_phrases = L.d(clause, 'phrase')
    cl_words = L.d(clause, 'word')
    cla_daughts = E.mother.t(clause_atom)
    cl_daughts = E.mother.t(clause)
    
    # add any of three different types of objects;
    # priority is determined by order, so e.g. a phrase takes
    # priority over a word in this case
    args = []
    covered_slots = set()
    slots_node = lambda node: L.d(node, 'word')
    slots_word = lambda word: (word,)
    for p in cl_phrases:
        args.extend(get_args(p, slots_node, tag_ph_arg, covered_slots, api, functions))
    for w in cl_words:
        args.extend(get_args(w, slots_word, tag_word_arg, covered_slots, api, functions))
    for d in cl_daughts:
        args.extend(get_args(d, slots_node, tag_cl_arg, covered_slots, api, functions))

    # iterate through the arguments, sorted by slots, and add the tags
    arg_strs = []
    for slots_arg in sorted(args):
        arg_strs.append(slots_arg[-1]) # last item is the argument tag 

    arg_str = '-'.join(arg_strs)

    return arg_str

def tag_position(clarg):
    """Tag time positions in the clause."""
    tags = {
        'medial': '.*V.*T.*[OCS]',
        'preposed': '.*T.*V',
    }
    if re.match(tags['medial'], clarg):
        return 'medial'
    elif re.match(tags['preposed'], clarg):
        return 'preposed'
    else:
        return 'final'
