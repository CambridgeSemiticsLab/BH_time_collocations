import tools.nav_tree as nt

def formal_tokens(nodes, tf, feature):
    """Make 'dumb' tokens"""
    Fs, L = tf.Fs, tf.L

    lexs = [
        Fs(feature).v(s) for node in nodes
            for s in L.d(node, 'word')
    ]
    tokens = '.'.join(lexs)

    return tokens

def tokenize_lexemes(parses, tf, 
        feature='lex_utf8', sep='.', 
        cardtok='מ׳', ordntok='ס׳', headtok='זמן' ,
        heads=[]):
    '''
    Return a surface token string of a BHSA node.
    The words are dot-separated and heh consonants
    are added if they are present in vocalized form. 
    '''
    Fs, L = tf.Fs, tf.L
    lextokens = []
    for ph in parses:
        
        if len(ph) == 1:
            head_path = [[0, ph[0], None]]
        else:
            head_path = list(nt.get_head_path(ph))

        for sp in head_path:
            src, tgt, rela = sp
            srcslot = nt.get_head(src)
            if Fs('lex').v(srcslot) == 'H':
                lextokens.append(
                    (srcslot, Fs(feature).v(6)) # first ה in HB; ensure H always shows
                )
            elif Fs('ls').v(srcslot) == 'card':
                if cardtok not in lextokens:
                    lextokens.append(
                        (srcslot, cardtok)
                    )
            elif Fs('ls').v(srcslot) == 'ordn':
                if ordntok not in lextokens:
                    lextokens.append(
                        (srcslot, ordntok)
                    )
            elif rela == 'PP':
                if type(src) == int:
                    lextokens.append(
                        (srcslot, Fs(feature).v(srcslot))       
                    )
                else:
                    slots = sorted(nt.get_slots(src))
                    pptext = ''.join(
                        Fs(feature).v(s) for s in slots
                    )
                    lextokens.append(
                        (scslot, pptext)
                    )

            elif Fs('pdp').v(srcslot) == 'prde':
                lextokens.append(
                    (srcslot, Fs(feature).v(srcslot))       
                )
        # executes only at end of head path
        if tgt in heads: 
            if Fs('nu').v(tgt) == 'du':
                lextokens.append(
                    (tgt, cardtok)
                )    
            if lextokens:
                lextokens.append((tgt, headtok))
            elif Fs('pdp').v(tgt) != 'advb':
                lextokens.append((tgt, headtok))
            else:
                lextokens.append(
                    (tgt, Fs(feature).v(tgt))
                )

    lextokens.sort()
    lextokens = [t[1] for t in lextokens]
    lex_token = f'{sep}'.join(lextokens)

    # add some sequence normalizations
    lex_token = lex_token.replace(
        'זמן.מ׳',
        'מ׳.זמן'
    )
    return lex_token
