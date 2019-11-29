'''
This module contains functions for tokenizing
the surface forms of time adverbials.
'''

def tokenize_surface(nodes, tf):
    '''
    Return a surface token string of a BHSA node.
    The words are dot-separated and heh consonants
    are added if they are present in vocalized form. 
    '''
    F, L = tf.F, tf.L
    subtokens = []
    slots = L.d(nodes,'word') if type(nodes) == int else nodes
    for s in slots:
        if F.lex.v(s) == 'H':
            subtokens.append('ה')
        else:
            subtokens.append(F.g_cons_utf8.v(s))
    return '.'.join(subtokens)
