'''
This module contains functions for tokenizing
the surface forms of time adverbials.
'''

def tokenize_surface(nodes, tf, feature='g_cons_utf8', sep='.'):
    '''
    Return a surface token string of a BHSA node.
    The words are dot-separated and heh consonants
    are added if they are present in vocalized form. 
    '''
    Fs, L = tf.Fs, tf.L
    subtokens = []
    slots = L.d(nodes,'word') if type(nodes) == int else nodes
    for s in slots:
        lex = L.u(s, 'lex')[0]
        subtokens.append(Fs(feature).v(lex))       
    return f'{sep}'.join(subtokens)
