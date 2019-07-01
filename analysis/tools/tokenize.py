'''
This module contains functions for tokenizing
the surface forms of time adverbials.
'''

def tokenize_surface(node, tf):
    '''
    Return a surface token string of a BHSA node.
    The words are dot-separated and heh consonants
    are added if they are present in vocalized form. 
    '''
    F, L = tf.F, tf.L
    subtokens = []
    for w in L.d(node, 'word'):
        if F.lex.v(w) == 'H':
            subtokens.append('ה')
        else:
            subtokens.append(F.g_cons_utf8.v(w))
    return '.'.join(subtokens)
