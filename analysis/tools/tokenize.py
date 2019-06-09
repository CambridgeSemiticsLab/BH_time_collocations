'''
This module contains functions for tokenizing
the surface forms of time adverbials.

This module assumes that an instance of BHSA
has been loaded with F class instantiated and the
features lex and g_cons_utf8 loaded.
'''

from __main__ import F, L

def tokenize_surface(node):
    '''
    Return a surface token string of a BHSA node.
    The words are dot-separated and heh consonants
    are added if they are present in vocalized form. 
    '''
    subtokens = []
    for w in L.d(node, 'word'):
        if F.lex.v(w) == 'H':
            subtokens.append('ה')
        else:
            subtokens.append(F.g_cons_utf8.v(w))
    return '.'.join(subtokens)