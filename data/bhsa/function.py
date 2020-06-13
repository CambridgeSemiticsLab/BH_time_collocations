"""
Re-map BHSA phrase functions based on
automatic and manual selections.
These are added to correct the data on
time phrases.
"""

from positions import PositionsTF
from tf.app import use

def modify_function(functions, api):
    """Apply automatic and manual function edits to BHSA

    Args:
        functions: a dict of node:function mappings to edit
        api: an instance of BHSA in Text-Fabric
    """

    F, L = api.F, api.L
    A = use('bhsa', api=api, silent='deep')

    # apply automatic mods
    for phrase in F.otype.s('phrase'):

        cases = [] 
        clause = L.u(phrase,'clause')[0]
        clause_speechs = set(F.pdp.v(w) for w in L.d(clause,'word'))
        phrase_funct = F.function.v(phrase)

        for word in L.d(phrase,'word'):
            P = PositionsTF(word, 'phrase', A).get
    
            # -- single-particle with verb required -- 
            single_lexset = {
                'TMJD/', '>Z', 'VRM/', '>ZJ'

            # put a few terms on ice for now
            # these have not seen good success -CK, 2020.06.13
            # '<WD', 'JWMM', 'CNJ', 'RG<', 'PT>M'
            }
            good_functs = {'Modi',} 
            cases.append({
                'funct': 'Time',
                'conds': (
                    not P(1) and not P(-1),
                    'verb' in clause_speechs, 
                    F.lex.v(word) in single_lexset,
                    phrase_funct in good_functs,
                )
            })

            # -- multiple word phrases with verb required -- 
            multi_lexset_modi = {
                # put these terms on ice for now
                # frequency terms are not yet in 
                # view in the project, durational terms
                # e.g. TMJD (above) ARE needed, but 
                # the terms below are frequentive, -CK 2020.06.13
                # 'P<M/', 'MHR[', 'MHRH/', 
                # 'PT<', 'RGL/', 
            }
            cases.append({
                'funct': 'Freq',
                'conds': (
                    'verb' in clause_speechs,
                    F.lex.v(word) in multi_lexset_modi,
                    phrase_funct == 'Modi',
                )
            })


            # see note above about frequentives -CK 2020.06.13
            # -- Special Cases -- 
#            cases.append(
#                { # plural/dual cardinals 
#                    'funct':'Time',
#                    'conds': (
#                        not P(1) and not P(-1),
#                        'verb' in clause_speechs,
#                        phrase_funct == 'Modi',
#                        F.ls.v(word) == 'card',
#                        F.nu.v(word) in {'pl','du'},
#                    )
#                },
#            )
                         

        # apply the mods
        # NB, in cases of overlap the last case prevails
        for case in cases:
            if all(case['conds']):
                functions[phrase] = case['funct']            

    # manually remap phrase functions
    functions.update({
        849296:'Loca', 
        825329:'Loca',
        828081:'Cmpl',
        774349:'Adju',
        774352:'Adju',
        775948:'Adju',
        775985:'Adju',
        876172:'Adju',
        881665:'Objc', # phrase belongs with previous as adjectival element
        870273:'Time', # prep and conj belong with time phrase
        870274:'Time', # modifier "KBR" belongs with time phrase
        731915:'Modi', # no reviewed sources took this as temporal

        # similar to the frequentive note above,
        # I no longer want to include incremental
        # phrases in the sample -CK 2020.06.13
        #677350:'Time', # use indicates incremental activity
        #677350:'Time', # ^
        #706975:'Time', # ^
        #706975:'Time', # ^
    })  
