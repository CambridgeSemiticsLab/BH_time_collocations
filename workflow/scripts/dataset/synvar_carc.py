# Calculate clause relation tags for BHSA data in Text-Fabric
# source: Martijn Naaijer and Marianne Kaajan
# https://github.com/MartijnNaaijer/phdthesis/blob/master/Various/main_subordinate_clauses.ipynb

"""
# Main and subordinate clauses
# Martijn Naaijer

For the Syntactic Variation project, Marianne Kaajan developed a procedure which decides whether a clause is a main clause or a subordinate clause. Further, it distinguishes three types of subordinate clauses. The function in_dep_calc() is a Python implementation of this procedure. It requires installation of Text-Fabric.

The function in_dep_calc() can be called with a clause node as argument.

It returns one of the following values:
    * Main (main clause)
    * SubAdv (subordinate adverbial clause)
    * SubArg (subordinate argument clause)
    * SubMod (subordinate attributive clause)
    * Undc (Undecided)

These functions given here are translated from a number of MQL queries made by Marianne Kaajan.
"""

def in_dep_calc(cl, api):  
    """
    Whether a clause is a main clause or one of the types of dependent clauses is found out in several steps.
    An important role is played by the features "code", which retrieves the clause atom relation code (carc), 
    and the feature "rela", which retrieves the clause constituent relation (ccr) of a clause.
    """
      
    F, E, L = api.F, api.E, api.L

    in_dep = ''        
    if F.rela.v(cl) == 'ReSu': # is the clause resumptive?
        moth_obj = E.mother.f(cl)[0]
        in_dep = rela_calc(moth_obj, api)
    else:
        in_dep = rela_calc(cl, api) # does the clause have a dependent ccr?

    # if the previous step did not give a result, check if there is a wayyiqtol verb in the clause
    if in_dep == '':
        words = L.d(cl, 'word') # is there a wayyiqtol?
        for word in words:
            if F.vt.v(word) == 'wayq':
                in_dep += 'Main'
                      
    # if everything else does not give a result, we look at the carc.
    if in_dep == '':  
        cl_atoms = L.d(cl, 'clause_atom')
        in_dep = carc_calc(cl_atoms, api)
        
    return(in_dep)


def carc_calc(cl_atoms, api):
    """
    The values of the feature carc are related to whether a clause is a main clause or a subordinate clause.
    """
    F, E = api.F, api.E

    in_dep_c = ''
    carc = F.code.v(cl_atoms[0])
    if 999 > int(carc) > 499:
        in_dep_c += 'SubAdv'
    elif int(carc) in {0, 999}:
        in_dep_c = 'Main'
    elif 17 > int(carc) > 9:
        in_dep_c += 'SubAdv'
    elif 75 > int(carc) > 50:
        in_dep_c += 'SubAdv'
    elif 168 > int(carc) > 99:
        in_dep_c += 'Main'
    elif 500 > int(carc) > 299:
        in_dep_c += 'Main'
        
    elif int(carc) in {200, 201}:         
        while F.code.v(cl_atoms[0]) in {200, 201}:
            cl_atoms = E.mother.f(cl_atoms[0])
        carc = F.code.v(cl_atoms[0])
        
        if 999 > int(carc) > 499:
            in_dep_c += 'SubAdv'
        elif int(carc) in {0, 999}:
            in_dep_c = 'Main'
        elif 17 > int(carc) > 9:
            in_dep_c += 'SubAdv'
        elif 75 > int(carc) > 50:
            in_dep_c += 'SubAdv'
        elif 168 > int(carc) > 99:
            in_dep_c += 'Main'
        elif 500 > int(carc) > 299:
            in_dep_c += 'Main'
        elif int(carc) in {220, 221, 222, 223}:
            in_dep_c += 'Undc'
        
    else:
        in_dep_c += 'Undc'
        
    return(in_dep_c)


def rela_calc(cl, api):
    """
    In this function the clause constituent relation is retrieved, to be able to distinguish between 
    different kinds of subordinate clauses. In the case of a coordinate clause, we move up to the mother of the
    clause recursively, until we find a clause that is not a coordinate clause.
    
    """

    F, E, L = api.F, api.E, api.L
    
    in_dep_r = ''
    ccr = F.rela.v(cl)
    
    # which function does the clause have in another clause?
    if ccr in {'Subj', 'Objc', 'Cmpl', 'PreC', 'Voct', 'Frnt'}:
        in_dep_r += 'SubArg'
        
    elif ccr in {'Attr', 'RgRc', 'Spec'}:
        in_dep_r += 'SubMod'
        
    elif ccr in {'Adju', 'PrAd'}:
        in_dep_r += 'SubAdv'
        
    # is the clause a coordinate clause?
    elif ccr == 'Coor':
        moth_obj = E.mother.f(cl)[0]
        
        # check if the mother is a word or a phrase
        if F.otype.v(moth_obj) in {'word', 'phrase'}:
            in_dep_r += 'SubMod'
            
        # move up to the mother clause
        else:
            while F.rela.v(moth_obj) == 'Coor':
                moth_obj = E.mother.f(moth_obj)[0]
            ccr = F.rela.v(cl)
            if ccr in {'Subj', 'Objc', 'Cmpl', 'PreC', 'Voct', 'Frnt'}:
                in_dep_r += 'SubArg'
            elif ccr in {'Attr', 'RgRc', 'Spec'}:
                in_dep_r += 'SubMod'
            elif ccr in {'Adju', 'PrAd'}:
                in_dep_r += 'SubAdv'
                
        if in_dep_r == '':
            if F.otype.v(moth_obj) != 'clause':
                in_dep_r += 'SubMod'
            else:
                cl_atoms = L.d(moth_obj, 'clause_atom')
                in_dep_r = carc_calc(cl_atoms, api)
                
    return in_dep_r
