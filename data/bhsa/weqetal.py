def get_grandma(clause_atom, api):
    """Return a wayyiqtol/yiqtol/imperative from a clause tree.

    Recursively climbs up a clause's ancestorial tree.
    Stops upon identifying either wayyiqtol
    or a yiqtol|impv grand(mother).

    Returns:
        a string of the ancestor's tense, or None.
    """

    F, E, L = api.F, api.E, api.L

    this_verb = next((F.vt.v(w) for w in L.d(clause_atom) if F.pdp.v(w)=='verb'), '') 
    mother = next((m for m in E.mother.f(clause_atom)), 0)
    mom_verb = next((F.vt.v(w) for w in L.d(mother) if F.pdp.v(w)=='verb'), '')    

    if mom_verb in {'wayq', 'impf', 'impv'}:
        return mom_verb
    elif not mother:
        return this_verb
    else:
        return get_grandma(mother, api)

def convert_tense(word, api):
    """Convert weqetal tenses where present.

    If not weqetal, return the tense.
    """
    
    F, L = api.F, api.L
    
    # check for weqetal
    if F.vt.v(word) == 'perf' and F.lex.v(word-1) == 'W':

        # get tense of the ancestor of the verb's clause
        clause = L.u(word, 'clause_atom')[0]
        qatal_ancestor = get_grandma(clause, api)

        # check for whether ancestor triggers weqatal analysis
        if qatal_ancestor in {'impf', 'impv'}:
            return 'weqt' # change tense to weqt
        else:
            return F.vt.v(word) # no change on tense

    # return tense as-is
    else:
        return F.vt.v(word)