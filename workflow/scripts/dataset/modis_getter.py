"""
Retrieve modifiers from a phrase.
"""

import tools.nav_tree as nt

def get_modis(ph_parse, API, boolean=True):
    """Iterate through a phrase and retrieve its modifiers."""

    F, L = API.F, API.L

    modifiers = {
    }

    # get modifier head words from the phrase parses
    if len(ph_parse) == 1:
        head_phrases = [[None, ph_parse[0], None]]
    elif len(ph_parse) == 3:
        head_phrases = list(nt.get_head_path(ph_parse))
    for sp in head_phrases:
        src, tgt, rela = sp
        if src:
            src_head = nt.get_head(src)
            modifiers[rela] = src_head

    # get modifiers on head-word
    if F.nu.v(tgt) == 'pl':
        modifiers['PL'] = tgt
    elif F.nu.v(tgt) == 'du':
        modifiers['DU'] = tgt
    if F.uvf.v(tgt) == 'H':
        modifiers['HLOC'] = tgt
    if F.prs.v(tgt) not in {'absent', 'n/a'}:
        person = F.prs_ps.v(tgt).replace('p','')
        gender = F.prs_gn.v(tgt).upper()
        number = F.prs_nu.v(tgt).upper()
        sfx = f'SFX{person}'
        modifiers[sfx] = tgt 
        modifiers['SFX'] = tgt

    # do modifiers from a semantic perspective
    if 'APPO' in modifiers:
        lexset = F.ls.v(tgt)
        pdp = F.pdp.v(tgt)
        if lexset in {'ordn'}:
            modifiers['ORDN'] = tgt
        elif pdp in {'prde', 'prps'}:
            modifiers['DEMON'] = tgt

    if 'PP' not in modifiers:
        modifiers['ØPP'] = tgt

    if boolean:
        for key, val in modifiers.items():
            modifiers[key] = 1 * bool(val)

    return modifiers
