"""
Retrieve modifiers from a phrase.
"""

import collections
import tools.nav_tree as nt

def get_modis(ph_parse, API, boolean=True):
    """Iterate through a phrase and retrieve its modifiers."""

    F, L = API.F, API.L

    modifiers = collections.defaultdict(list)

    # get modifier head words from the phrase parses
    if len(ph_parse) == 1:
        head_phrases = [[None, ph_parse[0], None]]
    elif len(ph_parse) == 3:
        head_phrases = list(nt.get_head_path(ph_parse))
    for sp in head_phrases:
        src, tgt, rela = sp
        if src:
            src_head = nt.get_head(src)
            modifiers[rela].append(src_head)

    # get modifiers on head-word
    if F.nu.v(tgt) == 'pl':
        modifiers['PL'].append(tgt)
    elif F.nu.v(tgt) == 'du':
        modifiers['DU'].append(tgt)
    if F.uvf.v(tgt) == 'H':
        modifiers['HLOC'].append(tgt)
    if F.prs.v(tgt) not in {'absent', 'n/a'}:
        person = F.prs_ps.v(tgt).replace('p','')
        gender = F.prs_gn.v(tgt).upper()
        number = F.prs_nu.v(tgt).upper()
        sfx = f'SFX{person}'
        modifiers[sfx].append(tgt) 
        modifiers['SFX'].append(tgt)

    # do modifiers from a semantic perspective
    if 'APPO' in modifiers:
        src = modifiers['APPO'][0]
        lexset = F.ls.v(src)
        pdp = F.pdp.v(src)
        if lexset in {'ordn'}:
            modifiers['ORDN'].append(src)
        elif pdp in {'prde', 'prps'}:
            modifiers['DEMON'].append(src)

    if 'PP' not in modifiers:
        modifiers['ØPP'].append(tgt)

    if boolean:
        for key, val in modifiers.items():
            modifiers[key] = 1 * bool(val)

    return modifiers
