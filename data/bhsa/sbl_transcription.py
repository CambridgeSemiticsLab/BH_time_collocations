"""
This module produces a consonantal
transcription feature which adheres to
the SBL handbook guidelines for 
transcribing Biblical Hebrew
"""

import collections
import re

def transcribe_lexemes(nodeFeatures, api):
    """Convert BHSA lexemes to SBL transcription.

    The BHSA lexeme transcription contains disambiguating
    symbols to distinguish nominal and verbal roots,
    namely "/" and "[". In addition, BHSA distinguishes
    homographs with any number of trailing "="s. How should
    this information be incorporated in a new transcription?
    The new transcription should be interpretable and easy on
    the eyes, but without losing any of the disambiguation
    info. 
    
    To achieve the above goals, this function will convert all
    BHSA symbols to interpretable latin characters. The "/"
    noun disambiguator will become a 'N'; the "[" will become 
    a "V". Homographs will be converted to Arabic numerals if
    present.
    """

    etcbc2sbl = { 
        '>': 'ʾ',
        'B': 'b',
        'G': 'g',
        'D': 'd',
        'H': 'h',
        'W': 'w',
        'Z': 'z',
        'X': 'ḥ',
        'V': 'ṭ',
        'J': 'y',
        'K': 'k',
        'L': 'l',
        'M': 'm',
        'N': 'n',
        'S': 's',
        '<': 'ʿ',
        'P': 'p',
        'Y': 'ṣ',
        'Q': 'q',
        'R': 'r',
        'F': 'ś',
        'C': 's̆',
        'T': 't',
        '_': ' ',
    }

    F, L = api.F, api.L
    nodeFeatures['lex_sbl'] = {}
    features = nodeFeatures['lex_sbl'] 

    for lnode in F.otype.s('lex'):
        bhsa_lex = F.lex.v(lnode)
        raw_lex = re.sub('\[|/|=', '', bhsa_lex)
        sbl_lex = ''.join(etcbc2sbl[c] for c in raw_lex)
        if '/' in bhsa_lex:
            pos = '-N'
        elif '[' in bhsa_lex:
            pos = '-V'
        else:
            pos = ''
        homo = bhsa_lex.count('=')
        homo_count = f'{homo+1}' if homo else ''
        sbl_lex_final = f'{sbl_lex}{pos}{homo_count}'
        features[lnode] = sbl_lex_final
        for w in L.d(lnode,'word'):
            features[w] = sbl_lex_final
