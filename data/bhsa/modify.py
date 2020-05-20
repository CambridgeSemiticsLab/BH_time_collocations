import collections
from tf.fabric import Fabric
from tf.compose import modify
from tenses import convert_tense
from function import modify_function
from funct_assoc import calculate
from transcriptions import transcribe_lexemes

def mod_features(locs, base_metadata, do_assoc=False):
    """Remap node features in BHSA

    BHSA contains features of  nodes that 
    are sometimes inaccurate or not quite
    fit for purpose for this study. In this
    module, we modify the dataset with Text-Fabric
    and export the new version to a new directory.
    """

    # setup data and methods
    bhsa, heads, output = locs['bhsa'], locs['heads'], locs['custom']
    TF = Fabric(locations=[bhsa,heads], silent='deep')
    api = TF.load('''
        function vt mother pdp lex
        nhead ls nu lex_utf8
    ''')
    F = api.F

    # features to be modded
    mod_features = {
        'function2': {n:F.function.v(n) for n in F.otype.s('phrase')},
        'vt': {},
    }

    # remap certain phrase functions
    modify_function(mod_features['function2'], api)

    # modify verb tenses: 
    #   - add weqatal
    #   - change tense codes   
    for verb in api.F.pdp.s('verb'):
        mod_features['vt'][verb] = convert_tense(verb, api)

    # add lex_sbl features
    transcribe_lexemes(mod_features, api)

    meta_data = {}

    # add statistical association features for head-words
    if do_assoc:
        assocs = calculate(mod_features['function2'], api)
        mod_features['top_assoc'] = assocs['top_assoc']
        mod_features['funct_assoc'] = assocs['funct_assoc']
        meta_data.update({
            'funct_assoc': {
            'description':'a feature on words that function as a head in their enclosing phrase; integer tells how attracted the head word is to its phrase\'s functions',
            'interpreting scores':'score > 1.3 is significantly attracted; score < -1.3 is significantly repelled',
            'valueType':'int',
        },
        'top_assoc': {
            'description':'top associated function to this word',
            'interpreting scores':'score > 1.3 is significantly attracted; score < -1.3 is significantly repelled',
            'valueType':'str',
        },
    })

    # data for new features
    meta_data.update({
        '': base_metadata,
        'function2': {
            'description': 'function of a phrase in a clause modified for this project',
            'valueType': 'str',
        },
        'vt': {
            'description': 'tense of a verb',
            'valueType': 'str',
        },
        'lex_sbl': {
            'description': 'A lexeme string for a word/lexeme node, converted from BHSA to SBL style',
            'valueType': 'str',
        },
        'lex_sbl_l': {
            'description': 'A light transcribed (no disambiguators) lexeme string for a word/lexeme node, converted from BHSA to SBL style',
            'valueType': 'str',
        },
        'lex_utf8': {
            'description': 'UTF8 consonantal text without accents; same as BHSA feature but with lexeme nodes',
            'valueType': 'str',
        },
    })

    # enact the changes
    TF = Fabric(locations=output, silent=True)
    TF.save(
        nodeFeatures=mod_features,
        metaData=meta_data,
    ) 
