import collections
from tf.fabric import Fabric
from tf.compose import modify
from weqatal import convert_tense
from funct_assoc import calculate

def mod_features(locs, base_metadata):
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
        nhead
    ''')
    F = api.F

    # features to be modded
    mod_features = {
        'function': {n:F.function.v(n) for n in F.otype.s('phrase')},
        'vt': {},
    }

    # manually remap phrase functions
    mod_features['function'].update({
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
    })

    # modify verb tenses to add weqatal
    for verb in api.F.pdp.s('verb'):
        mod_features['vt'][verb] = convert_tense(verb, api)

    # add statistical association features for head-words
    assocs = calculate(mod_features['function'], api)
    mod_features['top_assoc'] = assocs['top_assoc']
    mod_features['funct_assoc'] = assocs['funct_assoc']

    # data for new features
    meta_data = {
        '': base_metadata,
        'function': {
            'description': 'function of a phrase in a clause',
            'valueType': 'str',
        },
        'vt': {
            'description': 'tense of a verb',
            'valueType': 'str',
        },
        'funct_assoc': {
            'description':'a feature on words that function as a head in their enclosing phrase; integer tells how attracted the head word is to its phrase\'s functions',
            'interpreting scores':'score > 1.3 is significantly attracted; score < -1.3 is significantly repelled',
            'valueType':'int',
        },
        'top_assoc': {
            'description':'top associated function to this word',
            'interpreting scores':'score > 1.3 is significantly attracted; score < -1.3 is significantly repelled',
            'valueType':'str',
        }
    }

    # enact the changes
    TF = Fabric(locations=output, silent=True)
    TF.save(
        nodeFeatures=mod_features,
        metaData=meta_data,
    ) 
